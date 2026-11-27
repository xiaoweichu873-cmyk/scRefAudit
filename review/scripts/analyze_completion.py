"""Post-review analyses of frozen outputs; no model fitting or scenario selection."""
from pathlib import Path
import argparse, json, hashlib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def analyze(data, out):
    out.mkdir(parents=True, exist_ok=True)
    used = []
    rows = []
    for f in sorted((data/'range').glob('*/baseline.npz')):
        b = np.load(f); P = b['probabilities']; score = b['scores']; used.append(f)
        key = f.parent.name
        task, tail = key.split('_m', 1); refs, taxonomy = tail.split('_', 1)
        gap = np.sort(P, axis=1)[:, -1] - np.sort(P, axis=1)[:, -2]
        assert np.allclose(-gap, score[:, 2], atol=1e-12)
        assert len(P) == 1000
        for radius in [.03, .1, .3]:
            target = np.zeros(len(P), bool)
            for seed in [20261106, 20261107, 20261108]:
                ff=f.parent/f'{seed}_r{radius}_random.npz'; used.append(ff)
                z=np.load(ff); target |= z['switched'].any(axis=1)
                assert z['marginal_residual'].max() <= 1e-7
            for method, rank_score in [('confidence',score[:,0]),('entropy',score[:,1]),('margin',-gap),('linear_radius',score[:,6])]:
                idx=np.argsort(-rank_score,kind='stable')[:100]
                rows.append(dict(task=task,references=int(refs),taxonomy=taxonomy,radius=radius,method=method,reviewed=100,n_query=len(P),sensitive_total=int(target.sum()),found_sensitive=int(target[idx].sum())))
    review=pd.DataFrame(rows); review.to_csv(out/'review_margin_all.csv',index=False)
    primary=review[(review.taxonomy=='native')&(review.radius==.1)].pivot(index=['task','references'],columns='method',values='found_sensitive')
    primary['increment_over_margin']=primary.linear_radius-primary.margin
    primary.to_csv(out/'review_margin_primary.csv')
    review.groupby(['taxonomy','radius','method']).found_sensitive.mean().unstack().to_csv(out/'review_margin_by_scope.csv')
    # Read only the frozen held-out HCC probabilities. The kNN is the same 15-neighbor
    # baseline already used to initialize marginal class weights, not a tuned competitor.
    samples=[]; cells=[]; changes=[]; support=[]
    inp=pd.read_csv(data/'tissue'/'HCC_Myeloid'/'input_cells.csv.gz'); used.append(data/'tissue'/'HCC_Myeloid'/'input_cells.csv.gz')
    for f in sorted((data/'tissue').glob('*/*/*.npz')):
        b=np.load(f); meta=json.loads(f.with_suffix('.json').read_text())
        if str(b['role'])!='test': continue
        used += [f,f.with_suffix('.json')]
        group=f.parents[1].name; P=b['baseline_probabilities']; p0=P.mean(axis=0)
        ps=b['active_proportions']; blocks=meta['active_blocks']; B=len(ps)//len(blocks)
        for j,block in enumerate(blocks):
            pp=ps[j*B:(j+1)*B]; tv=.5*np.abs(pp-p0).sum(axis=1); largest=np.abs(pp-p0).max(axis=1)
            for k in range(B): changes.append(dict(group=group,subject=meta['subject'],sample=meta['sample'],radius=block['radius'],method=block['method'],scenario=k,total_variation=float(tv[k]),max_class_change=float(largest[k])))
        if group!='HCC_Myeloid':continue
        classes=b['classes']; truth=b['author_labels']; pred=classes[P.argmax(axis=1)]; knn=classes[b['knn_probabilities'].argmax(axis=1)]
        assert len(truth)==len(set(b['barcodes']))
        row={k:meta[k] for k in ['group','fold','subject','sample','condition']}
        row.update(n_cells=len(P),ot_agreement=float((pred==truth).mean()),knn_agreement=float((knn==truth).mean()),missing_reference_labels=int((~b['reference_supported']).sum()))
        samples.append(row)
        cells.append(pd.DataFrame(dict(subject=meta['subject'],sample=meta['sample'],condition=meta['condition'],author_label=truth,ot_label=pred,knn_label=knn)))
        ref=json.loads((f.parent/'reference.json').read_text()); used.append(f.parent/'reference.json')
        pool=inp[inp.subject.isin(ref['reference_donors'])]
        for label,n in pd.Series(truth).value_counts().items():
            r=pool[pool.author_fine==label]; ri=inp[inp['Cell'].isin(ref['reference_barcodes']) & (inp.author_fine==label)]
            support.append(dict(subject=meta['subject'],sample=meta['sample'],condition=meta['condition'],label=label,query_cells=n,reference_pool_cells=len(r),reference_pool_donors=r.subject.nunique(),sampled_reference_cells=len(ri)))
    sample=pd.DataFrame(samples); sample.to_csv(out/'hcc_sample_agreement.csv',index=False)
    donor=sample.groupby(['subject','condition'])[['ot_agreement','knn_agreement']].mean().reset_index()
    condition=donor.groupby('condition').agg(n_donors=('subject','nunique'),ot_agreement=('ot_agreement','mean'),knn_agreement=('knn_agreement','mean'))
    condition['n_cells']=sample.groupby('condition').n_cells.sum(); condition.to_csv(out/'hcc_condition_agreement.csv')
    cell=pd.concat(cells,ignore_index=True)
    for method in ['ot','knn']:
        cm=cell.groupby(['condition','author_label',method+'_label']).size().rename('n').reset_index().rename(columns={method+'_label':'predicted_label'})
        cm.to_csv(out/f'hcc_{method}_confusion.csv',index=False)
    pd.DataFrame(support).to_csv(out/'hcc_reference_support.csv',index=False)
    composition=pd.crosstab(inp.author_fine,inp.subject); concentration=pd.DataFrame({'total_cells':composition.sum(axis=1),'donors_present':(composition>0).sum(axis=1),'largest_donor_fraction':composition.max(axis=1)/composition.sum(axis=1)})
    concentration.to_csv(out/'hcc_label_donor_concentration.csv')
    pd.crosstab(inp.author_fine,inp.condition).to_csv(out/'hcc_label_condition_counts.csv')
    variation=pd.DataFrame(changes);variation.to_csv(out/'realized_composition_changes.csv',index=False)
    vs=variation.groupby(['group','radius','method']).agg(n_scenarios=('total_variation','size'),median_tv=('total_variation','median'),max_tv=('total_variation','max'),median_max_class_change=('max_class_change','median'),max_class_change=('max_class_change','max'))
    vs.to_csv(out/'realized_composition_summary.csv')
    # Empirical annotation-composition variation is context, not a probability model
    # for the perturbations. Each reference donor is compared with the mean of the
    # donors in that fold; proportions come from pre-cap QC-retained cells.
    emp=[]
    for ip in sorted((data/'tissue').glob('*/input_cells.csv.gz')):
        used.append(ip); ii=pd.read_csv(ip)
        for refp in sorted(ip.parent.glob('fold*/reference.json')):
            used.append(refp); rr=json.loads(refp.read_text()); sub=ii[ii.subject.isin(rr['reference_donors'])]
            fractions=pd.crosstab(sub.subject,sub.author_fine,normalize='index').reindex(columns=rr['reference_classes'],fill_value=0)
            center=fractions.mean(axis=0)
            for donor_id,v in fractions.iterrows():emp.append(dict(group=ip.parent.name,fold=rr['fold'],subject=donor_id,total_variation=float(abs(v-center).sum()/2)))
    ed=pd.DataFrame(emp);ed.to_csv(out/'empirical_reference_composition.csv',index=False)
    ed.groupby('group').total_variation.agg(['count','median','min','max']).to_csv(out/'empirical_reference_summary.csv')
    confusion=pd.read_csv(out/'hcc_ot_confusion.csv'); wrong=confusion[confusion.author_label!=confusion.predicted_label].sort_values('n',ascending=False)
    wrong.to_csv(out/'hcc_error_pairs.csv',index=False)
    summary={'analysis':'post-review diagnostic addition; frozen inputs, targets and scores; no fitting or selection of models',
      'primary_review_means':{c:float(primary[c].mean()) for c in ['confidence','entropy','margin','linear_radius']},
      'increment_over_margin':float(primary.increment_over_margin.mean()),'relative_gain_over_margin':float(primary.linear_radius.mean()/primary.margin.mean()-1),
      'configurations_better_than_margin':int((primary.increment_over_margin>0).sum()),'n_configurations':len(primary),
      'hcc_missing_reference_label_cells':int(sample.missing_reference_labels.sum()),'hcc_test_cells':len(cell),
      'hcc_condition_agreement':condition.reset_index().to_dict('records'),
      'hcc_top_confusions':wrong.head(10).to_dict('records'),
      'primary_radius':.1,'review_budget':100,'query_count':1000,'random_evaluation_scenarios':96,
      'no_new_independent_cohorts':True,'uncertainty':'configurations reuse query cells; no cell-level p values or cohort bootstrap of this eight-configuration table'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    (out/'input_checksums.json').write_text(json.dumps({str(f.relative_to(data)):digest(f) for f in sorted(set(used))},indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--data',type=Path,default=ROOT/'data');parser.add_argument('--out',type=Path,default=ROOT/'results');args=parser.parse_args();analyze(args.data,args.out)
