"""Frozen donor-separated validation of sensitivity readouts and active search."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,pathlib,json,time,hashlib,traceback,zlib,signal
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial.distance import cdist
from scipy.special import softmax
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
R=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'previous/source'))
from sopc_ot import CompositionAudit
CFG=json.load(open(R/'protocol.json'))['geometry'];SEED=CFG['seed']
def norm(X):
 z=X.astype(np.float32).tocsr();z.data*=np.repeat(1e4/np.maximum(np.asarray(z.sum(1)).ravel(),1),np.diff(z.indptr));np.log1p(z.data,out=z.data);return z
def log_jacobian(audit):
 D=audit._cache[(1,'chord')][0];q=len(audit.classes);H=(audit.G@((D['g'][:,:,None]*audit.Y[:,None,:]).reshape(len(audit.b),q*q))).reshape(len(audit.a),q,q)/audit.a[:,None,None];H+=D['f'][:,:,None]*audit.base_probabilities[:,None,:];Q=np.diag(audit.p0)-np.outer(audit.p0,audit.p0);return np.einsum('ijk,jl->ilk',H,Q)
def scores(P,J,C,local):
 n,q=P.shape;ii=np.arange(n);base=P.argmax(1);entropy=-(P*np.log(np.maximum(P,1e-300))).sum(1)/np.log(q);h=local.argmax(2);flip=(h[:,1:]!=base[:,None]).mean(1);competitor=local.copy();competitor[ii[:,None],np.arange(local.shape[1])[None,:],base[:,None]]=-np.inf
 worst=(local[ii[:,None],np.arange(local.shape[1])[None,:],base[:,None]]-competitor.max(2)).min(1);g=P[ii,base,None]-P;dg=J[ii,:,base][:,:,None]-J;gn=np.linalg.norm(dg,axis=1);rad=g/np.maximum(gn,1e-14);rad[ii,base]=np.inf;j=rad.argmin(1);radius=rad[ii,j];direction=-dg[ii,:,j]/np.maximum(gn[ii,j,None],1e-14);direction-=direction.mean(1,keepdims=True);direction/=np.maximum(np.linalg.norm(direction,axis=1,keepdims=True),1e-14)
 vals=np.column_stack([1-P.max(1),entropy,-np.diff(np.sort(P,axis=1)[:,-2:],axis=1).ravel(),flip,np.abs(local-P[:,None,:]).max((1,2)),-worst,-np.log(np.maximum(radius,1e-9)),np.partition(C,4,axis=1)[:,:5].mean(1)])
 return vals,['uncertainty','entropy','negative_margin','switch_fraction','max_probability_shift','negative_worst_margin','negative_log_linear_radius','nearest5_distance'],direction,radius
def init_audit(C,rl,props):
 classes,ci=np.unique(rl,return_inverse=True);freq=np.bincount(ci);nb=np.argpartition(C,min(14,len(rl)-1),axis=1)[:,:min(15,len(rl))];w=1/np.maximum(np.take_along_axis(C,nb,axis=1),.02);w*=(props/(freq/len(rl)))[ci[nb]];w/=w.sum(1,keepdims=True);knn=np.stack([(w*(ci[nb]==j)).sum(1) for j in range(len(classes))],axis=1);p0=(knn.sum(0)+.5)/(len(C)+.5*len(classes));audit=CompositionAudit(C,rl,b=p0[ci]/freq[ci],epsilon=CFG['epsilon']);return audit,knn
def make_folds(o,group):
 ids=sorted(o.subject.unique());rng=np.random.default_rng(SEED);out={}
 strata=[ids] if group.startswith('HCC') else [[s for s in ids if ((o.loc[o.subject==s,'condition'].iloc[0]=='Control')==is_control)] for is_control in [True,False]]
 for names in strata:
  for k,s in enumerate(rng.permutation(names)):out[str(s)]=int(k%5)
 return out
def one_group(group,only_fold=None):
 a=ad.read_h5ad(R/'previous/data'/f'{group}_annotated_counts.h5ad');keep=a.obs.passes_QC.to_numpy(copy=True)
 if group.startswith('HCC'):keep &= a.obs.condition.isin(['Normal','Tumor']).values
 a=a[keep].copy();o=a.obs.copy();o['subject']=o.subject.astype(str);o['label']=o.author_fine.astype(str);o['sample']=o['sample'].astype(str) if group.startswith('HCC') else o.subject.astype(str);X=a.X.tocsr();genes=a.var_names;folds=make_folds(o,group);root=R/'results/frozen'/group;root.mkdir(parents=True,exist_ok=True);(root/'donor_folds.json').write_text(json.dumps(folds,indent=2));o.to_csv(root/'input_cells.csv.gz')
 excluded=set(v for g in json.load(open(R/'previous/metadata/program_sets_human.json')).values() for v in g)|{'HAS1','PLIN2','KRT17','KRT5','VEGFA','PPARG','ACTA2'}
 for fold in range(5):
  if only_fold is not None and fold!=only_fold:continue
  fd=root/f'fold{fold}';fd.mkdir(exist_ok=True);test=[s for s,f in folds.items() if f==fold];cal=[s for s,f in folds.items() if f==(fold+1)%5];refs=[s for s,f in folds.items() if f not in [fold,(fold+1)%5]];assert not (set(refs)&set(cal) or set(refs)&set(test) or set(cal)&set(test))
  pool=np.flatnonzero(o.subject.isin(refs));classes=np.sort(o.label.iloc[pool].unique());ro=o.iloc[pool];props=pd.crosstab(ro.subject,ro.label,normalize='index').reindex(columns=classes,fill_value=0).mean().values;avail=np.array([(ro.label==c).sum() for c in classes]);nref=min(CFG['reference_cells'],len(pool));alloc=np.minimum(4+np.floor((nref-4*len(classes))*props).astype(int),avail)
  while alloc.sum()<nref:
   ix=np.flatnonzero(alloc<avail);alloc[ix[np.argmax((nref*props-alloc)[ix])]]+=1
  rng=np.random.default_rng(SEED+fold);ri=[]
  for cl,n in zip(classes,alloc):
   ix=pool[o.label.iloc[pool].values==cl];subjects=o.subject.iloc[ix];w=1/subjects.map(subjects.value_counts()).values;w=w/w.sum();ri.extend(rng.choice(ix,int(n),replace=False,p=w).tolist())
  ri=np.array(ri);rl=o.label.iloc[ri].values;Z=norm(X[ri]);mu=np.asarray(Z.mean(0)).ravel();va=np.asarray(Z.power(2).mean(0)).ravel()-mu**2;eligible=(np.asarray((Z>0).sum(0)).ravel()>=5)&(mu>.005)&~genes.str.match(r'(?i)^(MT-|RPL|RPS)')&~genes.isin(excluded);ix=np.flatnonzero(eligible);bins=pd.qcut(mu[ix],20,labels=False,duplicates='drop');disp=np.log(np.maximum(va,1e-8)/np.maximum(mu,1e-8));sc=np.full(len(genes),-np.inf)
  for b in np.unique(bins):
   j=ix[bins==b];sc[j]=(disp[j]-disp[j].mean())/max(disp[j].std(),1e-6)
  hvg=np.argsort(sc)[-min(CFG['features'],len(ix)):];mat=Z[:,hvg].toarray();mean=mat.mean(0);sd=np.maximum(mat.std(0),.05);pca=PCA(n_components=min(CFG['PCA'],len(ri)-1,len(hvg)-1),svd_solver='randomized',random_state=SEED+fold);rp=pca.fit_transform(np.clip((mat-mean)/sd,-10,10))
  support=ro.groupby('label').subject.nunique();eligible_missing=[c for c in classes if support[c]>=3 and (rl==c).sum()>=10];removed=max(eligible_missing,key=lambda c:int((ro.label==c).sum())) if eligible_missing else None
  plan={'group':group,'fold':fold,'reference_donors':refs,'calibration_donors':cal,'test_donors':test,'reference_barcodes':o.index[ri].tolist(),'features':genes[hvg].tolist(),'reference_classes':list(classes),'controlled_removed_class':removed,'selection_of_removed_class':'most abundant reference-pool class with>=3 reference donors and>=10 sampled references; no query labels used','geometry':CFG};(fd/'reference.json').write_text(json.dumps(plan,indent=2))
  samples=sorted(o.loc[o.subject.isin(cal+test),'sample'].unique())
  for sample in samples:
   role='test' if o[o['sample']==sample].subject.iloc[0] in test else 'calibration';dest=fd/(sample+'.npz');done=fd/(sample+'.json')
   if dest.exists() and done.exists():continue
   t=time.perf_counter();print('START',group,fold,role,sample,flush=True)
   try:
    qi=np.flatnonzero(o['sample'].values==sample);seed=SEED+fold*1000+zlib.crc32(sample.encode())%100000;rng=np.random.default_rng(seed)
    if len(qi)>CFG['query_cap_per_sample']:qi=np.sort(rng.choice(qi,CFG['query_cap_per_sample'],replace=False))
    qmat=norm(X[qi])[:,hvg].toarray();qp=pca.transform(np.clip((qmat-mean)/sd,-10,10));C=cdist(qp,rp,metric='cosine');scale=float(np.median(C));C/=scale;audit,knn=init_audit(C,rl,props);logdraw=np.zeros((17,len(classes)));prng=np.random.default_rng(SEED+fold)
    for k,sdlog in enumerate([.05,.1]):
     logdraw[1+8*k:9+8*k]=prng.normal(0,sdlog,(8,len(classes)))
    logdraw-=logdraw.mean(1,keepdims=True);ps=softmax(np.log(audit.p0)[None,:]+logdraw,axis=1);local=audit.run(ps,order=1,corrector='chord',tol=CFG['tol']);J=log_jacobian(audit);P=audit.base_probabilities;s,names,direction,radius=scores(P,J,C,local['probabilities']);base=P.argmax(1);truth=o.label.iloc[qi].values;supports=np.isin(truth,classes)
    result={'barcodes':np.asarray(o.index[qi].tolist(),str),'author_labels':truth.astype(str),'classes':classes.astype(str),'baseline_probabilities':P,'local_probabilities':local['probabilities'],'scores':s,'score_names':np.asarray(names,str),'knn_probabilities':knn,'reference_supported':supports,'linear_radius':radius,'subject':np.asarray(o.subject.iloc[qi[0]]),'sample':np.asarray(sample),'role':np.asarray(role)}
    meta={'group':group,'fold':fold,'sample':sample,'subject':o.subject.iloc[qi[0]],'role':role,'n_query':len(qi),'condition':str(o.condition.iloc[qi[0]]),'n_reference':len(ri),'baseline_author_agreement':float((classes[base]==truth).mean()),'missing_reference_classes':sorted(set(truth)-set(classes)),'local':local['diagnostics'],'max_residual':float(local['marginal_residual'].max()),'distance_scale':scale}
    if role=='test':
     anchors=rng.choice(len(qi),min(32,len(qi)),replace=False);random=rng.normal(size=(len(anchors),len(classes)));random-=random.mean(1,keepdims=True);random/=np.linalg.norm(random,axis=1,keepdims=True);ds=direction[anchors];assert np.allclose(np.linalg.norm(ds,axis=1),1,atol=1e-6)
     blocks=[];blockinfo=[]
     for rad in [.1,.2]:
      for name,dirs in [('directed',ds),('random',random)]:blocks.append(softmax(np.log(audit.p0)+rad*dirs,axis=1));blockinfo.append({'radius':rad,'method':name})
     active=audit.run(np.concatenate(blocks),order=1,corrector='chord',tol=CFG['tol']);ap=active['probabilities'];result.update(active_probabilities=ap,active_anchors=anchors,active_proportions=np.concatenate(blocks));meta['active_blocks']=blockinfo;meta['active_diagnostics']=active['diagnostics'];meta['active_max_residual']=float(active['marginal_residual'].max());meta['active_detection']=[]
     for k,info in enumerate(blockinfo):
      hp=ap[:,k*len(anchors):(k+1)*len(anchors)].argmax(2);meta['active_detection'].append({**info,'anchors':len(anchors),'anchor_hits':int((hp[anchors,np.arange(len(anchors))]!=base[anchors]).sum()),'all_query_hits':int(np.any(hp!=base[:,None],axis=1).sum())})
     # Independent numerical and derivative checks on the first test query of each group.
     if not (root/'numerical_checks.json').exists():
      import ot
      check=[];d=random[0];h=1e-4;pv=softmax(np.log(audit.p0)[None,:]+np.array([-h,h])[:,None]*d,axis=1);exact=[]
      for qmass in pv:
       coupling=ot.sinkhorn(audit.a,audit.W@qmass,C,CFG['epsilon'],numItermax=100000,stopThr=1e-12);res=max(abs(coupling.sum(0)-audit.W@qmass).sum(),abs(coupling.sum(1)-audit.a).sum());assert res<1e-8;exact.append(coupling@audit.Y/audit.a[:,None])
      fdif=(exact[1]-exact[0])/(2*h);analytic=np.einsum('ijk,j->ik',J,d);rel=np.linalg.norm(fdif-analytic)/max(np.linalg.norm(fdif),1e-12);assert rel<1e-4
      for si in [0,len(anchors),2*len(anchors)]:
       bb=audit.W@np.concatenate(blocks)[si];gg=ot.sinkhorn(audit.a,bb,C,CFG['epsilon'],numItermax=100000,stopThr=1e-12);pp=gg@audit.Y/audit.a[:,None];res=max(abs(gg.sum(0)-bb).sum(),abs(gg.sum(1)-audit.a).sum());dif=float(abs(pp-ap[:,si]).max());assert res<1e-8 and dif<1e-3;check.append({'scenario':int(si),'residual':float(res),'max_class_score_difference':dif,'label_disagreements':int((pp.argmax(1)!=ap[:,si].argmax(1)).sum())})
      (root/'numerical_checks.json').write_text(json.dumps({'sample':sample,'fold':fold,'central_difference_h':h,'relative_derivative_error':float(rel),'independent_POT':check},indent=2))
    if removed is not None:
     take=rl!=removed;rclasses=classes[classes!=removed];rprops=props[classes!=removed];rprops/=rprops.sum();ma,mknn=init_audit(C[:,take],rl[take],rprops);mp=softmax(np.log(ma.p0)+logdraw[:,classes!=removed],axis=1);ml=ma.run(mp,order=1,corrector='chord',tol=CFG['tol']);mj=log_jacobian(ma);ms,_,_,_=scores(ma.base_probabilities,mj,C[:,take],ml['probabilities']);result.update(missing_classes=rclasses.astype(str),missing_baseline_probabilities=ma.base_probabilities,missing_scores=ms,missing_removed_label=np.asarray(removed),controlled_missing_target=(truth==removed),missing_reference_supported=np.isin(truth,rclasses));meta['missing_class_diagnostics']=ml['diagnostics'];meta['missing_class_max_residual']=float(ml['marginal_residual'].max());meta['controlled_removed_class']=removed
    meta['seconds']=time.perf_counter()-t;np.savez_compressed(dest,**result);done.write_text(json.dumps(meta,indent=2));print('DONE',group,fold,role,sample,'n',len(qi),'agreement',round(meta['baseline_author_agreement'],3),'seconds',round(meta['seconds'],2),flush=True)
   except Exception as e:
    failure={'group':group,'fold':fold,'sample':sample,'role':role,'seconds':time.perf_counter()-t,'error':repr(e),'traceback':traceback.format_exc()};(fd/(sample+'.failure.json')).write_text(json.dumps(failure,indent=2));print('FAILED',group,fold,sample,repr(e),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=1):one_group(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else None)
