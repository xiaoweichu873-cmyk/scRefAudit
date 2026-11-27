"""Add genuine nested reference cells using the already frozen RNA embedding."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import json,hashlib,gc
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse
from prepare_specialty import norm
ROOT=Path(__file__).resolve().parent

def prepare(name,X,obs,genes):
    z=dict(np.load(ROOT/'data/prepared'/name/'inputs.npz'))
    donors=z['refdonors'].astype(str)
    eligible=obs.donor.isin(donors).to_numpy()
    if 'condition' in obs:eligible&=(obs.condition=='ctrl').to_numpy()
    assert not set(donors).intersection(z['q_donor'])
    old=set(z['r_id']); pool=np.flatnonzero(eligible&~obs.index.isin(old))
    assert len(old)==1200 and set(z['r_id']).issubset(set(obs.index[eligible]))
    rng=np.random.default_rng(20261102+4701+(0 if name.startswith('kang') else 1))
    extra=rng.choice(pool,3600,replace=False)
    genes=np.asarray(genes).astype(str)
    assert len(np.unique(genes))==len(genes)
    idx=pd.Index(genes).get_indexer(z['selected_genes']);assert np.all(idx>=0)
    normalized=norm(X[extra])[:,idx].toarray()
    x=np.clip((normalized-z['ref_gene_mean'])/z['ref_gene_sd'],-10,10)
    emb=(x-z['pca_mean'])@z['pca_components'].T
    labs=np.searchsorted(z['classes'],obs.iloc[extra].type.to_numpy(str))
    assert np.array_equal(z['classes'][labs],obs.iloc[extra].type.to_numpy(str))
    for m in [2400,4800]:
        k=m-1200;zz=dict(z)
        zz['zr']=np.vstack([z['zr'],emb[:k]])
        zz['yr']=np.r_[z['yr'],labs[:k]]
        zz['r_id']=np.r_[z['r_id'],obs.index[extra[:k]].to_numpy(str)]
        assert len(np.unique(zz['r_id']))==m
        assert not set(zz['r_id']).intersection(z['q_id'])
        zz['p0']=np.bincount(zz['yr'],minlength=len(z['classes']))/m
        zz['ps_mild']=.75*zz['p0']+.25*z['ps_full']
        dest=ROOT/'data/prepared'/f'{name}_m{m}';dest.mkdir(exist_ok=True)
        np.savez_compressed(dest/'inputs.npz',**zz)
        meta={'source_task':name,'n_reference':m,'n_query':len(z['zq']),
              'unique_references':len(np.unique(zz['r_id'])),'original_references_unchanged':True,
              'query_cells_unchanged':True,'new_cells_projected_with_original_reference_fit':True,
              'reference_donors':donors.tolist(),'reference_counts':np.bincount(zz['yr']).tolist(),
              'sha256':hashlib.sha256((dest/'inputs.npz').read_bytes()).hexdigest()}
        (dest/'meta.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)

def run():
    a=ad.read_h5ad(ROOT/'data/kang_2018.h5ad')
    obs=pd.DataFrame({'donor':a.obs.replicate.astype(str),'type':a.obs.cell_type.astype(str),
                      'condition':a.obs.label.astype(str)},index=a.obs_names)
    prepare('kang_f0',a.X.tocsr(),obs,a.var_names);del a,obs;gc.collect()
    ann=pd.read_csv(ROOT/'data/hao_annotations.csv.gz',index_col='barcode');xs=[];obs=[];genes=None
    for donor in ['P2','P3','P4']:
        with h5py.File(ROOT/'data'/f'hao_{donor}.h5') as f:
            a=f['matrix'];isrna=a['features/feature_type'][:].astype(str)=='Gene Expression'
            gn=a['features/name'][:].astype(str)[isrna]
            if genes is not None:assert np.array_equal(genes,gn)
            genes=gn
            x=sparse.csc_matrix((a['data'][:],a['indices'][:].astype(np.int32),a['indptr'][:].astype(np.int32)),shape=tuple(a['shape'][:])).T.tocsr()
            ob=ann.loc[a['barcodes'][:].astype(str)].copy()
            keep=(~ob['celltype.l1'].isin(['other','unknown']))&(ob['celltype.l2']!='Doublet')
            ob=ob.loc[keep]
            xs.append(x[keep.to_numpy()][:,isrna]);obs.append(pd.DataFrame({'donor':ob.donor.astype(str),'type':ob['celltype.l1'].astype(str)},index=ob.index))
    x=sparse.vstack(xs,format='csr');del xs;gc.collect()
    prepare('hao_f0',x,pd.concat(obs),genes)

if __name__=='__main__':run()
