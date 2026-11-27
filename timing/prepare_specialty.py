"""Author-processed public immune-cell matrices -> donor-held-out OT tasks.
No query annotation is used for feature fitting or composition scenarios.
"""
from pathlib import Path
import os,json,time,hashlib,gc
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse
from sklearn.decomposition import PCA
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data'/'prepared';OUT.mkdir(exist_ok=True)
SEED=20261102

def norm(X):
    X=X.astype(np.float64).tocsr();total=np.asarray(X.sum(1)).ravel()
    if np.any(total<=0):raise ValueError('empty cell')
    X=sparse.diags(1e4/total)@X;X.data=np.log1p(X.data)
    return X.tocsr()

def prepare(name,X,obs,genes,fold,refmask,qmask,proteins=None,pnames=None,nr=1200,nq=3000):
    t=time.perf_counter();rng=np.random.default_rng(SEED+fold)
    ri=np.flatnonzero(refmask);qi=np.flatnonzero(qmask)
    assert not set(obs.iloc[ri].donor).intersection(obs.iloc[qi].donor)
    classes=np.unique(obs.iloc[ri].type.to_numpy(str));rlabels=obs.iloc[ri].type.to_numpy(str)
    if not set(obs.iloc[qi].type).issubset(classes):raise ValueError('unrepresented query class')
    reserve=np.concatenate([rng.choice(ri[rlabels==c],min(5,np.sum(rlabels==c)),replace=False) for c in classes])
    rsel=np.sort(np.r_[reserve,rng.choice(np.setdiff1d(ri,reserve),min(nr-len(reserve),len(ri)-len(reserve)),replace=False)])
    qsel=np.sort(rng.choice(qi,min(nq,len(qi)),replace=False))
    Xref=norm(X[ri]);mean=np.asarray(Xref.mean(0)).ravel();variance=np.asarray(Xref.power(2).mean(0)).ravel()-mean**2
    detect=np.asarray((X[ri]>0).sum(0)).ravel();variance[detect<10]=-np.inf
    gi=np.argsort(variance)[-min(3000,np.sum(detect>=10)):]
    sd=np.maximum(np.sqrt(np.maximum(variance[gi],0)),.1)
    xr=np.clip((norm(X[rsel])[:,gi].toarray()-mean[gi])/sd,-10,10)
    xq=np.clip((norm(X[qsel])[:,gi].toarray()-mean[gi])/sd,-10,10)
    pca=PCA(n_components=30,svd_solver='randomized',random_state=SEED)
    zr=pca.fit_transform(xr);zq=pca.transform(xq)
    labels=obs.type.to_numpy(str);yr=np.searchsorted(classes,labels[rsel]);yq=np.searchsorted(classes,labels[qsel])
    p0=np.bincount(yr,minlength=len(classes))/len(yr)
    refdonors=sorted(obs.iloc[ri].donor.unique());fractions=[]
    for d in refdonors:
        y=np.searchsorted(classes,obs.loc[refmask & (obs.donor.to_numpy()==d),'type'].to_numpy(str))
        c=np.bincount(y,minlength=len(classes))+.5;fractions.append(c/c.sum())
    fractions=np.array(fractions)
    mix=rng.dirichlet(np.ones(len(refdonors)),256)@fractions
    ps_mild=.75*p0+.25*mix
    ps_full=mix
    dest=OUT/name;dest.mkdir(exist_ok=True)
    extra={}
    if proteins is not None:extra={'protein_q':proteins[qsel].astype(np.float32),'protein_names':pnames}
    np.savez_compressed(dest/'inputs.npz',zr=zr,zq=zq,yr=yr,yq=yq,classes=classes,
         p0=p0,ps_mild=ps_mild,ps_full=ps_full,donor_fractions=fractions,refdonors=np.array(refdonors),
         q_donor=obs.iloc[qsel].donor.to_numpy(str),q_id=obs.index[qsel].to_numpy(str),r_id=obs.index[rsel].to_numpy(str),
         selected_genes=np.array(genes)[gi].astype(str),pca_components=pca.components_,pca_mean=pca.mean_,ref_gene_mean=mean[gi],ref_gene_sd=sd,**extra)
    meta={'name':name,'fold':fold,'n_reference':len(rsel),'n_query':len(qsel),'classes':classes.tolist(),
          'reference_donors':refdonors,'query_donors':sorted(obs.iloc[qsel].donor.unique()),
          'reference_pool':len(ri),'query_pool':len(qi),'ref_counts':np.bincount(yr,minlength=len(classes)).tolist(),
          'query_counts':np.bincount(yq,minlength=len(classes)).tolist(),'selected_genes':len(gi),
          'preprocess_seconds':time.perf_counter()-t,'protein_channels':0 if proteins is None else len(pnames),
          'input_sha256':hashlib.sha256((dest/'inputs.npz').read_bytes()).hexdigest()}
    (dest/'meta.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
    del Xref,xr,xq;gc.collect()

def kang():
    a=ad.read_h5ad(ROOT/'data'/'kang_2018.h5ad');obs=pd.DataFrame({'donor':a.obs.replicate.astype(str),'type':a.obs.cell_type.astype(str),'condition':a.obs.label.astype(str)},index=a.obs_names)
    X=a.X.tocsr();donors=sorted(obs.donor.unique());halves=[donors[::2],donors[1::2]]
    for fold in range(2):
        rm=(obs.donor.isin(halves[1-fold]) & (obs.condition=='ctrl')).to_numpy()
        qm=(obs.donor.isin(halves[fold]) & (obs.condition=='stim')).to_numpy()
        prepare(f'kang_f{fold}',X,obs,a.var_names,fold,rm,qm)

def hao():
    ann=pd.read_csv(ROOT/'data'/'hao_annotations.csv.gz',index_col='barcode');Xs=[];os=[];Ps=[];gn=None;pn=None
    for donor in ['P1','P2','P3','P4']:
        with h5py.File(ROOT/'data'/f'hao_{donor}.h5') as f:
            m=f['matrix'];names=m['features/name'][:].astype(str);types=m['features/feature_type'][:].astype(str)
            X=sparse.csc_matrix((m['data'][:],m['indices'][:].astype(np.int32),m['indptr'][:].astype(np.int32)),shape=tuple(m['shape'][:])).T.tocsr()
            ob=ann.loc[m['barcodes'][:].astype(str)].copy();keep=(~ob['celltype.l1'].isin(['other','unknown'])) & (ob['celltype.l2']!='Doublet')
            ob=ob.loc[keep];X=X[keep.to_numpy()];isrna=types=='Gene Expression'
            if gn is not None:assert np.array_equal(gn,names[isrna]) and np.array_equal(pn,names[~isrna])
            gn=names[isrna];pn=names[~isrna]
            Xs.append(X[:,isrna]);Ps.append(X[:,~isrna].toarray());os.append(ob)
            print('loaded',donor,len(ob),'RNA genes',sum(isrna),'protein',sum(~isrna),flush=True)
            del X;gc.collect()
    X=sparse.vstack(Xs,format='csr');del Xs;gc.collect();P=np.concatenate(Ps);ann=pd.concat(os)
    obs=pd.DataFrame({'donor':ann.donor.astype(str),'type':ann['celltype.l1'].astype(str)},index=ann.index)
    for fold,donor in enumerate(['P1','P2']):
        prepare(f'hao_f{fold}',X,obs,gn,fold,(obs.donor!=donor).to_numpy(),(obs.donor==donor).to_numpy(),P,pn)

if __name__=='__main__':
    import sys
    if len(sys.argv)==1 or sys.argv[1]=='kang':kang()
    if len(sys.argv)==1 or sys.argv[1]=='hao':hao()
