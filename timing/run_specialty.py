"""Prespecified, single-thread CPU specialty benchmark with resumable artifacts."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,traceback,platform,importlib.metadata
ROOT=Path(__file__).resolve().parent
import numpy as np,pandas as pd
from scipy.spatial.distance import cdist
from sklearn.metrics import balanced_accuracy_score,roc_auc_score
from threadpoolctl import threadpool_limits,threadpool_info
import ot
from sopc_ot import derivatives,predict,solve,batch_solve,readout

OUT=ROOT/'results'/'specialty';OUT.mkdir(exist_ok=True)
METHODS=['reuse','linear','quadratic']

def make_init(u,D,delta):
    if D is None:return np.tile(u[:,None],(1,len(delta))),0
    U=[];bad=0
    for d in delta:
        v,flag=predict(np.log(u),D,d,quadratic='f2' in D);U.append(v);bad+=bool(flag)
    return np.column_stack(U),bad

def biological(z,base,P,work,eps,out):
    pbase=base.argmax(1);pred=P.argmax(2);y=z['yq'];err=pbase!=y
    flip=(pred!=pbase[:,None]).mean(1);span=(P.max(1)-P.min(1)).max(1)
    change=np.abs(P-base[:,None,:]).max((1,2));entropy=-(base*np.log(np.maximum(base,1e-300))).sum(1)
    top=max(1,int(np.ceil(.1*len(y))))
    def auc(x):return float(roc_auc_score(err,x)) if np.unique(err).size==2 else None
    def toprate(x):return float(err[np.argsort(-x,kind='stable')[:top]].mean())
    metrics={'eps':eps,'workload':work,'n':len(y),'base_accuracy':float(np.mean(~err)),
        'base_balanced_accuracy':float(balanced_accuracy_score(y,pbase)),
        'switch_any_fraction':float(np.mean(flip>0)),'mean_switch_frequency':float(flip.mean()),
        'mean_probability_span':float(span.mean()),'median_max_probability_change':float(np.median(change)),
        'max_probability_change':float(change.max()),'scenario_accuracy_min':float((pred==y[:,None]).mean(0).min()),
        'scenario_accuracy_max':float((pred==y[:,None]).mean(0).max()),
        'error_auc_flip':auc(flip),'error_auc_change':auc(change),'error_auc_entropy':auc(entropy),
        'error_rate_top10_flip':toprate(flip),'error_rate_top10_change':toprate(change),'error_rate_top10_entropy':toprate(entropy),
        'switcher_error_rate':float(err[flip>0].mean()) if (flip>0).any() else None,
        'stable_error_rate':float(err[flip==0].mean()) if (flip==0).any() else None}
    df=pd.DataFrame({'cell':z['q_id'],'donor':z['q_donor'],'author_type':z['classes'][y],
       'base_prediction':z['classes'][pbase],'base_error':err,'switch_frequency':flip,'probability_span':span,
       'max_probability_change':change,'entropy':entropy})
    df.to_csv(out/f'cells_{work}.csv',index=False)
    np.savez_compressed(out/f'biology_{work}.npz',base=base,scenario_accuracy=(pred==y[:,None]).mean(0),
        scenario_pred=pred.astype(np.int8),flip=flip,span=span,max_change=change)
    df.groupby('author_type').agg(n=('cell','size'),switch_any=('switch_frequency',lambda x:float(np.mean(x>0))),
       accuracy=('base_error',lambda x:1-float(x.mean())),mean_span=('probability_span','mean')).to_csv(out/f'types_{work}.csv')
    (out/f'biology_{work}.json').write_text(json.dumps(metrics,indent=2));return metrics

def run_task(name,eps,tol=1e-7,strict=False):
    tag=f'{name}_eps{eps:g}'+('_strict' if strict else '')
    out=OUT/tag;out.mkdir(exist_ok=True)
    if (out/'complete.json').exists():print('SKIP',tag,flush=True);return
    z=np.load(ROOT/'data'/'prepared'/name/'inputs.npz',allow_pickle=False)
    t=time.perf_counter();C=cdist(z['zq'],z['zr'],'sqeuclidean');C/=np.median(C);K=np.exp(-C/eps)
    n,m=K.shape;a=np.full(n,1/n);b=np.full(m,1/m);Y=np.eye(len(z['classes']))[z['yr']]
    p0=b@Y;W=Y*b[:,None]/p0;db=W-b[:,None]
    kernel_s=time.perf_counter()-t
    t=time.perf_counter();G,u,v,it,err=solve(K,a,b,tol=1e-11,memory=2);base_s=time.perf_counter()-t
    base=G@Y/a[:,None]
    meta={'name':name,'eps':eps,'tol':tol,'n':n,'m':m,'q':Y.shape[1],'kernel_seconds':kernel_s,
       'base_seconds':base_s,'base_iterations':it,'base_residual':err,'min_kernel':float(K.min()),
       'base_accuracy':float(np.mean(base.argmax(1)==z['yq']))}
    (out/'setup.json').write_text(json.dumps(meta,indent=2));print('SETUP',tag,json.dumps(meta),flush=True)
    rows=[];checks=[]
    settings=[('mild',16),('mild',64),('mild',256),('full',256)] if not strict else [('mild',256),('full',256)]
    for work,B in settings:
        key=f'{work}_{B}';tsfile=out/f'timing_{key}.csv'
        if tsfile.exists():
            rows+=pd.read_csv(tsfile).to_dict('records');print('SKIP CONFIG',tag,key,flush=True);continue
        ps=z[f'ps_{work}'][:B];bs=W@ps.T;delta=ps-p0
        cur=[];Pcheck=None;pref=None;lastP=None
        for rep in range(3):
            order=METHODS[rep:]+METHODS[:rep]
            for method in order:
                t0=time.perf_counter()
                D=None if method=='reuse' else derivatives(G,a,b,db,second=method=='quadratic')
                t1=time.perf_counter();U0,fallback=make_init(u,D,delta);t2=time.perf_counter()
                U,V,its,res=batch_solve(K,a,bs,U0,tol=tol,memory=2);t3=time.perf_counter()
                P=readout(K,a,Y,U,V);t4=time.perf_counter()
                row={'task':name,'eps':eps,'tol':tol,'workload':work,'B':B,'rep':rep,'method':method,
                     'setup_s':t1-t0,'predict_s':t2-t1,'correct_s':t3-t2,'readout_s':t4-t3,'total_s':t4-t0,
                     'max_residual':float(res.max()),'mean_iterations':float(its.mean()),'max_iterations':int(its.max()),'fallback_count':fallback}
                if pref is None:pref=P.copy()
                row['max_probability_difference_reuse']=float(np.max(np.abs(P-pref)))
                row['label_disagreement_reuse']=float(np.mean(P.argmax(2)!=pref.argmax(2)))
                row['row_probability_error']=float(np.max(np.abs(P.sum(2)-1)))
                if row['max_residual']>tol*1.01 or not np.isfinite(P).all():raise AssertionError(row)
                cur.append(row);print('TIMING',tag,key,rep,method,round(row['total_s'],4),round(row['mean_iterations'],2),flush=True)
                if method=='quadratic':lastP=P.copy()
        # Three independent POT solves validate absolute annotation probabilities.
        for idx in sorted(set([0,B//2,B-1])):
            ref,lg=ot.bregman.sinkhorn_knopp(a,bs[:,idx],C,eps,numItermax=50000,stopThr=1e-12/np.sqrt(m),log=True,warn=False)
            r=max(np.abs(ref.sum(1)-a).sum(),np.abs(ref.sum(0)-bs[:,idx]).sum())
            if r>1e-10:raise AssertionError(('POT residual',r))
            Q=ref@Y/a[:,None]
            ck={'task':name,'eps':eps,'tol':tol,'workload':work,'B':B,'scenario':idx,'pot_residual':float(r),
                'max_probability_difference':float(np.max(np.abs(Q-lastP[:,idx,:]))),
                'label_disagreement':float(np.mean(Q.argmax(1)!=lastP[:,idx,:].argmax(1)))}
            checks.append(ck)
        pd.DataFrame(cur).to_csv(tsfile,index=False);rows+=cur
        pd.DataFrame(checks).to_csv(out/'pot_validation.csv',index=False)
        if B==256:
            bio=biological(z,base,lastP,work,eps,out);print('BIOLOGY',tag,json.dumps(bio),flush=True)
        # Official shared-cost POT implementation: integration comparator, one timing per configuration.
        # Aggregate Frobenius tolerance <= tol/sqrt(m) guarantees each marginal L1 <= tol.
        # Its mandatory kernel construction and objective readout are retained and disclosed.
        if B==256 and not strict:
            potrows=[]
            for rep in range(3):
                start=time.perf_counter()
                _,lg=ot.bregman.sinkhorn_knopp(a,bs,C,eps,numItermax=10000,stopThr=tol/np.sqrt(m),log=True,warn=False,
                   warmstart=(np.tile(np.log(u)[:,None],(1,B)),np.tile(np.log(v)[:,None],(1,B))))
                P=readout(K,a,Y,lg['u'],lg['v']);end=time.perf_counter()
                residual=np.abs(lg['v']*(K.T@lg['u'])-bs).sum(0).max()
                if residual>tol*1.01:raise AssertionError(('POT batch',residual))
                potrows.append({'task':name,'eps':eps,'workload':work,'B':B,'rep':rep,'method':'POT_batch',
                    'total_s':end-start,'max_residual':float(residual),'iterations':lg['niter']+1,
                    'max_probability_difference':float(np.max(np.abs(P-lastP)))})
                print('POT',tag,key,rep,round(end-start,3),lg['niter']+1,flush=True)
            pd.DataFrame(potrows).to_csv(out/f'pot_batch_{key}.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'timings.csv',index=False)
    (out/'complete.json').write_text(json.dumps({'task':tag,'rows':len(rows),'utc':pd.Timestamp.now(tz='UTC').isoformat()},indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--task',nargs='+');parser.add_argument('--eps',nargs='+',type=float,default=[.05,.1]);parser.add_argument('--strict',action='store_true');args=parser.parse_args()
    env={'python':platform.python_version(),'platform':platform.platform(),'processor':platform.processor(),
         'packages':{k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','scikit-learn','POT','anndata','h5py']},'threadpools':threadpool_info()}
    (OUT/'environment.json').write_text(json.dumps(env,indent=2))
    tasks=args.task or ['kang_f0','kang_f1','hao_f0','hao_f1'];fails=[]
    with threadpool_limits(limits=1):
        for task in tasks:
            for eps in args.eps:
                try:run_task(task,eps,tol=1e-9 if args.strict else 1e-7,strict=args.strict)
                except Exception as e:
                    rec={'task':task,'eps':eps,'strict':args.strict,'error':str(e),'traceback':traceback.format_exc()};fails.append(rec);print('FAILED',json.dumps(rec),flush=True)
                    (OUT/('failures_strict.json' if args.strict else 'failures.json')).write_text(json.dumps(fails,indent=2))
