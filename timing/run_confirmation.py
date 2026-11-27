"""Independent-donor confirmation of fixed-Schur correction.

Protocol 2 was frozen after the Kang-fold0 B16 pilot. All outcomes retained.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,traceback,signal
import numpy as np,pandas as pd,ot
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits
from sopc_ot import derivatives,predict,solve,batch_solve,readout
from sopc_ot.chord import factor_target,chord_solve
from run_specialty import biological
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'/'confirmation';OUT.mkdir(exist_ok=True)
METHODS=['reuse_AA2','quadratic_AA2','reuse_chord','linear_chord','quadratic_chord']

def timeout(signum,frame):raise TimeoutError('exceeded 180 s per solver call')
signal.signal(signal.SIGALRM,timeout)

def execute(method,K,a,b,bs,G,u,v,db,delta,Y,tol):
    t0=time.perf_counter();D=None;ch=None
    if method=='reuse_chord':ch=factor_target(G,a,b)
    elif method!='reuse_AA2':
        D=derivatives(G,a,b,db,side='target',second='quadratic' in method)
        ch=D['factor']
    t1=time.perf_counter();B=bs.shape[1]
    if 'chord' in method:
        if D is None:V0=np.tile(v[:,None],(1,B))
        else:
            predD=dict(D,f=D['g'])
            if 'f2' in D:predD['f2']=D['g2']
            V0=np.column_stack([predict(np.log(v),predD,d,quadratic='f2' in D)[0] for d in delta])
        t2=time.perf_counter();U,V,its,res,info=chord_solve(K,a,bs,V0,ch,tol=tol)
    else:
        if D is None:U0=np.tile(u[:,None],(1,B))
        else:U0=np.column_stack([predict(np.log(u),D,d)[0] for d in delta])
        t2=time.perf_counter();U,V,its,res=batch_solve(K,a,bs,U0,tol=tol,memory=2);info={'fallback_count':0,'evaluations':[],'backtracks':[]}
    t3=time.perf_counter();P=readout(K,a,Y,U,V);t4=time.perf_counter()
    row={'method':method,'setup_s':t1-t0,'predict_s':t2-t1,'correct_s':t3-t2,'readout_s':t4-t3,'total_s':t4-t0,
         'max_residual':float(res.max()),'mean_iterations':float(its.mean()),'max_iterations':int(its.max()),
         'fallback_count':info['fallback_count'],'mean_chord_evaluations':float(np.mean(info['evaluations'])) if info['evaluations'] else 0.,
         'row_probability_error':float(np.max(np.abs(P.sum(2)-1))),'status':'ok'}
    if res.max()>tol*1.01 or not np.isfinite(P).all():raise AssertionError('invalid output')
    return row,P

def task(name,eps,strict=False,scale=False):
    tol=1e-9 if strict else 1e-7
    tag=f'{name}_eps{eps:g}'+('_strict' if strict else '')
    out=OUT/tag;out.mkdir(exist_ok=True)
    if (out/'complete.json').exists():return
    z=np.load(ROOT/'data'/'prepared'/name/'inputs.npz');t=time.perf_counter()
    C=cdist(z['zq'],z['zr'],'sqeuclidean');C/=np.median(C);K=np.exp(-C/eps);n,m=K.shape
    a=np.full(n,1/n);b=np.full(m,1/m);Y=np.eye(len(z['classes']))[z['yr']];p0=b@Y;W=Y*b[:,None]/p0;db=W-b[:,None]
    kernel_s=time.perf_counter()-t;t=time.perf_counter();G,u,v,it,err=solve(K,a,b,tol=1e-11);base_s=time.perf_counter()-t
    base=G@Y/a[:,None];meta={'name':name,'eps':eps,'tol':tol,'n':n,'m':m,'q':Y.shape[1],'kernel_s':kernel_s,'base_s':base_s,
         'base_iterations':it,'base_residual':err,'base_accuracy':float(np.mean(base.argmax(1)==z['yq']))}
    (out/'setup.json').write_text(json.dumps(meta,indent=2));print('SETUP',tag,json.dumps(meta),flush=True)
    settings=[('mild',16),('mild',64),('mild',256),('full',256)]
    if strict:settings=[('mild',64),('full',64)]
    if scale:settings=[('mild',256)]
    rows=[];valid=pd.read_csv(out/'validation.csv').to_dict('records') if (out/'validation.csv').exists() else [];fail=[]
    for work,B in settings:
        key=f'{work}_{B}';dest=out/f'timings_{key}.csv'
        if (out/f'done_{key}.json').exists():continue
        ps=z[f'ps_{work}'][:B];bs=W@ps.T;delta=ps-p0;pref=None;last=None
        config=pd.read_csv(dest).to_dict('records') if dest.exists() else []
        for rep in range(3):
            order=METHODS[rep:]+METHODS[:rep]
            for method in order:
                common={'task':name,'eps':eps,'tol':tol,'workload':work,'B':B,'rep':rep}
                if any(r['rep']==rep and r['method']==method for r in config):continue
                if any(r['method']==method and r['status']=='failed' for r in config):
                    config.append({**common,'method':method,'status':'not_run_after_failure','error':'Deterministic batch previously exceeded time budget or failed; repeat omitted.'})
                    pd.DataFrame(config).to_csv(dest,index=False);continue
                try:
                    signal.alarm(180);row,P=execute(method,K,a,b,bs,G,u,v,db,delta,Y,tol);signal.alarm(0)
                    if pref is None:pref=P.copy()
                    row['max_probability_difference']=float(np.max(np.abs(P-pref)))
                    row['label_disagreement']=float(np.mean(P.argmax(2)!=pref.argmax(2)))
                    if method=='quadratic_chord':last=P.copy()
                    row.update(common);config.append(row)
                    print('TIME',tag,key,rep,method,round(row['total_s'],3),'fallback',row['fallback_count'],flush=True)
                except Exception as e:
                    signal.alarm(0);row={**common,'method':method,'status':'failed','error':str(e)};config.append(row);fail.append(row)
                    print('FAIL',json.dumps(row),flush=True)
                pd.DataFrame(config).to_csv(dest,index=False)
        if last is None:last=pref
        if last is None:
            # A fully censored stress configuration is a reported outcome.
            (out/f'done_{key}.json').write_text(json.dumps({'rows':len(config),'successful':0,'audit_available':False,'reason':'All attempted methods failed or exceeded the budget; no biological readout inferred.'},indent=2))
            rows+=config;continue
        for idx in sorted(set([0,B//2,B-1])):
            ref,lg=ot.bregman.sinkhorn_knopp(a,bs[:,idx],C,eps,numItermax=50000,stopThr=1e-12/np.sqrt(m),log=True,warn=False)
            err=max(np.abs(ref.sum(1)-a).sum(),np.abs(ref.sum(0)-bs[:,idx]).sum())
            ck={'task':name,'eps':eps,'tol':tol,'workload':work,'B':B,'scenario':idx,'pot_residual':float(err),
               'max_probability_difference':float(np.max(np.abs(ref@Y/a[:,None]-last[:,idx,:]))),
               'label_disagreement':float(np.mean((ref@Y).argmax(1)!=last[:,idx,:].argmax(1)))}
            if err>1e-10:raise AssertionError(ck)
            valid.append(ck)
        pd.DataFrame(valid).to_csv(out/'validation.csv',index=False)
        if B==256 or strict:
            bio=biological(z,base,last,work,eps,out);print('BIO',tag,json.dumps(bio),flush=True)
        # Official POT integration timing is bounded; timeouts are censored, never called speedups.
        if B==16 and not strict and not scale:
            pots=[]
            for rep in range(3):
                common={'task':name,'eps':eps,'workload':work,'B':B,'rep':rep,'method':'POT_batch'}
                try:
                    signal.alarm(180);t=time.perf_counter()
                    _,lg=ot.bregman.sinkhorn_knopp(a,bs,C,eps,numItermax=50000,stopThr=tol/np.sqrt(m),log=True,warn=False,
                         warmstart=(np.tile(np.log(u)[:,None],(1,B)),np.tile(np.log(v)[:,None],(1,B))))
                    P=readout(K,a,Y,lg['u'],lg['v']);sec=time.perf_counter()-t;signal.alarm(0)
                    res=np.abs(lg['v']*(K.T@lg['u'])-bs).sum(0).max()
                    if res>tol*1.01:raise AssertionError('POT residual failure')
                    pots.append({**common,'total_s':sec,'iterations':lg['niter']+1,'max_residual':res,'status':'ok'})
                    print('POT',tag,rep,round(sec,3),flush=True)
                except Exception as e:
                    signal.alarm(0);pots.append({**common,'status':'failed','error':str(e)})
                    print('POT FAILED',tag,str(e),flush=True)
                    if isinstance(e,TimeoutError):break
                pd.DataFrame(pots).to_csv(out/'POT_batch.csv',index=False)
            pd.DataFrame(pots).to_csv(out/'POT_batch.csv',index=False)
        (out/f'done_{key}.json').write_text(json.dumps({'rows':len(config),'failures':sum(r['status']!='ok' for r in config)},indent=2))
        rows+=config
    (out/'complete.json').write_text(json.dumps({'task':tag,'utc':pd.Timestamp.now(tz='UTC').isoformat(),'failures':fail},indent=2))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--tasks',nargs='+',default=['kang_f0','kang_f1','hao_f0','hao_f1']);p.add_argument('--eps',type=float,nargs='+',default=[.05,.1]);p.add_argument('--strict',action='store_true');p.add_argument('--scale',action='store_true');p.add_argument('--threads',type=int,default=1);args=p.parse_args()
    from threadpoolctl import threadpool_info
    with threadpool_limits(limits=args.threads):
        (OUT/'thread_environment.json').write_text(json.dumps({'threads':args.threads,'threadpools':threadpool_info(),'cpu_affinity':sorted(os.sched_getaffinity(0))},indent=2))
        for name in args.tasks:
            for eps in args.eps:
                try:task(name,eps,args.strict,args.scale)
                except Exception as e:
                    signal.alarm(0);rec={'task':name,'eps':eps,'error':str(e),'traceback':traceback.format_exc()};print('TASK_FAILED',json.dumps(rec),flush=True)
                    (OUT/f'failure_{name}_{eps}.json').write_text(json.dumps(rec,indent=2))
