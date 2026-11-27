"""Sequential isolated-process published-method and genuine-cell scale benchmark."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import sys,time,json,subprocess,traceback,resource,hashlib
import numpy as np,pandas as pd,psutil
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits,threadpool_info
ROOT=Path(__file__).resolve().parent
OUT=ROOT/os.environ.get('SOPC_FOLLOWUP_OUT','results/followup');OUT.mkdir(parents=True,exist_ok=True)
CACHE=ROOT/'.cache_followup';CACHE.mkdir(exist_ok=True)
PROTOCOL=json.loads((ROOT/'protocol_followup.json').read_text())
if (ROOT/'protocol_followup_amendment.json').exists():
    PROTOCOL.update(json.loads((ROOT/'protocol_followup_amendment.json').read_text())['overrides'])

def key(c):return f"{c['task']}_m{c['m']}_e{c['eps']:g}_B{c['B']}"
def ck(c):return f"{c['task']}_m{c['m']}_e{c['eps']:g}"
def read_trial(path):
    row=json.loads(path.read_text())
    monitor_path=path.with_name(path.stem+'_monitor.json')
    if monitor_path.exists():row.update(json.loads(monitor_path.read_text()))
    return row
def prepare(c):
    dest=CACHE/ck(c)
    if (dest/'ready.json').exists():return
    from sopc_ot import solve
    dest.mkdir(exist_ok=True)
    name=c['task']+(f"_m{c['m']}" if c['m']!=1200 else '')
    z=np.load(ROOT/'data/prepared'/name/'inputs.npz')
    orig=np.load(ROOT/'data/prepared'/c['task']/'inputs.npz')
    t=time.perf_counter()
    scale=float(np.median(cdist(orig['zq'],orig['zr'],'sqeuclidean')))
    C=cdist(z['zq'],z['zr'],'sqeuclidean')/scale
    K=np.exp(-C/c['eps']);n,m=K.shape;a=np.full(n,1/n);b=np.full(m,1/m)
    Y=np.eye(len(z['classes']))[z['yr']];p0=b@Y;W=Y*b[:,None]/p0
    preprocessing=time.perf_counter()-t;t=time.perf_counter()
    G,u,v,it,err=solve(K,a,b,tol=PROTOCOL['common_base_tolerance'])
    base_s=time.perf_counter()-t
    arrays=dict(C=np.asfortranarray(C),K=K,G=G,u=u,v=v,a=a,b=b,Y=Y,p0=p0,W=W,ps=z['ps_mild'])
    for name,array in arrays.items():np.save(dest/(name+'.npy'),array)
    info={**c,'n':n,'cost_median_original':scale,'common_preprocess_s':preprocessing,'common_base_s':base_s,'base_iterations':it,'base_residual':err}
    (dest/'ready.json').write_text(json.dumps(info,indent=2));(OUT/(ck(c)+'_setup.json')).write_text(json.dumps(info,indent=2))
    print('PREPARED',ck(c),round(base_s,3),it,flush=True)

def ordering(ps,p0):
    # Greedy nearest-next sequence; each start is the closest solved scenario or base.
    remaining=set(range(len(ps)));done=[];curr=p0;order=[]
    while remaining:
        j=min(remaining,key=lambda i:(float(np.abs(ps[i]-curr).sum()),i))
        candidates=[p0]+[ps[i] for i in done]
        ix=int(np.argmin([np.abs(ps[j]-p).sum() for p in candidates]))
        start=-1 if ix==0 else done[ix-1]
        order.append((j,start));done.append(j);remaining.remove(j);curr=ps[j]
    return order

def published(method,C,K,a,bs,eps,u,v,ps,p0,Y,tol):
    import regot_followup
    from sopc_ot import batch_solve,readout
    t=time.perf_counter();f=eps*np.log(u);g=eps*np.log(v)
    base=np.r_[f+g[-1],g[:-1]-g[-1]];duals={-1:base};order=ordering(ps,p0)
    setup=time.perf_counter()-t;correct=0.;read=0.;resids=[];maxrowerr=0.;pending=[];starts=[];native_resids=[];nonfinite_count=0
    n=len(a);B=len(ps);P=np.empty((n,B,Y.shape[1]))
    its=np.zeros(B,dtype=int)
    for j,start in order:
        tt=time.perf_counter()
        result=regot_followup.solve(C,a,bs[:,j],eps,duals[start],tol=tol,max_iter=40,method=0,splr=method.startswith('SPLR'))
        duals[j]=np.asarray(result.dual).copy();correct+=time.perf_counter()-tt
        tt=time.perf_counter();G=np.asarray(result.plan)
        residual=max(float(np.abs(G.sum(1)-a).sum()),float(np.abs(G.sum(0)-bs[:,j]).sum()))
        if np.isfinite(residual):native_resids.append(residual)
        with np.errstate(invalid='ignore',over='ignore'):Q=G@Y/a[:,None]
        P[:,j,:]=Q
        finite=np.isfinite(Q).all() and np.isfinite(duals[j]).all()
        if not finite:
            nonfinite_count+=1;duals[j]=base.copy()
        if finite and residual<=tol*1.001:
            maxrowerr=max(maxrowerr,float(np.max(np.abs(Q.sum(1)-1))));resids.append(residual)
        else:
            lu=duals[j][:n]/eps
            if np.isfinite(lu).all() and np.ptp(lu)<1300:
                starts.append(np.exp(lu-.5*(lu.max()+lu.min())))
            else:starts.append(u.copy())
            pending.append(j)
        read+=time.perf_counter()-tt;its[j]=int(result.niter)
        del G,result,Q
    if pending:
        tt=time.perf_counter();U,V,nit,err=batch_solve(K,a,bs[:,pending],np.column_stack(starts),tol=tol,memory=2)
        correct+=time.perf_counter()-tt;tt=time.perf_counter();PP=readout(K,a,Y,U,V)
        P[:,pending,:]=PP;read+=time.perf_counter()-tt;resids+=err.tolist()
        maxrowerr=max(maxrowerr,float(np.max(np.abs(PP.sum(2)-1))))
        its[pending]+=nit
    if max(resids)>tol*1.001:raise RuntimeError('Published hybrid failed common marginal tolerance')
    return {'setup_s':setup,'predict_s':0.,'correct_s':correct,'readout_s':read,
       'total_s':time.perf_counter()-t,'max_residual':max(resids),'mean_iterations':float(np.mean(its)),
       'max_iterations':int(max(its)),'fallback_count':len(pending),'native_max_residual':max(native_resids) if native_resids else None,'native_nonfinite_count':nonfinite_count,'linear_solver':'CG','native_iteration_cap':40,'row_probability_error':maxrowerr,
       'status':'ok','published_initializations':order},P

def worker(c,method,rep):
    if (ROOT/'.pause_followup_after_trial').exists() and (ROOT/'.pause_followup_after_trial').read_text().strip()=='pause':
        import signal
        os.kill(os.getppid(),signal.SIGINT)
        return
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    sys.path.insert(0,str(ROOT/'source/regot_followup'))
    import regot_followup
    from run_confirmation import execute
    dest=OUT/key(c);dest.mkdir(exist_ok=True)
    tag=f'{method}_r{rep}';load=lambda name:np.load(CACHE/ck(c)/(name+'.npy'))
    common={x:load(x) for x in ['a','b','u','v','Y','p0','W','ps']}
    a,b,u,v,Y,p0,W,ps=[common[x] for x in ['a','b','u','v','Y','p0','W','ps']]
    ps=ps[:c['B']];bs=W@ps.T;db=W-b[:,None];delta=ps-p0
    if method.startswith(('SSNS','SPLR')):C=load('C');K=load('K')
    else:
        K=load('K');G=load('G') if method!='reuse_AA2' else None
    host_pid=int(os.readlink('/proc/self'))
    row={**c,'method':method,'rep':rep,'pid':os.getpid(),'procfs_pid':host_pid,'baseline_rss_mib':psutil.Process(host_pid).memory_info().rss/2**20,
         'threads':threadpool_info(),'eigen_threads':regot_followup.eigen_threads(),'affinity':sorted(os.sched_getaffinity(0))}
    row['adapter_binary_sha256']=hashlib.sha256(Path(regot_followup.__file__).read_bytes()).hexdigest()
    (dest/(tag+'_phase.json')).write_text(json.dumps({'phase':'running','start':time.time(),'procfs_pid':host_pid}))
    t=time.perf_counter();cpu_start=time.process_time()
    try:
        if method.startswith(('SSNS','SPLR')):
            result,P=published(method,C,K,a,bs,c['eps'],u,v,ps,p0,Y,1e-7)
        else:result,P=execute(method,K,a,b,bs,G,u,v,db,delta,Y,1e-7)
        row.update(result)
    except Exception as e:row.update(status='failed',error=str(e),traceback=traceback.format_exc(),elapsed_s=time.perf_counter()-t)
    row['process_cpu_s']=time.process_time()-cpu_start
    row['process_peak_rss_mib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    (dest/(tag+'_phase.json')).write_text(json.dumps({'phase':'done'}))
    if row['status']=='ok':np.save(dest/(tag+'_probabilities.npy'),P)
    (dest/(tag+'.json')).write_text(json.dumps(row,indent=2))
    print('RESULT',key(c),method,rep,row['status'],round(row.get('total_s',row.get('elapsed_s',0)),3),flush=True)

def monitor(c,method,rep):
    dest=OUT/key(c);dest.mkdir(exist_ok=True);tag=f'{method}_r{rep}'
    path=dest/(tag+'.json')
    if path.exists():return read_trial(path)
    if PROTOCOL.get('skip_repeated_failure'):
        for oldrep in range(rep):
            old=dest/(f'{method}_r{oldrep}.json')
            if old.exists() and json.loads(old.read_text())['status'] in ['failed','timeout']:
                row={**c,'method':method,'rep':rep,'status':'not_run_after_failure','error':'Earlier deterministic trial failed or exhausted the fixed time budget; retained as planned but not run.'}
                path.write_text(json.dumps(row,indent=2));return row
    phase=dest/(tag+'_phase.json')
    if phase.exists():phase.unlink()
    with (dest/(tag+'.log')).open('w') as log:
        proc=subprocess.Popen([sys.executable,__file__,'worker',json.dumps(c),method,str(rep)],stdout=log,stderr=subprocess.STDOUT)
        p=None;peak=0;timed=False;start=time.monotonic();running_start=None;samplecount=0
        while proc.poll() is None:
            try:
                state=json.loads(phase.read_text()) if phase.exists() else {}
                if state.get('phase')=='running':
                    if running_start is None:running_start=time.monotonic()
                    if p is None:p=psutil.Process(state['procfs_pid'])
                    peak=max(peak,p.memory_info().rss/2**20);samplecount+=1
                if (running_start is not None and time.monotonic()-running_start>PROTOCOL['method_timeout_s']) or time.monotonic()-start>PROTOCOL['method_timeout_s']+90:
                    proc.kill();proc.wait();timed=True;break
            except (psutil.NoSuchProcess,json.JSONDecodeError):pass
            time.sleep(.01)
    if path.exists():row=json.loads(path.read_text())
    else:row={**c,'method':method,'rep':rep,'status':'timeout' if timed else 'failed','error':'180 s method budget' if timed else f'exit {proc.returncode}; see log'}
    row['sampled_audit_peak_rss_mib']=peak;row['rss_sample_count']=samplecount
    if 'baseline_rss_mib' in row:row['sampled_audit_increment_mib']=max(0,peak-row['baseline_rss_mib'])
    # Keep parent-owned sampling metadata in its own file, separately from the
    # worker-owned result, so resuming never silently loses those observations.
    metadata={k:row[k] for k in ['sampled_audit_peak_rss_mib','rss_sample_count','sampled_audit_increment_mib'] if k in row}
    path.with_name(path.stem+'_monitor.json').write_text(json.dumps(metadata,indent=2))
    path.write_text(json.dumps(row,indent=2))
    print('COMPLETE',key(c),method,rep,row['status'],round(row.get('total_s',0),3),'peak MiB',round(peak,1),flush=True)
    return row

def main():
    configs=PROTOCOL['main_configs']+PROTOCOL['scale_configs']
    for c in configs:prepare(c)
    env={'python':sys.version,'threads':threadpool_info(),'cpu':Path('/proc/cpuinfo').read_text().split('model name')[1].split('\n')[0],
         'memory_limit':Path('/sys/fs/cgroup/memory.max').read_text().strip(),'protocol_sha256':hashlib.sha256((ROOT/'protocol_followup.json').read_bytes()).hexdigest()}
    if (ROOT/'protocol_followup_amendment.json').exists():env['amendment_sha256']=hashlib.sha256((ROOT/'protocol_followup_amendment.json').read_bytes()).hexdigest()
    env['adapter_provenance_sha256']=hashlib.sha256((ROOT/'source/regot_followup/provenance.json').read_bytes()).hexdigest()
    (OUT/'environment.json').write_text(json.dumps(env,indent=2));rows=[]
    for ci,c in enumerate(configs):
        methods=PROTOCOL['methods']
        for rep in range(PROTOCOL['repeats']):
            shift=(ci+rep)%len(methods)
            for method in methods[shift:]+methods[:shift]:
                row=monitor(c,method,rep);rows.append(row)
                pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(list,dict))} for r in rows]).to_csv(OUT/'timings.csv',index=False)
    (OUT/'complete.json').write_text(json.dumps({'trials':len(rows),'status_counts':pd.Series([r['status'] for r in rows]).value_counts().to_dict()},indent=2))

if __name__=='__main__':
    with threadpool_limits(limits=1):
        if len(sys.argv)>1 and sys.argv[1]=='worker':worker(json.loads(sys.argv[2]),sys.argv[3],int(sys.argv[4]))
        else:main()
