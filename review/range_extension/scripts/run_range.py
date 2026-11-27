import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,pathlib,json,time,traceback,signal,importlib.util
import numpy as np,pandas as pd
from scipy.spatial.distance import cdist
from scipy.special import softmax
from threadpoolctl import threadpool_limits
R=pathlib.Path(__file__).resolve().parents[1];P=R.parent
spec=importlib.util.spec_from_file_location('v',P/'scripts/02_validation.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
CFG=json.load(open(R/'protocol.json'))
def timeout(*args):raise TimeoutError('240s configuration budget exceeded')
signal.signal(signal.SIGALRM,timeout)
def mapping(a,mode):
 t={'CD4 T','CD8 T','other T','CD4 T cells','CD8 T cells'}
 if mode=='binary_T':return np.array(['T' if x in t else 'non-T' for x in a])
 if mode=='merged':return np.array(['T' if x in t else 'Mono' if x in ['CD14+ Monocytes','FCGR3A+ Monocytes','CD14 Mono','FCGR3A Mono','CD14 mono','FCGR3A mono'] else x for x in a])
 return a
def main(name):
 z=np.load(P/'inputs/original_tasks'/name/'inputs.npz');yr=z['classes'][z['yr']] if np.issubdtype(z['yr'].dtype,np.integer) else z['yr'];yq=z['classes'][z['yq']] if np.issubdtype(z['yq'].dtype,np.integer) else z['yq'];rng=np.random.default_rng(20261106);qi=np.sort(rng.choice(len(yq),1000,replace=False));Cfull=cdist(z['zq'],z['zr'],'sqeuclidean');scale=np.median(Cfull);Cfull=Cfull[qi]/scale
 small=[]
 for cl in np.unique(yr):small.extend(rng.choice(np.flatnonzero(yr==cl),min(2,(yr==cl).sum()),replace=False))
 small=np.sort(np.r_[small,rng.choice(np.setdiff1d(np.arange(len(yr)),small),300-len(small),replace=False)])
 for m,mode in CFG['configurations']:
  key=f'{name}_m{m}_{mode}';folder=R/'results'/key;folder.mkdir(parents=True,exist_ok=True)
  if (folder/'complete.json').exists():continue
  print('START',key,flush=True);t=time.perf_counter();rows=[];signal.alarm(240)
  try:
   ri=small if m==300 else np.arange(len(yr));rl=mapping(yr[ri],mode);truth=mapping(yq[qi],mode);C=Cfull[:,ri];audit=v.CompositionAudit(C,rl,epsilon=.05);audit.run(audit.p0[None,:],order=1,corrector='chord',tol=1e-7);J=v.log_jacobian(audit);P0=audit.base_probabilities;s,_,direction,radius=v.scores(P0,J,C,np.repeat(P0[:,None,:],2,axis=1));base=P0.argmax(1);all_changes=[];checks=[]
   np.savez_compressed(folder/'baseline.npz',classes=audit.classes,probabilities=P0,scores=s,radius=radius,direction=direction,q_id=z['q_id'][qi],q_donor=z['q_donor'][qi],q_index=qi,r_index=ri,author_labels=truth)
   for seed in CFG['seeds']:
    rng=np.random.default_rng(seed);anchors=rng.choice(len(P0),32,replace=False);rd=rng.normal(size=(32,len(audit.classes)));rd-=rd.mean(1,keepdims=True);rd/=np.linalg.norm(rd,axis=1,keepdims=True)
    for rad in CFG['radii']:
     record={}
     for method,dirs in [('directed',direction[anchors]),('random',rd)]:
      ps=softmax(np.log(audit.p0)+rad*dirs,axis=1);tt=time.perf_counter();o=audit.run(ps,order=1,corrector='chord',tol=1e-7);elapsed=time.perf_counter()-tt;pp=o['probabilities'];labels=pp.argmax(2);sw=labels!=base[:,None];record[method]=sw;all_changes.append(sw);maxerr=float(o['marginal_residual'].max());assert maxerr<=1e-7
      np.savez_compressed(folder/f'{seed}_r{rad}_{method}.npz',switched=sw,scenario_labels=labels,proportions=ps,anchors=anchors,marginal_residual=o['marginal_residual'])
      for B in CFG['budgets']:
       hit=sw[:,:B].any(1);rows.append({'task':name,'references':m,'taxonomy':mode,'n_classes':len(audit.classes),'seed':seed,'radius':rad,'budget':B,'method':method,'n_query':len(P0),'hits':int(hit.sum()),'fraction':float(hit.mean()),'anchor_hits':int(sw[anchors[:B],np.arange(B)].sum()),'unique_directions':len(np.unique(np.round(dirs[:B],8),axis=0)),'max_residual':maxerr,'full32_seconds':elapsed,'fallbacks_full32':o['diagnostics']['fallback_count']})
      # Independent solution for one largest-radius directed scenario per configuration.
      if seed==CFG['seeds'][0] and rad==.3 and method=='directed':
       import ot
       gg=ot.sinkhorn(audit.a,audit.W@ps[0],C,.05,numItermax=100000,stopThr=1e-12);pr=gg@audit.Y/audit.a[:,None];res=max(np.abs(gg.sum(0)-audit.W@ps[0]).sum(),np.abs(gg.sum(1)-audit.a).sum());dif=float(abs(pr-pp[:,0]).max());assert res<1e-8 and dif<1e-3;checks.append({'max_probability_difference':dif,'residual':float(res),'label_disagreements':int((pr.argmax(1)!=labels[:,0]).sum())})
     # Both hit sets remain available even when one is worse.
   df=pd.DataFrame(rows);df.to_csv(folder/'counts.csv',index=False)
   review=[]
   for rad in CFG['radii']:
    union=np.zeros(len(P0),bool)
    for seed in CFG['seeds']:
     for method in ['directed','random']:union|=np.load(folder/f'{seed}_r{rad}_{method}.npz')['switched'].any(1)
    for method,score in [('confidence',s[:,0]),('entropy',s[:,1]),('linear_radius',s[:,6])]:
     idx=np.argsort(-score,kind='stable')[:100];review.append({'task':name,'references':m,'taxonomy':mode,'radius':rad,'method':method,'reviewed':100,'sampled_union_sensitive':int(union.sum()),'found_sensitive':int(union[idx].sum()),'fraction_of_sampled_union_found':float(union[idx].sum()/union.sum()) if union.any() else None})
   pd.DataFrame(review).to_csv(folder/'review_priority.csv',index=False)
   meta={'key':key,'seconds':time.perf_counter()-t,'classes':audit.classes.tolist(),'reference_native_classes':np.unique(yr[ri]).tolist(),'base_author_agreement':float((audit.classes[base]==truth).mean()),'POT_checks':checks,'seed_repeats_are_not_independent_biological_replicates':True};(folder/'complete.json').write_text(json.dumps(meta,indent=2));print('DONE',key,round(meta['seconds'],2),'s',flush=True)
  except Exception as e:
   (folder/'failure.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc(),'seconds':time.perf_counter()-t},indent=2));pd.DataFrame(rows).to_csv(folder/'partial_counts.csv',index=False);print('FAILED',key,repr(e),flush=True)
  finally:signal.alarm(0)
if __name__=='__main__':
 import fcntl
 with open(R/'logs'/f'{sys.argv[1]}.lock','w') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  with threadpool_limits(limits=1):main(sys.argv[1])
