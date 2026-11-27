"""Shared-kernel multi-RHS Sinkhorn / AA2 baseline and Taylor initializer.
Batching is established numerical engineering, not a novelty claim.
No full per-scenario coupling tensor is materialized for class readouts.
"""
import numpy as np

def batch_solve(K,a,bs,U0,tol=1e-9,memory=2,maxiter=10000):
 # bs shape m x B; source scaling U shape n x B.
 U=U0.copy();n,B=U.shape;active=np.arange(B);fu=np.empty_like(U);fv=np.empty_like(bs);its=np.zeros(B,int);hx=[];hf=[]
 for it in range(1,maxiter+1):
  if memory:
   x=np.log(U)
   if it-1>memory:
    X,F=np.stack(hx,axis=2),np.stack(hf,axis=2);R=F-X
    gram=np.einsum('nbk,nbl->bkl',R,R);norm=np.linalg.norm(gram,axis=(1,2));norm=np.maximum(norm,1e-300)
    co=np.linalg.solve(gram/norm[:,None,None]+.01*np.eye(memory)[None,:,:],np.ones((len(active),memory,1)))[...,0];co/=co.sum(1)[:,None]
    x=np.einsum('nbk,bk->nb',F,co);U=np.exp(x)
   hx.append(x.copy());hx=hx[-memory:]
  V=bs[:,active]/(K.T@U);U=a[:,None]/(K@V)
  if memory:hf.append(np.log(U));hf=hf[-memory:]
  if it%5==0:
   err=np.abs(V*(K.T@U)-bs[:,active]).sum(0)
   if not np.isfinite(err).all():raise FloatingPointError('nonfinite batch iterate')
   done=err<=tol
   if done.any():
    fu[:,active[done]]=U[:,done];fv[:,active[done]]=V[:,done];its[active[done]]=it
    keep=~done;active=active[keep]
    if not len(active):break
    U=U[:,keep];hx=[x[:,keep] for x in hx];hf=[x[:,keep] for x in hf]
 else:raise RuntimeError('batch did not converge')
 residual=np.abs(fv*(K.T@fu)-bs).sum(0)
 return fu,fv,its,residual

def readout(K,a,Y,U,V):
 # m x B x q -> a single GEMM; output n x B x q.
 m,B=V.shape;q=Y.shape[1]
 rhs=(V[:,:,None]*Y[:,None,:]).reshape(m,B*q)
 return (K@rhs).reshape(len(a),B,q)*(U/a[:,None])[:,:,None]
