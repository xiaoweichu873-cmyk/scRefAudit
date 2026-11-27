"""Exploratory fixed-Schur chord corrector for repeated target marginals.

Classical modified Newton/chord iteration, with shared base factorization.
Not a new mathematical method. A residual-based line search and AA2 fallback
protect the output accuracy; the speed depends on distance from the base.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from .batch import batch_solve

def factor_target(G,a,b):
    return cho_factor(np.diag(b)-(G.T/a)@G+np.outer(b,b),lower=True,check_finite=False)

def chord_solve(K,a,bs,V0,chol,tol=1e-7,maxiter=40,max_backtrack=12):
    B=bs.shape[1];V=V0.copy();active=np.arange(B);fu=np.empty((len(a),B));fv=np.empty_like(bs)
    its=np.zeros(B,int);fallback=[];evaluations=np.zeros(B,int);backtracks=np.zeros(B,int)
    def evaluate(v,idx):
        u=a[:,None]/(K@v);r=bs[:,idx]-v*(K.T@u)
        return u,r,np.abs(r).sum(0)
    U,R,err=evaluate(V,active);evaluations+=1
    for it in range(maxiter+1):
        done=err<=tol
        if done.any():
            fu[:,active[done]]=U[:,done];fv[:,active[done]]=V[:,done];its[active[done]]=it
            keep=~done;active=active[keep];V=V[:,keep];U=U[:,keep];R=R[:,keep];err=err[keep]
        if not len(active):break
        if it==maxiter:break
        dg=cho_solve(chol,R,check_finite=False)
        dg-=np.mean(dg,axis=0,keepdims=True)
        # Backtracking begins with a <= 4 log-unit update. No clipping of dg.
        step=np.minimum(1.,4./np.maximum(np.max(np.abs(dg),axis=0),1e-300))
        good=np.zeros(len(active),bool);newV=V.copy();newU=U.copy();newR=R.copy();newerr=err.copy()
        for bt in range(max_backtrack):
            ix=np.flatnonzero(~good)
            if not len(ix):break
            v=V[:,ix]*np.exp(dg[:,ix]*step[ix])
            # Multiplicative gauge normalization keeps dynamic range safe.
            lv=np.log(v);v=np.exp(lv-np.mean(lv,axis=0,keepdims=True))
            u,r,e=evaluate(v,active[ix]);evaluations[active[ix]]+=1
            accept=np.isfinite(e)&(e<=err[ix]*(1-1e-4*step[ix]))
            accepted=ix[accept];good[accepted]=True
            newV[:,accepted]=v[:,accept];newU[:,accepted]=u[:,accept];newR[:,accepted]=r[:,accept];newerr[accepted]=e[accept]
            rejected=ix[~accept];step[rejected]*=.5;backtracks[active[rejected]]+=1
        if not good.all():
            ix=np.flatnonzero(~good);ids=active[ix]
            u,v,nit,_=batch_solve(K,a,bs[:,ids],U[:,ix],tol=tol,memory=2)
            fu[:,ids]=u;fv[:,ids]=v;its[ids]=nit+it;fallback+=ids.tolist()
            active=active[good];newV=newV[:,good];newU=newU[:,good];newR=newR[:,good];newerr=newerr[good]
        V,U,R,err=newV,newU,newR,newerr
        if not len(active):break
    if len(active):
        u,v,nit,_=batch_solve(K,a,bs[:,active],U,tol=tol,memory=2)
        fu[:,active]=u;fv[:,active]=v;its[active]=nit+maxiter;fallback+=active.tolist()
    residual=np.abs(fv*(K.T@fu)-bs).sum(0)
    if not np.isfinite(residual).all() or residual.max()>tol*1.01:raise RuntimeError('failed chord residual')
    return fu,fv,its,residual,{'fallback_count':len(fallback),'fallback_indices':fallback,'evaluations':evaluations.tolist(),'backtracks':backtracks.tolist()}
