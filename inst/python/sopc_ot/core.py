"""Repeated-marginal balanced entropic OT. Known implicit/Taylor algebra.
Source/target Schur choice and mixed second moments are implementation work;
this file does not assert mathematical novelty. Positive fixed kernel only.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve

def derivatives(G,a,b,db,side='auto',second=True):
    n,m=G.shape
    if side=='auto': side='source' if n<=m else 'target'
    if side=='source':
        A=np.diag(a)-(G/b)@G.T+np.outer(a,a)
        chol=cho_factor(A,lower=True,check_finite=False)
        f=cho_solve(chol,-G@(db/b[:,None]),check_finite=False)
        g=(db-G.T@f)/b[:,None]
    else:
        A=np.diag(b)-(G.T/a)@G+np.outer(b,b)
        chol=cho_factor(A,lower=True,check_finite=False)
        g=cho_solve(chol,db,check_finite=False)
        f=-(G@g)/a[:,None]
    # Keep only source Taylor coefficients for a v-first Sinkhorn corrector.
    # Canonical source gauge also allows comparison between Schur choices.
    gauge=a@f; f=f-gauge; g=g+gauge
    result={'f':f,'g':g,'side':side,'factor_dim':len(A),'factor':chol}
    if not second: return result
    q=db.shape[1]; ij=np.array([(i,j) for i in range(q) for j in range(i,q)])
    i,j=ij.T
    # Mixed second moments; no n*m*q^2 tensor.
    row=G@(g[:,i]*g[:,j])-a[:,None]*f[:,i]*f[:,j]
    col=G.T@(f[:,i]*f[:,j])+g[:,i]*db[:,j]+g[:,j]*db[:,i]-b[:,None]*g[:,i]*g[:,j]
    if side=='source':
        f2=cho_solve(chol,-row+G@(col/b[:,None]),check_finite=False)
        g2=(-col-G.T@f2)/b[:,None]
    else:
        g2=cho_solve(chol,-col+G.T@(row/a[:,None]),check_finite=False)
        f2=(-row-G@g2)/a[:,None]
    gauge=a@f2;f2-=gauge;g2+=gauge
    result.update(f2=f2,g2=g2,ij=ij)
    return result

def predict(logu,D,delta,quadratic=True):
    x=logu+D['f']@delta
    if quadratic:
        i,j=D['ij'].T
        weights=delta[i]*delta[j]*np.where(i==j,.5,1.)
        x+=D['f2']@weights
    if not np.isfinite(x).all() or np.ptp(x)>1300:
        return np.exp(logu),'invalid_predictor'
    # Multiplicative gauge does not affect final coupling.
    x-=.5*(x.max()+x.min())
    return np.exp(x),''

def solve(K,a,b,u0=None,tol=1e-9,memory=2,maxiter=10000):
    u=np.ones(len(a)) if u0 is None else u0.copy()
    hx,hf=[],[]
    for it in range(1,maxiter+1):
        if memory:
            x=np.log(u)
            if it-1>memory:
                X,F=np.column_stack(hx),np.column_stack(hf)
                H=F-X;gram=H.T@H;norm=np.linalg.norm(gram)
                if norm>0:
                    coeff=np.linalg.solve(gram/norm+.01*np.eye(memory),np.ones(memory))
                    coeff/=coeff.sum();x=F@coeff;u=np.exp(x)
            hx.append(x.copy());hx=hx[-memory:]
        v=b/(K.T@u);u=a/(K@v)
        if memory:hf.append(np.log(u));hf=hf[-memory:]
        if it%5==0:
            err=float(np.abs(v*(K.T@u)-b).sum())
            if not np.isfinite(err) or not np.isfinite(u).all():raise FloatingPointError('nonfinite solve')
            if err<=tol:break
    else:raise RuntimeError('no convergence')
    G=u[:,None]*K*v[None,:]
    err=max(float(np.abs(G.sum(1)-a).sum()),float(np.abs(G.sum(0)-b).sum()))
    if not np.isfinite(G).all() or err>tol*1.001:raise RuntimeError(f'bad residual {err}')
    return G,u,v,it,err
