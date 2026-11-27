"""Small reproducible API for fixed-cost OT composition audits."""
import time
import numpy as np
from .core import derivatives,predict,solve
from .batch import batch_solve,readout
from .chord import factor_target,chord_solve

class CompositionAudit:
    """Fit a balanced OT base plan and audit reference composition changes.

    Parameters
    ----------
    cost : (n_query, n_reference) finite array
        Fixed distance matrix; choose normalization before constructing this object.
    labels : length n_reference
        Existing reference annotations. No query labels are required.
    epsilon : float
        Fixed strictly positive entropic regularization.
    a, b : optional positive marginals summing to one
        Defaults are uniform. Zero weights and unbalanced OT are unsupported.

    Notes
    -----
    This estimates sensitivity to specified class proportions, not confidence
    intervals or the true class composition. It does not improve the accuracy
    of the underlying balanced-OT model. Dense memory is O(n*m + m*m).
    """
    def __init__(self,cost,labels,epsilon=.05,a=None,b=None,base_tol=1e-11):
        t=time.perf_counter();C=np.asarray(cost,dtype=np.float64)
        if C.ndim!=2 or not np.isfinite(C).all() or C.size==0:raise ValueError('cost must be a finite nonempty matrix')
        if not np.isfinite(epsilon) or epsilon<=0:raise ValueError('epsilon must be positive')
        if not np.isfinite(base_tol) or base_tol<=0:raise ValueError('base_tol must be positive')
        n,m=C.shape;labels=np.asarray(labels)
        if labels.shape!=(m,):raise ValueError('labels must align with cost columns')
        self.classes,inverse=np.unique(labels.astype(str),return_inverse=True)
        self.Y=np.eye(len(self.classes))[inverse]
        self.a=self._marginal(a,n);self.b=self._marginal(b,m)
        with np.errstate(over='ignore',under='ignore'):
            self.K=np.exp(-C/epsilon)
        if not np.isfinite(self.K).all() or np.any(self.K<=0):raise ValueError('kernel overflow or underflow: rescale cost or use a stabilized OT solver')
        self.p0=self.b@self.Y;self.W=self.Y*self.b[:,None]/self.p0
        self.db=self.W-self.b[:,None]
        self.base_solver='AA2'
        self.base_fallback_reason=''
        try:
            self.G,self.u,self.v,_,self.base_residual=solve(self.K,self.a,self.b,tol=base_tol)
        except RuntimeError as err:
            self.base_fallback_reason=str(err)
            self.base_solver='plain Sinkhorn after AA2 nonconvergence'
            self.G,self.u,self.v,_,self.base_residual=solve(self.K,self.a,self.b,tol=base_tol,memory=0,maxiter=100000)
        self.base_probabilities=self.G@self.Y/self.a[:,None]
        self.base_seconds=time.perf_counter()-t;self._cache={}

    @staticmethod
    def _marginal(x,n):
        x=np.full(n,1/n) if x is None else np.asarray(x,dtype=float)
        if x.shape!=(n,) or not np.isfinite(x).all() or np.any(x<=0) or not np.isclose(x.sum(),1,rtol=0,atol=1e-12):
            raise ValueError('marginals must be positive, aligned, and sum to one')
        return x

    def run(self,proportions,order=0,corrector='chord',tol=1e-7):
        """Return exact-to-tolerance class readouts and sensitivity summaries.

        proportions is (B, number of classes), following ``self.classes``.
        order=0 reuses the base scaling; 1/2 use Taylor prediction.
        corrector='chord' caches the base target Schur factor; 'aa2' uses
        the common batched Sinkhorn corrector with Anderson memory two.
        Each returned class readout is from a numerically corrected plan.
        """
        ps=np.asarray(proportions,dtype=float)
        if ps.ndim!=2 or ps.shape[1]!=len(self.classes) or ps.shape[0]==0:raise ValueError('invalid composition matrix shape')
        if not np.isfinite(ps).all() or np.any(ps<=0) or not np.allclose(ps.sum(1),1,rtol=0,atol=1e-12):raise ValueError('each composition must be positive and sum to one')
        if order not in [0,1,2] or corrector not in ['chord','aa2']:raise ValueError('unsupported order/corrector')
        if not np.isfinite(tol) or tol<=0:raise ValueError('tol must be positive')
        t=time.perf_counter();cached=(order,corrector) in self._cache
        if not cached:
            D=None if order==0 else derivatives(self.G,self.a,self.b,self.db,side='target',second=order==2)
            ch=(factor_target(self.G,self.a,self.b) if D is None else D['factor']) if corrector=='chord' else None
            self._cache[(order,corrector)]=(D,ch)
        D,ch=self._cache[(order,corrector)];setup=time.perf_counter()-t
        bs=self.W@ps.T;delta=ps-self.p0;B=len(ps);fallback_predictor=0
        scaling=self.v if corrector=='chord' else self.u
        if D is None:init=np.tile(scaling[:,None],(1,B))
        else:
            E=dict(D)
            if corrector=='chord':
                E['f']=E['g']
                if order==2:E['f2']=E['g2']
            vals=[predict(np.log(scaling),E,d,quadratic=order==2) for d in delta]
            init=np.column_stack([x[0] for x in vals]);fallback_predictor=sum(bool(x[1]) for x in vals)
        if corrector=='chord':U,V,its,res,info=chord_solve(self.K,self.a,bs,init,ch,tol=tol)
        else:U,V,its,res=batch_solve(self.K,self.a,bs,init,tol=tol);info={'fallback_count':0}
        P=readout(self.K,self.a,self.Y,U,V)
        if not np.isfinite(P).all() or res.max()>tol*1.01:raise RuntimeError('invalid corrected result')
        base_labels=self.base_probabilities.argmax(1);scenario_labels=P.argmax(2)
        return {'classes':self.classes.copy(),'proportions':ps.copy(),'probabilities':P,
            'base_probabilities':self.base_probabilities.copy(),'switch_frequency':(scenario_labels!=base_labels[:,None]).mean(1),
            'max_probability_change':np.abs(P-self.base_probabilities[:,None,:]).max((1,2)),
            'probability_span':(P.max(1)-P.min(1)).max(1),'marginal_residual':res,'iterations':its,
            'diagnostics':{**info,'predictor_fallback_count':fallback_predictor,'cached_setup':cached,
                          'setup_seconds':setup,'total_seconds':time.perf_counter()-t,'base_seconds':self.base_seconds}}

    def local_sensitivity(self):
        """Return dP/dt along p(t)=(1-t)*p0+t*e_j for each reference class j.

        Shape is (query cell, perturbed reference class, predicted class).
        This is a local derivative, not a corrected finite-perturbation result.
        To report change per unit increase of p_j with proportional reduction
        in other classes, divide axis j by (1-p0[j]).
        """
        D=derivatives(self.G,self.a,self.b,self.db,side='target',second=False)
        q=len(self.classes);rhs=(D['g'][:,:,None]*self.Y[:,None,:]).reshape(len(self.b),q*q)
        out=(self.G@rhs).reshape(len(self.a),q,q)/self.a[:,None,None]
        out+=D['f'][:,:,None]*self.base_probabilities[:,None,:]
        return out
