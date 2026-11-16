"""Stage 4A distinct implementation: Gaussian marginal matrix solve, no import of Stage 4.
Reviewer had read Stage 4 and its code. Not a blind or separate-human review.
"""
import itertools as it
from collections import Counter
from math import sqrt, pi
import numpy as np
import sympy as sp

def primitive(L,kappa):
    N=len(L); a=np.array([1.]*(N//2)+[-1.]*(N//2)); z=np.array(L)
    C=np.eye(N)+(z[:,None]==z[None,:])/kappa
    b=np.linalg.solve(C,a)
    q=a@b; x=a*b
    assert q>0 and np.all(x>0)
    return q,x,x/sqrt(q)

def exact_matrix(L,kappa):
    N=len(L); a=sp.Matrix([1]*(N//2)+[-1]*(N//2))
    C=sp.Matrix(N,N,lambda i,j: int(i==j)+sp.Rational(int(L[i]==L[j]),1)/kappa)
    b=C.inv()*a
    return (a.T*b)[0], [a[i]*b[i] for i in range(N)]

def symbolic():
    k=sp.symbols('k',positive=True); d=(k+1)*(k+3)
    assert sp.simplify((1-2/d)**2-k/(k+2)-2*(k**4+6*k**3+10*k**2+4*k+1)/((k+1)**2*(k+2)*(k+3)**2))==0
    print('threshold-order squared gap:',sp.factor((1-2/d)**2-k/(k+2)))
    n,m,s=sp.symbols('n m s',positive=True)
    print('majority departure numerator:',sp.expand(s*s*(k+n-1)-(s-1)**2*(k+n)))
    L=(0,1,2,0,1,1); R=(1,1,2,0,1,1)
    q,x=exact_matrix(L,sp.Rational(5,4));Q,X=exact_matrix(R,sp.Rational(5,4))
    assert Q<q and X[0]**2*q>x[0]**2*Q
    print('exact harmful deviation:',q,Q,x[0],X[0], 'concentration',sum(L.count(g)**2 for g in set(L)),sum(R.count(g)**2 for g in set(R)))
    print('V, B:',*[float(sp.sqrt(v*sp.Rational(4,5)/(2*sp.pi))) for v in (q,Q)],*[float(sp.Rational(1,2)*y*sp.sqrt(sp.Rational(4,5)/(2*sp.pi*v))) for v,y in ((q,x[0]),(Q,X[0]))])

def base():
    profiles_count=deviations=ties=harm=0
    for N,K in [(2,2),(4,3),(6,3),(8,2)]:
        for k in [.001,.1,1.25,10,1000]:
            a=np.array([1]*(N//2)+[-1]*(N//2)); P=list(it.product(range(K),repeat=N)); data={L:primitive(L,k) for L in P}; c=Counter()
            for L,(q,x,u) in data.items():
                S=np.bincount(L,weights=a,minlength=K); bal=np.all(S==0); gains=[]
                for i in range(N):
                    for g in range(K):
                        if g==L[i]:continue
                        R=L[:i]+(g,)+L[i+1:]; Q,X,U=data[R]; gain=U[i]-u[i]; gains.append(gain); deviations+=1
                        if abs(gain)<1e-10:ties+=1
                        if gain>1e-10 and Q<q-1e-10:
                            harm+=1; assert R.count(g)>L.count(L[i])
                assert (max(gains)<=1e-10)==bal,(N,K,k,L,max(gains))
                if bal:assert max(gains)<-1e-10;c['NE']+=1
                else:
                    i=int(np.argmin(x));g=next(g for g in range(K) if a[i]*S[g]<0)
                    R=L[:i]+(g,)+L[i+1:];Q,X,U=data[R]
                    assert U[i]>u[i]+1e-10 and Q>q+1e-10
                    assert np.abs(np.bincount(R,weights=a,minlength=K)).sum()==np.abs(S).sum()-2
                profiles_count+=1
            print('base',N,K,k,'profiles',len(P),'NE',c['NE'])
    print('base totals profiles',profiles_count,'deviations',deviations,'ties within tolerance',ties,'harmful',harm)

def market():
    # Independent exact location payoff from primitive normal covariance and demand continuation.
    P=list(it.product(range(2),repeat=4)); rivals=(2,3,0,1); calls=0; cases=0
    for k in [.0001,.01,.1,1.25,5,100,10000.]:
        data={L:primitive(L,k)[2] for L in P}
        # compute switching margins from actual profiles, not Stage 4 formulas
        p0=(0,0,0,0);p1=(0,1,0,1);p2=(0,0,1,1)
        flip=lambda L,i: L[:i]+(1-L[i],)+L[i+1:]
        alpha=data[p0][0]-data[flip(p0,0)][0];delta=data[p1][0]-data[flip(p1,0)][0];gamma=data[flip(p2,0)][0]-data[p2][0]
        assert 0<delta<alpha<gamma
        for lam in [0,delta/2,delta,(alpha+delta)/2,alpha,(alpha+gamma)/2,gamma,2*gamma]:
            ne=[]
            for L in P:
                u=data[L]+lam*np.array([L[i]!=L[rivals[i]] for i in range(4)])
                d=[]
                for i in range(4):
                    R=flip(L,i);d.append(data[R][i]+lam*(R[i]!=R[rivals[i]])-u[i]);calls+=1
                if max(d)<=1e-10:ne.append(L)
            expected=2+2*(lam<=alpha)+2*(lam<=delta)+2*(lam>=gamma)
            assert len(ne)==expected,(k,lam,len(ne),expected)
            cases+=1
        if k==1.25:
            scale=.5/sqrt(2*pi*1.25)
            print('original t thresholds (delta,alpha,gamma):',*[2*scale*z for z in (delta,alpha,gamma)])
            a=primitive((0,1,1,0),k)[0]; b=primitive((0,0,1,1),k)[0]; V=lambda q:sqrt(q/1.25/(2*pi));print('V efficient, inefficient, loss fraction',V(a),V(b),1-V(b)/V(a))
    print('product',cases,'parameter cases',calls,'full finite deviations')
    # Direct piecewise demand with full coverage, p>=0, tie at a point irrelevant.
    price_calls=0
    for t in [.001,.1,1,100]:
        for p in np.r_[np.linspace(0,3*t,301),100*t]:
            share=min(1,max(0,(2*t-p)/(2*t)))
            assert p*share<=t/2+1e-10*max(1,t);price_calls+=1
    print('endpoint price deviations',price_calls)

if __name__=='__main__':
    print('numpy',np.__version__,'sympy',sp.__version__)
    symbolic();base();market(); print('PASS on enumerated tests. General proofs in certification memo; no FIP theorem attempted.')
