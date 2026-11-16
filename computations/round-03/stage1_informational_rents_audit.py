"""Stage 1 audit, 2026-11-16. Independent closed-form evaluator and finite deviations.
Requires numpy and sympy. Firms A1..Am,B1..Bm; two equal fixed investment projects.
All payoffs are signed, risk-neutral, commitment assumed; full market coverage.
For k>2 the product market is a declared extension: equispaced points on [0,1].
Finite enumerations are computational evidence, not independent certification.
"""
from itertools import product
from collections import Counter
import numpy as np
import sympy as sp


def evaluate(m,k,T=1.,s=.25,r=1.,beta=.5):
    a=np.r_[np.ones(m),-np.ones(m)]; N=2*m; h=T+s; A=T*T/h
    L=np.array(list(product(range(k),repeat=N)),dtype=int)
    z=np.zeros_like(L,dtype=float)
    for g in range(k):
        M=L==g; n=M.sum(1); S=M@a
        z+=M*A*(1-a[None,:]*(r*S/(h+r*n))[:,None])
    Q=z.sum(1); V=np.sqrt(Q/(2*np.pi)); B=beta*z/np.sqrt(2*np.pi*Q[:,None])
    for j,l in enumerate(L):
        J=T*T*np.linalg.inv(h*np.eye(N)+r*(l[:,None]==l[None,:]))
        assert np.allclose(z[j],a*(J@a),rtol=1e-11,atol=1e-12)
    x=L/(k-1); rivals=np.r_[np.arange(m,N),np.arange(m)]
    xx=x[:,rivals]; d=abs(x-xx); sm=x+xx
    P=d*np.where(x<xx,(sm+2)**2,(4-sm)**2)/18
    lo=np.minimum(x[:,:m],x[:,m:]); hi=np.maximum(x[:,:m],x[:,m:]);cut=(lo+hi+2)/6
    C=(((cut-lo)**3+lo**3+(1-hi)**3-(cut-hi)**3)/3).sum(1)
    pw=k**np.arange(N-1,-1,-1);ids=L@pw
    destinations=np.stack([np.stack([ids+(g-L[:,i])*pw[i] for g in range(k)],axis=1) for i in range(N)],axis=1)
    cols=np.arange(N)[None,:,None]
    dB=B[destinations,cols]-B[:,:,None];dP=P[destinations,cols]-P[:,:,None]
    dV=V[destinations]-V[:,None,None];dC=C[destinations]-C[:,None,None]
    mask=np.arange(k)[None,None,:]!=L[:,:,None]
    return dict(m=m,k=k,L=L,z=z,Q=Q,V=V,B=B,P=P,C=C,dB=dB,dP=dP,dV=dV,dC=dC,mask=mask,dest=destinations)


def summary(d,t):
    gain=np.where(d['mask'],d['dB']+t*d['dP'],-np.inf)
    maxgain=gain.max((1,2));ne=maxgain<=1e-10;strict=maxgain<-1e-10
    bad=(d['dB']>1e-10)&(d['dV']<-1e-10)
    harmful=bad&(gain>1e-10)&(d['dV']-t*d['dC']<-1e-10)
    W=d['V']-t*d['C']
    print('m,k,t',d['m'],d['k'],t,'profiles',len(ne),'NE',ne.sum(),'strict',strict.sum(),
          'NE V values',np.unique(np.round(d['V'][ne],9)),
          'NE W range', (W[ne].min(),W[ne].max()) if ne.any() else None,'Wmax',W.max(),
          'information-harmful moves',bad.sum(),'profitable & welfare-harmful moves',harmful.sum(),
          'latter directly entering NE', (harmful&ne[d['dest']]).sum())
    for v in np.unique(np.round(d['V'][ne],9)):
        j=np.where(ne&np.isclose(d['V'],v,rtol=0,atol=2e-9))[0][0]
        print(' representative NE',d['L'][j],'V',d['V'][j],'W',W[j],'largest actual deviation gain',maxgain[j])
    return ne


def move(d,before,after,i,t):
    for L in [before,after]:
        j=np.flatnonzero((d['L']==L).all(1))[0];l=d['L'][j]
        J=np.linalg.inv(1.25*np.eye(len(l))+(l[:,None]==l[None,:]))
        print('move',l,'z',d['z'][j,i],'Q',d['Q'][j],'V',d['V'][j],'B',d['B'][j,i],
              'product profit',t*d['P'][j,i],'total profit',d['B'][j,i]+t*d['P'][j,i],
              'travel coefficient',d['C'][j],'W',d['V'][j]-t*d['C'][j],
              'individual posterior variance',1-J[i,i])


def minimal_symbolic():
    # All partitions, one possible unoccupied destination; identical region labels irrelevant.
    k=sp.symbols('kappa',positive=True)
    def configs(n):
        def rec(p):
            if len(p)==n:yield p;return
            for g in range(max(p)+2):yield from rec(p+(g,))
        yield from rec((0,))
    def pars(L,a,i):
        groups=[(L.count(g),sum(a[j] for j,x in enumerate(L) if x==g)) for g in set(L)]
        n=L.count(L[i]);S=sum(a[j] for j,x in enumerate(L) if x==L[i])
        return sp.factor(len(a)-sum(S*S/(k+m) for m,S in groups)),sp.factor(1-a[i]*S/(k+n))
    def sign(f):
        num,den=sp.factor(sp.together(f)).as_numer_denom()
        c=sp.Poly(num,k).all_coeffs();b=sp.Poly(den,k).all_coeffs()
        # Denominator sign verified explicitly, avoiding accidental sign flip.
        if all(v<=0 for v in b):c=[-v for v in c];b=[-v for v in b]
        assert all(v>=0 for v in b) and any(v>0 for v in b)
        if all(v==0 for v in c):return 0
        if all(v>=0 for v in c):return 1
        if all(v<=0 for v in c):return -1
        return None
    for N in [2,4]:
        a=(1,)*(N//2)+(-1,)*(N//2);pairs=set()
        for L in configs(N):
            for i in range(N):
                q,z=pars(L,a,i)
                for g in range(max(L)+2):
                    if g==L[i]:continue
                    qq,zz=pars(L[:i]+(g,)+L[i+1:],a,i);pairs.add((q,z,qq,zz))
        result=Counter((sign(qq-q),sign(zz**2*q-z**2*qq)) for q,z,qq,zz in pairs)
        assert all(x is not None and y is not None and not(x<0 and y>0) for x,y in result)
        print('minimal-N symbolic sign pairs (delta q, delta squared payoff), N=',N,dict(result))
    u,S=sp.symbols('u S',positive=True)
    assert sp.simplify(1-S*S/u+(S-1)**2/(u-1)-(1-S/u)**2*u/(u-1))==0

if __name__=='__main__':
    minimal_symbolic()
    d=evaluate(3,3)
    for t in [0,.01,.05,.075,.1,.2]:summary(d,t)
    print('Original six-firm move, line product-market extension at t=.075:')
    move(d,(0,1,2,0,1,1),(1,1,2,0,1,1),0,.075)
    print('Eight-firm move, constant product-market profit and transport cost:')
    d8=evaluate(4,3);move(d8,(0,2,1,1,1,2,2,0),(2,2,1,1,1,2,2,0),0,.075)
    print('Public information comparative static, original two-site four-firm game:')
    for r in [1,.5,.25,.1,0]:
        print('residual r=',r);summary(evaluate(2,2,r=r),.075)
