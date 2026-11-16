'''Stage 4: exact finite deviations, proof inequalities, two-site threshold screen.
Stage 4A independent certification not performed. Python >=3.10; numpy and sympy.
'''
from fractions import Fraction as F
from itertools import product
from math import sqrt, pi
from collections import Counter, deque
import numpy as np
import sympy as sp

def exact(L,m,k,kap):
    types=[1]*m+[-1]*m
    n=[L.count(g) for g in range(k)]
    S=[sum(a for a,l in zip(types,L) if l==g) for g in range(k)]
    x=[1-a*F(S[l],kap+n[l]) for a,l in zip(types,L)]
    q=sum(x)
    assert q==2*m-sum(F(s*s,kap+ng) for s,ng in zip(S,n))
    return q,x,S,n

def direct(L,m,k,kap):
    # Independent primitive Gaussian covariance solve, h=1, tau^2=1, r=1/kap.
    a=np.r_[np.ones(m),-np.ones(m)]
    D=np.eye(k)[np.asarray(L)]
    J=np.linalg.solve(np.eye(2*m)+D@D.T/float(kap),a)
    return float(a@J),a*J

def audit(m,k,kap):
    states=list(product(range(k),repeat=2*m)); d={L:exact(L,m,k,kap) for L in states}
    edges=0;harm=0;ne=0;maxpath=0; badne=0; ties=0; adj={L:[] for L in states};ind={L:0 for L in states}
    for L,(q,x,S,n) in d.items():
        # Primitive comparison and exact rational cross multiplication.
        q0,x0=direct(L,m,k,kap);assert np.allclose([float(q),*map(float,x)],[q0,*x0],atol=1e-10,rtol=1e-10)
        for i in range(2*m):
            for g in range(k):
                if L[i]==g:continue
                P=L[:i]+(g,)+L[i+1:];p,y,_,nn=d[P]
                diff=y[i]**2*q-x[i]**2*p
                if diff==0:ties+=1
                if diff>0:
                    edges+=1;adj[L].append(P);ind[P]+=1
                    if p<q:
                        harm+=1;assert nn[g]>n[L[i]]
        if not adj[L]:ne+=1;badne+=int(any(S))
        assert bool(adj[L])==bool(any(S))
        if any(S):
            # constructive majority-to-minority path, starting at globally minimum x
            P=L;steps=0
            while any(d[P][2]):
                p,y,T,_=d[P];i=min(range(2*m),key=lambda i:y[i]);a=1 if i<m else -1
                assert a*T[P[i]]>0
                g=next(g for g in range(k) if a*T[g]<0)
                R=P[:i]+(g,)+P[i+1:];rr,yy,U,_=d[R]
                assert yy[i]**2*p>y[i]**2*rr and rr>p
                assert sum(map(abs,U))==sum(map(abs,T))-2
                P=R;steps+=1
            assert steps==sum(map(abs,S))//2<=m
            maxpath=max(maxpath,steps)
    todo=deque(L for L in states if ind[L]==0);visited=0;longest={L:0 for L in states}
    while todo:
        L=todo.popleft();visited+=1
        for R in adj[L]:
            longest[R]=max(longest[R],longest[L]+1);ind[R]-=1
            if not ind[R]:todo.append(R)
    print('exact',m,k,str(kap),'states',len(states),'strict edges',edges,'ties',ties,'harmful',harm,'NE',ne,'inefficient NE',badne,'constructive max length',maxpath,'DAG',visited==len(states),'longest DAG path',max(longest.values()) if visited==len(states) else 'UNRESOLVED')

def thresholds(kap):
    q=4-1/(kap+1)-1/(kap+3)
    uS=kap/(kap+1)/sqrt(q);uL=(kap+2)/(kap+3)/sqrt(q);uH=(kap+4)/(kap+3)/sqrt(q);uP=.5*sqrt(kap/(kap+2))
    return .5-uS,.5-uL,uH-uP

def label(L):
    if len(set(L))==1:return 'all'
    if sum(L) in (1,3):return '3+1'
    if L[0]==L[1]:return 'portfolio'
    if L[0]==L[2]:return 'paired'
    return 'crossed'

def product_audit():
    # u=x/sqrt(q); lambda = (t/2) / [beta*tau^2/sqrt(2*pi*h)]
    profiles=list(product(range(2),repeat=4));R=[2,3,0,1]
    calls=0
    for kap in [1e-4,.01,.1,1.25,5,100,1e4]:
        A,D,G=thresholds(kap);assert 0<D<A and D<G
        for lam in sorted(set([0,A,D,G,1,2*A,2*G]+[max(0,t-eps) for t in [A,D,G] for eps in [1e-8,1e-6]]+[t+eps for t in [A,D,G] for eps in [1e-8,1e-6]])):
            def payoff(L,i):
                q,x,S,n=exact(L,2,2,F(str(kap)))
                return float(x[i])/sqrt(float(q))+lam*(L[i]!=L[R[i]])
            kinds=Counter()
            for L in profiles:
                gain=max(payoff(L[:i]+(1-L[i],)+L[i+1:],i)-payoff(L,i) for i in range(4));calls+=4
                expected=label(L)=='crossed' or (label(L)=='all' and lam<=A+1e-12) or (label(L)=='paired' and lam<=D+1e-12) or (label(L)=='portfolio' and lam>=G-1e-12)
                assert (gain<=1e-12)==expected,(kap,lam,L,gain,A,D,G)
                if gain<=1e-12:kinds[label(L)]+=1
            assert kinds['crossed']==2
        print('thresholds kappa',kap,'all ceiling',A,'paired ceiling',D,'portfolio floor',G)
    print('product full finite deviation calls',calls)
    kap=1.25;A,D,G=thresholds(kap);c=.5/sqrt(2*pi*1.25)
    print('original t thresholds',*[2*c*x for x in (A,D,G)],'scale',c)
    print('V optimal and portfolio',sqrt(4/1.25/(2*pi)),sqrt(4*kap/(kap+2)/1.25/(2*pi)))
    # Exact independent endpoint demand & travel check. Price t at opposite sites, epsilon deviations and corners.
    for t in [.001,.075,1,100]:
        for p in [0,t*.2,t*.99,t,t*1.01,2*t,3*t,100*t]:
            share=np.clip((t+t-p)/(2*t),0,1)
            assert p*share <= t/2+1e-12*max(1,t)
        # balanced split of consumers, integrate squared distances
        x=sp.symbols('x'); assert sp.integrate(x*x,(x,0,sp.Rational(1,2)))+sp.integrate((1-x)**2,(x,sp.Rational(1,2),1))==sp.Rational(1,12)
    print('endpoint demand deviations and transportation cost checked')

if __name__=='__main__':
    print('numpy',np.__version__,'sympy',sp.__version__)
    d,S=sp.symbols('d S',positive=True)
    assert sp.factor(1-S**2/d+(S-1)**2/(d-1)-d/(d-1)*(1-S/d)**2)==0
    # Reduction of each regional squared imbalance by majority exit/minority entry.
    print('majority-out numerator',sp.expand(S**2*(d-1)-(S-1)**2*d))
    # Strict payoff loss at balanced state analytic numerator written in report.
    for m,k in [(1,2),(2,3),(3,3),(4,3)]:
        for kap in [F(1,100),F(5,4),F(100)]:audit(m,k,kap)
    L=(0,1,2,0,1,1);R=(1,1,2,0,1,1);q,x,*_=exact(L,3,3,F(5,4));p,y,*_=exact(R,3,3,F(5,4));assert p<q and x[0]==y[0] and y[0]**2*q>x[0]**2*p
    print('exact harmful regression',q,p,x[0],y[0], 'V',sqrt(float(q)*.8/(2*pi)),sqrt(float(p)*.8/(2*pi)))
    product_audit()
