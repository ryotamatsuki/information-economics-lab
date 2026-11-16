"""Stage 7: welfare accounting using certified four-firm endpoint model.
No new strategic primitives. Prices are equilibrium continuations; full coverage,
risk neutrality, binding linear sharing, no extra project costs. W excludes 2*v.
Outputs are diagnostics, not novelty or formal certification.
"""
from itertools import product
from math import sqrt, pi, acos
from collections import Counter
import numpy as np
import sympy as sp
A=np.array([1,1,-1,-1]); R=[2,3,0,1]; Ls=list(product(range(2),repeat=4))

def evaluate(L,t,r,beta=.5,tau2=1.,s=.25):
    # reconstruct Gaussian posterior and integrate equilibrium market prices/travel
    D=np.eye(2)[list(L)]; J=tau2**2*np.linalg.solve((tau2+s)*np.eye(4)+r*D@D.T,A)
    Q=A@J; V=sqrt(Q/(2*pi)); B=beta*A*J/sqrt(2*pi*Q)
    sep=np.array([L[i]!=L[R[i]] for i in range(4)]); d=int(sep.sum()/2)
    p=t*sep/2; travel=t*(d/12+(2-d)/3); CS=-d*t-travel; I=(1-beta)*V
    W=V-travel
    assert np.isclose(W,CS+(p+B).sum()+I)
    return dict(Q=Q,V=V,B=B,profit=p+B,I=I,CS=CS,travel=travel,W=W,d=d)

def label(L):
    if len(set(L))==1:return 'all'
    if sum(L) in (1,3):return '3+1'
    if L[0]==L[1]:return 'portfolio'
    if L[0]==L[2]:return 'paired'
    return 'crossed'

def audit(t,r,beta=.5,tau2=1.,s=.25):
    E={L:evaluate(L,t,r,beta,tau2,s) for L in Ls}
    eq=[L for L in Ls if all(E[L[:i]+(1-L[i],)+L[i+1:]]['profit'][i]<=E[L]['profit'][i]+1e-11 for i in range(4))]
    Wmax=max(e['W'] for e in E.values()); best=max(E[L]['W'] for L in eq); worst=min(E[L]['W'] for L in eq)
    C=E[(0,1,1,0)]; P=E[(0,0,1,1)]; gap=C['V']-P['V']
    assert abs(Wmax-C['W'])<1e-10 and abs(best-Wmax)<1e-10
    assert np.isclose(C['CS'],P['CS']) and np.allclose(C['profit']-P['profit'],beta*gap/4) and np.isclose(C['I']-P['I'],(1-beta)*gap)
    if r==0 and t==0:
        assert len(eq)==16 and np.isclose(Wmax,worst)
    else:
        assert set(label(L) for L in eq)<={'all','paired','crossed','portfolio'}, (t,r,beta,tau2,s,eq)
    predicted=max([0]+[t/2 for L in eq if label(L) in ('all','paired')]+[gap for L in eq if label(L)=='portfolio'])
    assert np.isclose(Wmax-worst,predicted)
    return E,eq,Wmax-worst

if __name__=='__main__':
    print('numpy',np.__version__,'sympy',sp.__version__)
    x=sp.symbols('x'); print('travel separated, coincident', sp.integrate(x*x,(x,0,sp.Rational(1,2)))+sp.integrate((1-x)**2,(x,sp.Rational(1,2),1)),sp.integrate(x*x,(x,0,1)))
    for t in [0,.01,.05,.075,.1312552263136127,.15,.2]:
        E,eq,gap=audit(t,1); print('t',t,'NE',dict(Counter(map(label,eq))),'best W',max(E[L]['W'] for L in eq),'worst W',min(E[L]['W'] for L in eq),'gap',gap)
    E,eq,gap=audit(.2,1)
    for L in [(0,1,1,0),(0,0,1,1),(0,0,0,0),(0,1,0,1),(0,0,0,1)]:
        e=E[L]; print(label(L),'Q',e['Q'],'V',e['V'],'I',e['I'],'firm profit',e['profit'],'CS-2v',e['CS'],'travel',e['travel'],'W-2v',e['W'],'wrong choice prob',acos(sqrt(e['Q']/4))/pi)
    print('absolute gap',gap,'loss share investment',gap/E[(0,1,1,0)]['V'])
    # product competition and public residual noise, fixed other primitives
    for r in [1,.5,.25,.1,.01,0]:
        E,eq,gap=audit(.075,r);print('r',r,'NE',dict(Counter(map(label,eq))),'worst W',min(E[L]['W'] for L in eq),'loss',gap)
    count=0
    for tau2 in [.1,1,10]:
      for s in [0,.25,10]:
       for r in [0,.001,.1,1,100]:
        for beta in [.01,.5,.99]:
         for t in [0,.001,.075,.2,10]:
          audit(t,r,beta,tau2,s);count+=1
    print('grid parameter cases',count,'profiles',count*16,'deviation checks',count*64,'all passed (tolerance 1e-11)')
