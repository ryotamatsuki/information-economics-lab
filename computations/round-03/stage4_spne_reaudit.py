"""Stage 4/4A reopening: analytic-support checks, independent demand and covariance evaluator.
Maintained domain: compulsory unit demand, prices >=0, fixed two endpoint sites,
public location/price/signal histories; binding linear sharing, mandatory investment.
Numerical checks supplement the analytic proof; no earlier solver imported.
"""
import itertools as it
from collections import Counter
import sympy as sp
import mpmath as mp
mp.mp.dps=70
k=sp.symbols('k', positive=True)
q=4-1/(k+1)-1/(k+3)
print('delta positivity numerator:',sp.factor(q-4*((k+2)/(k+3))**2))
print('alpha>delta numerator:',sp.factor((k+2)/(k+3)-k/(k+1)))
D=(k+1)*(k+3)
print('gamma>alpha sufficient squared difference:',sp.factor((1-2/D)**2-k/(k+2)))
# demand evaluated by integrating the choice rule (affine cost difference), not clipped formula
F=sp.Rational
def demand(l1,l2,p1,p2,t,tie=F(1,2)):
    f0=p1+t*l1*l1-p2-t*l2*l2
    f1=p1+t*(1-l1)**2-p2-t*(1-l2)**2
    if f0==f1:return sp.Integer(1) if f0<0 else sp.Integer(0) if f0>0 else tie
    cuts=[sp.Integer(0),sp.Integer(1)]
    root=-f0/(f1-f0)
    if 0<root<1:cuts.insert(1,root)
    return sum(b-a for a,b in zip(cuts,cuts[1:]) if f0+(f1-f0)*(a+b)/2<0)
# exhaustive finite grid as attack, exact arithmetic; continuous price proof in report.
checks=0
for t in [F(0),F(1,1000),F(1),F(1000)]:
  scale=t if t>0 else F(1)
  for l1,l2 in it.product(range(2),repeat=2):
    p=t if l1!=l2 else F(0)
    d=demand(l1,l2,p,p,t);eq=p*d
    for tie in [F(0),F(1,2),F(1)]:
      for z in [F(0),F(1,10000),F(1,4),F(1,2),F(999,1000),F(1),F(1001,1000),F(3,2),F(2),F(3),F(1000000)]:
        P=z*scale;assert P*demand(l1,l2,P,p,t,tie)<=eq;checks+=1
print('exact off-path individual price deviations:',checks)
# For arbitrary rival prices, reconstruct global candidate deviations from endpoints of piecewise pricing.
# Test every candidate on the exact demand evaluator (including zero/full share and extreme prices).
checks=0; other_eq=[]
for t in [F(1,1000),F(1),F(1000)]:
  grid=[v*t for v in [F(0),F(1,1000),F(1,2),F(999,1000),F(1),F(1001,1000),F(2),F(3),F(1000000)]]
  for p1,p2 in it.product(grid,repeat=2):
    is_ne=True
    for own,rival,l1,l2 in [(p1,p2,0,1),(p2,p1,1,0)]:
      # candidate maximizer on both pieces and boundaries using exact demand evaluator
      cand={F(0),max(F(0),rival-t),rival+t,max(F(0),(rival+t)/2)}
      value=own*demand(l1,l2,own,rival,t)
      best=max(z*demand(l1,l2,z,rival,t) for z in cand)
      is_ne &= value==best;checks+=len(cand)
    if is_ne and (p1,p2)!=(t,t):other_eq.append((t,p1,p2))
assert not other_eq
print('arbitrary-price pair global piecewise BR candidate comparisons:',checks)
# independent covariance evaluation across every location (r=1/kappa,h=1), exact rationals.
P=list(it.product(range(2),repeat=4)); a=sp.Matrix([1,1,-1,-1]); rivals=[2,3,0,1]
def flip(L,i):return L[:i]+(1-L[i],)+L[i+1:]
def label(L):
    if len(set(L))==1:return 'all'
    if sum(L) in (1,3):return '3+1'
    if L[0]==L[1]:return 'portfolio'
    if L[0]==L[2]:return 'paired'
    return 'crossed'
def ev(L,kap):
    C=sp.Matrix(4,4,lambda i,j: int(i==j)+(int(L[i]==L[j])/kap if kap is not None else 0))
    b=C.inv()*a; Q=(a.T*b)[0]
    x=[a[i]*b[i] for i in range(4)]
    return [z/sp.sqrt(Q) for z in x],Q
m=lambda z:mp.mpf(str(sp.N(z,75)))
calls=cases=continuations=0
for kap in [F(1,10**6),F(1,100),F(1,10),F(5,4),F(5),F(100),F(10**6)]:
    E={L:ev(L,kap) for L in P}; U={L:list(map(m,E[L][0])) for L in P};continuations+=len(E)
    allL=(0,0,0,0); paired=(0,1,0,1); port=(0,0,1,1)
    alpha=U[allL][0]-U[flip(allL,0)][0]
    delta=U[paired][0]-U[flip(paired,0)][0]
    gamma=U[flip(port,0)][0]-U[port][0]
    assert 0<delta<alpha<gamma
    for lam in [0,delta/2,delta,(delta+alpha)/2,alpha,(alpha+gamma)/2,gamma,gamma*2]:
      EQ=[]
      for L in P:
        gains=[]
        for i in range(4):
          R=flip(L,i)
          gains.append(U[R][i]-U[L][i]+lam*(int(R[i]!=R[rivals[i]])-int(L[i]!=L[rivals[i]])));calls+=1
        if max(gains)<=mp.mpf('1e-60'):EQ.append(L)
      expected={L for L in P if label(L)=='crossed' or label(L)=='all' and lam<=alpha or label(L)=='paired' and lam<=delta or label(L)=='portfolio' and lam>=gamma}
      assert set(EQ)==expected
      if kap==F(5,4):print('benchmark lambda',mp.nstr(lam,12),'NE',dict(Counter(map(label,EQ))))
      cases+=1
    if kap==F(5,4):
      c=mp.mpf('.5')/mp.sqrt(2*mp.pi*mp.mpf('1.25'))
      print('benchmark t thresholds',*[mp.nstr(2*c*x,17) for x in (delta,alpha,gamma)])
# r=0 boundary (independent signals): all information payoffs invariant, price separation payoff only
for lam in [0,mp.mpf(1)]:
    E={L:ev(L,None) for L in P}
    assert all(z==F(1,2) for U,_ in E.values() for z in U)
    eq=[L for L in P if all(lam*(int(flip(L,i)[i]!=flip(L,i)[rivals[i]])-int(L[i]!=L[rivals[i]]))<=0 for i in range(4))]
    print('r=0 lambda',lam,'location NE',len(eq))
print('location parameter cases',cases,'covariance continuations',continuations,'unilateral deviation checks',calls)
print('PASS; no solver failure; analytical proofs control full continuous domains. High-precision tie checks are diagnostics only.')
