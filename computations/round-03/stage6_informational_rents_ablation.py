"""Stage 6: matched-level ablation, exact harmful-deviation regression.
No external source model implemented here; ablation is internal, not evidence of priority.
"""
from fractions import Fraction as F
from itertools import product
from math import sqrt
import sympy as s

k=s.symbols('k',positive=True)
a=1/(2*(k+1)); d=1/(2*(k+3)); g=(1/(k+3)+2/(k+2))/2
print('local alpha-delta:',s.factor(a-d))
print('local gamma-alpha:',s.factor(g-a))

def counts(A,B,k):
    n=[a+b for a,b in zip(A,B)]; S=[a-b for a,b in zip(A,B)]
    q=sum(n)-sum(F(z*z)/(k+m) for z,m in zip(S,n))
    x=[1-F(z)/(k+m) for z,m in zip(S,n)]
    return q,x
q,x=counts([0,1,1,6],[0,3,5,0],F(5,4))
Q,X=counts([0,0,2,6],[0,3,5,0],F(5,4))
assert Q<q and X[2]<x[1] and X[2]**2*q>x[1]**2*Q
print('strict decline in local score and aggregate information:',q,Q,x[1],X[2])
print('full payoff before/after:',float(x[1])/sqrt(q),float(X[2])/sqrt(Q))
print('local payoff before/after:',float(x[1])/4,float(X[2])/4)

def pay(L,k,theta):
    A=[1,1,-1,-1]; ng=[L.count(g) for g in range(2)]; sg=[sum(A[i] for i in range(4) if L[i]==g) for g in range(2)]
    q=4-sum(F(S*S)/(k+n) for n,S in zip(ng,sg))
    x=[1-F(A[i]*sg[L[i]])/(k+ng[L[i]]) for i in range(4)]
    return q,[float(v)/(sqrt(q) if theta else 2) for v in x]
profiles=list(product(range(2),repeat=4)); rivals=[2,3,0,1]
flip=lambda L,i:L[:i]+(1-L[i],)+L[i+1:]
for k in [F(1,1000),F(1,100),F(1,10),F(5,4),F(5),F(100),F(10000)]:
    for theta in [False,True]:
        p={L:pay(L,k,theta)[1] for L in profiles}
        all0=(0,0,0,0); pair=(0,1,0,1); port=(0,0,1,1)
        alpha=p[all0][0]-p[flip(all0,0)][0];delta=p[pair][0]-p[flip(pair,0)][0];gamma=p[flip(port,0)][0]-p[port][0]
        assert 0<delta<alpha<gamma
        for lam in [0,delta/2,delta,(delta+alpha)/2,alpha,(alpha+gamma)/2,gamma,2*gamma]:
            ne=[]
            for L in profiles:
                gains=[]
                for i in range(4):
                    R=flip(L,i);gains.append(p[R][i]-p[L][i]+lam*(int(R[i]!=R[rivals[i]])-int(L[i]!=L[rivals[i]])))
                if max(gains)<=1e-10:ne.append(L)
            assert len(ne)==2+2*(lam<=delta)+2*(lam<=alpha)+2*(lam>=gamma),(k,theta,lam,ne)
        print('kappa',k,'full' if theta else 'local','delta alpha gamma',delta,alpha,gamma)
print('All 112 parameter/regime cases; 7168 finite deviations; matched phase topology.')
print('Note: local model is a payoff diagnostic, not a microfounded alternative contract or prior theorem.')
