"""Stage 1: finite location-game reproduction, 2026-11-16.
Dependencies: numpy. Signed project cashflows, risk-neutral committed participants.
Order: A1,A2,B1,B2. Same-industry rivals: (A1,B1),(A2,B2).
"""
from itertools import product
from math import sqrt, pi
import numpy as np

def solve(a=(1,1,-1,-1), t=.075, beta=.5, T=1., r=1., e=.25, outside=False):
    a=np.array(a)
    rival=[2,3,0,1]
    rows={}
    for L in product(range(2), repeat=4):
        l=np.array(L)
        J=T**2*np.linalg.inv((T+e)*np.eye(4)+r*(l[:,None]==l[None,:]))
        Q=float(a@J@a)
        b=a*(J@a)/sqrt(2*pi*Q)
        if outside:
            # Gaussian integration by parts; project chosen iff it beats both rival and zero.
            old=b.copy()
            for sign in [1,-1]:
                w=(a==sign).astype(float)
                mask=a==sign
                b[mask]=.5*(old[mask]+(J@w)[mask]/sqrt(2*pi*float(w@J@w)))
        u=t/2*(l!=l[rival])+beta*b
        rows[L]=(Q,sqrt(Q/(2*pi)),u)
    for L,(Q,V,u) in rows.items():
        gain=np.array([rows[L[:i]+(1-L[i],)+L[i+1:]][2][i]-u[i] for i in range(4)])
        yield L,Q,V,u,gain,bool(max(gain)<=1e-12)

def thresholds(T=1.,r=1.,e=.25,beta=.5):
    h=T+e; A=T*T/h
    x=r/(h+3*r); y=r/(h+r)
    q0=4*A; q31=A*(4-x-y); qp=4*A*(1-2*r/(h+2*r))
    b0=A/sqrt(2*pi*q0)
    bm=A*(1-x)/sqrt(2*pi*q31)
    bn=A*(1+x)/sqrt(2*pi*q31)
    bs=A*(1-y)/sqrt(2*pi*q31)
    bp=A*(1-2*r/(h+2*r))/sqrt(2*pi*qp)
    return 2*beta*(b0-bm),2*beta*(b0-bs),2*beta*(bn-bp)

if __name__=="__main__":
    print("thresholds (same-industry, all-together, project sorting):", thresholds())
    for name,a,outside,e in [
       ("baseline",(1,1,-1,-1),False,.25),
       ("reassigned",(1,-1,-1,1),False,.25),
       ("higher iid noise",(1,1,-1,-1),False,1.),
       ("outside option",(1,1,-1,-1),True,.25),
       ("outside option reassigned",(1,-1,-1,1),True,.25)]:
        print("\n"+name)
        rows=list(solve(a=a,outside=outside,e=e))
        for L,Q,V,u,g,ne in rows:
            print("".join(map(str,L)), f"Q={Q:.8f}", f"forced-choice value={V:.8f}",
                  "u=",np.round(u,8), "deviation gains=",np.round(g,8), "NE=",ne)
        print("Equilibria:",[L for L,Q,V,u,g,ne in rows if ne])
