"""Stage 3 screening; aggregated count states; a cycle in labeled states projects to a cycle in count states, and a count cycle can be repeated to restore labels. Numerical screening only."""
import numpy as np
from itertools import product
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from math import comb

def compositions(m,k):
    if k==1:
        yield (m,);return
    for i in range(m+1):
        for t in compositions(m-i,k-1):yield (i,)+t

def graph(m,k):
    C=np.array(list(compositions(m,k)),dtype=np.int32);c=len(C)
    code={tuple(x):i for i,x in enumerate(C)}
    A=np.repeat(C,c,axis=0);B=np.tile(C,(c,1)); S=A-B; n=A+B
    trans=[]
    for i,x in enumerate(C):
        for g in range(k):
            if x[g]==0:continue
            for h in range(k):
                if g==h:continue
                y=x.copy();y[g]-=1;y[h]+=1
                trans.append((i,code[tuple(y)],g,h))
    u,v,g,h=np.array(trans).T; a=np.arange(c)
    uA=np.repeat(u,c)*c+np.tile(a,len(u));vA=np.repeat(v,c)*c+np.tile(a,len(u))
    uB=np.tile(u,c)+np.repeat(a,len(u))*c;vB=np.tile(v,c)+np.repeat(a,len(u))*c
    src=np.r_[uA,uB];dst=np.r_[vA,vB];g=np.r_[np.repeat(g,c),np.tile(g,c)];h=np.r_[np.repeat(h,c),np.tile(h,c)]
    sign=np.r_[np.ones(len(uA)), -np.ones(len(uB))]
    return S,n,src,dst,g,h,sign

def scan(m,k,kappas):
    S,n,u,v,g,h,a=graph(m,k);N=len(n);balanced=(S==0).all(1)
    print('m,k,states,feasible transitions',m,k,N,len(u),flush=True)
    for kap in kappas:
        q=2*m-(S*S/(kap+n)).sum(axis=1)
        x=1-a*S[u,g]/(kap+n[u,g]);y=1-a*S[v,h]/(kap+n[v,h]);diff=y*y*q[u]-x*x*q[v]
        good=diff>1e-10*np.maximum(q[u],q[v])
        i,j=u[good],v[good]
        graphmat=coo_matrix((np.ones(len(i),dtype=np.int8),(i,j)),shape=(N,N)).tocsr()
        ncomp,labels=connected_components(graphmat,directed=True,connection='strong');sizes=np.bincount(labels)
        sinks=np.bincount(i,minlength=N)==0
        harm=(q[j]<q[i]-1e-10)
        print('kap',kap,'edges',len(i),'cycle Scc',int((sizes>1).sum()),'max_scc',sizes.max(),'sinks',sinks.sum(),'bad_sinks',sum(sinks&~balanced),'balanced_non_sinks',sum(balanced&~sinks),'harm',harm.sum(),flush=True)
        if sizes.max()>1:
            c=labels[i]==labels[j];print('example cycle SCC nodes',[(S[d].tolist(),n[d].tolist()) for d in np.flatnonzero(labels==sizes.argmax())[:10]],flush=True)
            return False
    return True
if __name__=='__main__':
    for m,k in [(3,3),(4,3),(5,3),(6,3),(8,3),(10,3),(15,3),(20,3),(4,4),(6,4),(8,4)]:
        if not scan(m,k,[.01,.1,.5,1.25,5.,20.,100.]):break
    # Institutional falsifier: remove the global normalizer, retaining the local contribution.
    for m,k,kap in [(3,3,1.25),(8,4,1.25),(20,3,1.25)]:
        S,n,u,v,g,h,a=graph(m,k);q=2*m-(S*S/(kap+n)).sum(1)
        x=1-a*S[u,g]/(kap+n[u,g]);y=1-a*S[v,h]/(kap+n[v,h])
        full=y*y*q[u]-x*x*q[v]>1e-10*np.maximum(q[u],q[v]);local=y>x+1e-10
        for name,good in [('full',full),('local',local)]:
            N=len(q);i,j=u[good],v[good]
            z=coo_matrix((np.ones(len(i),dtype=np.int8),(i,j)),shape=(N,N)).tocsr()
            _,labels=connected_components(z,directed=True,connection='strong')
            sinks=np.bincount(i,minlength=N)==0
            print('ablation',m,k,kap,name,'edges',sum(good),'harmful',sum(good&(q[v]<q[u]-1e-10)),'cyclic SCC',sum(np.bincount(labels)>1),'sinks',sum(sinks),'bad sinks',sum(sinks&~(S==0).all(1)))
        print('full only',sum(full&~local),'local only',sum(local&~full))
        assert not np.any(full&(q[v]<q[u]-1e-10)&(n[v,h]<=n[u,g]))
    # Demonstrate the discrete size factor identity symbolically.
    import sympy as sp
    d,Sg=sp.symbols('d Sg', positive=True)
    assert sp.factor(1-Sg**2/d+(Sg-1)**2/(d-1)-d/(d-1)*(1-Sg/d)**2)==0
    print('leave-one-out identity: symbolic check passed')
    # A strict payoff improvement with lower individual contribution and total information value.
    m,k,kap=8,4,1.25;S,n,u,v,g,h,a=graph(m,k);q=2*m-(S*S/(kap+n)).sum(1)
    x=1-a*S[u,g]/(kap+n[u,g]);y=1-a*S[v,h]/(kap+n[v,h]);full=y*y*q[u]-x*x*q[v]>1e-10*np.maximum(q[u],q[v]);j=np.flatnonzero(full&(y<x-1e-10))[0];i=u[j];o=v[j]
    print('declining local and total information example: type',a[j],'from',g[j],'to',h[j],'A before',((n[i]+S[i])/2).astype(int),'B before',((n[i]-S[i])/2).astype(int),'A after',((n[o]+S[o])/2).astype(int),'B after',((n[o]-S[o])/2).astype(int),'q',q[i],q[o],'x',x[j],y[j],'normalized payoff',x[j]/np.sqrt(q[i]),y[j]/np.sqrt(q[o]))
