import mpmath as mp
mp.mp.dps=40
def P_exact(a,b): return 4*a*mp.ellipe(1-(b/a)**2) if b>0 else 4*a
QC=[1,2,15,150,1707,20910]          # q = sum QC[i] lam^(4i+1)
SIG=[1,3,4,7,6,12,8,15]              # sigma_1(n)
def nome(lam,nq): return sum(QC[i]*lam**(4*i+1) for i in range(nq))
def th3sq(q,nt): return (1+2*sum(q**(j*j) for j in range(1,nt+1)))**2
def E2(Q,ne): return 1-24*sum(SIG[n-1]*Q**n for n in range(1,ne+1))
def P_mod(a,b,nq,nt,ne):
    a=mp.mpf(a);b=mp.mpf(b); r=b/a
    if r*r>=mp.mpf(1)/2:   # ramo A
        lam=(1-mp.sqrt(r))/(2*(1+mp.sqrt(r))); q=nome(lam,nq); T=th3sq(q,nt)
        return 2*mp.pi*a/3*(T*(1+r*r)+E2(q*q,ne)/T)
    else:                  # ramo B
        e=mp.sqrt(1-r*r)
        mu=r*r/(2*(1+e)*(1+mp.sqrt(e))**2); p=nome(mu,nq); T=th3sq(p,nt)
        return 4*a*(1/T+mp.log(1/p)/6*(T*(1+r*r)-E2(p*p,ne)/T))
xs=[mp.mpf(i)/4000 for i in range(1,4001)]+[mp.mpf(10)**-k for k in range(4,40,3)]
E={x:P_exact(1,x) for x in xs}
for nq,nt,ne in [(1,1,1),(1,1,2),(2,1,2),(2,1,3),(2,2,4),(3,2,5),(3,3,8)]:
    e=[(abs(P_mod(1,x,nq,nt,ne)/E[x]-1),x) for x in xs]
    print(nq,nt,ne,'max',mp.nstr(max(e)[0],3),'at',mp.nstr(max(e)[1],3))
print('exact check (many terms):',mp.nstr(P_mod(1,mp.mpf('0.3'),6,4,8)/P_exact(1,mp.mpf('0.3'))-1,5))
R2=lambda a,b: mp.pi*(a+b)*(1+3*((a-b)/(a+b))**2/(10+mp.sqrt(4-3*((a-b)/(a+b))**2)))
print('b/a  simples(1,1,2)  precisa(2,2,4)  Ram2')
for x in ['0.99','0.9','0.8','0.71','0.7','0.5','0.3','0.2','0.1','0.01','0.001','1e-6','1e-12']:
    x=mp.mpf(x); ex=P_exact(1,x)
    print(mp.nstr(x,3), mp.nstr(P_mod(1,x,1,1,2)/ex-1,3), mp.nstr(P_mod(1,x,2,2,4)/ex-1,3), mp.nstr(R2(1,x)/ex-1,3))
