# Verificação: perímetro da elipse via polilogaritmos (assintótica dos números de Catalan)
import mpmath as mp
mp.mp.dps=40
def f_exact(h):
    if h==1: return 4/mp.pi
    t=mp.sqrt(h); b=(1-t)/(1+t)
    return 4*mp.ellipe(1-b*b)/(mp.pi*(1+b))
a=lambda n: mp.binomial(mp.mpf(1)/2,n)**2
def coefs(M,K):
    ks=list(range(3,3+K+1)); rows=[];rhs=[]
    for n in range(M+1,M+K+1):
        rows.append([mp.mpf(1)/mp.mpf(n)**k for k in ks]); rhs.append(a(n))
    rows.append([mp.zeta(k)-sum(mp.mpf(1)/mp.mpf(n)**k for n in range(1,M+1)) for k in ks])
    rhs.append(4/mp.pi-sum(a(n) for n in range(M+1)))
    return ks,list(mp.lu_solve(mp.matrix(rows),mp.matrix(rhs)))
def make(M,ks,g):
    def f(h):
        s=sum(a(n)*h**n for n in range(M+1))
        for gi,k in zip(g,ks): s+=gi*(mp.polylog(k,h)-sum(h**n/mp.mpf(n)**k for n in range(1,M+1)))
        return s
    return f
R2=lambda h: 1+3*h/(10+mp.sqrt(4-3*h))
hs=[mp.mpf(i)/2000 for i in range(0,2001)]+[1-mp.mpf(10)**-k for k in range(4,18)]
E={h:f_exact(h) for h in hs}
for nome,M,K in [('simples',0,2),('principal',2,3)]:
    ks,g=coefs(M,K)
    g=[mp.mpf(mp.nstr(x,16)) for x in g]   # constantes com 16 dígitos
    F=make(M,ks,g)
    print('==',nome,'M=',M,' k=',ks)
    for k,x in zip(ks,g): print('  g%d = %s'%(k,mp.nstr(x,16)))
    e=[(abs(F(h)/E[h]-1),h) for h in hs]
    print('  max geral %s ; max b/a>=0.2 %s'%(mp.nstr(max(e)[0],3),mp.nstr(max(t[0] for t in e if t[1]<=mp.mpf(16)/36),3)))
    print('  série:',[mp.nstr(v,10) for v in mp.taylor(F,0,7)])
print('exata  :',[mp.nstr(a(n),10) for n in range(8)])
ks,g=coefs(2,3); F=make(2,ks,g); ks0,g0=coefs(0,2); F0=make(0,ks0,g0)
print('b/a  principal  simples  Ram2')
for x in ['0.9','0.7','0.5','0.3','0.2','0.1','0.05','0.01','0.001','0']:
    x=mp.mpf(x); h=((1-x)/(1+x))**2; ex=f_exact(h)
    print(mp.nstr(x,3),*[mp.nstr(G(h)/ex-1,3) for G in (F,F0,R2)])
