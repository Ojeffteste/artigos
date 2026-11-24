import mpmath as mp
mp.mp.dps=30
def f_exact(h):
    if h==1: return 4/mp.pi
    t=mp.sqrt(h); b=(1-t)/(1+t)
    return 4*mp.ellipe(1-b*b)/(mp.pi*(1+b))
hs=[mp.mpf(i)/5000 for i in range(0,5001)]+[1-mp.mpf(10)**-k for k in range(4,25)]
E={h:f_exact(h) for h in hs}
def mx(F):
    e=[(F(h)/E[h]-1,h) for h in hs]; return max(e,key=lambda t:abs(t[0]))
P4=4/mp.pi; K=P4-mp.mpf(81)/64
def Lg(h): return mp.log(1-h) if h<1 else mp.mpf(0)
J=lambda h: 1+h/4+h**2/64+K*h**3+(P4-mp.mpf(325)/256)*h**2*(1-h)*Lg(h)
def E2(a3,g1,g2):
    return lambda h: 1+h/4+h**2/64+a3*h**3+(K-a3)*h**4+h**2*(1-h)*Lg(h)*(g1+g2*(1-h))
R2=lambda h: 1+3*h/(10+mp.sqrt(4-3*h))
Can=lambda h: R2(h)+(P4-mp.mpf(14)/11)*h**12
for nm,F in [('J',J),('E2 lp',E2(0.00071675,0.00150822,-0.00530526)),('E2 round',E2(mp.mpf('0.000717'),mp.mpf('0.001508'),mp.mpf('-0.005305'))),('R2',R2),('Cantrell',Can)]:
    m=mx(F); print(nm, mp.nstr(m[0],4),'at h=',mp.nstr(m[1],5))
F=E2(mp.mpf('0.000717'),mp.mpf('0.001508'),mp.mpf('-0.005305'))
print('h-series E2',[mp.nstr(x,8) for x in mp.taylor(F,0,5)])
print('exact     ',[mp.nstr(mp.binomial(0.5,n)**2,8) for n in range(6)])
print('b/a | J | E2 | Ram2 | Cantrell')
for ba in ['0.99','0.9','0.7','0.5','0.3','0.2','0.1','0.05','0.01','0.001','0']:
    ba=mp.mpf(ba); h=((1-ba)/(1+ba))**2; ex=f_exact(h)
    print(mp.nstr(ba,3),*[mp.nstr(G(h)/ex-1,3) for G in (J,F,R2,Can)])
