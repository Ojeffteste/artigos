# Verificação: P ≈ π(a+b)·(2/√π)·Γ(s+3/2)/Γ(s+1),  s = h·N(h)/D(h)
import mpmath as mp
mp.mp.dps=40
N=[1/(8*(1-mp.log(2)))]+[mp.mpf(v) for v in ['-0.5017006478790753','0.1280463626575522','0.006629020704070353']]
D=[mp.mpf(v) for v in ['1','-1.404769353861877','0.5261406749235058','-0.04069901177880171']]
def P_nova(a,b):
    a=mp.mpf(a);b=mp.mpf(b); h=((a-b)/(a+b))**2
    s=h*mp.polyval(N[::-1],h)/mp.polyval(D[::-1],h)
    return mp.pi*(a+b)*2/mp.sqrt(mp.pi)*mp.gamma(s+mp.mpf(3)/2)/mp.gamma(s+1)
def P_exact(a,b): return 4*a*mp.ellipe(1-(b/a)**2) if b>0 else 4*a
R2=lambda a,b: mp.pi*(a+b)*(1+3*((a-b)/(a+b))**2/(10+mp.sqrt(4-3*((a-b)/(a+b))**2)))
print('c0 exato 1/(8(1-ln2)) =',mp.nstr(1/(8*(1-mp.log(2))),12))
xs=[mp.mpf(i)/4000 for i in range(0,4001)]+[mp.mpf(10)**-k for k in range(4,16)]
e=[(abs(P_nova(1,x)/P_exact(1,x)-1),x) for x in xs]
print('max geral',mp.nstr(max(e)[0],3),'em b/a',mp.nstr(max(e)[1],3))
print('max b/a>=0.2',mp.nstr(max(t for t in e if t[1]>=0.2)[0],3),' Ram2:',mp.nstr(max(abs(R2(1,x)/P_exact(1,x)-1) for x in xs if x>=0.2),3))
print('max b/a>=0.1',mp.nstr(max(t for t in e if t[1]>=0.1)[0],3),' Ram2:',mp.nstr(max(abs(R2(1,x)/P_exact(1,x)-1) for x in xs if x>=0.1),3))
print('série de f(h) vs exata:')
fn=lambda h: 2/mp.sqrt(mp.pi)*mp.gamma(h*mp.polyval(N[::-1],h)/mp.polyval(D[::-1],h)+mp.mpf(3)/2)/mp.gamma(h*mp.polyval(N[::-1],h)/mp.polyval(D[::-1],h)+1)
print([mp.nstr(v,10) for v in mp.taylor(fn,0,8)])
print([mp.nstr(mp.binomial(mp.mpf(1)/2,n)**2,10) for n in range(9)])
for x in ['0.9','0.7','0.5','0.4','0.3','0.2','0.1','0.05','0.01','0']:
    x=mp.mpf(x); ex=P_exact(1,x); print(mp.nstr(x,3), mp.nstr(P_nova(1,x)/ex-1,3), mp.nstr(R2(1,x)/ex-1,3))
