# Verificação: perímetro da elipse via harmônicos toroidais e funções de Bessel (aproximação uniforme de Olver)
import mpmath as mp
mp.mp.dps=40
P_=[mp.mpf(0),mp.mpf(7)/320]+[mp.mpf(v) for v in ['0.0027421486237688702','8.0863833392703183e-5','1.4084949212793479e-6','2.9916307576727871e-9']]
Q_=[mp.mpf(v) for v in ['1','0.25525899377591932','0.024338672837890872','0.0010607656517810025','2.010001314444072e-5','1.1066450226912341e-7']]
def P_bessel(a,b,corr=True):
    a=mp.mpf(a); b=mp.mpf(b)
    if a==b: return 2*mp.pi*a
    xi=mp.log(a/b); t=xi*xi
    B=(xi*mp.coth(xi)-1)/(8*t)
    phi=mp.polyval(P_[::-1],t)/mp.polyval(Q_[::-1],t) if corr else 0
    return 2*mp.pi*mp.sqrt(a*b)*mp.sqrt(xi/mp.sinh(xi))*(mp.besseli(0,xi)+xi*mp.besseli(1,xi)*B*(1+phi))
def P_exact(a,b): return 4*a*mp.ellipe(1-(b/a)**2) if b>0 else 4*a
R2=lambda a,b: mp.pi*(a+b)*(1+3*((a-b)/(a+b))**2/(10+mp.sqrt(4-3*((a-b)/(a+b))**2)))
print('identidade toroidal:',mp.nstr(2*mp.pi*mp.sqrt(0.3)*mp.legenp(0.5,0,(1+0.09)/(0.6),type=3),20),mp.nstr(P_exact(1,mp.mpf('0.3')),20))
xs=[mp.mpf(i)/2000 for i in range(1,2000)]+[mp.mpf(10)**-k for k in range(4,30)]
E={x:P_exact(1,x) for x in xs}
def rng(F,lo,hi): return mp.nstr(max(abs(F(1,x)/E[x]-1) for x in xs if lo<=x<=hi),3)
for nm,F in [('nova',P_bessel),('Olver 2 termos (sem correção)',lambda a,b:P_bessel(a,b,False)),('Ramanujan II',R2)]:
    print(nm,'| b/a>=0.2:',rng(F,0.2,1),'| 0.05-0.2:',rng(F,0.05,0.2),'| 0.01-0.05:',rng(F,0.01,0.05),'| 1e-3-1e-2:',rng(F,1e-3,1e-2),'| 1e-6-1e-3:',rng(F,1e-6,1e-3),'| <1e-6:',rng(F,0,1e-6))
print('b/a  nova  Olver2  Ram2')
for x in ['0.9','0.7','0.5','0.3','0.2','0.1','0.05','0.01','0.001','1e-6']:
    x=mp.mpf(x); ex=P_exact(1,x)
    print(mp.nstr(x,3),*[mp.nstr(F(1,x)/ex-1,3) for F in (P_bessel,lambda a,b:P_bessel(a,b,False),R2)])
