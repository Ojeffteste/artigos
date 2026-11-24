# Verificação da fórmula "qui com graus de liberdade efetivos" para o perímetro da elipse
import mpmath as mp
mp.mp.dps=30
def P_exact(a,b): return 4*a*mp.ellipe(1-(b/a)**2) if b>0 else 4*a
al=1/(6-8*mp.log(2))
p=mp.mpf('1.613'); r=mp.mpf('-0.107'); q=al/p-1
def P_chi(a,b):
    a=mp.mpf(a); b=mp.mpf(b)
    W=(2*a*b/(a*a+b*b))**p; u=1-W
    nu=1+W/(1+q*u+r*u*u)
    return 4*mp.sqrt(mp.pi*(a*a+b*b)/nu)*mp.gamma((nu+1)/2)/mp.gamma(nu/2)
# identidade estatística: P = sqrt(8 pi) E||(aX,bY)||
a,b=mp.mpf(1),mp.mpf('0.4')
EN=mp.quad(lambda x,y: mp.sqrt(a*a*x*x+b*b*y*y)*mp.npdf(x)*mp.npdf(y),[-12,0,12],[-12,0,12])
print('identidade:',mp.nstr(mp.sqrt(8*mp.pi)*EN,15),mp.nstr(P_exact(a,b),15))
print('q =',mp.nstr(q,10),' alpha =',mp.nstr(al,10))
R2=lambda a,b: mp.pi*(a+b)*(1+3*((a-b)/(a+b))**2/(10+mp.sqrt(4-3*((a-b)/(a+b))**2)))
xs=[mp.mpf(i)/5000 for i in range(0,5001)]+[mp.mpf(10)**-k for k in range(4,16)]
e=[(P_chi(1,x)/P_exact(1,x)-1,x) for x in xs]; m=max(e,key=lambda t:abs(t[0]))
print('erro relativo máximo:',mp.nstr(m[0],4),'em b/a =',mp.nstr(m[1],4))
print('b/a   nova   Ramanujan II')
for x in ['1','0.99','0.9','0.7','0.5','0.3','0.2','0.1','0.05','0.01','0.001','0']:
    x=mp.mpf(x); ex=P_exact(1,x)
    print(mp.nstr(x,3), mp.nstr(P_chi(1,x)/ex-1,3), mp.nstr(R2(1,x)/ex-1,3))
