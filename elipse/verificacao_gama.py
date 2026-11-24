import mpmath as mp
mp.mp.dps=40
def f_exact(h):
    if h==1: return 4/mp.pi
    t=mp.sqrt(h); b=(1-t)/(1+t)
    return 4*mp.ellipe(1-b*b)/(mp.pi*(1+b))
Q=lambda s: 2/mp.sqrt(mp.pi)*mp.gamma(s+mp.mpf(3)/2)/mp.gamma(s+1)
# exact series of u=s*/h by series reversion with mpmath
fser=lambda h: mp.hyp2f1(-0.5,-0.5,1,h)
mp.mp.dps=60
def sst(h): return mp.findroot(lambda s: Q(s)-fser(h), h*mp.mpf('0.4'))
us=mp.taylor(lambda h: sst(h)/h, mp.mpf('1e-25'),4, method='step') if False else None
# do reversion symbolically: Q(s)=sum q_k s^k, f=sum f_k h^k
q=mp.taylor(Q,0,6); fk=[mp.binomial(mp.mpf(1)/2,n)**2 for n in range(7)]
# unknown s = sum a_k h^k, solve order by order
a=[mp.mpf(0)]*6
for n in range(1,6):
    def coef(an):
        a[n]=an
        return mp.taylor(lambda h: Q(sum(a[k]*h**k for k in range(1,n+1))),0,n)[n]-fk[n]
    a[n]=mp.findroot(coef,mp.mpf('0.01'))
print('s series a1..a5:',[mp.nstr(x,12) for x in a[1:]])
print('s0 closed', mp.nstr(1/(8*(1-mp.log(2))),12))
mp.mp.dps=30
s0,u1,u2,u3=[mp.mpf(mp.nstr(x,10)) for x in a[1:5]]
c4=mp.mpf('0.012317'); c5=mp.mpf('-0.026923'); c6=mp.mpf(1)/2-(s0+u1+u2+u3+c4+c5)
print('constants',s0,u1,u2,u3,c4,c5,'c6',mp.nstr(c6,8))
co=[s0,u1,u2,u3,c4,c5,c6]
F=lambda h: Q(h*sum(co[k]*h**k for k in range(7)))
hs=[mp.mpf(i)/5000 for i in range(0,5001)]+[1-mp.mpf(10)**-k for k in range(4,22)]
e=[(F(h)/f_exact(h)-1,h) for h in hs]; m=max(e,key=lambda t:abs(t[0]))
print('MAX rel err',mp.nstr(m[0],4),'at h',mp.nstr(m[1],5))
R2=lambda h: 1+3*h/(10+mp.sqrt(4-3*h)); Can=lambda h:R2(h)+(4/mp.pi-mp.mpf(14)/11)*h**12
print('b/a  nova  Ram2  Cantrell')
for ba in ['0.9','0.7','0.5','0.3','0.2','0.1','0.05','0.01','0.001','0']:
    ba=mp.mpf(ba); h=((1-ba)/(1+ba))**2; ex=f_exact(h)
    print(mp.nstr(ba,3),*[mp.nstr(G(h)/ex-1,3) for G in (F,R2,Can)])
print('series F',[mp.nstr(x,10) for x in mp.taylor(F,0,6)])
print('exact   ',[mp.nstr(x,10) for x in fk])
