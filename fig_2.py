# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as si
from scipy.optimize import minimize
from tabulate import tabulate


def dCdt(C,t,k,n):
    return k*C**n

def r2(data):
    Ct=data[0]
    cm=data[1]
    xm=sum(Ct)/len(Ct)
    ym=sum(Cm)/len(Cm)
    Sxy=0
    Sx=0
    Sy=0
    for i in range(len(Ct)):
        Sxy+=(Ct[i]-xm)*(Cm[i]-ym)
        Sx+=(Ct[i]-xm)**2
        Sy+=(Cm[i]-ym)**2
    r2=Sxy/(Sx**0.5*Sy**.5)
    return r2[0]
#
#    Gerar ponto experimentais
#
k=-1;Co=1;n=3;erro=0.3
t=np.arange(0.,10.,0.1)
Ct=si.odeint(dCdt,Co,t,args=(k,n))
Ct=[(Ct[i]*(1+erro* np.random.random()))for i in range(len(Ct))]
##
##    Representar os pontos experimentais
##
fig=plt.figure(figsize=(5.5,5))
plt.title(u'ODE     $\\frac{dC}{dt}=b\\cdot C^c$ com $C(0)=a$')
plt.ylabel('$[M]$')
plt.xlabel('t (s)')
plt.plot(t,Ct,'r.', label='$C(t)_{exp}$')
plt.grid('on')
#
#    Modelo
#
#var('a b c')
def model(parm):
    a,b,c=parm
    y=si.odeint(dCdt,a,t,args=(b,c))
    srq=sum((y[i]-Ct[i])**2 for i in range(len(Ct)))
    return srq
#
#    optimização com minimize
# 
print '############################################'
sol=minimize(model,[2,-0.1,1],method='Nelder-Mead')
Com=sol.x[0]
km=sol.x[1]
nm=sol.x[2]
Cm=si.odeint(dCdt,Com,t,args=(km,nm))
print '############################################'
table=[['Dados',str(Co),' ',str(k),' ',str(n),'  ','%.2f'%(erro*100)],
['Obtidos','%.4f'%Com,' ','%.4f'%km,' ','%.4f'%nm,' ','r=%.4f'%r2([Ct,Cm])]]
print tabulate(table,headers = [' ','Co',' ','k',' ','n','  ','%erro'])
#
#    representação do ajuste
# 
plt.plot(t,Cm,'g',lw=2,label='$C(t)_{mod}$')
plt.legend(loc=0)
plt.savefig('fig.png')
plt.show()
plt.close()