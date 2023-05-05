import numpy as np 
import math as m  
import matplotlib.pyplot as plt 

MeV=1.6*10**-13
fm=10**-15
c=3*10**8

V0=60*MeV
a=1.45*fm
h=197.92*MeV*fm/c
m=938*MeV/c**2


def V(r):
    if r<=a:
        return -V0 
    else: 
        return 0

def f(r,u,u1,E,V):
    s=2*m/h**2*(V(r)-E)*u
    return s


def RK4(f,r,u,u1,h,E,V):
    k1=h*u1
    j1=h*f(r,u,u1,E,V)
    k2=h*(u1+j1/2)
    j2=h*f(r+h/2,u+k1/2,u1+j1/2,E,V)
    k3=h*(u1+j2/2)
    j3=h*f(r+h/2,u+k2/2,u1+j2/2,E,V)
    k4=h*(u1+j3)
    j4=h*f(r+h,u+k3,u1+j3,E,V)
    u1+=1/6*(j1+2*j2+2*j3+j4)
    u+=1/6*(k1+2*k2+2*k3+k4)
    r+=h
    return r,u,u1 


N=1000
def U(r,r0,u0,u10,E,V):
    h=(r-r0)/N
    R,WF,WFS=[r0],[u0],[u10]
    for i in range(N):
        oput=RK4(f,r0,u0,u10,h,E,V)
        r0,u0,u10=oput[0],oput[1],oput[2]
        R.append(r)
        WF.append(u0)
        WFS.append(u10)
    return R,WF,WFS

R,WF=U(5*a,0,0,10**15,-V0*0.2,V)[0],U(0.5*a,0,0,10**15,-V0*0.2,V)[1]
plt.plot(R,WF,'r')
plt.ylim(0,10**-3)
plt.show()



   

        


