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
        R.append(r0)
        WF.append(u0)
        WFS.append(u10)
    return R,WF,WFS


A=1#This is the amplitude of the wavefunction the value of which doesn't 
# matter. Derivative of u at r=0 is, u'(0)=Ak.
M=1000
E0=-V0*0.5
e=(0-E0)/M
r=5*a
for i in range(M):
    E=E0+e*i
    k=np.sqrt(2*m*(V0+E))/h
    u10=A*k 
    WF_r=U(r,0,0,u10,E,V)[2][-1]
    if abs(WF_r<10**(-18)):
        print(E/MeV)
        break
    else:
        continue 


def V1(r):
    a,b,c=1,4,7
    Va,Vb,Vc=-10.463*MeV,-1650.6*MeV,6484.3*MeV
    x=0.7/fm*r
    return (Va*np.exp(-a*x)+Vb*np.exp(-b*x)+Vc*np.exp(-c*x))/x

A=1
E0=-10*MeV
V10=60*MeV
# k=np.sqrt(2*m*(V10+E0))/h
# L=U(0.2*a,5*a,0,k*A,E0,V1)
# R,WF=L[0],L[1]
# plt.plot(R,WF,'r')
# plt.show()

    
M=1000
E0=-V10*0.5
e=(0-E0)/M
for i in range(M):
    E=E0+e*i
    k=np.sqrt(2*m*(V0+E))/h
    u10=A*k 
    WF_r=U(0.2*a,5*a,0,u10,E,V1)[2][-1]
    if abs(WF_r<10**(0)):
        print(E/MeV)
        break
    else:
        continue 






   

        


