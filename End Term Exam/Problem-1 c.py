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
# M=1000
# E0=-V0*0.5
# e=(0-E0)/M
# for i in range(M):
#     E=E0+e*i
#     k=np.sqrt(2*m*(V0+E))/h
#     u10=A*k 
#     WF_r=U(r,0,0,u10,E,V)[2][-1]
#     if WF_r<10**(-19):
#         print("Eigenvalue using shooting method is",E/MeV)
#         break
#     else:
#         continue 

# r=10*a
# def u(E):
#     k=np.sqrt(2*m*(V0+E))/h
#     u10=A*k 
#     s=U(r,0,0,u10,E,V)[2][-1]
#     return s

# e0=-V0*0.4 
# h1=V0*10**-4
# for i in range(50):
#     u1=(u(e0+h1)-u(e0))/h1
#     e0-=u(e0)/u1

# print("The energy eigenvalue is", e0/MeV)


# E=np.linspace(-V0*0.3,0,100)
# plt.plot(E,u(E),'r',label="u(10a)")
# plt.plot(E,E*0,'b')
# plt.xlabel("Energy--->")
# plt.ylabel("u(10a)--->")
# plt.legend()
# plt.show()


# E0=-15.75*MeV
# k=np.sqrt(2*m*(V0+E0))/h
# L=U(5*a,0,0,k*A,E0,V)
# R,WF=L[0],L[1]
# plt.plot(R,WF,'b',label="Energy=-15.75 MeV")
# plt.xlabel("r--->")
# plt.ylabel("$u(r)$--->")
# plt.legend()
# plt.show()


#The new potential 

def V1(r):
    a,b,c=1,4,7
    Va,Vb,Vc=-10.463*MeV,-1650.6*MeV,6484.3*MeV
    x=0.7/fm*r
    return (Va*np.exp(-a*x)+Vb*np.exp(-b*x)+Vc*np.exp(-c*x))/x


A=1
V10=60*MeV
r=10*a
def u1(E):
    k=np.sqrt(2*m*(V0+E))/h
    u10=A*k 
    s=U(r,0.2*a,0,u10,E,V1)[2][-1]
    return s

e0=-V0*0.5
h1=V0*10**-4
for i in range(50):
    u10=(u1(e0+h1)-u1(e0))/h1
    e0-=u1(e0)/u10

print("The energy eigenvalue is", e0/MeV)

E=np.linspace(-V0*0.12,0,100)
plt.plot(E,u1(E),'r',label="u(10a)")
plt.plot(E,E*0,'b')
plt.xlabel("Energy--->")
plt.ylabel("u(10a)--->")
plt.legend()
plt.show()

# E10=-5.94*MeV
# k1=np.sqrt(2*m*(-E10))/h
# L=U(0.2*a,5*a,0,-k1*A,E10,V1)
# R,WF=L[0],L[1]
# plt.plot(R,WF,'b',label="Energy=-5.939.. MeV")
# plt.xlabel("Energy--->")
# plt.ylabel("$u(r)$--->")
# plt.legend()
# plt.show()

    
# M=1000
# E0=-V10*0.5
# e=(0-E0)/M
# for i in range(M):
#     E=E0+e*i
#     k1=np.sqrt(2*m*(-E))/h
#     u10=A*k1
#     WF_r=U(0.1*a,10*a,0,u10,E,V1)[2][-1]
#     if WF_r<1:
#         print("Eigenvalue using shooting method is",E/MeV)
#         break
#     else:
#         continue 






   

        


