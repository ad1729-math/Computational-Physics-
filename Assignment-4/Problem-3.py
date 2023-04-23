import numpy as np 
import math as m
import matplotlib.pyplot as plt

k,M=1,1
t0,x0,x10=0,1,0
tf=12*m.pi
E0=0.5*(M*x10**2+k*x0**2)
N=10**4

def f(t,x,x1):
    return -k/M*x

def RK4(f,t,x,x1,h):
    k1=h*x1
    j1=h*f(t,x,x1)
    k2=h*(x1+j1/2)
    j2=h*f(t+h/2,x+k1/2,x1+j1/2)
    k3=h*(x1+j2/2)
    j3=h*f(t+h/2,x+k2/2,x1+j2/2)
    k4=h*(x1+j3)
    j4=h*f(t+h,x+k3,x1+j3)
    x1+=1/6*(j1+2*j2+2*j3+j4)
    x+=1/6*(k1+2*k2+2*k3+k4)
    t+=h
    return t,x,x1

def Deq(f,t0,tf,x0,x10): #This function returns the values of x(t) upto t_f.
    T,X,X1=[],[],[]
    h=(tf-t0)/N
    for i in range(N):
        L=RK4(f,t0,x0,x10,h)
        T.append(L[0])
        X.append(L[1])
        X1.append(L[2])
        t0,x0,x10=RK4(f,t0,x0,x10,h)
    return T,X,X1

def E(T): #This returns the value/ listed value of energy and delta E as a list.
    if isinstance (T,list):
        X,X1=Deq(f,t0,tf,x0,x10)[1],Deq(f,t0,tf,x0,x10)[2]
        L=[]
        dL=[]
        for i in range(len(T)):
            x,x1=X[i],X1[i]
            L.append(0.5*(M*x1**2+k*x**2))
            dL.append(E0-0.5*(M*x1**2+k*x**2))
        return L,dL
    else:
        i=int((T-t0)/(tf-t0)*N)
        x,x1=Deq(f,t0,tf,x0,x10)[1][i],Deq(f,t0,tf,x0,x10)[2][i]
        E=0.5*(M*x1**2+k*x**2)
        return E,E0-E
              

# plt.subplot(2,1,1)
plt.plot(Deq(f,t0,tf,x0,x10)[0],Deq(f,t0,tf,x0,x10)[1],'r',label="$x(t)$") #Plots x(t) versus t
plt.xlabel("$t$ (Time) --->")
plt.ylabel("$x(t)$--->")
plt.legend()
# # plt.subplot(2,1,2)
# plt.plot(Deq(f,t0,tf,x0,x10)[0],E(Deq(f,t0,tf,x0,x10)[0])[1],'g',label="$\Delta(t)$") #Plots Delta E versus t
# plt.xlabel("$t$ (Time) --->")
# plt.ylabel("$\Delta(t)$--->")
# plt.legend()
plt.show()


