import numpy as np 
import math 
import matplotlib.pyplot as plt 

m,g,l=1,1,1
Mesh=[100,1000,100000]
v=5
A=1.2
w,v1=2/3,1/2
Tp=2*math.pi/w

def F(t):
    return 0

w0=math.sqrt(g/l)

def f(t,x,x1,F,v):
    Phi=v*g/(m*l)
    return F(t)/(m*g)-Phi*x1-np.sin(x)  #This is the dimensionless form of the equation

def RK4(f,t,x,x1,h,F,v):
    k1=h*x1
    j1=h*f(t,x,x1,F,v)
    k2=h*(x1+j1/2)
    j2=h*f(t+h/2,x+k1/2,x1+j1/2,F,v)
    k3=h*(x1+j2/2)
    j3=h*f(t+h/2,x+k2/2,x1+j2/2,F,v)
    k4=h*(x1+j3)
    j4=h*f(t+h,x+k3,x1+j3,F,v)
    x1+=1/6*(j1+2*j2+2*j3+j4)
    x+=1/6*(k1+2*k2+2*k3+k4)
    t+=h
    return t,x,x1

def Theta(f,t0,tf,x0,x10,N,F,v): #This function returns the values of x(t) upto t_f.
    T,X,X1=[],[],[]
    h=(tf-t0)/N
    for i in range(N):
        L=RK4(f,t0,x0,x10,h,F,v)
        T.append(L[0])
        X.append(L[1])
        X1.append(L[2])
        t0,x0,x10=RK4(f,t0,x0,x10,h,F,v)
    return T,X,X1

for N in Mesh:
    L=Theta(f,0,20*math.pi,0.2,0,N,F,v)
    T,X=L[0],L[1]
    T1=[x/w0 for x in T]
    plt.plot(T1,X,'r',label="Viscocity="+str(v)+"")
    plt.xlabel("Time--->"+str(N)+"Mesh points")
    plt.ylabel("$Theta (t)--->$")
    plt.legend()
    plt.show()

def F1(t):
    return A*np.cos(w*t)

for N in Mesh:
    L=Theta(f,0,50*Tp,0,0,N,F1,v1)
    T,X=L[0],L[1]
    T1=[x/w0 for x in T]
    plt.plot(T1,X,'g',label="Driven oscillation ; Viscocity="+str(v1)+"")
    plt.xlabel("Time--->"+str(N)+"Mesh points")
    plt.ylabel("$Theta(t)$--->")
    plt.legend()
    plt.show()