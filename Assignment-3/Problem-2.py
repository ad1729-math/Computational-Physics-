#Integration using differernt methods
import numpy as np 
import math as m
from Gauss_Legendre import GL

W10=np.loadtxt("W_10.dat",unpack=True)
R10=np.loadtxt("R_10.dat",unpack=True)
W40=np.loadtxt("W_40.dat",unpack=True)
R40=np.loadtxt("R_40.dat",unpack=True)
W100=np.loadtxt("W_100.dat",unpack=True)
R100=np.loadtxt("R_100.dat",unpack=True)
W1000=np.loadtxt("W_1000.dat",unpack=True)
R1000=np.loadtxt("R_1000.dat",unpack=True)

N=[10,40,100,1000]
W=[W10,W40,W100,W1000]
R=[R10,R40,R100,R1000]

def Trap(f,a,b,n):
    h=abs(b-a)/n
    sum=(f(a)+f(b))*h/2
    for k in range(1,n):
        sum+=f(a+k*h)*h 
    return sum 


def Simp(f,a,b,n):  #n is even
    h=abs(b-a)/n
    sum=(f(a)+f(b))*h/3
    for k in range(1,int(n/2)):
        sum+=(4*f(a+(2*k-1)*h)+2*f(a+2*k*h))*h/3

    sum+=f(b-h)*4*h/3
    return sum

def GL_(f, a, b, n): 
    if n in N:
       def F(y):
           x=0.5*((b-a)*y+(b+a))
           return f(x)*(b-a)/2
       s=0
       v=N.index(n)
       for k in range(n):
           s+=W[v][k]*F(R[v][k])
       return s
           

#First integration using three methods

a,b=0,3

def F1(x):
    return 1/(2+x**2)

I1_Trap, I1_Simp,I1_GL=[],[],[]

for n in N:
    I1_Trap.append(Trap(F1,a,b,n))
    I1_Simp.append(Simp(F1,a,b,n))
    I1_GL.append(GL_(F1,a,b,n))


print("The first integration by three different methods are")
print("Trapezoida method",I1_Trap)
print("Simpson's method",I1_Simp)
print("Gauss-Legendre method",I1_GL)

#Second integration using three methods

e=10**-6 #Offset term to remove singular point
a,b=0,1-e


def F2(x,y):
    return 1/(1-x*y)

def Px_trap(x,n):
    def f(y):
        return F2(x,y)
    return Trap(f,0,1,n)

def Px_simp(x,n):
    def f(y):
        return F2(x,y)
    return Simp(f,0,1,n)

def Px_GL(x,n):
    def f(y):
        return F2(x,y)
    return GL_(f,0,1,n)

I2_Trap,I2_Simp,I2_GL=[],[],[]

for n in N:
    def PxT(x):
        return Px_trap(x,n)
    def PxS(x):
        return Px_simp(x,n)
    def PxGL(x):
        return Px_GL(x,n)
    
    I2_Trap.append(Trap(PxT,a,b,n))
    I2_Simp.append(Simp(PxS,a,b,n))
    I2_GL.append(GL_(PxGL,a,b,n))

print("The second integration found by three different methods")
print("Trapezoidal method",I2_Trap)
print("Simpson method",I2_Simp)
print("Gauss-Legendre method",I2_GL)
