#Integration using differernt method
import numpy as np 
import math as m
from scipy.special import roots_legendre, eval_legendre

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

def GL(f, a, b, n):  
    def F(y):
        x=(b-a)/2*y+(b+a)/2 
        return f(x)*(b-a)/2

    T=roots_legendre(n)
    sum=0
    for i in range(0,n):
        sum+=F(T[0][i])*T[1][i]
    return sum



N=[10,40,100,1000]

#First integration using three methods

a,b=0,3

def F1(x):
    return 1/(2+x**2)

I1_Trap, I1_Simp, I1_GL=[],[],[]

for n in N:
    I1_Trap.append(Trap(F1,a,b,n))
    I1_Simp.append(Simp(F1,a,b,n))
    I1_GL.append(GL(F1,a,b,n))


#print(I1_Trap,I1_Simp,I1_GL)

#Second integration using three methods

e=10**-5 #Offset term to remove singular point
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
    return GL(f,0,1,n)

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
    I2_GL.append(GL(PxGL,a,b,n))

print("Trapezoidal method",I2_Trap)
print("Simpson method",I2_Simp)
print("Gauss-Legendre method",I2_GL)
