import numpy as np 
import math as m 
import numpy.polynomial.legendre as legendre
import numpy.random as random

N=100
# R=np.loadtxt("R_100.dat", unpack=True)
# W=np.loadtxt("W_100.dat", unpack=True)

R,W=legendre.leggauss(N)

def GL(f,a,b): 
    def F(y):
        x=0.5*((b-a)*y+(b+a))
        return f(x)*(b-a)/2
    s=0
    for k in range(N):
        s+=W[k]*F(R[k])
    return s

ar,br=0,10*m.sqrt(3)
ax,bx=0,m.pi

#Using the symmetry of the integrand we reduce the 6-dimensional integration
#to 3-dimensional. 

def f(x,r1,r2):
    return (r1*r2)**2*np.exp(-2*(r1+r2))/np.sqrt(r1**2+r2**2-2*r1*r2*np.cos(x))*np.sin(x)


def f1(r1,r2):
    def Temp(x):
        return f(x,r1,r2)
    return GL(Temp,ax,bx)

def f2(r2):
    def Temp(x):
        return f1(x,r2)
    return GL(Temp,ar,br)

I=8*m.pi**2*GL(f2,ar,br) #This 8pi^2 factor comes from phi1, phi2 and theta2 intgeral
print("Value of the integral in question using Gauss quadrature is", I)

#Using Brute force Monte-Carlo Method: In this method we again use the 
#function that is simplified using spherical symmetry

N1,N2,Nx=200,200,200

def G(x,r1,r2):
    X=ax+(bx-ax)*x
    R1=ar+(br-ar)*r1
    R2=ar+(br-ar)*r2 
    return (br-ar)**2*(bx-ax)*f(X,R1,R2)
    
S=0
for i in range(N1):
    for j in range(N2):
        for k in range(Nx):
            r10,r20,x0=random.rand(), random.rand(), random.rand()
            S+=G(x0,r10,r20)

I=8*m.pi**2*S/(N1*N2*Nx)
print("The value of the integration using Brute force Monte caarlo method is",I)


#Importance Sampling

