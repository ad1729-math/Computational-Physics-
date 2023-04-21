import numpy as np 

#Brute Force Monte Carlo Method for evaluating integration
N=10**7

def f(x):
    return np.exp(-x**2)

def BF_M(f,a,b,N):
    def F(y):
        return (b-a)*f(a+(b-a)*y)     
    s=0
    for i in range(N):
        y=np.random.rand()
        s+=F(y)
    return s/N

print("The Integration of the funtion found using Brute-Force Monte Carlo method is",BF_M(f,0,1,N))

#Improved Monte Carlo Integration 
N=[10**3,10**4,10**5]

A=(1-1/np.exp(1))**(-1)

def p(x):
    return A*np.exp(-x)

def y(x):
    return -np.log(1-x/A)

I=[]
for n in N:
    s=0
    for j in range(n):
        x=np.random.rand()
        s+=f(y(x))/p(y(x)) 
    I.append(s/n)

print("The Integration using Importance sampling is",I)