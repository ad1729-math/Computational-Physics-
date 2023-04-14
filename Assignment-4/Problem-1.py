import numpy as np 
from numpy import random 
import matplotlib.pyplot as plt
import math as m

#Poisson distribution of random decay number

def fact(n):
    if n==0:
        return 1
    else:
        p=1
        for i in range(1,n+1):
            p=p*i
        return p

def Poisson(n,l0):
    return l0**n*m.exp(-l0)/fact(n)

def Cu_Po(n,l):
    if n==-1:
        return 0
    else:  
       s=0
       for i in range(0,n+1):
            s+=Poisson(i,l)
       return s

def P_b(x,l):
    for i in range(-1,10):
        if Cu_Po(i,l)<x<Cu_Po(i+1,l):
            return i+1
        else:
            continue

def PP(l0):
    x=np.random.rand()
    return P_b(x,l0)

L=[]
for i in range(0,200):
    L.append(PP(1))

X=[]
Y=[]
for k in range(0,10):
    X.append(k)
    Y.append(L.count(k))

#Decay
N=1000
dt=0.1
l=0.3

def Decay(n,dt):
    l0=n*dt*l
    return PP(l0)

T=np.arange(0,10,dt)
Dec=[]
for t in T:
    s=Decay(N,dt)
    Dec.append(s)
    N=N-s

plt.plot(T,Dec,'ro')
plt.show()

     
