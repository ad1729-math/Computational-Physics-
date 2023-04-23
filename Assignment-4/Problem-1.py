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
    for i in range(-1,20):
        if Cu_Po(i,l)<x<Cu_Po(i+1,l):
            return i+1
        else:
            continue

def PP(l0):
    x=np.random.rand()
    return P_b(x,l0)

#Decay
N=1000
dt=0.01
l=0.3

def Decay(n,dt):
    l0=n*dt*l
    return PP(l0)

def Exp(N,dt):
    i=0
    T=[]
    Dec=[]
    Res=[]
    while N>0:
        T.append(i*dt)
        s=Decay(N,dt)
        if s<=N:
           Dec.append(s)
           N=N-s
           Res.append(N)
        else:
            continue
        i+=1
    return T,Dec,Res


T,Dec,Res=Exp(N,dt)

plt.subplot(2,1,1)
plt.plot(T,Dec,'ro',label="Number of decays")
plt.xlabel("$t$(Time)--->")
plt.ylabel("Number of decays")
plt.legend()
plt.subplot(2,1,2)
plt.plot(T,Res,'bo',label="Number of esidues")
plt.xlabel("$t$(Time)--->")
plt.ylabel("Number of particles left--->")
plt.legend()


#Ensemble avergae at time t
M=1000 #Number of ensembles
N0=10 #Number of initital nucleuses
En=[]
for i in range(M):
    En.append(Exp(N0,dt)[2])

En_av=[]
N_ln=[]
Sl=[]
T=np.arange(0,2,dt)
for i in range(len(T)):
    s=0
    for k in range(M):
        s+=En[k][i]
    av=s/M
    En_av.append(av)
    N_ln.append(np.log(av))

for j in range(len(T)-1):
    v=(En_av[j+1]-En_av[j])/(dt*N0)
    Sl.append(v)
Sl.append(0)

#Determining the slope
s=0
for i in range(len(T)):
    s+=(N_ln[i]-N_ln[-1])*dt 

m=-2*s/(T[-1]**2)
print("The slope of the curve is",m)

plt.subplot(2,1,1)
plt.plot(T,N_ln,"ro",label="Log Ensemble average, $N_0$=1000")
plt.xlabel("$t$(Time)--->")
plt.ylabel("$\ln(N(t))$ (Log Ensemble average)--->")
plt.legend()
plt.subplot(2,1,2)
plt.plot(T,Sl,'g')
plt.xlabel("$t$ (Time)--->")
plt.ylabel("Slope--->")
plt.show()


