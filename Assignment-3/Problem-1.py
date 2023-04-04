import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

BM=[0.938,1.440,1.535,1.650,1.710] #These are Baryon masses 
MM=[0.782,1.170,1.282,1.420,1.512] #These are Meson masses 

m0=0.5

def Rho(m,A,TH):
    return A/(m**2+m0**2)**(5/4)*np.exp(m/TH) 

gB,gM=4,3 

def Theta(x):
    if x>=0:
       v=1
    else:
        v=0
    return v 

def Ne_s(m):
    S=0
    for M in BM:
        S+=gB*Theta(m-M)
    for M in MM:
        S+=gM*Theta(m-M)
    return S 

def Ne(m):
    if isinstance(m,list):
        S=[]
        for x in m:
            S.append(Ne_s(x))
        return S
    else:
        return Ne_s(m)


def Ne(m):
    L=[]
    for x in m:
        L.append(Ne_s(x))
    return L 

def GL(f, a, b, n):
    x,w = np.polynomial.legendre.leggauss(n)
    xp=0.5*(b-a)*x+0.5*(b+a)
    wp=0.5*(b-a)*w
    I=np.sum(wp*f(xp))
    return I

def Nth(m, A, TH):
    if isinstance(m, list):
        s = []
        for i in m:
            def f(x):
                return Rho(x, A, TH)
            s.append(GL(f, 0, i, 100))
        return s
    else:
        def f(x):
            return Rho(x, A, TH)
        s = GL(f, 0, m, 100)
        return s

def Chi_sq(m,A,TH):
    return 5*(np.log(Ne(m))-np.log(Nth(m,A,TH)))**2/(np.log(10))**2

x=np.linspace(0.7,2,100) 
Y=Ne(x)

popt, pcov=curve_fit(Nth, x, Y)

A,TH=popt[0],popt[1]
print(popt)

#plt.plot(x,Chi_sq(x,1,1),'b')
plt.plot(x,Nth(x.tolist(),A,TH),'g')
plt.plot(x,Y,'r')
plt.show()