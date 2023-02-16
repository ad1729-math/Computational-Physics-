import numpy as np 
import math as m 
import scipy.special as sp 
import matplotlib.pyplot as plt 

#Up method 

def j0(x):
    return m.sin(x)/x 

def j1(x):
    return (m.sin(x)-x*m.cos(x))/x**2

X=np.arange(0.1,50,0.2)
J10=[]

for x in X:
    a=j0(x)
    b=j1(x)
    for i in range(1,10):
        temp=b
        b=(2*i+1)*b/x-a
        a=temp 
    J10.append(b)


#Down method: Here we start from J50, J49 to which we assign some random value 

#First we take a guess and scale as per the output  

def j100(x):
    return 1

def j99(x):
    return 1

J1_test=[]
for x in X: 
    a=j100(x)
    b=j99(x)
    for i in range(99,1,-1):
        temp=b
        b=(2*i+1)*b/x-a
        a=temp 
    J1_test.append(b)

S=((np.sin(X)-X*np.cos(X))/X**2)/J1_test

#After we find the required scaling we use this as the funtion j50 and j49
l=10 #l is the order of the Bessel funtion that we want to find, here l=10

def j100(x):
    i=int((x-0.1)/0.2)
    return S[i]

def j99(x):
    return j100(x)

Jl=[]
for x in X: 
    a=j100(x)
    b=j99(x)
    for i in range(99,l,-1):
        temp=b
        b=(2*i+1)*b/x-a
        a=temp 
    Jl.append(b)


e10=[]
for x in X:
    i=int((x-0.1)/0.2)
    v=(abs(J10[i]-Jl[i]))/(abs(J10[i])+abs(Jl[i]))
    e10.append(v)

plt.subplot(2,2,1)
plt.plot(X,J10,'g',label="$j_{10}$ by up method")
plt.plot(X,X*0,'r')
#plt.plot(X,sp.spherical_jn(10,X)+0.01,'b', label="Actual $j_{10}$, shifted little bit for distinguishability")
plt.xlim(1,50)
plt.ylim(-0.5,0.5)
plt.legend()
plt.subplot(2,2,2)
plt.plot(X,Jl,'r',label="$j_1$ by down method, starting from $j_{50}, j_{49}$")
#plt.plot(X,sp.spherical_jn(1,X),'g+',label="Actual $j_1$, shifted little bit for distinguishability")
plt.plot(X,X*0,'k')
plt.xlim(0.1,50)
plt.ylim(-0.3,0.3)
plt.legend()
plt.subplot(2,2,3)
plt.loglog(X,J10,Jl)
plt.legend(["log-log plot for up $j_{10}$   ","log-log plot for down $j_{10}$"])
plt.subplot(2,2,4)
plt.loglog(X,e10,'b',label="log-log plot of the error")
#plt.ylim(-1,1)
plt.legend()
plt.show() 


