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

def j50(x):
    return 1

def j49(x):
    return 1

J1_test=[]
for x in X: 
    a=j50(x)
    b=j49(x)
    for i in range(49,1,-1):
        temp=b
        b=(2*i+1)*b/x-a
        a=temp 
    J1_test.append(b)

S=((np.sin(X)-X*np.cos(X))/X**2)/J1_test

def j50(x):
    i=int((x-0.1)/0.2)
    return S[i]

def j49(x):
    return j50(x)

J1=[]
for x in X: 
    a=j50(x)
    b=j49(x)
    for i in range(49,1,-1):
        temp=b
        b=(2*i+1)*b/x-a
        a=temp 
    J1.append(b)

plt.subplot(1,2,1)
plt.plot(X,J10,'g',label="$j_{10}$ by up method")
plt.plot(X,X*0,'r')
plt.plot(X,sp.spherical_jn(10,X)+0.01,'b', label="Actual $j_{10}$, shifted little bit for distinguishability")
plt.xlim(1,50)
plt.ylim(-0.5,0.5)
plt.legend()
plt.subplot(1,2,2)
plt.plot(X,J1,'r',label="$j_1$ by down method, starting from $j_{50}, j_{49}$")
plt.plot(X,sp.spherical_jn(1,X),'g+',label="Actual $j_1$, shifted little bit for distinguishability")
plt.plot(X,X*0,'k')
plt.xlim(0.1,50)
plt.ylim(-10,10)
plt.legend()
plt.show() 


