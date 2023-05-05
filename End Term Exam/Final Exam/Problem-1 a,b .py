import numpy as np 
import matplotlib.pyplot as plt 
import math as m 

MeV=1.6*10**-13
fm=10**-15
c=3*10**8

V0=60*MeV
a=1.45*fm
h=197.92*MeV*fm/c
m=938*MeV/c**2

def V(r):
       if r<=a:
          return -V0 
       else: 
          return 0
       
def Eigen(V0):
    
    def fl(E):
        k=np.sqrt(2*m*(V0-abs(E)))/h
        return k/np.tan(k*a)

    def fr(E):
        k1=np.sqrt(2*m*(abs(E)))/h 
        return -k1


    # E=np.linspace(-V0,0,1000)
    # plt.plot(E,fl(E),'r')
    # plt.plot(E,fr(E),'g')
    # plt.show()

#Finding the eigenvalue: From the plot we can see that there is only ine eigenvalue for this particular case which we find using Newton-Rhapson method
    b=np.sqrt(2*m*V0)/h

    def f(k):
        return k/np.tan(k*a)+np.sqrt(b**2-k**2)

    def f1(k):
        return (1/np.tan(k*a)-(k*a)/np.sin(k*a)**2)-k/np.sqrt(b**2-k**2)

    k0=b*0.8
    for i in range(100):
        k0-=(f(k0))/f1(k0)

    E0=-(h)**2/(2*m)*(b**2-k0**2)
    # K=np.linspace(0,b,100)
    # plt.plot(K,f(K),'r')
    # plt.plot(K,K*0,'g')
    plt.show()
    return E0/MeV

print("The eigneenergy of V0=60MeV is"+str(Eigen(60*MeV))+"MeV")

V=np.arange(60,10,-0.25)
plt.plot(V,Eigen(V*MeV),'ro')
plt.show()
for v in V:
    if Eigen(v*MeV)<0:
        continue
    else:
        print("The min V0 for which an energy eigenvalue exists is"+str(v)+"MeV")
        break

# So, the minimum V0 for which there is a bound energy value is 39.0 MeV

#Solving this differential equation using Runge-Kutta-4 method using 
#the eigenvalue E for V0=60 MeV
