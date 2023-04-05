import numpy as np 
import math as m 
import matplotlib.pyplot as plt 
from scipy.special import eval_legendre

def Leg(x,n):
   A,B=1,x
   for i in range(1,n):
      C=B
      B=((2*i+1)*x*B-(i)*A)/(i+1)
      A=C
   return B 


def Leg_diff(x,n):
   return n*(x*Leg(x,n)-Leg(x,n-1))/(x**2-1)


def Root_Asymp(k,n):
   return (1-1/(8*n**2)+1/(8*n**3))*np.cos(m.pi*(4*k-1)/(4*n+2))

def Roots(N):
   S=[]
   for i in range(1,N+1):
       v=Root_Asymp(i,N)
       S.append(v)

   S_pol=[]
   for k in range(0,N):
      g=S[k]
      for i in range(0,10):
         g-=Leg(g,N)/Leg_diff(g,N)
      S_pol.append(g)
   return S_pol 

def W(N):
   W=[]
   for k in range(0,N):
      r=Roots(N)[k]
      v=2/((1-r**2)*Leg_diff(r,N)**2)
      W.append(v)
   return W 


def GL(f, a, b, n):  
    def F(y):
        x=(b-a)/2*y+(b+a)/2 
        return f(x)*(b-a)/2
    
    sum=0
    for i in range(0,n):
        sum+=W(n)[i]*F(Roots(n)[i])
    return sum


#W_10,R_10=W(10),Roots(10)
#$W_40,R_40=W(40),Roots(40)
#W_100,R_100=W(100),Roots(100)
W_1000=W(1000) #Roots(1000)

with open('W_1000.dat', 'w') as file1:
    for item in W_1000:
        file1.write("%s\n" % item)

