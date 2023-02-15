import numpy as np 
import math as m



def f(x):
    return m.cos(x)-0.5 

#Finding root using bisection method 

def rootBS(n): #Here, n is the index of the period 2npi to 2(n+1)pi
    x0,x1=n*m.pi*2, n*m.pi*2+m.pi
    y0,y1=n*m.pi*2+m.pi, (n+1)*m.pi*2
    N=60
    for i in range(0,N):
       v1=(x0+x1)/2
       v2=(y0+y1)/2
       if f(v1)>0:
           x0=v1
       else:
           x1=v1

       if f(v2)>0:
            y1=v2
       else:
            y0=v2
    return [v1/(m.pi/3),v2/(m.pi/3)]

BS=[] #Now we list roots in between the range -10*(2pi) to +10*(2pi) as multiplees of pi/3

for n in range(-10,10):
    BS.append(rootBS(n)[0])
    BS.append(rootBS(n)[1])

#print(BS) #we can see that the roots/multiples of pi/3 are accurate upto 13 places after decimal


#Finding roots of the function using Newton-Rhapson method  


N1=5
def f1(x):
    return -m.sin(x)

def rootNR(n):
    x0=n*2*m.pi+0.5
    x1=(n+1)*2*m.pi-0.5 #Symmetric choice
    for i in range(0,N1):
       x0=x0-f(x0)/f1(x0)
       x1=x1-f(x1)/f1(x1)
    return [x0/(m.pi/3),x1/(m.pi/3)]


#At last we use Secant's method to find the roots
 
def rootSec(n):
    N2=8
    #For the two halfs we choose two pairs of close values of initial seeds (x0,x1) and (y0,y1)
    x0,x1=n*m.pi*2+0.5 ,n*m.pi*2+0.6
    y0,y1=(n+1)*m.pi*2-0.5,(n+1)*m.pi*2-0.6
    for i in range(0,N2):
        a=x1
        x1=(f(x1)*x0-f(x0)*x1)/(f(x1)-f(x0))
        x0=a
        #And for the other half
        b=y1
        y1=(f(y1)*y0-f(y0)*y1)/(f(y1)-f(y0))
        y0=b
    return [x1/(m.pi/3),y1/(m.pi/3)]

n=int(input("Enter the number n to get the roots between 2npi and 2(n+1)pi"))

print("Here roots are printed as multiple of pi/3")
print(rootBS(n))
print(rootNR(n))
print(rootSec(n))

#We see by varying N,N1,N2 that Newton-Rhapson method is th fastest followed by Secant mehtod and Bisection method. The Secant method 
#fails at a point becuase x_{n+1} and x_{n} becomes close to the extent that the difference f(x_{n+1})-f(x_n) fall out of machine 
#precision, to the computer reads it as 0 and this gives the divison by zero error.


