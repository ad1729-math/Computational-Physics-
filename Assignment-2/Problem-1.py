import parser
import numpy as np 
import matplotlib.pyplot as plt 

#Defining the Matrix 
N=1000
h=1/(N-1)

A1=[1]
for j in range(1,N):
    A1.append(0)

L=[A1]

for i in range(1,N-1):
    S=[]
    for j in range(0,N):
        if abs(i-j)==1:
            S.append(1)
        else: 
            if j==i:
                S.append(-2) 
            else: 
                S.append(0)
    L.append(S)

AN=[]
for i in range(0,N-1):
    AN.append(0)
AN.append(1) 

L.append(AN)

#Defining the vector 

def f(x):
    return (3*x+x**2)*np.exp(x)

U=[]
for i in range(0,N-1):
    x=(i)*h
    U.append(-(h**2)*f(x))
U.append(0)

#Gauss Elimination
M=L
Y=U
n=len(M)

#Forward substitution

for i in range(0,n):
    v=M[i][i]
    M[i]=[x/v for x in M[i]] 
    Y[i]=Y[i]/v 
    for j in range(i+1, n):
        a=M[j][i]
        Temp=[]
        for k in range(0,n):
            t=M[j][k]-M[i][k]*a
            Temp.append(t)
        M[j]=Temp
        Y[j]=Y[j]-Y[i]*a


#Backward substitution 

for i in range(n-1,-1,-1):
    for j in range(i-1,-1,-1):
        a=M[j][i]
        Temp=[]
        for k in range(0,n):
            t=M[j][k]-M[i][k]*a
            Temp.append(t)
        M[j]=Temp
        Y[j]=Y[j]-Y[i]*a


#print(np.array(M),Y,S)

def Uact(x):
    return x*(1-x)*np.exp(x)
    
Err=[]
for k in range(len(M)):
    c=abs((Y[k]-Uact(k*h))/Uact(k*h))
    Err.append(c)    

X=np.linspace(0, 1, N)
plt.plot(X, Y,'b',label="Numerically found value")
plt.plot(X, Uact(X),'r',label="$u^{exact}$")
plt.legend()
plt.show()