import parser
import numpy as np 
import matplotlib.pyplot as plt 

#Defining the Matrix 
N=10**4
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

for i in range(0,n-1):
    v=M[i][i]
    M[i]=[x/v for x in M[i]] 
    Y[i]=Y[i]/v 
    for j in range(i+1, n-1):
        a=M[j][i]
        for k in range(j-1,j+2):
            M[j][k]=M[j][k]-M[i][k]*a
        Y[j]=Y[j]-Y[i]*a


#Backward substitution 

for i in range(n-1,-1,-1):
    for j in range(i-1,-1,-1):
        for k in range(j-1,j+2):
            a=M[j][i]
            M[j][k]=M[j][k]-M[i][k]*a
        Y[j]=Y[j]-Y[i]*a

#print(np.array(M),Y)

def Uact(x):
    return x*(1-x)*np.exp(x)
    
Err=[]
for k in range(len(M)):
    c=abs((Y[k]-Uact(k*h))/Uact(k*h))
    Err.append(c)    

#Err_sort=sorted(Err) 
#print(Err_sort[-2])   

X=np.linspace(0, 1, N)
plt.plot(X, Y,'k+',label="Numerically found value")
plt.plot(X, Uact(X),'r',label="$u^{exact}$")
plt.legend()
plt.show()