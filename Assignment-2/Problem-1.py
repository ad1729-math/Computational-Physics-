import numpy as np 
import matplotlib.pyplot as plt 

#Defining the Matrix 
N=5000
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
    a=M[i+1][i]
    for k in range(i,i+2):
        M[i+1][k]=M[i+1][k]-M[i][k]*a
    Y[i+1]=Y[i+1]-Y[i]*a


#Backward substitution 

for i in range(n-1,0,-1):
    a=M[i-1][i]
    for k in range(i-1,i+1):
        M[i-1][k]=M[i-1][k]-M[i][k]*a
    Y[i-1]=Y[i-1]-Y[i]*a

#print(np.array(M),Y)

def Uact(x):
    return x*(1-x)*np.exp(x)
    
Err=[]
for k in range(len(M)):
    c=abs((Y[k]-Uact(k*h))/Uact(k*h))
    Err.append(c)    

Err_sort=sorted(Err) 
print(Err_sort[-2])   

X=np.linspace(0, 1, N)
plt.plot(X, Y,'k+',label="Numerically found value")
plt.plot(X, Uact(X),'r',label="$u^{exact}$")
plt.legend()

#
N=[100,250,500,750,1000,2000,2500,3000,4000,5000,7500,10000]
E=[6.659835336247023e-05,1.0560407239650785e-05,2.6322226044365292e-06,1.1687127686445837e-06,6.570740200539755e-07,
1.6414560411029998e-07,1.0503666113608006e-07,7.293369284391954e-08,4.1017049273798214e-08,2.6246275045706256e-08,
1.1659389737943674e-08,6.5483221261819044e-09]

#plt.plot(-np.log(N), np.log(E),'bo',label="Log plot of maximum errors")
#plt.xlabel("$\log(h)$--->")
#plt.ylabel("$\log(\epsilon)$--->")
#plt.legend()
plt.show()