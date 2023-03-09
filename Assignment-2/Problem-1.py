import numpy as np 
import matplotlib.pyplot as plt 

#Defining the Matrix 
N=10**5
h=1/(N-1)

A1=[1]
for j in range(1,N):
    A1.append(0)

L=[[1,0,0]]

for i in range(1,N-1):
    S=[1,-2,1]
    L.append(S)

L.append([0,0,1])


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
M[1][0]-=M[0][0]
v=M[1][1]
M[1]=[x/v for x in M[1]]
U[1]=U[1]/v
for i in range(1,n-1):
    v=M[i][1]
    M[i]=[x/v for x in M[i]]
    U[i]=U[i]/v
    a=M[i+1][0]
    M[i+1][0]-=M[i][1]*a
    M[i+1][1]-=M[i][2]*a
    Y[i+1]-=Y[i]*a


#Backward substitution 

M[n-2][2]=0

for i in range(n-2,-1,-1):
    a=M[i-1][2]
    M[i-1][2]-=M[i][1]*a
    Y[i-1]-=Y[i]*a



#print(np.array(M),Y)

def Uact(x):
    return x*(1-x)*np.exp(x)
    
Err=[]
for k in range(len(M)):
    c=abs((Y[k]-Uact(k*h))/Uact(k*h))
    Err.append(c)    

Err_sort=sorted(Err) 
#print(Err_sort[-2])   

X=np.linspace(0, 1, N)
plt.plot(X, Y,'k+',label="Numerically found value")
plt.plot(X, Uact(X),'r',label="$u^{exact}$")
plt.legend()


#These values of maximum error have been found by inputting different values of N

N=[10,100,250,500,750,1000,2000,2500,3000,4000,5000,7500,10**4,5*10**4,10**5,5*10**5,10**6]
E=[0.0076482162546865845,6.659835336247023e-05,1.0560407239650785e-05,2.6322226044365292e-06,1.1687127686445837e-06,6.570740200539755e-07,
1.6414560411029998e-07,1.0503666113608006e-07,7.293369284391954e-08,4.1017049273798214e-08,2.6246275045706256e-08,
1.1659389737943674e-08,6.5483221261819044e-09,5.73739581940873e-10,6.349200159782446e-10,3.234784674211385e-07,
8.300185081108075e-07]

N1=[10,100,1000,10000,100000,1000000]
E1=[E[N.index(x)] for x in N1]

a=np.log(10)

#plt.plot(-np.log(N1)/a, np.log(E1)/a,'bo',label="Log plot of maximum errors")
#plt.xlabel("$\log(h)$--->")
#plt.ylabel("$\log(\epsilon)$--->")
#plt.legend()
plt.show()