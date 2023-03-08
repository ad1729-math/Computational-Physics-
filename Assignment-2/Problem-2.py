import numpy as np 

I0=1
R1,R2,R3,R4,R5=2,4,6,8,10

M=[[1,0,0,1,0],
   [0,1,1,0,0],
   [0,R2,-R3,0,-R5],
   [R1,0,0,-R4,R5],
   [1,-1,0,0,-1],
   ]

R=[I0,I0,0,0,0]
n=len(M)

U=M
LT=[]
for i in range(len(U)):
   S=[]

   for k in range(0,i):
      S.append(0)
   S.append(1)

   for j in range(i+1,len(U)):
      v=U[j][i]/U[i][i]

      P=[]
      for k in range(len(U)):
         a=U[j][k]-U[i][k]*v
         P.append(a)
      U[j]=P

      S.append(v)
   LT.append(S)


L=[]
for i in range(len(LT)):
   S=[]
   for j in range(len(LT)):
      S.append(LT[j][i])
   L.append(S)


Y=[R[0]]
for i in range(1,len(U)):
   s=R[i]
   for j in range(0,i):
      s-=L[i][j]*Y[j]
   Y.append(s)


S=[Y[n-1]/U[n-1][n-1]]

for k in range(n-2,-1,-1):
   s=Y[k]
   for j in range(n-1,k,-1):
      s-=S[n-1-j]*U[k][j]
   S.append(s/U[k][k]) 

S.reverse()
print("The LU decompoaistion of the matrix M is")
print("L=",np.array(L))
print("U=",np.array(U))
print("Solution of LY=R is",Y)
print("Solution of US=Y is",S)
print("Hence I_1="+str(S[0]),"I2="+str(S[1]),"I3="+str(S[2]),"I4="+str(S[3]),"I5="+str(S[4]))






