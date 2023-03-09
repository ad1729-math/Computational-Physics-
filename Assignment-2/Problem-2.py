import numpy as np 

I0=1
R1,R2,R3,R4,R5=2,4,6,8,10

Z=[[1,0,0,1,0],
   [0,1,1,0,0],
   [1,-1,0,0,-1],
   [R1,0,0,-R4,R5],
   [0,R2,-R3,0,-R5],
   ]

R=[I0,I0,0,0,0]
n=len(Z)

L=[]
for i in range(len(Z)):
   S=[]
   for j in range(len(Z)):
      if i==j:
         S.append(1)
      else:
         S.append(0)
   L.append(S)


for i in range(len(Z)):

   if Z[i][i]==0:
      for k in range(i+1,len(Z)):
         if Z[k][i]!=0:
            Z[i], Z[k]=Z[k], Z[i]
            R[i],R[k]=R[k],R[i]
            for j in range(0,i):
               L[i][j], L[k][j]= L[k][j], L[i][j]
            break

   for j in range(i+1,len(Z)):
      v=Z[j][i]/Z[i][i]
      for k in range(i,len(Z)):
         Z[j][k]-=Z[i][k]*v 
      L[j][i]=v
         


Y=[R[0]]
for i in range(1,len(Z)):
   s=R[i]
   for j in range(0,i):
      s-=L[i][j]*Y[j]
   Y.append(s)


S=[Y[n-1]/Z[n-1][n-1]]

for k in range(n-2,-1,-1):
   s=Y[k]
   for j in range(n-1,k,-1):
      s-=S[n-1-j]*Z[k][j]
   S.append(s/Z[k][k]) 

S.reverse()
print("The LU decompoaistion of the matrix M is")
print("L=",np.array(L))
print("U=",np.array(Z))
print("Solution of LY=R is",Y)
print("Solution of US=Y is",S)
print("Hence I_1="+str(S[0]),"I2="+str(S[1]),"I3="+str(S[2]),"I4="+str(S[3]),"I5="+str(S[4]))






