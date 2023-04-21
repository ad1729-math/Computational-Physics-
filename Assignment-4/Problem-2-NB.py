import numpy as np 
import matplotlib.pyplot as plt 

def f(x):
    return np.exp(-x**2)

#Here we restirct the f to monotonically increasing or decreasing functions
M=1
def Monte(f,a,b,N):
    L_in=[[],[]]
    L_out=[[],[]]
    count=0
    for k in range(N):
        x,y=np.random.rand(),np.random.rand()
        X,Y=x*(b-a),y*M
        if Y<=f(X):
            count+=1
            L_in[0].append(X)
            L_in[1].append(Y)
        else:
            count+=0
            L_out[0].append(X)
            L_out[1].append(Y)
    plt.show()
    return count*M*(b-a)/N, L_in, L_out

#The mean of the result, the distribution of which is Gaussian for la
#large m 

a,b,n=0,1,10**4
m=1000
S=[]
for k in range(0,m):
    S.append(Monte(f,a,b,n)[0])
print(np.sum(S)/m)

