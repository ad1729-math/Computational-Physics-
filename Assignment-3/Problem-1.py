import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
#from Gauss_Legendre import GL

W1000=np.loadtxt("W_1000.dat",unpack=True)
R1000=np.loadtxt("R_1000.dat",unpack=True)

BM=[0.938,1.440,1.535,1.650,1.710] #These are Baryon masses 
MM=[0.782,1.170,1.282,1.420,1.512] #These are Meson masses 

m0=0.5

def Rho(m,A,TH):
    return A/(m**2+m0**2)**(5/4)*np.exp(m/TH) 

gB,gM=4,3 

def Theta(x):
    if x>=0:
        v=1
    else:
        v=0
    return v 

def Ne_Bs(m):
    S=0
    for M in BM:
        S+=gB*Theta(m-M)
    return S 

def Ne_B(m):
    if isinstance(m,list):
        S=[]
        for x in m:
            S.append(Ne_Bs(x))
        return S
    else:
        return Ne_Bs(m)
    
def Ne_Ms(m):
    S=0
    for M in MM:
        S+=gM*Theta(m-M)
    return S 

def Ne_M(m):
    if isinstance(m,list):
        S=[]
        for x in m:
            S.append(Ne_Ms(x))
        return S
    else:
        return Ne_Ms(m)

def GL(f,a,b): #Here we use 1000 datapoints
    def F(y):
        x=0.5*((b-a)*y+(b+a))
        return f(x)*(b-a)/2
    s=0
    for k in range(1000):
        s+=W1000[k]*F(R1000[k])
    return s

def Nth(m, A, TH):
    if isinstance(m, list):
        s = []
        for i in m:
            def f(x):
                return Rho(x, A, TH)
            s.append(GL(f,0,i))
        return s
    else:
        def f(x):
            return Rho(x,A,TH)
        s=GL(f,0,m)
        return s
    

def Chi_B(A,TH):
    s=0
    for m in BM:
        s+=(np.log(Ne_B(m))-np.log(Nth(m,A,TH)))**2/(np.log(10))**2
    return s 

def Chi_M(A,TH):
    s=0
    for m in MM:
        s+=(np.log(Ne_M(m))-np.log(Nth(m,A,TH)))**2/(np.log(10))**2
    return s 


####
# T=100
# a1b,a2b=0.05,0.3
# b1b,b2b=0.05,0.3

# a_b=np.linspace(a1b,a2b,T)
# th_b=np.linspace(b1b,b2b,T)

# ChiB=[]
# for A in a_b:
#     ChiB_A=[]
#     for TH in th_b:
#         ChiB_A.append(Chi_B(A,TH))
#     ChiB.append(ChiB_A)

# LM=[]
# IM=[]
# for i in range(T):
#     s=np.sort(ChiB[i])[0]
#     i=ChiB[i].index(s)
#     IM.append(i)
#     LM.append(s)

# L_min=np.sort(LM)[0]
# Im=LM.index(L_min)
# Jm=IM[Im]

# Am_B=a1b+(a2b-a1b)/T*Im
# Th_B=b1b+(b2b-b1b)/T*Jm


####

T=100
a1m,a2m=0.05,0.3
b1m,b2m=0.05,0.3

a_m=np.linspace(a1m,a2m,T)
th_m=np.linspace(b1m,b2m,T)

ChiM=[]
for A in a_m:
    ChiM_A=[]
    for TH in th_m:
        ChiM_A.append(Chi_M(A,TH))
    ChiM.append(ChiM_A)

LM=[]
IM=[]
for i in range(T):
    s=np.sort(ChiM[i])[0]
    i=ChiM[i].index(s)
    IM.append(i)
    LM.append(s)

L_min=np.sort(LM)[0]
Im=LM.index(L_min)
Jm=IM[Im]

Am_M=a1m+(a2m-a1m)/T*Im
Th_M=b1m+(b2m-b1m)/T*Jm

# print("The fitted values for Baryon and Mesons are respectively",[Am_B,Th_B])

x=np.linspace(0.7,2,100)
# #plt.subplot(2,1,1)
# plt.plot(x,Nth(x,0.2,0.25),'g',label="Theoretical,Baryons")
# plt.plot(BM,Ne_B(BM),'ko',label="Data points")
# plt.legend()

# plt.subplot(2,1,2)
plt.plot(x ,Nth(x,Am_M,Th_M),'g',label="Theoretical,Meson")
plt.plot(MM,Ne_M(MM),'ko',label="Data points")
plt.legend()
plt.show()
