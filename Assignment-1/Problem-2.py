#Derivative of the arctan function 

import math as m 
import numpy as np 
import matplotlib.pyplot as plt 

h=np.linspace(10**-2,1,1000) #This is the step length 
a=m.sqrt(2)

def f2c(h):
    return (np.arctan(a+h)-np.arctan(a))/h

def f3c(h):
     
    return (np.arctan(a+h)-np.arctan(a-h))/(2*h)
    


#Comments: We can see that the 3 point method approaches the actual value much faster than 2 point method when we decrease the 
#value of h. Also we can see from the plot that the deviation of 3 point method from the actual value 1/3 is quadratic (O(h^2))
#while that the deviation of 2 point method from the actual value is apporximately linear (O(h)). 



#Logarithmic scale error plot:

def E1(h):
    v=abs(3*(1/3-(np.arctan(a+h)-np.arctan(h))/h))
    return np.log10(v)

def E2(h):
    v=abs(3*(1/3-(np.arctan(a+h)-np.arctan(a-h))/(2*h)))
    return np.log10(v)

plt.subplot(1,2,1)
plt.plot(h,f2c(h),'r',label="$f'_{2c}$")
plt.plot(h,f3c(h),'g',label="$f'_{3c}$")
plt.plot(h,h/h*1/3,'k',label="$f'=1/3$")
plt.legend()
plt.subplot(1,2,2)
plt.plot(np.log10(h),E1(h),'r',label="Logscale deviation for $f'_{2c}$")
plt.plot(np.log10(h),E2(h),'g',label="Logscale deviation for $f'_{3c}$")
plt.legend()
plt.show() 




