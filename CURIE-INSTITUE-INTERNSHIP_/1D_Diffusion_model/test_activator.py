import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def Activator(A1,A2,f23,f85,delta_phi,t0,t1,N):
    t =  np.linspace(t0,t1,N,True)
    w23 = 2*np.pi*f23
    w85 = 2*np.pi*f85
    s = (A1*A2/4)*(1+np.sin(w23*t-np.pi/2))*(1+np.sin(w85*t-np.pi/2+delta_phi*np.pi/180))
    return s

A1 = 2
A2 = 2
f23 = 1
f85 = 1
w23 = 2*np.pi*f23
w85 = 2*np.pi*f85
delta_phi = -180
t0 = 0
t1 = 1
N = 1000

t =  np.linspace(t0,t1,N,True)

n23 = A1*0.5*(1+np.sin(w23*t-np.pi/2))
n85 = A2*0.5*(1+np.sin(w85*t-np.pi/2+delta_phi*np.pi/180))

Activator = Activator(A1,A2,f23,f85,delta_phi,t0,t1,N)

plt.plot(t,n23,c='b', label = 'NHR23')
plt.plot(t,n85,c='k', label = 'NHR85')
plt.plot(t,Activator, c='r', label = 'Activator')
plt.legend()
plt.grid(True)
plt.show()
