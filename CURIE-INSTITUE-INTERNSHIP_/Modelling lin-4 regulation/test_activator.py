import math
import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import signal

def Activator(A1,A2,f23,f85,delta_phi,t0,t1,N):
    t =  np.linspace(t0,t1,N,True)
    w23 = 2*np.pi*f23
    T23 = 1/f23
    t0_23 = 0
    t1_23 = t0_23 + 2*T23
    print(t1_23)
    w85 = 2*np.pi*f85
    T85 = 1/f85
    t0_85 = T85*((3/4)-phi85-delta_phi)
    t1_85 = t0_85 + T85
    s = (A1*A2/4)*(1+np.sin(w23*t-np.pi/2))*(1+np.sin(w85*t-np.pi/2+delta_phi*np.pi/180))
    return s

def compute_truncated_sin(s,f,phi,nT):
    T = 1/f
    dt = T/100
    w = 2*np.pi*f
    t0 = (asin(s)-phi)/w
    if t0<0:
        while(t0<0):
            t0+=T
    t1 = t0+nT*T
    N = 1+int(round((t1-t0)/dt))
    L_t = np.linspace(t0,t1,N,True)
    sin = np.sin(w*L_t+phi)
    return [L_t,sin]

f = 2
phi = -np.pi/2
dt = 0.1
nT = 2
[L_t1,sin1] = compute_truncated_sin(-1,f,phi,nT)
[L_t2,sin2] = compute_truncated_sin(-1,f,phi-np.pi/2,nT-1)
##print(L_t)
##print(sin)
plt.plot(L_t1,sin1,c='b', label = 'sin1')
plt.plot(L_t2,sin2,c='b', label = 'sin2')
plt.legend()
plt.grid(True)
plt.show()
##
##A1 = 2
##A2 = 2
##f23 = 1
##f85 = 1
##w23 = 2*np.pi*f23
##w85 = 2*np.pi*f85
##phi85 = -np.pi/2
##delta_phi = -90
##T23 = 1/f23
##t0_23 = 0
##t1_23 = t0_23 + T23
##t0 = 0.75
##t1 = 2
##N = 1000
##T85 = 1/f85
##t0_85 = T85*((3/4)-phi85-delta_phi)
##t1_85 = t0_85 + 85
##
##t =  np.linspace(t0,t1,N,True)
##
##n23 = A1*0.5*(1+np.sin(w23*t))
##n85 = A2*0.5*(1+np.sin(w85*t+delta_phi*np.pi/180))
##
##Activator = Activator(A1,A2,f23,f85,delta_phi,t0,t1,N)
##
##plt.title("NHR23 - NHR85 combined activator")
##plt.plot(t,n23,c='b', label = 'NHR23')
##plt.plot(t,n85,c='k', label = 'NHR85')
##plt.plot(t,Activator, c='r', label = 'Activator')
##plt.legend()
##plt.grid(True)
##plt.show()
