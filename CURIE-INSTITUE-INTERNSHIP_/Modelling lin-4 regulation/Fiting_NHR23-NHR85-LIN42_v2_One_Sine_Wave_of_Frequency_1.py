# This code tries to find the best parameters for sin waves to fit the experimental curves of {NHR23;NHR85;LIN42} DNA-binding over time.
# The found parameters will serve in the other codes where these three species are taken into account.
import math
import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import signal
import numpy, scipy.optimize

def constant(y1,N):
    y = np.empty(0)
    for i in range(N):
        y = np.append(y,y1)
    return y

def sinus(y0,y1,f,phi,t0,t1,dt):
    N = 1+int(round((t1-t0)/dt))
    phi = phi*np.pi/180 # conversion préliminaire de la phase en radians

    A = (y1-y0)/2
    B = A+y0
    t = np.linspace(t0,t1,N,True)
    return A*np.sin(2*np.pi*f*t+phi)+B

def fit_sin(tt, yy):
    # Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"
    tt = numpy.array(tt)
    yy = numpy.array(yy)
    ff = numpy.fft.fftfreq(len(tt),(tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(numpy.fft.fft(yy))
    guess_freq = abs(ff[numpy.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*guess_freq, 0., guess_offset])
    return guess

def sinfunc(t,A,w,p,c):
    teta = w*t+p
    s = np.sin(teta)
    #s = (1 + s)/2
    s = A*s
    s = s+c
    # w = 2πf
    # sin = A*(1+sin(2πft+p))/2+c
    return s

def LeastSquaresMethod(Ytheor,Yexp):
    N = len(Ytheor)
    eps = 0
    for i in range(N):
        dy = Ytheor[i]-Yexp[i]
        eps+=dy*dy
    return eps

Ltime = np.linspace(0,0.9,10,True)
# 0 to 8.8cm = 0 to 1000(-) DNA binding relative intensity
NHR23_experimental = np.array([0    , 1.5   , 4.0   , 4.6   , 3.3   , 0     , 0     , 0     , 0     , 0])
NHR23_experimental*=1000/8.8
NHR23_exp_first_indice = 1
NHR23_exp_last_indice = 6
NHR23_experimental_for_fit = NHR23_experimental[NHR23_exp_first_indice:NHR23_exp_last_indice]

# 0 to 11.7cm = 0 to 4000(-) DNA binding relative intensity
NHR85_experimental = np.array([5.2 , 3.3    , 1.5   , 0.2   , 0.6   , 0.3   , 1.8   , 4.8   , 7.9   , 7.2])
NHR85_experimental*=4000/11.7
NHR85_exp_first_indice = 0
NHR85_exp_last_indice = 10
NHR85_experimental_for_fit = NHR85_experimental[NHR85_exp_first_indice:NHR85_exp_last_indice]

# 0 to 10.2cm = 0 to 2500(-) DNA binding relative intensity
LIN42_experimental = np.array([0    , 0.2   , 1.5   , 2.9   , 4.2   , 2.8   , 0.2   , 0.1   , 0.1   , 0.1])
LIN42_experimental*=2500/10.2
LIN42_exp_first_indice = 1
LIN42_exp_last_indice = 7
LIN42_experimental_for_fit = LIN42_experimental[LIN42_exp_first_indice:LIN42_exp_last_indice]

NAmplitude = 200
NPhase = 200
NFrequency = 1
NConstant = 5000
delta_tmin_tmax = 0

f = plt.figure(1)

# NHR23
Error = [[[0 for k in range(NFrequency)] for u in range(NPhase)] for i in range(NAmplitude)]
Yexp = NHR23_experimental_for_fit
Amax = 2*max(Yexp)
Amplitudes = np.linspace(0,Amax,NAmplitude,True)
Phases = np.linspace(0,2*np.pi,NPhase,True)
Frequencies = np.array([1])
Cmin = -max(Yexp)
Cmax = max(Yexp)
Constants = np.linspace(Cmin,Cmax,NConstant,True)
t = Ltime[NHR23_exp_first_indice:NHR23_exp_last_indice]
for i in range(NAmplitude):
    A = Amplitudes[i]
    for j in range(NPhase):
        phase = Phases[j]
        for k in range(NFrequency):
            frequency = Frequencies[k]
            Ytheor = sinfunc(t,A,2*np.pi*frequency,phase,0)
            error = LeastSquaresMethod(Ytheor,Yexp)
            Error[i][j][k] = error
eps = Error[0][0][0]
for i in range(NAmplitude):
    for j in range(NPhase):
        for k in range(NFrequency):
            error = Error[i][j][k]
            if error<eps:
                eps = error
                A = Amplitudes[i]
                phase = Phases[j]
                frequency = Frequencies[k]

A23 = A
Phase23 = phase
Frequency23 = frequency
C23 = 0
Error23 = error
Error_Constant = [0 for i in range(NConstant)]

for i in range(NConstant):
    c = Constants[i]
    Ytheor = sinfunc(t,A23,2*np.pi*Frequency23,Phase23,c)
    error_constant = LeastSquaresMethod(Ytheor,Yexp)
    Error_Constant[i] = error_constant
eps = Error_Constant[0]
for i in range(NConstant):
    error = Error_Constant[i]
    if error<eps:
        eps = error
        C23 = Constants[i]

print("A23 =",A23)
print("Phase23 =",Phase23)
print("Freq23 =",Frequency23)

t_fit = np.linspace(0,1,101,True)
NHR23_fit = sinfunc(t_fit,A23,2*np.pi*Frequency23,Phase23,C23) #pi is in radians and phase too
NHR23_fit_scaled = NHR23_fit + A23 - C23
NHR23_fit_scaled/=max(NHR23_fit_scaled)
NHR23_fit_scaled*=2*A23
NHR23_fit_scaled+=1000-A23
print(" ")

plt.subplot(311)
title = "LIN23 fit to DNA-binding results"
label_str = "LIN23_fit = "+str(round(A23))+"sin(2π"+str(round(Frequency23,2))+"t+"+str(round(Phase23,2))+")+"+str(round(C23))
plt.title(title)
plt.plot(Ltime,NHR23_experimental   , color = 'r', marker = '+' , label = 'NHR23_exp')
plt.plot(t_fit,NHR23_fit            , color = 'b', label = label_str)
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

# NHR85
Error = [[[0 for k in range(NFrequency)] for u in range(NPhase)] for i in range(NAmplitude)]
Yexp = NHR85_experimental_for_fit
Amax = 2*max(Yexp)
Amplitudes = np.linspace(0,Amax,NAmplitude,True)
Phases = np.linspace(0,2*np.pi,NPhase,True)
Frequencies = np.array([1])
Cmax = max(Yexp)
Constants = np.linspace(0,Cmax,NConstant,True)
t = Ltime[NHR85_exp_first_indice:NHR85_exp_last_indice]
for i in range(NAmplitude):
    A = Amplitudes[i]
    for j in range(NPhase):
        phase = Phases[j]
        for k in range(NFrequency):
            frequency = Frequencies[k]
            Ytheor = sinfunc(t,A,2*np.pi*frequency,phase,0)
            error = LeastSquaresMethod(Ytheor,Yexp)
            Error[i][j][k] = error
eps = Error[0][0][0]
for i in range(NAmplitude):
    for j in range(NPhase):
        for k in range(NFrequency):
            error = Error[i][j][k]
            if error<eps:
                eps = error
                A = Amplitudes[i]
                phase = Phases[j]
                frequency = Frequencies[k]

A85 = A
Phase85 = phase
Frequency85 = frequency
C85 = 0
Error85 = error
Error_Constant = [0 for i in range(NConstant)]

for i in range(NConstant):
    c = Constants[i]
    Ytheor = sinfunc(t,A85,2*np.pi*Frequency85,Phase85,c)
    error_constant = LeastSquaresMethod(Ytheor,Yexp)
    Error_Constant[i] = error_constant
eps = Error_Constant[0]
for i in range(NConstant):
    error = Error_Constant[i]
    if error<eps:
        eps = error
        C85 = Constants[i]

print("A85 =",A85)
print("Phase85 =",Phase85)
print("Freq85 =",Frequency85)

t_fit = np.linspace(0,1,101,True)
NHR85_fit = sinfunc(t_fit,A85,2*np.pi*Frequency85,Phase85,C85) #pi is in radians and phase too
NHR85_fit_scaled = NHR85_fit + 1000 # A85 - C85
#NHR85_fit_scaled/=max(NHR85_fit_scaled)
#NHR85_fit_scaled*=2*A85
#NHR85_fit_scaled+=1000-A85
print(" ")

plt.subplot(312)
title = "LIN85 fit to DNA-binding results"
label_str = "LIN85_fit = "+str(round(A85))+"sin(2π"+str(round(Frequency85,2))+"t+"+str(round(Phase85,2))+")+"+str(round(C85))

plt.title(title)
plt.plot(Ltime,NHR85_experimental   , color = 'k', marker = '+' , label = 'NHR85_exp')
plt.plot(t_fit,NHR85_fit            , color = 'b', label = label_str)
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

# LIN42
Error = [[[0 for k in range(NFrequency)] for u in range(NPhase)] for i in range(NAmplitude)]
Yexp = LIN42_experimental_for_fit
Amax = 2*max(Yexp)
Amplitudes = np.linspace(0,Amax,NAmplitude,True)
Phases = np.linspace(0,2*np.pi,NPhase,True)
Frequencies = np.array([1])
Cmin = -max(Yexp)
Cmax = max(Yexp)
Constants = np.linspace(Cmin,Cmax,NConstant,True)
t = Ltime[LIN42_exp_first_indice:LIN42_exp_last_indice]
for i in range(NAmplitude):
    A = Amplitudes[i]
    for j in range(NPhase):
        phase = Phases[j]
        for k in range(NFrequency):
            frequency = Frequencies[k]
            Ytheor = sinfunc(t,A,2*np.pi*frequency,phase,0)
            error = LeastSquaresMethod(Ytheor,Yexp)
            Error[i][j][k] = error
eps = Error[0][0][0]
for i in range(NAmplitude):
    for j in range(NPhase):
        for k in range(NFrequency):
            error = Error[i][j][k]
            if error<eps:
                eps = error
                A = Amplitudes[i]
                phase = Phases[j]
                frequency = Frequencies[k]

A42 = A
Phase42 = phase
Frequency42 = frequency
C42 = 0
Error42 = error
Error_Constant = [0 for i in range(NConstant)]

for i in range(NConstant):
    c = Constants[i]
    Ytheor = sinfunc(t,A42,2*np.pi*Frequency42,Phase42,c)
    error_constant = LeastSquaresMethod(Ytheor,Yexp)
    Error_Constant[i] = error_constant
eps = Error_Constant[0]
for i in range(NConstant):
    error = Error_Constant[i]
    if error<eps:
        eps = error
        C42 = Constants[i]

print("A42 =",A42)
print("Phase42 =",Phase42)
print("Freq42 =",Frequency42)

t_fit = np.linspace(0,1,101,True)
LIN42_fit = sinfunc(t_fit,A42,2*np.pi*Frequency42,Phase42,C42) #pi is in radians and phase too
LIN42_fit_scaled = LIN42_fit + 1000 #+ A42 - C42

plt.subplot(313)
title = "LIN42 fit to DNA-binding results"
plt.title(title)
label_str = "LIN42_fit = "+str(round(A42))+"sin(2π"+str(round(Frequency42,2))+"t+"+str(round(Phase42,2))+")+"+str(round(C42))

plt.plot(Ltime,LIN42_experimental   , color = 'g', marker = '+' , label = 'LIN42_exp')
plt.plot(t_fit,LIN42_fit            , color = 'b', label = label_str)
plt.xlabel('')
plt.ylabel('')
plt.legend()
plt.grid(True)

f.show()

g = plt.figure(2)

title = "NHR23 and NHR85 conc. during one LS"
plt.title(title)
#label_str = "LIN42_fit = "+str(round(A42))+"sin(2π"+str(round(Frequency42,2))+"t+"+str(round(Phase42,2))+")+"+str(round(C42))

plt.plot(t_fit,NHR23_fit_scaled , color = 'r' , label = "NHR23")
plt.plot(t_fit,NHR85_fit_scaled , color = 'k' , label = "NHR85")
#plt.plot(t_fit,LIN42_fit_scaled , color = 'g' , label = "LIN42")

NHR23_85_fit_scaled = np.multiply(NHR23_fit_scaled,NHR85_fit_scaled)
NHR23_85_fit_scaled = np.multiply(NHR23_85_fit_scaled,1/max(NHR23_85_fit_scaled))
index = list(NHR23_85_fit_scaled).index(max(NHR23_85_fit_scaled))
NHR23_85_fit_scaled = np.multiply(NHR23_85_fit_scaled,NHR23_fit_scaled[index])

#plt.plot(t_fit,NHR23_85_fit_scaled, color = 'b' , label = "NHR23-NHR85")
plt.xlabel("Scaled Larval Stage")
plt.ylabel("Fitted+scaled fluorescence intensity")
plt.ylim(0,3500)
plt.legend()
plt.grid(True)

g.show()
