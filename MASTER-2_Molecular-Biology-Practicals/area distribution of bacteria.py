#tracé d'histogramme
from math import *
import hist
import statistics
import numpy as np
import matplotlib.pyplot as plt
import random
from random import gauss

def fit_gauss(x,y):
    N = len(x)
    e = 0

    a_inf = statistics.mean(y)/200
    a_sup = statistics.mean(y)/150
    a_step = (a_sup-a_inf)/10

    m_inf = 0.5*statistics.mean(x)
    m_sup = 1.5*statistics.mean(x)
    m_step = (m_sup-m_inf)/10

    s_inf = (max(x)-min(x))/(N)
    s_sup = max(x)-min(x)
    s_step = (s_sup-s_inf)/10
    
    Ea = []
    Em = []
    Em_copie = []
    Es = []
    Es_copie = []
    a_i = a_inf
    while(a_i<a_sup):
        m_i = m_inf
        while(m_i<m_sup):
            s_i = s_inf
            while(s_i<s_sup):
                e_i = 0
                for k in range(N):
                    x_k = x[k]
                    y_theor_k = y[k]
                    b_k = ((x_k-m_i)**2)/(2*s_i*s_i)
                    y_exp_k = a_i*exp(-b_k)
                    e_i = e_i + (y_theor_k-y_exp_k)**2
                Es.append(e_i)
                s_i = s_i + s_step
            Es_copie = list(Es)
            Em.append(Es_copie)
            Es.clear()
            m_i = m_i + m_step
        Em_copie = list(Em)
        Ea.append(Em_copie)
        Em.clear()
        a_i = a_i + a_step
    print(Ea)
    emin = Ea[0][0][0]
    A = a_inf
    m = m_inf
    s = s_inf
    for i in range(len(Ea)):
        for k in range(len(Ea[0])):
            for j in range(len(Ea[0][0])):
                if (Ea[i][j][k]<emin):
                    emin = Ea[i][j][k]
                    A = a_inf + i*a_step
                    m = m_inf + j*m_step
                    s = s_inf + k*s_step
    return A,m,s

def histog(X,L):
    [A,moy,sigma] = fit_gauss(X,L)
    # A = 
    #sigma = 
    var = sigma**2
    xgauss = []
    ygauss = []
    maxi = 0
    l = len(L)
    for i in range(l):
        if L[i]>maxi:
            maxi = L[i]
    for i in range(round(min(L)-10),round(maxi+10),1):
        b = ((i-moy)**2)/(2*sigma*sigma)
        y = A*exp(-b)
        xgauss.append(i)
        ygauss.append(y)
    plt.plot(xgauss,ygauss,"r")
    plt.hist(L,bins = [min(L)+pas-1 for pas in range(0,round(max(L)-min(L))+1,5000)], range=(round(min(L)-1),round(max(L)+1)),color='blue',edgecolor='yellow')
    plt.xlabel("Intensity/Area (µm-2)")
    plt.ylabel("Number of occurrences")
    moy = round(moy,1)
    sigma = round(sigma,1)
    titre="'Intensity/Area' distribution in bacteria, Mean_I/A= ",moy,", Sigma= ",sigma
    plt.title(titre)
    plt.show()
    return moy,sigma

f = np. genfromtxt(r'C:\Users\lilia\Desktop\all_areas_bacteria_txt.txt')
Area = []
Mean_Intensity = []
Mean_I_over_Area = []
for i in range(len(f)):
    A = f[i][1]
    I = f[i][2]
    Area.append(A)
    Mean_Intensity.append(I)
    Mean_I_over_Area.append(I/A)

h = histog(Area,Mean_I_over_Area)
