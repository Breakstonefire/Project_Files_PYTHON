from math import *
import matplotlib.pyplot as plt
import numpy as np
import statistics

def affichage(L,consigne):
    l=len(L)
    c=342
    #Lv=[]
    Ldt=[]
    Ltemps=[]
    Ldtof=[]
    t0=L[0][0]
    #Lpaquet=[]
    #Lprime=[]
    
    for i in range(l-1):
        dt=abs(L[i+1][0]-L[i][0])
        dtof=(abs(L[i+1][1]-L[i][1]))
        #v=abs(c*dtof*100/(2*dt)) #v=vitesse voiture en cm/s
        #Lv.append(v)
        if dt<consigne*1000+5000:
            Ldt.append(dt)#ajout du délai calculé en micros dans Ldt
            Ltemps.append(L[i][0]-t0)
            if dtof>(2*0.5*consigne*1000/343)*110/100:
                Ldtof.append(2*0.5*consigne*1000/343)
            else:
                Ldtof.append(dtof)
        #Lpaquet.append(2*v*dt/c)
    ldt=len(Ldt)
    moy=int(statistics.mean(Ldt))
    sigma=int(statistics.pstdev(Ldt))
    Lmoy=[]
    Lconsig=[]
    for i in range(ldt):
        Lmoy.append(moy)
        Lconsig.append(consigne*1000)
    #le paquet 2vT/c vaut dans cette expérience (mur en face) la différence "tof2-tof1" on trouve ainsi u(paquet) grâce à u(tof2-tof1)        
    plt.figure(1)
    plt.scatter(Ltemps,Ldt)
    plt.plot(Ltemps,Lmoy,label="Espérance")
    plt.plot(Ltemps,Lconsig,label="Valeur consigne")
    plt.legend()
    titre="Délai en fonction du temps pour ",consigne*1000,"μs, T(μs)= ",moy," +/- ",sigma
    plt.title(titre)
    plt.xlabel('Temps')
    plt.ylabel('Délai')
    plt.grid()
    plt.show()

    plt.figure(2)
    plt.scatter(Ltemps,Ldtof)
    plt.legend()
    titre="Tau en fonction du temps pour ",consigne*1000,"μs"
    plt.title(titre)
    plt.xlabel('Temps')
    plt.ylabel('Tau')
    plt.grid()
    plt.show()

    Ltempscorrig=[]
    Ldtofcorrig=[]
    for i in range(ldt):
        if Ltemps[i]>=2500000:
            Ltempscorrig.append(Ltemps[i])
            Ldtofcorrig.append(Ldtof[i])
    moydtofcorrig=int(statistics.mean(Ldtofcorrig))
    sigmadtofcorrig=int(statistics.pstdev(Ldtofcorrig))
    plt.figure(3)
    plt.scatter(Ltempscorrig,Ldtofcorrig)
    plt.legend()
    titre="Tau stationnaire en fonction du temps pour ",consigne*1000,"μs, Tau(μs)= ",moydtofcorrig," +/- ",sigmadtofcorrig
    plt.title(titre)
    plt.xlabel('Temps')
    plt.ylabel('Tau')
    plt.grid()
    plt.show()
    return None
