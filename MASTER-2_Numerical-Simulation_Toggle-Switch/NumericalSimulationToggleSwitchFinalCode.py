import math
from math import *
#import numpy as np
import matplotlib.pyplot as plt

n = 100
dt = 0.1
tmax = n*dt

kil = 1
klt = 1
ktl = 1
kat = 1
kpL = 1
kpT = 1
kdL = 0.9
kdT = 0.8
K = 1
K1 = 1
K2 = 1
kR1 = 100
kiR1 = 1
kR2 = 100
kiR2 = 1
B = 0.7
n_hills = 2

CGFP0 = 0
CRFP0 = 0
CLacIt = 0.2
CLacIa = 0
CTetRt = 0
CTetRa = 0
CIPTG = 0
CATC = 1

#Création d'une fonction créneau qui permet de générer un créneau de période T, de rapport cyclique PWM, sur n points à intervalle de temps régulier dt
def creneau(y0,y1,n,PWM,T,dt):  #Np: nombre de périodes, n: nombre de points sur tout le signal
    y = []
    t = 0
    cpt = 0 #la variable cpt permet de dénombrer le nombre de points qui ont déjà été calculés et stockés dans y au cours du programme
    Th = PWM*T #temps passé à l'état haut = PWM * temps d'une période
    if(y0<y1):
        #ON EST DANS LA CONFIGURATION OU ON COMMENCE A L'ETAT BAS
        while(cpt!=n):
            while(t<T-Th):
                y.append(y0)
                cpt = cpt + 1
                t = t + dt
            while(t<T):
                y.append(y1)
                cpt = cpt + 1
                t = t + dt
            t = 0 #réinitialisation de la variable t qui permet simplement de générer les points sur une période
    else:
        #ON EST DANS LA CONFIGURATION OU ON COMMENCE A L'ETAT HAUT
        while(cpt!=n):
            while(t<Th):
                y.append(y0)
                cpt = cpt + 1
                t = t + dt
            while(t<T):
                y.append(y1)
                cpt = cpt + 1
                t = t + dt
            t = 0 #réinitialisation de la variable t qui permet simplement de générer les points sur une période
    return y

def constante(y0,n):
    y = []
    for i in range(n):
        y.append(y0)
    return y

def Prod_LacI(CTetRt,CTetRa,CATC,K1,kR1,kiR1,kpL):
    return (kpL*kiR1*CTetRt/(kR1*CATC+kiR1))

def Degrad_LacI(CLacI,kdL):
    return (CLacI*kdL)

def Prod_TetR(CLacIt,CLacIa,CIPTG,K2,kR2,kiR2,kpT):
    return (kpT*kiR2*CLacIt/(kR2*CIPTG+kiR2))

def Degrad_TetR(CTetR,kdT):
    return (CTetR*kdT)

def Prod_GFP(CLacIt,CLacIa,CIPTG,K2,kR2,kiR2,coef_prod_T):
    return (coef_prod_T*kiR2*CLacIt/(kR2*CIPTG+kiR2))

def Prod_RFP(CTetRt,CTetRa,CATC,K1,kR1,kiR1,coef_prod_L):
    return (coef_prod_L*kiR1*CTetRt/(kR1*CATC+kiR1))

def H(B,y,K,n_hills): ## cette fonction de Hills traduit l'inhibition d'une molécule sur une autre directement
## en inversant le rapport y/K dans la formule, on traduirait cette fois l'inhibition de l'inhibition d'une molécule par une autre (avec un intermédiaire entre les deux nécessairement)
    x = y/K
    denom = 1+ (x**n_hills)
    return (B/denom)

def IPTG_LacI(CLaci,B,CIPTG,K,n_hills):
    return CLaci*H(B,CIPTG*CLacI,K,n_hills)

def ATC_TetR(CTetR,B,CATC,K,n_hills):
    return (CTetR*H(B,CATC*CTetR,K,n_hills))

def LacI_tetR(CtetR,B,CLacI,K,n_hills):
    return (CtetR*H(B,CLacI,K,n_hills))

def TetR_lacI(ClacI,B,CTetR,K,n_hills):
    return (ClacI*H(B,CTetR,K,n_hills))

##def Process(Concentrations,Constants,dt):
##    kil = Constants[0]
##    klt = Constants[1]
##    ktl = Constants[2]
##    kat = Consta[0.182, 0.16563980198nts[3]
##    kpl = Constants[4]
##    kpt = Constants[5]
##    K = Constants[6]
##
##    CLacI = Concentrations[0]
##    CTetR = Concentrations[1]
##    CIPTG = Concentrations[2]
##    CATC = Concentrations[3]
##    ClacI = Concentrations[4]
##    CtetR = Concentrations[5]
##    
##    CLacIf = CLacI + dt*(Prod_LacI(CLacI,kpl) + Degrad_LacI(CLacI,kdl) + ATC_TetR(CTetR,CATC,K,n))
##    CTetRf = CTetR + dt*(Prod_TetR(CTetR,kpt) + Degrad_TetR(CTetR,kdt) + IPTG_LacI(CLaci,B,CIPTG,K,n))
##    ClacIf = ClacI + dt*(TetR_lacI(ClacI,B,CTetR,K,n))
##    CtetRf = Cte[0.182, 0.16563980198tR + dt*(LacI_tetR(CtetR,B,CLacI,K,n))
##    Concentrations = [CLacIf,CTetRf,CIPTG,CATC,ClacIf,CtetRf]
##    return Concentrations


def Process(CLacIa,CLacIt,CTetRa,CTetRt,CGFP,CRFP,CIPTG,CATC,kR1,kR2,kiR1,kiR2,K1,K2,B,n,n_hills,tmax,dt):
    CTetRaf = kiR2*CTetRt/(kR2*CATC+kiR2)
    CLacIaf = kiR1*CLacIt/(kR1*CIPTG+kiR1)
    CGFPf = CGFP + Prod_GFP(CLacIt,CLacIa,CIPTG,K2,kR2,kiR2,kpT)
    CRFPf = CRFP + Prod_RFP(CTetRt,CTetRa,CATC,K1,kR1,kiR1,kpL)
    CLacItf = CLacIt + dt*(Prod_LacI(CTetRt,CTetRa,CATC,K1,kR1,kiR1,kpL) - Degrad_LacI(CLacIt,kdL))
    CTetRtf = CTetRt + dt*(Prod_TetR(CLacIt,CLacIa,CIPTG,K2,kR2,kiR2,kpT) - Degrad_TetR(CTetRt,kdT))
    Concentrations = [CTetRaf,CLacIaf,CGFPf,CRFPf,CLacItf,CTetRtf]
    return Concentrations


def programme(CLacIa,CLacIt,CTetRa,CTetRt,CGFP0,CRFP0,CIPTG,CATC,All_CIPTG,All_CATC,kR1,kR2,kiR1,kiR2,K1,K2,B,n,n_hills,tmax,dt):
    All_C = []
    Concentrations = []
    t=0
    cpt = 1
    CGFP = CGFP0
    CRFP = CRFP0
    C_Initiales = [CTetRa,CLacIa,CGFP,CRFP,CLacIt,CTetRt]
    All_C.append(C_Initiales)
    while (t<tmax and cpt<n):
        print(cpt)
        CIPTG = All_CIPTG[cpt]
        CATC = All_CATC[cpt]
        Concentrations = Process(CLacIa,CLacIt,CTetRa,CTetRt,CGFP,CRFP,CIPTG,CATC,kR1,kR2,kiR1,kiR2,K1,K2,B,n,n_hills,tmax,dt)
        All_C.append(Concentrations)
        ##mise à jour des valeurs des concentrations
        CTetRa = Concentrations[0]
        CLacIa = Concentrations[1]
        CGFP = Concentrations[2]
        CRFP = Concentrations[3]
        CLacIt = Concentrations[4] #même commentaire que pour LacI que pour TetR
        CTetRt = Concentrations[5] #Concentration de protéine TetR totale dans le milieu (={TetR actifs + TetR greffés sur les sites du gène lacI})
        
        t=t+dt
        cpt+=1
    return (All_C)

All_CIPTG = constante(CIPTG,n)
All_CATC = constante(CATC,n)

print(All_CIPTG)
print(len(All_CIPTG))


Time = [k for k in range (n)]
All_C = programme(CLacIa,CLacIt,CTetRa,CTetRt,CGFP0,CRFP0,CIPTG,CATC,All_CIPTG,All_CATC,kR1,kR2,kiR1,kiR2,K1,K2,B,n,n_hills,tmax,dt)

CTetRa=[]
CLacIa=[]
CGFP=[]
CRFP=[]
CLacIt=[]
CTetRt=[]


for k in range (n):
    CTetRa.append(All_C[k][0])
    CLacIa.append(All_C[k][1])
    CGFP.append(All_C[k][2])
    CRFP.append(All_C[k][3])
    CLacIt.append(All_C[k][4])
    CTetRt.append(All_C[k][5])


print(CLacIt)
print(CTetRt)

plt.subplot(211)
plt.plot (Time,CLacIt, c = 'r', label = 'LacI')
plt.plot (Time,CTetRt, c = 'g', label = 'TetR')
plt.legend()
plt.title("Concentrations des protéines")
plt.grid(True)

plt.subplot(212)
plt.plot (Time,All_CIPTG, c = 'b', label = 'IPTG')
plt.plot (Time,All_CATC, c = 'y', label = 'ATC')
plt.legend()
plt.title("Concentrations des molécules de contrôle")
plt.grid(True)
plt.show()

#CODE COULEUR:
#IPTG en bleu
#ATC en jaune
#LACI en rouge (RFP)
#TETR en vert (GFP)

#NOTE IMPORTANTE:
#les protéines TetR et LacI qui sont greffées sur les gènes (resp.) lacI et tetR, peuvent être dé-greffées si on ajoute (resp.) ATC et IPTG.
#Cela signifie qu'elles ne disparaissent pas complètement lorsqu'elles sont greffées, dans le code, on doit donc faire en sorte de considérer les intéractions entre "ATC et TetR qui donne le complexe ATC-TetR" ainsi que pour Laci et IPTG.



##     iptgdelay1: 4.3722e+03
##            k_l: 221.4087
##            k_t: 113.6997
##           katc: 36.8292
##          kiptg: 0.3090
##           matc: 2
##          miptg: 2
##             nl: 2
##             nt: 2.1403
##            crl: 0.3337
##            crt: 0.7350
##         deltal: 0.0165
##         deltat: 0.0165
##    delta_mrnal: 0.1386
##    delta_mrnat: 0.1386
##            cil: 10.4000
##            cit: 8.3200
##             cl: 0.8000
##             ct: 0.3000
##      atcdelay1: 0.1000
##      atcdelay2: 0.1000
