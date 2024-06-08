# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from math import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Definition of the functions # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def constant(y1,N):
    y = np.empty(0)
    for i in range(N):
        y = np.append(y,y1)
    return y

def normalize(Y):
    Y = np.array(Y)
    maxi= max(Y)
    if(maxi==0):
        print(" ")
        print("WARNING - the function cannot be normalized.")
        print(" ")
    else:
        Y = Y/maxi
    return Y
    
def sinus(y0,y1,f,phi,t0,t1,dt):
    N = 1+int(round((t1-t0)/dt))
    phi = phi*np.pi/180 # conversion of the phase in radians

    A = (y1-y0)/2
    B = A+y0
    t = np.linspace(t0,t1,N,True)
    return A*np.sin(2*np.pi*f*t+phi)+B

def compute_truncated_sin(s,f,phi,nT,dt):
    T = 1/f
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

def Occupancy(Kr,Complex_DNA_TF_concentration,DNA_site_concentration):
    if(Complex_DNA_TF_concentration==0):
        print("Be careful, the concentration of Protein-DNA complex is null.")
        return 0
    else:
        ratio = DNA_site_concentration/(Kr*Complex_DNA_TF_concentration)
        s = pow(1+ratio,-1)
        return s

def Initiation_rate(Imax,kefficiency,Occupancy,t,Tresponse_delay):
    Ir = Imax*(1-np.exp(-kefficiency*Occupancy*(t-Tresponse_delay)/Imax))
    return Ir

def Time_delay_in_response(Translocation_rate,gene_size):
    Tresponse_delay = Translocation_rate*gene_size
    return Tresponse_delay

def Translocation_rate(Trate_T0,Q10_law_Temperature_increase,Temperature,Temperature_T0):
    # This function is a consequence from the Q10 law which says that we ASSUME the translocation rate as a power function of the Q10 value.
    # Q10 value corresponds to the increase of the translocation rate at Temperature = T based on its value at Temperature = T0
    # For mammalian cells and organisms, it has been shown that the Q10 value is framed between 2 and 3.
    # Thus, this law is based on the fact that the Translocation rate is "Q10" time higher at Temperature=T0+10 than at T=T0.
    power = (Temperature-Temperature_T0)/10
    Trate = Trate_T0*pow(Q10_law_Temperature_increase,power)
    return Trate

def Irate_NHR23_NHR85_lin4 (Imax,Dn,ka_specific_common_value,Kn,Kr,KCn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_site1,efficiency23_site1,efficiency_85_site1,efficiency23_85_site1,gene_size_site2,efficiency23_site2,efficiency23_85_site2,efficiency85_site2):

    if(efficiency23_site1 + efficiency_85_site1 + efficiency23_85_site1>1.1): # This condition ensures the good conversion of the efficiencies from percentages ( going in [0:100] ) to values going in [0:1].
                                                                                # This is needed to be sure that the maximal Irate doesn't exceed Imax, the maximal Initiation rate given in parameters (instead of 100*Imax if we don't do the conversion here).
        efficiency23_site1/=100
        efficiency_85_site1/=100
        efficiency23_85_site1/=100

    if(efficiency23_site2 + efficiency85_site2 + efficiency23_85_site2>1.1): # This condition ensures the good conversion of the efficiencies from percentages ( going in [0:100] ) to values going in [0:1].
                                                                                # This is needed to be sure that the maximal Irate doesn't exceed Imax, the maximal Initiation rate given in parameters (instead of 100*Imax if we don't do the conversion here).
        efficiency23_site2/=100
        efficiency85_site2/=100
        efficiency23_85_site2/=100
    
    Trate = Translocation_rate(Trate_T0,Q10_law_Temperature_increase,Temperature,Temperature_T0)
    
    ### Computing the inherent time delay in response ###
    Tdelay_in_response23_site1 = Time_delay_in_response(Trate,gene_size_site1) # Inherent delay in the response due to time it takes the RNA polymerase to transcribe the gene. The average delay is the product of the gene size and the RNA translocation rate.
    
    #########################################
    ### Block for the Irates computations ###
    Occupancy23 = Occupancy(Kr,Cn23,Dn)
    Occupancy85 = Occupancy(Kr,Cn85,Dn)
    OccupancyHm = Occupancy(Kr,Hm,Dn)

    #exp_term = -kb23   *(t-Tdelay_in_response23_site1)*Occupancy23*OccupancyHm              /Imax
    #if(exp_term>1.44*pow(10,3)):
        #print("WARNING on the exponential term in the Irate_NHR23_alone_site1")
        #print(exp_term)

    #t_term = t-Tdelay_in_response23_site1
    #if(t_term<0):
        #print("t=",t)
        #print("Tdelay=",Tdelay_in_response23_site1)


    
    #Irate_NHR23_alone_site_1    = efficiency23_site1        *(1-np.exp(-kb23   *(t-Tdelay_in_response23_site1)*Occupancy23*OccupancyHm              /Imax))
    #Irate_NHR23_NHR85_site_1    = efficiency23_85_site1     *(1-np.exp(-kb23_85*(t-Tdelay_in_response23_site1)*Occupancy23*Occupancy85*OccupancyHm  /Imax))
    #Irate_NHR85_alone_site_1    = efficiency85_site1        *(1-np.exp(-kb85   *(t-Tdelay_in_response23_site1)*Occupancy85*OccupancyHm              /Imax))
    #Irate_lin4                  = Irate_NHR23_alone_site_1  + Irate_NHR23_NHR85_site_1 + Irate_NHR85_alone_site_1

    Irate_NHR23_alone_site_1    = efficiency23_site1        *(1-np.exp(-kb23   *Occupancy23*OccupancyHm              /Imax))
    Irate_NHR23_NHR85_site_1    = efficiency23_85_site1     *(1-np.exp(-kb23_85*Occupancy23*Occupancy85*OccupancyHm  /Imax))
    Irate_NHR85_alone_site_1    = efficiency85_site1        *(1-np.exp(-kb85   *Occupancy85*OccupancyHm              /Imax))
    Irate_lin4                  = Irate_NHR23_alone_site_1  + Irate_NHR23_NHR85_site_1 + Irate_NHR85_alone_site_1

    Irate_lin4*= Imax
    if(Irate_lin4<0):
        return 0
    return Irate_lin4 

def Irate_NHR23_NHR85_lin14(Imax,Dn,ka_specific_common_value,Kn,Kr,Cn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_lin14,efficiency23_lin14,efficiency_85_lin14,efficiency23_85_lin14):

    if(efficiency23_lin14 + efficiency_85_lin14 + efficiency23_85_lin14>1.1): # This condition ensures the good conversion of the efficiencies from percentages ( going in [0:100] ) to values going in [0:1].
                                                                                # This is needed to be sure that the maximal Irate doesn't exceed Imax, the maximal Initiation rate given in parameters (instead of 100*Imax if we don't do the conversion here).
        efficiency23_lin14/=100
        efficiency_85_lin14/=100
        efficiency23_85_lin14/=100

    
    Trate = Translocation_rate(Trate_T0,Q10_law_Temperature_increase,Temperature,Temperature_T0)
    
    Tdelay_in_response23_site1 = Time_delay_in_response(Trate,gene_size_lin14) # Inherent delay in the response due to time it takes the RNA polymerase to transcribe the gene. The average delay is the product of the gene size and the RNA translocation rate.

    Occupancy23 = Occupancy(Kr,Cn23,Dn)
    Occupancy85 = Occupancy(Kr,Cn85,Dn)
    OccupancyHm = Occupancy(Kr,Hm,Dn)
    
    Irate_NHR23_alone_site_1    = efficiency23_site1        *(1-np.exp(-kb23   *Occupancy23*OccupancyHm              /Imax))
    #print("Irate_NHR23_alone_site_1",Irate_NHR23_alone_site_1)
    Irate_NHR23_NHR85_site_1    = efficiency23_85_site1     *(1-np.exp(-kb23_85*Occupancy23*Occupancy85*OccupancyHm  /Imax))
    #print("Irate_NHR23_NHR85_site_1",Irate_NHR23_NHR85_site_1)
    Irate_NHR85_alone_site_1    = efficiency85_site1        *(1-np.exp(-kb85   *Occupancy85*OccupancyHm              /Imax))
    #print("Irate_NHR85_alone_site_1",Irate_NHR85_alone_site_1)
    Irate_lin14                 = Irate_NHR23_alone_site_1  + Irate_NHR23_NHR85_site_1 + Irate_NHR85_alone_site_1
    Irate_lin14*= Imax
    if(Irate_lin14<0):
        return 0
    return Irate_lin14

def Activator(A1,A2,f23,f85,delta_phi,t0,t1,dt):
    N = 1+int(round((t1-t0)/dt))
    t =  np.linspace(t0,t1,N,True)
    w23 = 2*np.pi*f23
    w85 = 2*np.pi*f85
    s = (A1*A2/4)*(1+np.sin(w23*t))*(1+np.sin(w85*t+delta_phi))
    return s

def invert_color(Color):
    Color_inverted = np.array(list(Color))
    if(max(Color_inverted)>1):
        for i in range(len(Color_inverted)):
            Color_inverted/=255
    for i in range(len(Color_inverted)):
        x = Color_inverted[i]
        if(x<0.5):
            delta = 0.5-x
            x = 0.5+delta
        elif(x>0.5):
            delta = x-0.5
            x = 0.5-delta
        Color_inverted[i] = x
    return Color_inverted

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Computing the analytical solution for the miRNA lin4 evolution during the non-null activator signal # # # # # # # #

## DURING THE ACTIVATION PEAK:
## CONDITIONS to use this solution:
# The lin4 anaytical solution is written for this chemical equation : dlin4/dt = alpha_lin4*Activator(t) - beta_lin4*lin4
# Where Activator(t) = 0.5*A1*(sin(2pi*f*t-pi/2)+1)*0.5*A2*(sin(2pi*f*t-pi/2+delta_phi)+1)
# Compute the result only during the moment NHR23 and NHR85 sin waves are superposing and only beetween twonull values of the resulting solution
# Between two null values of the analytical solution, it corresponds to the peak of {NHR23 ; NHR85} we are looking to implement.

def C_Processing(L_alpha,L_beta,L_gamma,L_Irate,L_Params_NHR85,C_t,t,dt):
    
    [alpha_LIN14]                           = L_alpha
    [beta_lin4 , beta_lin14 , beta_LIN14]   = L_beta
    [gamma_lin4_lin14 , gamma_NHR42_NHR85]  = L_gamma
    [Irate_lin4 , Irate_lin14]              = L_Irate
    [lin4 , lin14 , lin14_numeric_Mutant_Type , LIN14 , NHR85 , LIN42]  = C_t
    [A85 , phi85] = L_Params_NHR85
    
    lin4_dt = lin4  *(1- beta_lin4 *dt)                         + Irate_lin4 *dt
    lin14_dt= lin14 *(1-(beta_lin14+gamma_lin4_lin14*lin4)*dt)  + Irate_lin14*dt
    lin14_numeric_Mutant_Type_dt = lin14_numeric_Mutant_Type *(1-beta_lin14*dt)     + Irate_lin14*dt
    LIN14_dt= LIN14 *(1- beta_LIN14*dt)                         + alpha_LIN14*lin14*dt
    NHR85_dt= NHR85 + A85*np.cos(2*np.pi*t + phi85 + np.pi*dt)*np.sin(np.pi*dt) - gamma_NHR42_NHR85*LIN42*dt

##    print(lin14_dt)
##    print(1-(beta_lin14+gamma_lin4_lin14*lin4)*dt)
##    print(Irate_lin14)
##    print("STOP")
    if(lin4_dt<0):
        lin4_dt=0
    if(lin14_dt<0):
        lin14_dt=0
    if(LIN14_dt<0):
        LIN14_dt=0
    if(NHR85_dt<0):
        NHR85_dt=0
    C_t_dt   = [lin4_dt , lin14_dt , lin14_numeric_Mutant_Type_dt , LIN14_dt , NHR85_dt , LIN42]
    
    return C_t_dt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Setting the global parameters for the Activator oscillations  # # # # # # # # # # # # #

########################################################################################
# D = kB*T/(6*pi*eta*R)
# kB = 1.380649 × 10-23 m2 kg s-2 K-1
# T belongs to this range of temperature: [15:25]°C => [288.15 : 303.15]°K
# eta is approximated by the viscosity of the cytoplasm, since cytoplasm is composed mainly of water, the dynamic viscosity of water can be taken for consideration which is, 8.90 × 10−4 Pa*s at about 25°C.
# R is the radius of a hormone wich can be approximated by a value belonging to [1.8 : 3.5]nm of diameter so ~1nm radius.
# Dtheoretical = 1.9 * 10^-10 ~= 2e-10 m²/s = 2e8 nm²/s
D_cytoplasm = 1.9e-10 # Diffusion constant in m²/s
D_nucleus = 2*D_cytoplasm
Radius_protein_hormone = pow(10,-9)
########################################################################################


########################################################################################
### NHR23 parameters for the activator ###
NHR23_min = 0                  # minimum value of the NHR23 signal
NHR23_max = 100               # maximum value of the NHR23 signal
                                # This value comes from the experiments from Chris Hammel's lab about the NHRs/LIN42 to DNA binding.
A23 = NHR23_max-NHR23_min
f23 = 1                         # frequency of the NHR23 signal during each larval stage
T23 = 1/f23
t0_23 = 0
t1_23 = t0_23 + T23
phi23 = 6.14                    # phase of the NHR23 signal in RADIANS
                                # This value comes from the experiments from Chris Hammel's lab about the NHRs/LIN42 to DNA binding.
delta_phi = 0                 # enter the phase shift in DEGREES in this line (it will be converted in radians after)
delta_phi*=-np.pi/180
########################################################################################


########################################################################################
### NHR85 parameters for the activator ###
NHR85_min = 0                   # minimum value of the NHR85 signal
NHR85_max = NHR23_max*1122/505  # maximum value of the NHR85 signal
                                # This value comes from the experiments from Chris Hammel's lab about the NHRs/LIN42 to DNA binding.
A85 = NHR85_max-NHR85_min
phi85 = 2.38                    # phase of the NHR85 signal in RADIANS
                                # This value comes from the experiments from Chris Hammel's lab about the NHRs/LIN42 to DNA binding.
L_Params_NHR85 = [A85 , phi85]
f85 = 1                         # frequency of the NHR85 signal during each larval stage
T85 = 1/f85
t0_85 = T85*((3/4)-phi85-delta_phi)
t1_85 = t0_85 + T85
########################################################################################


########################################################################################
### LIN42 parameters for the activator ###
LIN42_min = 0                   # minimum value of the NHR23 signal
LIN42_max = NHR23_max*916/505   # maximum value of the NHR23 signal
                                # This value comes from the experiments from Chris Hammel's lab about the NHRs/LIN42 to DNA binding.
A42 = LIN42_max-LIN42_min
phi42 = 5.54                    # phase of the NHR23 signal in RADIANS
                                # This value comes from the experiments from Chris Hammel's lab about the NHRs/LIN42 to DNA binding.
f42 = 1                         # frequency of the NHR23 signal during each larval stage
T42 = 1/f42
t0_42 = 0
t1_42 = t0_42 + T42
########################################################################################


########################################################################################
### Initiation rate common parameters ###
Imax = 1/6                              # UNKNOWN RANGE = ]0:1/6] molecule/sec
Dn = 3*pow(10,5)
ka_specific_common_value = D_nucleus*pow(Radius_protein_hormone,-2) # UNKNOWN but function of diffusion -> f(Diff constant) = Diff/Surface of the DNA site ~= Diff/(Radius_hormone_or_protein)².
Kr = pow(10,4)                          # UNKNOWN RANGE = [10^4:10^6]
Kn = 1                                  # UNKNOWN but can be taken as the same common value for all the transcription factors.
# Cn23 =                                # VARIABLE Concentration of non-specific-DNA-protein complex ~ Concentration of NHR23 proteins in the nucleus.
# Hm =                                  # VARIABLE Concentration of hormones in the nucleus given by a diffusion model in 1D,2D or 3D.
# Cn85 =                                  # VARIABLE Concentration of non-specific-DNA-protein complex ~ Concentration of NHR85 proteins in the nucleus.
Q10_law_Temperature_increase = 2        # RANGE = [2:3]
Temperature = 25
Temperature_T0 = 15                     # The article on the logic gates describes that for sea urchin embryos it is 6 to 9 Translocation rate at 15 degrees.
Trate_T0 = 6                            # UNKNOWN for C.elegans but the article on the logic gates describes that for sea urchin embryos it is 6 to 9 at 15 degrees so we will set a standard value between 6 and 9 after.
ka_specificHm = ka_specific_common_value
kd_specificHm = ka_specificHm/(Kr*Kn)   # UNKNOWN BUT DETERMINED BECAUSE ka is given and Kr RANGE = [10^4:10^6]
ka_specific23 = ka_specific_common_value
kd_specific23 = ka_specific23/(Kr*Kn)
ka_specific85 = ka_specific_common_value
kd_specific85 = ka_specific85/(Kr*Kn)
kb23 = pow(10,-6)                                # UNKNOWN unit is "1/s" or "amount of transcription initiation rate at a given occupancy"
kb23_85 = pow(10,-6)                             # UNKNOWN unit is "1/s" or "amount of transcription initiation rate at a given occupancy"
kb85 = pow(10,-6)                                # UNKNOWN unit is "1/s" or "amount of transcription initiation rate at a given occupancy"
########################################################################################
# Irate(Imax,Dn,ka_specific_common_value,Kn,Kr,KCn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,...


########################################################################################
### lin4 parameters for the activator ###
gene_size_site1         = 94 # See this link to have more information -> https://www.ncbi.nlm.nih.gov/gene/266860 From this we know that two microRNA fragments are formed, one of 21nt and one of 62 nt.
efficiency23_site1      = 1                   # UNKNOWN RANGE = ]0:7.5]%
efficiency85_site1      = efficiency23_site1/3
efficiency23_85_site1   = 100-efficiency23_site1-efficiency85_site1
gene_size_site2         = pow(10,8)
efficiency23_site2      = efficiency23_site1
efficiency23_85_site2   = efficiency23_85_site1
efficiency85_site2      = efficiency85_site1
########################################################################################
# Irate_NHR23_NHR85_lin4(...,gene_size_site1,efficiency23_site1,efficiency_85_site1,efficiency23_85_site1,gene_size_site2,efficiency23_site2,efficiency23_85_site2,efficiency85_site2)
# => Irate_NHR23_NHR85_lin4(Imax,Dn,ka_specific_common_value,Kn,Kr,KCn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_site1,efficiency23_site1,efficiency_85_site1,efficiency23_85_site1,gene_size_site2,efficiency23_site2,efficiency23_85_site2,efficiency85_site2)


########################################################################################
### lin14 parameters for the activator ###
gene_size_lin14 = 0#20312                     # This is the number of nucleotides of the gene coding sequence in the genome of the C.elegans.
                                            # This value can be found on the official website NCBI : https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=181337
efficiency23_lin14      = efficiency23_site1     # UNKNOWN RANGE = ]0:7.5]%
efficiency85_lin14      = efficiency23_site1/3
efficiency23_85_lin14   = 100-efficiency23_site1-efficiency85_site1
########################################################################################
# Irate_NHR23_NHR85_lin14(...,gene_size_lin14,efficiency23_lin14,efficiency_85_lin14,efficiency23_85_lin14)
# => Irate_NHR23_NHR85_lin14(Imax,Dn,ka_specific_common_value,Kn,Kr,KCn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_lin14,efficiency23_lin14,efficiency_85_lin14,efficiency23_85_lin14)

# => Irate_lin4  = Irate_NHR23_NHR85_lin4 (Imax,Dn,ka_specific_common_value,Kn,Kr,KCn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_site1,efficiency23_site1,efficiency_85_site1,efficiency23_85_site1,gene_size_site2,efficiency23_site2,efficiency23_85_site2,efficiency85_site2)
# => Irate_lin14 = Irate_NHR23_NHR85_lin14(Imax,Dn,ka_specific_common_value,Kn,Kr,KCn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_lin14,efficiency23_lin14,efficiency_85_lin14,efficiency23_85_lin14)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Declaring the chemical reaction constants # # # # # # # # # # # # # # # # # # # # # # #

# ASSUMPTIONS :
### (1)
# LIN14 protein is 539 units (=3 nucleotides) long so the mRNA lin-14 is approximately 3*539 = 1617 nucleotides long : https://www.uniprot.org/uniprot/Q21446
# lin-4 miRNA is only 22 nucleotides long
#   => one miRNA lin4 is read 1617/22= (~) 80 times faster than mRNA lin14
#   => the natural production rate of lin4 and lin14 are related the same way:
#   => alpha_lin4 = 1/tau_lin4 = 80/tau_lin14 = 80*alpha_lin14
#   => alpha_lin4 = 80*alpha_lin14

### (2)
# tau = half-life time of lin4 mRNA =~ 20.7h
#   => dlin14/dt =~ (after the activator has disappeared) -beta_lin14*lin14
#   => beta_lin14 = ln(2)/tau
#   => beta_lin14 =~ 10^-5 sec^-1

### (3)
# beta_lin4 << beta_lin14

# Known parameters due to measurements - reasonable assumptions
tau_lin4 = 21*3600                      # corresponds to 21 hours in sec
tau_lin14 = 2*3600                      # corresponds to 2 hours in sec

#tau_lin4/=10*3600                       # here we divide by the duration in seconds of one larval stage (We do the assumption that a larval stage at 25°C has always the same duration which is 10*3600 sec = 10 hours)
#tau_lin14/=10*3600                      # here we divide by the duration in seconds of one larval stage (We do the assumption that a larval stage at 25°C has always the same duration which is 10*3600 sec = 10 hours)

beta_lin4           = log(2)/tau_lin4   # (2) mRNAlin4 has a half-time period equals to ~21h => beta_lin4 =~ 10^-5 sec^-1 (at 25°C)
beta_lin14          = log(2)/tau_lin14  # (3) the lin4 miRNA is relatively less degraded compared to lin14 mRNA

# Unknown parameters
gamma_lin4_lin14    = 0.000001
alpha_LIN14         = 1 #at beta_LIN14 constant, diminishing alpha_LIN14 causes a fall of the amplitude bumps, it dampens the oscillations of LIN-14
beta_LIN14          = beta_lin14
gamma_NHR42_NHR85   = 0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Declaring the inital values of concentrations # # # # # # # #

# ASSUMPTIONS :
### (1)
# N(lin14 ; t=0) =~0
# N(lin4  ; t=0) =~0

### (2)
# LIN14(t=0)>>1 , very bright when the embryo hatches and then decreasing very fast at 9h in L1
# Number of cells in the C. elegans just after hatching is 558
# N(LIN-14 ; t=0) = ?
# (see this article) -> Caenorhabditis Elegans: Development from the Perspective of the Individual Cell - https://www.ncbi.nlm.nih.gov/books/NBK26861/#:~:text=Almost%20Perfectly%20Predictable-,C.,worm%20inside%20the%20egg%20shell.

n_lin4  = 0         # (1)
n_lin14 = 100      # (1)
n_LIN14 = 1000      # (2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Generating all the activator concentration values in time # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # #  Computing the Activator signal # # # # # # # # # # # # # # #

n_T = 4                         # number of periods of the whole Activator signal to compute and show (oscillating part + constant part = 1 period)
dt = 1
t_one_LS = 10*3600              # This corresponds to the supposed averaged duration (10h) of a larval stage at a given temperature. ENTER THE VALUE IN SECONDS.
N_one_LS = 1+int(t_one_LS/dt)

t_all_LS = n_T*t_one_LS
N_all_LS = n_T*N_one_LS
 
n_post_oscillation = 0

L_time_one_Larval_Stage = np.linspace(0,t_one_LS,N_one_LS,True)
Lt_whole_signal = np.linspace(0,t_all_LS,N_all_LS,True)

L_Activator = np.empty(0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Registering the first lists # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

f_effective = 1/t_one_LS
constant_non_measurable = 20 # This value corresponds to the lower value of concentration that is not known for instance but that can be imagine (for exemple, we can set 50 for 50 molecules at least whoever the species we consider).

L_alpha = [alpha_LIN14]
L_beta  = [beta_lin4 , beta_lin14 , beta_LIN14]
L_gamma = [gamma_lin4_lin14 , gamma_NHR42_NHR85]
L_NHR23 = A23*(1+np.sin(2*np.pi*f_effective*Lt_whole_signal+phi23))+constant_non_measurable
L_NHR85 = A85*(1+np.sin(2*np.pi*f_effective*Lt_whole_signal+phi85))+constant_non_measurable
L_LIN42 = A42*(1+np.sin(2*np.pi*f_effective*Lt_whole_signal+phi42))+constant_non_measurable

n_Hormones = 50
Hm = n_Hormones
Cn23 = L_NHR23[0]
Cn85 = L_NHR85[0]
n_LIN42 = L_LIN42[0]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Creating the result lists # # # # # # # # # # # # # # # # # # #

C_t     = [n_lin4 , n_lin14 , n_lin14 , n_LIN14 , Cn85 , n_LIN42]

L_Occupancy23 = [Occupancy(Kr,Cn23,Dn)]
L_Occupancy85 = [Occupancy(Kr,Cn85,Dn)]
L_OccupancyHm = [Occupancy(Kr,Hm,Dn)]
lin4_numeric    = [n_lin4]
lin14_numeric   = [n_lin14]
lin14_numeric_Mutant_Type = [n_lin14]
LIN14_numeric   = [n_LIN14]
L_Irate_lin4 = []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Launching the main computations # # # # # # # # # # # # # # # #

boolean = 0 # This variable ensures that we won't print the warning message below in the while loop more than once.
for i in range(0,len(Lt_whole_signal),1):
    t = Lt_whole_signal[i]
    Cn23 = L_NHR23[i] # The computations begin with the first value of the NHR23
    Cn85 = L_NHR85[i]
    LIN42 = L_LIN42[i] # The computations begin with the first value of the LIN42

    Occupancy23 = Occupancy(Kr,Cn23,Dn)
    Occupancy85 = Occupancy(Kr,Cn85,Dn)
    OccupancyHm = Occupancy(Kr,Hm,Dn)
    
    Irate_lin4 = Irate_NHR23_NHR85_lin4(Imax,Dn,ka_specific_common_value,Kn,Kr,Cn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_site1,efficiency23_site1,efficiency85_site1,efficiency23_85_site1,gene_size_site2,efficiency23_site2,efficiency23_85_site2,efficiency85_site2)
    Irate_lin14 = Irate_NHR23_NHR85_lin14(Imax,Dn,ka_specific_common_value,Kn,Kr,Cn23,Hm,Cn85,Q10_law_Temperature_increase,Temperature,Temperature_T0,Trate_T0,ka_specificHm,kd_specificHm,ka_specific23,kd_specific23,ka_specific85,kd_specific85,kb23,kb23_85,kb85,gene_size_lin14,efficiency23_lin14,efficiency85_lin14,efficiency23_85_lin14)

    L_Irate = [Irate_lin4 , Irate_lin14]
    
    C_t = C_Processing(L_alpha,L_beta,L_gamma,L_Irate,L_Params_NHR85,C_t,t,dt)
    [lin4_dt , lin14_dt , lin14_numeric_Mutant_Type_dt , LIN14_dt , Cn85_dt , LIN42] = C_t # refreshing the current concentration values

    L_Occupancy23.append(Occupancy23)
    L_Occupancy85.append(Occupancy85)
    L_OccupancyHm.append(OccupancyHm)
    lin4_numeric.append(lin4_dt)
    lin14_numeric.append(lin14_dt)
    LIN14_numeric.append(LIN14_dt)
    lin14_numeric_Mutant_Type.append(lin14_numeric_Mutant_Type_dt)
    L_Irate_lin4.append(Irate_lin4)

    if(1-beta_lin14*dt - gamma_lin4_lin14*lin4_dt*dt<0 and boolean == 0):
        print(" ")
        print("WARNING")
        print("Be careful, beta_lin14 and gamma_lin4_lin14 are such that lin14 can be totally consumed.")
        print(" ")
        boolean = 1

#################################################################################
##Parameters for the plots

#Enter here the colors you want to have with or without any 'color invert' tool after to put the figures in some slides and/or presentations
color_fig1_1    = [0.6,0,1] # magenta
color_fig1_2    = [1,0,0.6] # purple
color_fig1_3    = [1,1,0]   # yellow
color_fig2      = [1,1,0]   # yellow
color_fig3_1    = [0.5,0.5,0.5] # default color - grey
color_fig3_2    = [1,0,0]   # red
color_fig4      = [0,0,1]   # blue
color_fig5      = [1,0.5,0] # orange

will_the_colors_be_inverted_after = 1 #enter 0 for no or 1 for yes if you plan to invert the colors of the figures after

###############################################################################
##Eventually Normalizing / Scaling / Changing colors / 

if (will_the_colors_be_inverted_after==1):
    print(" ")
    print("4 - Inverting the colors of the plots.")
    print(" ")
    color_fig1_1    = invert_color(color_fig1_1)
    color_fig1_2    = invert_color(color_fig1_2)
    color_fig1_3    = invert_color(color_fig1_3)
    color_fig2      = invert_color(color_fig2)
    color_fig3_1    = invert_color(color_fig3_1)
    color_fig3_2    = invert_color(color_fig3_2)
    color_fig4      = invert_color(color_fig4)
    color_fig5      = invert_color(color_fig5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Color code for the plotting part  # # # # # # # # # # # # # # #

# NHR23         - cyan - c
# NHR85         - magenta - m
# Activator     - green - g
# lin4 numeric  - red - r
# lin4 analytic - blue - b
# lin14 numeric - orange = [255,165,0] (HAS TO BE SCALED BETWEEN 0 AND 1 => DIVIDE BY 255: [1,165/255,0])
# LIN14 numeric - yellow = y

# r: red    [0,0,1]
# orange    []
# y: yellow [0,0,1]
# g: green  [0,0,1]
# c: cyan   [0,0,1]
# b: blue   [0,0,1]
# m: magenta[0,0,1]

# k: black [1,1,1]
# w: white [0,0,0]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Last parameters for the plots in themselves # # # # # # # # # #

## Normalizing the graphs :
boolean_normalize = 1 # This boolean says whether the plots are normalized or not (0 or 1)

## Rescaling the time axis :
boolean_rescaling_time = 1 # This boolean says whether the time axis is rescaled from hours to just unit interval [0:1] that corresponds to a one larval stage

if(boolean_normalize==1):
    maxi                        = max(max(L_LIN42),max(max(L_NHR23),max(L_NHR85)))
    L_NHR23                     = np.multiply(L_NHR23,1/maxi)
    L_NHR85                     = np.multiply(L_NHR85,1/maxi)
    L_LIN42                     = np.multiply(L_LIN42,1/maxi)
    lin4_numeric                = normalize(lin4_numeric)
    #lin14_numeric               = normalize(lin14_numeric)
    #lin14_numeric_Mutant_Type   = normalize(lin14_numeric_Mutant_Type)
    LIN14_numeric               = normalize(LIN14_numeric)
    L_Irate_lin4                = normalize(L_Irate_lin4)

maxi = max(lin14_numeric)
maxi2 = max(lin14_numeric_Mutant_Type)
maxi3 = max(maxi,maxi2)

lin14_numeric = np.array(lin14_numeric)/maxi3
lin14_numeric_Mutant_Type = np.array(lin14_numeric_Mutant_Type)/maxi3

if(boolean_rescaling_time==1):
    rescaling_factor = n_T/((1-pow(len(L_NHR23),-1))*(Lt_whole_signal[-1]-Lt_whole_signal[0]))
    Lt_whole_signal = (np.zeros(len(Lt_whole_signal))+rescaling_factor)*Lt_whole_signal

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Plotting the results  # # # # # # # # # # # # # # # # # # # # #

#graphe 1
plt.subplot(231)
indice_one_LS = int(len(Lt_whole_signal)/4)
plt.plot(Lt_whole_signal[0:indice_one_LS],L_NHR23[0:indice_one_LS] , color = color_fig1_1, label = 'NHR23 Pattern')
plt.plot(Lt_whole_signal[0:indice_one_LS],L_NHR85[0:indice_one_LS] , color = color_fig1_2, label = 'NHR85 Pattern')
plt.plot(Lt_whole_signal[0:indice_one_LS],L_LIN42[0:indice_one_LS] , color = color_fig1_3, label = 'L_LIN42 Pattern')
plt.xlabel('Scaled larval stage')
plt.ylabel('A(t)')
plt.legend()
plt.grid(True)

#graphe 2
plt.subplot(232)
L_Occupancy23 = L_Occupancy23[0:len(L_Occupancy23)-1]
plt.plot(Lt_whole_signal,L_Occupancy23, color = 'r', label = 'Y_NHR23')
L_Occupancy85 = L_Occupancy85[0:len(L_Occupancy85)-1]
plt.plot(Lt_whole_signal,L_Occupancy85, color = 'k', label = 'Y_NHR85')
L_OccupancyHm = L_OccupancyHm[0:len(L_OccupancyHm)-1]
plt.plot(Lt_whole_signal,L_OccupancyHm, color = 'y', label = 'Y_Hormone')
plt.xlabel('Scaled larval stages')
plt.ylabel('Occupancies')
plt.ylim(0,1)
plt.legend()
plt.grid(True)

#graphe 3
plt.subplot(233)
L_product_occupancies = np.multiply(L_Occupancy23,L_Occupancy85)
L_product_occupancies = np.multiply(L_product_occupancies,L_OccupancyHm)
#L_product_occupancies = L_product_occupancies[0:len(L_product_occupancies)-1]
plt.plot(Lt_whole_signal,L_product_occupancies, color = 'm', label = 'Y_NHR23-NHR85*Y_Hormone')
plt.xlabel('Scaled larval stages')
plt.ylabel('Resulting occupancy')
plt.ylim(0,1)
plt.legend()
plt.grid(True)

#graphe 4
#lin4 numeric et analytic sur toutes les périodes
plt.subplot(234)
#plt.title("Normalized lin4 miRNA concentration")
lin4_numeric = lin4_numeric[0:len(lin4_numeric)-1]
plt.plot(Lt_whole_signal,lin4_numeric , color = 'b',label = 'numerical solution')#color_fig3_2, label = 'numerical solution') #/max(lin4_numeric)
plt.plot(Lt_whole_signal,L_Irate_lin4, color = [0.3,0.3,0.3], label = 'Irate of lin4 gene')
plt.xlabel('Scaled larval stages')
plt.ylabel('lin4(t)')
plt.legend()
plt.grid(True)

#graphe 5
#lin 14 numeric
plt.subplot(235)
#plt.title("Normalized lin14 concentration in time")
lin14_numeric = lin14_numeric[0:len(lin14_numeric)-1]
lin14_numeric_Mutant_Type = lin14_numeric_Mutant_Type[0:len(lin14_numeric_Mutant_Type)-1]
plt.plot(Lt_whole_signal,lin14_numeric , color = 'g', label = 'lin14 num. sol. (WT)')
plt.plot(Lt_whole_signal,lin14_numeric_Mutant_Type , color = 'r', label = 'lin14 num. sol. (lin4 NULL Mt)')
plt.xlabel('Scaled larval stages')
plt.ylabel('lin14(t)')
plt.legend()
plt.grid(True)

#graphe 6
#LIN14 numeric
plt.subplot(236)
#plt.title("Normalized LIN14 concentration in time")
LIN14_numeric = LIN14_numeric[0:len(LIN14_numeric)-1]
plt.plot (Lt_whole_signal,LIN14_numeric , color = color_fig5, label = 'numerical solution')
plt.xlabel('Scaled larval stages')
plt.ylabel('LIN14(t)')
plt.legend()
plt.grid(True)

plt.show()
