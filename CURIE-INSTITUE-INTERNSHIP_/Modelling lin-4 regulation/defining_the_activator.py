import math
import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import signal

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
    
def Irate_NHR23_NHR85_one_enhancer_site(ka_specific,Kn,Cn23,Hm,Cn85,Temperature_T0,Trate_T0):
#######################################
### Block for the COMMON PARAMETERS ###
    # PARAMETER THAT MIGHT BE PRECISED
    Imax = 1.6  # This is the maximal rate of initiation of transcription that we set for the simulation.
                # According to (Smadar Ben-Tabou de-Leon, Eric H. Davidson paper review 2008), the initiation rate doesn't exceed 10 molecules per minute.
                # Thus, we can set the Imax at 10/60 = 1/6 rate (1/6 molecule/s)

    # PARAMETER TO SPECIFY # MIGHT BE PRECISED
    Dn = 3*pow(10,5) # Dn is the number of unoccupied nonspecific sites. It can be approximated from the amount of open chromatin (non-occluded) which is about 90% of the total genome.
                    # In C. elegans, the genome size is about 101.169 Mb (haploid) =~ 101 million baisepairs so =~ 202 million nucleotides
                    # A typical enhancer site length is framed by 100 <-> 300 baisepairs so Dn ~=~ 101 million / 100~300  ~=~ 3*10^5~1*10^6 sites

    # PARAMETER TO SPECIFY # MIGHT BE PRECISED with a model on diffusion
    ka_specific = ??? # Here we define a common value for the association rate of proteins on specific binding sites according to (Smadar Ben-Tabou de-Leon, Eric H. Davidson paper review 2008)
                    #Considering the processes of formation of a transcription factor–DNA complex:
                    # its decay, the rate of change in the amount of the factor–DNA complex can be seen to depend on complex formation and dissociation rates, "ka_spe" and "kd_spe", respectively.
                    # Here the association rate, kaS, is in terms of "number/time unit".

    # PARAMETER TO SPECIFY #
    Kn = ??? # The non-specific DNA–protein interactions can be described similarly to the specific interactions (Emerson et al. 1985).
            # However, unlike Ks, Kn, the non-specific equilibrium constant, is almost the same for every type of factor.
            # So Kn is taken as the same value for every protein in the model.

    # PARAMETER TO SPECIFY #
    Cn23 = ??? # Cn23 is the NHR23_protein-nonspecific-DNA complex concentration.
            # One can suppose that the number of free proteins + the number of protein-specific-DNA complex in the nucleus are negligible compared to the number of non-specific-DNA complex.
            # Thus, the number of non-spe-DNA complex can be approximated by the number of the proteins we choose to set in a cell.

    # PARAMETER TO SPECIFY #
    Hm = ??? # This is the relative amount of hormones in the nucleus given by the 1D/2D/3D Diffusion Model

    # PARAMETER TO SPECIFY #
    Cn85 = ??? # As described for the NHR23, one can assume that Cn85 is approximately the total number of proteins in a nucleus.

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ### Block of parameters for computing the translocation rate ###

    # PARAMETER TO SPECIFY #
    Q10_law_Temperature_increase = 2 # This parameter is the factor that describes how high is increased the translocation rate at a temperature increased by 10 degrees (Kelvin)
                                    # This parameter for all the mammalians appears to be framed by 2 <-> 3
                                    # Thus, 2 <= Q10 <= 3

    # PARAMETER TO SPECIFY #
    Temperature = 25 # This parameter is the temperature at which the experiments were done.
                    # The value has to be in Celcius degrees, the conversion in Kelvin will be made automatically.

    # PARAMETER TO SPECIFY #
    Temperature_T0 = ??? # This is the temperature at which the basal translocation rate is known.

    # PARAMETER TO SPECIFY #
    Trate_T0 = ??? # This is the value of a basal translocation rate in the nucleus for a given T0 temperature.
                    # This value is needed to compute a different translocation rate at a different temperature (according to the Q10 law).

    Translocation_rate = Translocation_rate(Trate_T0,Q10_law_Temperature_increase,Temperature,Temperature_T0)

##################################################################  
### Block for the first enhancer site for the 21nt lin-4 gene ###
    gene_size_site1 = 21 # This is the length of the gene to express in terms of nucleotid length.

    # PARAMETER TO SPECIFY #
    efficiency23_site1 = # This number in percentage represents how (IN AVERAGE) high the protein23-DNA complex is present compared to the protein23-85-DNA complex
                            # Because of what follows, this value has to be lower than 7.5%:
                            # According to Wolfgang, we have a 3 time higher activation when we only have NHR23 compared to a hypothetically NHR23 NULL mutant with only NHR85:
                            #   => this is the direct ratio of my two parameters of efficiency of the NHR23 alone and NHR85 (without quasi any NHR23) cases: efficiency23/efficiency85 = 3.
                            # Moreover, we have to keep this relation between the thre efficiencies: efficiency23 + efficiency85 + efficiency23_85 = 1 or 100%
                            # We know that when we only have NHR23 OR NHR85 (quasi NHR23 NULL mutant) the activation is lower than 10% (10% is an arbitrary value but we assume that the simultaneous effect of NHR23-NHR85 binding is highly stronger than the two separated bindings).
                            # So we can consider this relationships between the three efficiencies:
                            # (1) efficiency23/efficiency85 = 3
                            # (2) eff23 + eff85 <= 10%
                            # (3) eff23 + eff85 + eff23_85 = 100%
                            # (1) and (2) in (3) gives:
                            #   eff23_85 >= 90%
                            # from (1) and (2) we also have:
                            #   eff23 + eff23/3 = 4*eff23/3 <= 10%
                            # ⇔eff23 <= 3*10/4 %
                            # ⇔eff23 <= 7.5%

    efficiency_85_site1 = efficiency23_site1/3
    efficiency23_85_site1 = 100 - efficiency23_site1 - efficiency85_site1

    # PARAMETER TO SPECIFY #
    kb23 = ??? # kb represents the efficiency with which a given degree of occupancy causes a given amount of transcription initiation and is a measure of the activator strength (number/minute; Bolouri and Davidson, 2003).
            # Here it is the occupancy-efficiency for NHR23

    # PARAMETER TO SPECIFY #
    kb23_85 = ??? #

    # PARAMETER TO SPECIFY #
    kb85 = ??? #
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ### Block for computing the inherent time delay in response ###
    Tdelay_in_response23 = Time_delay_in_response(Translocation_rate,gene_size_site1) # Inherent delay in the response due to time it takes the RNA polymerase to transcribe the gene. The average delay is the product of the gene size and the RNA translocation rate.

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ### Block for computing the different reaction constants and the concentrations ###
    
        #########################################
        ### Sub-block for Hormones parameters ###
    ka_specificHm = ka_specific # Considering the processes of formation of a transcription factor–DNA complex:
                                # its decay, the rate of change in the amount of the factor–DNA complex can be seen to depend on complex formation and dissociation rates, "ka_spe" and "kd_spe", respectively.
                                # Here the association rate, kaS, is in terms of "number/time unit".

    # PARAMETER TO SPECIFY #
    kd_specificHm = ??? #

    KsHm = ka_specificHm/kd_specificHm
    KrHm = KsHm/Kn
        ### End of sub-block for NHR23 parameters ###
        #############################################


    
        ######################################
        ### Sub-block for NHR23 parameters ###
    ka_specific23 = ka_specific # Considering the processes of formation of a transcription factor–DNA complex:
                                # its decay, the rate of change in the amount of the factor–DNA complex can be seen to depend on complex formation and dissociation rates, "ka_spe" and "kd_spe", respectively.
                                # Here the association rate, kaS, is in terms of "number/time unit".

    # PARAMETER TO SPECIFY #
    kd_specific23 = ??? # The dissociation rate, kdS, is in "per time unit".

    Ks23 = ka_specific23/kd_specific23  # Ks depends on the chemistry of the factor-to-DNA intercation.
                                        # Thus it basically depends on the primary sequence of the transcription factor.
                                        # It is an intrisic character of the protein that reflects the "affinity" of the transcription factor for the specific site to which it binds.
                                        # More correctly, Ks indicates the stability of the site-specific DNA-protein complex.
                                        # Thus in comparing diverse interactions that display widely different values of Ks, the differences are seen to depend almost entirely on the differents values of the dissociation rate "kd_specific".
    Kr23 = Ks23/Kn  # Kr is the ratio between the specific and the non-specific equilibrium constants (Ks/Kn).
                    # Kr, the relative equilibrium constant, is a quantitative measure of specific site binding in the presence of non-specific sites, the actual case in the nucleus.
                    # Kr is usually 10^4 <–> 10^6 (Calzone et al., 1988).
                    # This parameter value is framed by Kr_min = 10^4 and Kr_max = 10^6.
        ### End of sub-block for NHR23 parameters ###
        #############################################



        ######################################
        ### Sub-block for NHR85 parameters ###
    ka_specific85 = ka_specific # The association rate of NHR85 with DNA on specific sites of activation is in terms of "number/time unit".
                                # It is the same value for every transcription factor.
    
    # PARAMETER TO SPECIFY #
    kd_specific85 = #
    
    Ks85 = ka_specific85/kd_specific85
    Kr85 = Ks85/Kn
        ### End of sub-block for NHR23 parameters ###
        #############################################

##################################################################
### Block for the second enhancer site for the 62nt lin-4 gene ###
    gene_size_site2 = 62 #
    efficiency23_site2 = efficiency23_site1 #
    efficiency23_85_site2 = efficiency23_85_site1 #
    efficiency85_site2 = efficiency85_site1 #

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ### Block for computing the inherent time delay in response ###
    Tdelay_in_response23_site1 = Time_delay_in_response(Translocation_rate,gene_size_site1) # Inherent delay in the response due to time it takes the RNA polymerase to transcribe the gene. The average delay is the product of the gene size and the RNA translocation rate.
    Tdelay_in_response23_site2 = Time_delay_in_response(Translocation_rate,gene_size_site2) # Inherent delay in the response due to time it takes the RNA polymerase to transcribe the gene. The average delay is the product of the gene size and the RNA translocation rate.

#########################################
### Block for the Irates computations ###
    Occupancy23 = Occupancy(Kr23,Cn23,Dn)
    Occupancy85 = Occupancy(Kr85,Cn85,Dn)
    OccupancyHm = Occupancy(KrHm,Hm,Dn)
    
    Irate_NHR23_alone_site_1 = efficiency23_site1   *(1-np.exp(-kb23   *(t-Tdelay_in_response23_site1)*Occupancy23*OccupancyHm              /Imax))
    Irate_NHR23_NHR85_site_1 = efficiency23_85_site1*(1-np.exp(-kb23_85*(t-Tdelay_in_response23_site1)*Occupancy23*Occupancy85*OccupancyHm  /Imax))
    Irate_NHR85_alone_site_1 = efficiency85_site1   *(1-np.exp(-kb85   *(t-Tdelay_in_response23_site1)*Occupancy85*OccupancyHm              /Imax))
    Irate_site_1 = Irate_NHR23_alone_site_1 + Irate_NHR23_NHR85_site_1 + Irate_NHR85_alone_site_1

    Irate_NHR23_alone_site_2 = efficiency23_site2   *(1-np.exp(-kb23   *(t-Tdelay_in_response23_site2)*Occupancy23*OccupancyHm              /Imax))
    Irate_NHR23_NHR85_site_2 = efficiency23_85_site2*(1-np.exp(-kb23_85*(t-Tdelay_in_response23_site2)*Occupancy23*Occupancy85*OccupancyHm  /Imax))
    Irate_NHR85_alone_site_2 = efficiency85_site2   *(1-np.exp(-kb85   *(t-Tdelay_in_response23_site2)*Occupancy85*OccupancyHm              /Imax))
    Irate_site_2 = Irate_NHR23_alone_site_2 + Irate_NHR23_NHR85_site_2 + Irate_NHR85_alone_site_2
    
    return np.array([Irate_site_1,Irate_site_2])
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Setting the global parameters for the Activator oscillations  # # # # # # # # # # # # #

NHR23_min = 0                   # minimum value of the NHR23 signal
NHR23_max = 10000               # maximum value of the NHR23 signal
A23 = NHR23_max-NHR23_min
f23 = 1                         # frequency of the NHR23 signal during each larval stage
T23 = 1/f23
t0_23 = 0
t1_23 = t0_23 + T23
phi23 = (-np.pi/2)              # phase of the NHR23 signal in RADIANS

delta_phi = 100                # enter the phase shift in DEGREES in this line (it will be converted in radians after)
delta_phi*=-np.pi/180

NHR85_min = 0                   # minimum value of the NHR85 signal
NHR85_max = 10000                # maximum value of the NHR85 signal
A85 = NHR85_max-NHR85_min
phi85 = (-np.pi/2)              # phase of the NHR85 signal in RADIANS
f85 = 1                         # frequency of the NHR85 signal during each larval stage
T85 = 1/f85
t0_85 = T85*((3/4)-phi85-delta_phi)
t1_85 = t0_85 + 85

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
tau_lin14 = 5*3600                      # corresponds to 2 hours in sec

tau_lin4/=10*3600                       # here we divide by the duration in seconds of one larval stage (We do the assumption that a larval stage at 25°C has always the same duration which is 10*3600 sec = 10 hours)
tau_lin14/=10*3600                      # here we divide by the duration in seconds of one larval stage (We do the assumption that a larval stage at 25°C has always the same duration which is 10*3600 sec = 10 hours)

beta_lin4           = log(2)/tau_lin4   # (2) mRNAlin4 has a half-time period equals to ~21h => beta_lin14 =~ 10^-5 sec^-1 (at 25°C)
beta_lin14          = log(2)/tau_lin14  # (3) the lin4 miRNA is relatively less degraded compared to lin14 mRNA

# Deduced parameters due to relations with other parameters
alpha_lin14         = 0.2*beta_lin14      # (1) this assumption comes from the ratio of gene's length of lin4 and lin14
gamma_lin4_lin14    = 0#beta_lin4*beta_lin14/700

# Unknown parameters
alpha_lin4          = alpha_lin14
alpha_LIN14         = alpha_lin14/1000 #at beta_LIN14 constant, diminishing alpha_LIN14 causes a fall of the amplitude bumps, it dampens the oscillations of LIN-14
beta_LIN14          = beta_lin14
gamma_hormone_NHRs  = 
gamma_Aactivated_unbinding =

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
n_lin14 = 3000       # (1)
n_LIN14 = 1000       # (2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Generating all the activator concentration values in time # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # #  Computing the Activator signal # # # # # # # # # # # # # # #

n_T = 4                     # number of periods of the whole Activator signal to compute and show (oscillating part + constant part = 1 period)
n_post_oscillation = 0

dt = 1/(10000*min(f23,f85))

[t23,sin23] = compute_truncated_sin(-1,f23,phi23,1,dt)
[t85,sin85] = compute_truncated_sin(-1,f85,phi85+delta_phi,1,dt)

sin23 = A23*(1+sin23)/2
t0_23 = t23[0]
t1_23 = t23[-1]

sin85 = A85*(1+sin85)/2
t0_85 = t85[0]
t1_85 = t85[-1]

t85_i = t85[-1]
t23_i = t23[0]

N_flat23 = 0
while(t85_i>t1_23 and 0<=N_flat23<=len(t85)):
    N_flat23+=1
    t85_i = t85[len(t85)-N_flat23]
N_flat85 = 0
while(t23_i<t0_85 and 0<=N_flat85<=len(t23)):
    N_flat85+=1
    t23_i = t23[N_flat85]
if(N_flat23!=N_flat85):
    diff = abs(N_flat23-N_flat85)
    print("WARNING in 'Computing the Activator signal' part")
    print("Be careful the flat parts of the two NHRs haven't the same length, {",N_flat23,";",N_flat85,"}")
    if(N_flat23<N_flat85):
        message = "NHR-23"
        N_flat23 = N_flat85
    else:
        message = "NHR-85"
        N_flat85=N_flat23
    print("WARNING in 'Computing the Activator signal' part")
    print("The flat part of "+message+" has been completed with",diff,"point(s).")
Flat_signal_23 = constant(0,N_flat23)
Pattern_function_23 = np.append(sin23,Flat_signal_23)
Flat_signal_85 = constant(0,N_flat85)
Pattern_function_85 = np.append(Flat_signal_85,sin85)

N_pattern = len(Pattern_function_23)
t_pattern = np.linspace(t0_23,t1_85,N_pattern,True)
Pattern_function = Pattern_function_23 * Pattern_function_85

N_whole_signal = n_T*N_pattern + n_post_oscillation
Lt_whole_signal = np.linspace(t0_23,n_T*t1_85,N_whole_signal,True)

L_Activator = np.empty(0)
Amplitude_correction_term = 1/np.sqrt(A23*A85)
Pattern_function*=Amplitude_correction_term

for i in range(n_T):
    # This 'for' loop adds n_T time the pattern function to the whole Activator list
    L_Activator = np.append(L_Activator,Pattern_function)

# This line creates the constant value that will be taken by the activator after the oscillations
L_Activator_post_oscillation = constant(L_Activator[len(L_Activator)-1],n_post_oscillation)

# This line adds the constant function at the end of the Activator list
L_Activator_with_post_oscillation = np.append(L_Activator,L_Activator_post_oscillation)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Registering the first lists # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

L_alpha = [alpha_lin4   , alpha_lin14   , alpha_LIN14   ]
L_beta  = [beta_lin4    , beta_lin14    , beta_LIN14    ]
L_gamma = [gamma_hormone_NHRs , gamma_Aactivated_unbinding , gamma_lin4_lin14]
C_t     = [n_lin4       , n_lin14       , n_LIN14       , L_Activator[0]]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Creating the result lists # # # # # # # # # # # # # # # # # # #
nA              = L_Activator[0]
lin4_t0         = n_lin4
lin4_numeric    = [lin4_t0]
lin14_numeric   = [n_lin14]
LIN14_numeric   = [n_LIN14]
lin4_analytic   = np.empty(0)


##plt.title("NHR23 - NHR85 combined activator")
##plt.plot(t,n23,c='b', label = 'NHR23')
##plt.plot(t,n85,c='k', label = 'NHR85')
##plt.plot(t,Activator, c='r', label = 'Activator')
##plt.legend()
##plt.grid(True)
##plt.show()
