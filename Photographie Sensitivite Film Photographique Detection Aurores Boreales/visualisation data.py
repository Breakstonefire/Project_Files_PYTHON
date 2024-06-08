# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # Importing libraries # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import numpy as np
import matplotlib.pyplot as plt

L_ISO 				= [100,200,400,800,1600]
#L_ISO 				= [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600]
L_Texposition 			= [1/3200,1/2500,1/2000,1/1600,1/1250,1/1000,1/800,1/640,1/500,1/400,1/320,1/250,1/200,1/160,1/125,1/100,1/80,1/60,1/50,1/40,1/30,1/25,1/20,1/15,1/13,1/10,1/8,1/6,1/5,1/4,0.3,0.4,0.5,0.6,0.8,1,1.3,1.6,2,2.5,3.2,4,5,6,8,10,13,15]
L_Diametre_ouverture_diaphragme = [1/3.5,1/4,1/4.5,1/5,1/5.6,1/6.3,1/7.1,1/8]

I_Signal_Pur_0 = 0.10
alpha = 0.00001
L = 100*pow(10,3)
I_BPS_moyen = I_Signal_Pur_0
I_BPS_moyen = 0
I_BPL_moyen = 0.05
I_BPL_moyen = 0
Imoyen = I_Signal_Pur_0*np.exp(-alpha*L)+I_BPS_moyen+I_BPL_moyen
f = 80/1000
T = 1/3200

Hc = 10000000
for i in range(len(L_Diametre_ouverture_diaphragme)):
    Dod = L_Diametre_ouverture_diaphragme[i]
    ISO = Imoyen*T*Dod*Hc/f
    print(ISO)

L_ax = [-3.325252960978709e-08 , -2.3994580864466774e-08 , -1.1463423727570333e-08 , 6.268426303518674e-08 , 3.324785723594549e-08]
L_ay = [-1.5611511641258898e-05 , -1.0621340507513914e-05 , -8.969812476009555e-06 , 2.8228030019601423e-05 , 1.6246433061170508e-05]
L_b = [0.2451534531098893 , 0.24515603406806272 , 0.24515859957716968 , 0.2451486761240967 , 0.24514482580651636]

plt.subplot(221)
plt.plot(L_ISO,L_ax, "-^b", label = 'ax = f(ISO)')
plt.xlabel('ISO')
plt.ylabel('ax')
plt.legend()
plt.grid(True)

plt.subplot(222)
plt.plot(L_ISO,L_ay, "-^b", label = 'ay = f(ISO)')
plt.xlabel('ISO')
plt.ylabel('ay')
plt.legend()
plt.grid(True)

plt.subplot(223)
plt.plot(L_ISO,L_b, "-^b", label = 'b = f(ISO)')
plt.xlabel('ISO')
plt.ylabel('b')
plt.legend()
plt.grid(True)

plt.show()
