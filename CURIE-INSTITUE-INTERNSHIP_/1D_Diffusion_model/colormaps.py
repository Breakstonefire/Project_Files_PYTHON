# to see different color maps : https://matplotlib.org/3.3.1/tutorials/colors/colormaps.html

import matplotlib.pyplot as plt
import numpy as np

K = 2000    # value corresponding to the number of 'y(x)' curves to plot
N = 101     # value corresponding to the number of points to plot for each curve (= len(x))
xmin = 0
xmax = 100

x = np.linspace(xmin,xmax,N,True)
Y = [[]]
for k in range(K):
    Y.append([0 for u in range(N)])
for k in range(K):
    yk = np.linspace(0,100+10*k,N,True)
    Y[k] = yk

###################################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K curves composed of N VALUES ####
#RED TO RED

# K =     # value corresponding to the number of 'y(x)' curves to plot
# N =     # value corresponding to the number of points to plot for each curve (= len(x))

Lcolors = [[]]
for i in range(1531):
    Lcolors.append([0,0,0])

# Theoretically the number of triplet-color-values in Lcolors is equal to 6*(2^8-1) = 1530
# DECREASING THE VALUE
index = 0
for green in range(0,256,1):
    rgb = [255,green,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for red in range(254,-1,-1):
    rgb = [red,255,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for blue in range(1,256,1):
    rgb = [0,255,blue]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for green in range(254,-1,-1):
    rgb = [0,green,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for red in range(1,256,1):
    rgb = [red,0,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for blue in range(254,0,-1):
    rgb = [255,0,blue]
    Lcolors[index] = rgb
    index+=1
# DELETING THE DOUBLE VALUES IN THE END
for i in range(len(Lcolors)-1,0,-1):
    l1 = Lcolors[i]
    l2 = Lcolors[i-1]
    if l1==l2:
        print("doublon at index",i)
        del Lcolors[i]
del Lcolors[-1]
len_col = len(Lcolors)
# Theoretically the number of trilpet-color-values in Lcolors is equal to 6*(2^8-1) = 1530

# Then we normalize the values of Lcolors from [0:255] to [0:1]
for i in range(len(Lcolors)):
    for j in range(len(Lcolors[0])):
        Lcolors[i][j]/=255

#### END OF THE GENERIC BLOCK ####
##################################

###################################################################################################
#### THE FOLLOWING BLOCK CAN BE COPY-PASTE IN EVERY CODE TO PLOT K curves composed of N VALUES ####
#RED TO BLACK

# K =     # value corresponding to the number of 'y(x)' curves to plot
# N =     # value corresponding to the number of points to plot for each curve (= len(x))

Lcolors = [[]]
for i in range(1276):
    Lcolors.append([0,0,0])

# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275
# INCREASING THE VALUE
index = 0
for green in range(0,256,1):
    rgb = [255,green,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for red in range(254,-1,-1):
    rgb = [red,255,0]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# INCREASING THE VALUE
for blue in range(1,256,1):
    rgb = [0,255,blue]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for green in range(254,-1,-1):
    rgb = [0,green,255]
    Lcolors[index] = rgb
    index+=1
len_col = len(Lcolors)
# DECREASING THE VALUE
for blue in range(254,0,-1):
    rgb = [0,0,blue]
    Lcolors[index] = rgb
    index+=1
# DELETING THE DOUBLE VALUES IN THE END
for i in range(len(Lcolors)-1,0,-1):
    l1 = Lcolors[i]
    l2 = Lcolors[i-1]
    if l1==l2:
        print("doublon at index",i)
        del Lcolors[i]
del Lcolors[-1]
len_col = len(Lcolors)
# Theoretically the number of triplet-color-values in Lcolors is equal to 5*(2^8-1) = 1275

# Then we normalize the values of Lcolors from [0:255] to [0:1]
for i in range(len(Lcolors)):
    for j in range(len(Lcolors[0])):
        Lcolors[i][j]/=255

#### END OF THE GENERIC BLOCK ####
##################################

title = "title"
plt.title(title)
for k in range(K):
    yk = Y[k]
    if(K<len_col):
        ktild = int(k*len_col/K)
    elif(K==len_col):
        ktild = k
    elif(K>len_col):
        ktild = int(k*len_col/K)
    ck = Lcolors[ktild]
    plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
plt.legend()
plt.grid(True)
plt.show()
