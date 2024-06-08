# to see different color maps : https://matplotlib.org/3.3.1/tutorials/colors/colormaps.html

import matplotlib.pyplot as plt
import numpy as np

K = 1529 # value corresponding to the number of 'y(x)' curves to plot
N = 100 # value corresponding to the number of points to plot for each curve (= len(x))
xmin = 0
xmax = 100

x = np.linspace(xmin,xmax,N,True)
Y = [[]]
for k in range(K):
    Y.append([0 for u in range(N)])
for k in range(K):
    yk = np.linspace(0,100+10*k,N,True)
    Y[k] = yk

Lcolors = [[]]
for i in range(1535):
    Lcolors.append([0,0,0])

# Theoretically the number of trilpet-color-values in Lcolors is equal to 6*(2^8-1) = 1530
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
print(len_col)
# Theoretically the number of trilpet-color-values in Lcolors is equal to 6*(2^8-1) = 1530

# Then we normalize the values of Lcolors from [0:255] to [0:1]
for i in range(len(Lcolors)):
    for j in range(len(Lcolors[0])):
        Lcolors[i][j]/=255

title = "title"
plt.title(title)
for k in range(K):
    yk = Y[k]
    if(K<len_col):
        ktild = int(k*len_col/K)
        print("Kinf",ktild)
    elif(K==len_col):
        ktild = k
        print("K==",ktild)
    elif(K>len_col):
        if(k<len_col):
            ktild = k
            #print("Ksup ; k<=",ktild)
        if(k>=len_col):
            ktild = int(10*((k/len_col)-int(k/len_col)))
            #print("Ksup ; k>",ktild)
    ck = Lcolors[ktild]
    plt.plot(x,yk, color = ck , marker = '',linewidth = 1,  label = '')
plt.legend()
plt.grid(True)
plt.show()
