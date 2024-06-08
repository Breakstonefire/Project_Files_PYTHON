import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = [7.5, 3.5]
plt.rcParams["figure.autolayout"] = True

def f(x, y):
   return np.array([i * i + j * j for j in y for i in x]).reshape(5,5)

x = y = np.linspace(-3, 3, 5)
z = f(x, y)
titre = "titre"
plt.title(titre)
plt.imshow(z, interpolation='bilinear',cmap = 'magma')

plt.show()
