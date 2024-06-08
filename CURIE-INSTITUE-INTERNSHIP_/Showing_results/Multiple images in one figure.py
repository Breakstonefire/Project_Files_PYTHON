#add_subplot(rows, columns, i)
import numpy as np
import matplotlib.pyplot as plt

code = 4

width=5
height=5
rows = 2
cols = 2
axes=[]

#CODE 1
if (code==1):
    fig=plt.figure()

    for i in range(rows*cols):
        b = np.random.randint(7, size=(height,width))
        axes.append( fig.add_subplot(rows, cols, i+1))
        subplot_title=("Subplot"+str(i))
        axes[-1].set_title(subplot_title)
        plt.imshow(b)
    fig.tight_layout()
    plt.show()

#CODE 2
if (code==2):
    fig=plt.figure()

    x=np.linspace(-3,3,100)
    y1=np.sin(x)
    y2=1/(1+np.exp(-x))

    axes = []

    for i in range(cols*rows):
        b = np.random.randint(10, size=(height,width))
        axes.append(fig.add_subplot(rows, cols, i+1))
        subplot_title=("Subplot"+str(i))
        axes[-1].set_title(subplot_title)
        plt.imshow(b)

    axes[1].plot(x,y1)
    axes[3].plot(x,y2)
    fig.tight_layout()
    plt.show()

#CODE 3
if (code==3):
    x=np.linspace(0,3,100)
    y1=np.sin(x)
    y2=1/(1+np.exp(-x))
         
    figure, axes = plt.subplots(nrows=rows, ncols=cols)

    for a, b in enumerate(axes.flat):
        image = np.random.randint(7, size=(height,width))
        b.imshow(image, alpha=0.25)
        r = a // cols
        c = a % cols
        subtitle="Row:"+str(r)+", Col:"+str(c)
        b.set_title(subtitle)

    axes[0][1].plot(x, y1)
    axes[1][1].plot(x,y2)

    figure.tight_layout()
    plt.show()


#CODE 4
if (code==4):
    def display_multiple_img(images, rows = 1, cols=1):
        figure, ax = plt.subplots(nrows=rows,ncols=cols)
        for ind,title in enumerate(images):
            ax.ravel()[ind].imshow(images[title])
            ax.ravel()[ind].set_title(title)
            ax.ravel()[ind].set_axis_off()
        plt.tight_layout()
        plt.show()

    total_images = 4
    images = {'Image'+str(i): np.random.rand(100, 100) for i in range(total_images)}

    display_multiple_img(images, 2, 2)
