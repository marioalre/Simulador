import numpy as np
from matplotlib.widgets import Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')
        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices//2
        self.im = ax.imshow(self.X[:, :, self.ind], cmap='gray')
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices        
        self.update()

    def contrast(self, event):
        print('Changing contrast')
        print(smax.val)
        self.im.set_clim([0,smax.val])        
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


import tkinter as tk
root = tk.Tk()
fig = Figure()
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(fill="both", expand=True)

im = np.array(np.random.rand(10,10,10))

ax = fig.subplots(1,1)

axmax  = fig.add_axes([0.25, 0.01, 0.65, 0.03])
smax = Slider(axmax, 'Max', 0, np.max(im), valinit=50)
tracker = IndexTracker(ax, im)
canvas.mpl_connect('scroll_event', tracker.onscroll)
canvas.mpl_connect('button_release_event', tracker.contrast) #add this for contrast change
root.mainloop()
