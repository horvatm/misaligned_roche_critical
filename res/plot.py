"""
  Interative plotting of data: critical points and their properties at parameters
  
    theta in [-pi/2, pi/2]   Nth = 400 steps
    q in [0.1, 10]        Nq = 50 steps
    F in [0.1, 10]        NF = 50 steps
  
  Author: Martin Horvat, Jan 2018
"""

import numpy as np
import bzip2_pickle as compress

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

from matplotlib.widgets import Slider, Button, RadioButtons
from collections import defaultdict

def fn(a, a0):
  return a[np.abs(a - a0).argmin()]
     
def color_type(x):
  if x == 0:
    return "b"
  elif x == 2:
    return "r"
  elif x == -2:
    return "k"
  return "g"
    
def label_type(x):
  if x == 0:
    return "saddle"
  elif x == 2:
    return "min"
  elif x == -2:
    return "max"
  return "etc"
  
def plot_theta(p):
  global ptype
  
  ptype = 0
      
  T =  p[:,0]
  X =  p[:,1]
  Y =  p[:,2]
  
  ax.clear()
  dX = X.max() - X.min()
  ax.set_xlim([X.min()-0.1*dX, X.max()+0.1*dX])
  
  dY = Y.max() - Y.min()
  ax.set_ylim([Y.min()-0.1*dY, Y.max()+0.1*dY])
  im = ax.scatter(X, Y, c = T, cmap = cm.jet, lw = 0, s = 1)
  
  cax.clear()
  cax.get_xaxis().set_visible(True)
  cax.get_yaxis().set_visible(True)
  fig.colorbar(im, cax=cax)
 


def plot_type(p):
  global ptype
  
  ptype = 1
      
  T =  p[:,-1]
  X =  p[:,1]
  Y =  p[:,2]
  
  Tu = np.unique(T)
  
  ax.clear()
  dX = X.max() - X.min()
  ax.set_xlim([X.min()-0.1*dX, X.max()+0.1*dX])
  
  dY = Y.max() - Y.min()
  ax.set_ylim([Y.min()-0.1*dY, Y.max()+0.1*dY])
 
  for x in Tu:  
    sel = T == x
    ax.scatter(X[sel], Y[sel], c = color_type(x), lw = 0, s = 1, label = label_type(x))
  
  cax.clear()  
  cax.get_xaxis().set_visible(False)
  cax.get_yaxis().set_visible(False)
  
  ax.legend(frameon=False, loc=7, bbox_to_anchor=(1.2, 0.5), borderpad=0, handletextpad=0.1, markerscale=2)


def update_qF(val):
    global q, F, tp, p
    
    q1 = fn(qs, sq.val)
    F1 = fn(Fs, sF.val)
    
    if q1 != q or F1 != F:
      q = q1
      F = F1
            
      p = ch[(q,F)]
      
      if ptype == 0:
        plot_theta(p)
      else:
        plot_type(p)
        
      fig.canvas.draw_idle()


def update_ptype(label):
  global q, F, tp, p
  
  if label == 'theta':
    plot_theta(p)
  else:
    plot_type(p)
    
  fig.canvas.draw_idle()

# converting raw result to dictionary and 
# store it in compressed form via pickle
#ch = load2dict('res.dat')
#save_dict(ch, 'res.pkl')

# loading data 
ch = compress.load('res.pkl')

keys = np.array(ch.keys())
qs = np.unique(keys[:,0])
Fs = np.unique(keys[:,1])

# setup plots
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25, right=0.85)

# define colorbar
cax = fig.add_axes([0.9, 0.22, 0.02, 0.68])

# define first plot to be shown
ptype = 0
q = fn(qs, (qs.min()+qs.max())/2)
F = fn(Fs, (Fs.min()+Fs.max())/2)
p = ch[(q,F)]

plot_theta(p)

# define sliders
axq = plt.axes([0.25, 0.1, 0.65, 0.03])
axF = plt.axes([0.25, 0.15, 0.65, 0.03])

sq = Slider(axq, 'q', qs.min(), qs.max(), valinit = q)
sF = Slider(axF, 'F', Fs.min(), Fs.max(), valinit = F)

# define "callback functions"
sF.on_changed(update_qF)
sq.on_changed(update_qF)

# define radio button
rax = plt.axes([0.025, 0.5, 0.15, 0.15])
rax.set_title("Color based on:")
radio = RadioButtons(rax, ('theta', 'type'), active=ptype)

radio.on_clicked(update_ptype)

plt.show()
