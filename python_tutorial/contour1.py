#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys,os

figdir   = 'out'
if not os.path.exists(figdir):
   os.mkdir(figdir)

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

delta = 0.025
x = np.arange(-3.0, 3.0, delta)
y = np.arange(-2.0, 2.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = 10.0 * (Z2 - Z1)



# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
fig   = plt.figure()
cs    = plt.contour(X, Y, Z)
# plt.clabel(cs, inline=1, fontsize=10)
plt.title('Simplest default with labels')
plt.savefig(figdir+'/test1.png')
plt.close()
fig.clf()

coll  = cs.collections
fig   = plt.figure()
for nl in range(len(coll)):
   # loop over levels
   p     = coll[nl].get_paths()
   print(nl,cs.levels[nl],len(p))
   for ns in range(len(p)):
      # loop over segments
      v     = p[ns].vertices
      x     = v[:,0]
      y     = v[:,1]
      plt.plot(x,y)

plt.savefig(figdir+'/test1a.png')
plt.close()
fig.clf()
