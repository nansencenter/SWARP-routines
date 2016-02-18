# Calling in all the needed packages
from netCDF4 import Dataset
import sys,os
import time
import datetime
import glob
import numpy as np
import subprocess
import shutil
import matplotlib
#NOTE to be used only on servers (i.e. Hexagon)
#matplotlib.use('Agg')
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from matplotlib import pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from skimage import measure as msr
from skimage import morphology as morph
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

sys.path.append('./utilities')
import mod_reading as Mrdg
import Laplace_eqn_solution_2tw as Leqs



def plotter(dataset,name=None,color=None):
    DM = np.copy(linDmax)
    plt.plot(linDmax,dataset,str(color)+'o')
    plt.ylabel(str(name))
    plt.xlabel('Dmax')
    plt.show()

def linearizer(arr):
    linarr = np.reshape(arr,arr.size)
    return(linarr)

