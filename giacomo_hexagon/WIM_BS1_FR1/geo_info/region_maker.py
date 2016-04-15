# Sript to make region.npy file to be used with OSISAF grid
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

sys.path.append('../utilities')
import mod_reading as Mrdg 
import Laplace_eqn_solution_2tw as Leqs 

# Muting error coming from invalid values of binary_mod and binary_diff
np.seterr(invalid='ignore')


iqm = Basemap(width=7600000,height=11200000,resolution='i',rsphere=(6378273,6356889.44891), \
        projection='stere',lat_ts=70,lat_0=90,lon_0=-45)


def read_txt_file_polys(fname):

    print('Returning LONGITUDES and LATITUDES')
    # read in text file:
    # each line is:
    # polygon number, lon, lat, [function value]
    fid   = open(fname,'r')
    lins  = fid.readlines()    # 1 line
    fid.close()
    
    lon = []
    lat = []

    for lin in lins[1:]:
        ss   = lin.split(';')
        lon.append(float(ss[1]))
        lat.append(float(ss[2]))
            

    return(lon,lat) 

