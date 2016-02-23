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

lqm = Basemap(width=7600000,height=11200000,resolution='l',rsphere=(6378273,6356889.44891),projection='stere',lat_ts=70,lat_0=90,lon_0=-45)
