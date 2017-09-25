import sys,os
import numpy as np
import pyproj
import struct

from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.font_manager import FontProperties
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1 import make_axes_locatable

import time
from datetime import date, timedelta
import subprocess
import shapely.geometry as shgeom
import shapefile #(pyshp package)
import glob
import smtplib
from email.mime.text import MIMEText
from scipy.interpolate import griddata as grd
import scipy.ndimage as NDI

from operator import itemgetter as itg
from getopt import getopt
import shutil

# errors in import (2.7.2)
from skimage import measure as msr
from skimage.morphology import opening, closing
# 
# # errors in import (2.7.9)
from netCDF4 import Dataset
from netcdftime import netcdftime as NCT

# # not present (2.7.9)
import rtree.index      as Rindex         # not originally requested 


from distutils.util import strtobool
if 0:
   import ESMF # ESMF regridding
   from mpi4py import MPI

# pygrib doesn't work in conda
if 0:
   import pygrib
   from ncepgrib2 import Grib2Encode as g2e  # part of pygrib package
   from ncepgrib2 import Grib2Decode as g2d

if 1:
   # linking to C++
   from distutils.core import setup, Extension
   from Cython.Build import cythonize

