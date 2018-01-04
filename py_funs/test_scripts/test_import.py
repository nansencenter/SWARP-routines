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

from skimage import measure as msr
from skimage.morphology import opening, closing

from netCDF4 import Dataset
from netCDF4 import netcdftime as NCT

import rtree


from distutils.util import strtobool
from distutils.core import setup, Extension


try:
    # pygrib often doesn't work in conda
    import pygrib
    from ncepgrib2 import Grib2Encode as g2e  # part of pygrib package
    from ncepgrib2 import Grib2Decode as g2d
except:
    print('import pygrib failed')


try:
    from ncepgrib2 import Grib2Encode as g2e  # part of pygrib package
    from ncepgrib2 import Grib2Decode as g2d
except:
    print('import ncepgrib2 failed')


try:
    # linking to C++
    from Cython.Build import cythonize
except:
    print('import cython failed')


try:
   import gdal
except:
    print('import gdal failed')


try:
    import ESMF # ESMF regridding
except:
    print('import pygrib failed')


try:
    from mpi4py import MPI
except:
    print('import mpi4py failed')
