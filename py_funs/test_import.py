import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time
from datetime import date, timedelta
import subprocess
import shapely.geometry as shgeom
import shapefile #(pyshp package)
import glob
import smtplib
from email.mime.text import MIMEText
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

# errors in import (2.7.2)
from skimage import measure as msr
# 
# # errors in import (2.7.9)
from netCDF4 import Dataset

# # not present (2.7.9)
import rtree.index      as Rindex         # not originally requested 

# pygrib doesn't work in conda
if 0:
   import pygrib
   from ncepgrib2 import Grib2Encode as g2e  # part of pygrib package
