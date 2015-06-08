from netCDF4 import Dataset
import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time
from datetime import date, timedelta
import subprocess
import shapely.geometry as shgeom
import shapefile
import glob
import smtplib
from email.mime.text import MIMEText
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg

# # not present
# import pygrib
# from ncepgrib2 import Grib2Encode as g2e # part of pygrib
# import rtree.index      as Rindex
# 
# # errors in import
# from skimage import measure as msr # no measure in skimage?
