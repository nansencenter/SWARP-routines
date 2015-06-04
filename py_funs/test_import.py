from netCDF4 import Dataset
import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time
from datetime import date, timedelta
import subprocess
import pygrib
from ncepgrib2 import Grib2Encode as g2e
import shapely.geometry as shgeom # http://toblerity.org/shapely/manual.html
import rtree.index      as Rindex
from skimage import measure as msr
import shapefile
import glob
import smtplib
from email.mime.text import MIMEText
from scipy.interpolate import griddata as grd
from operator import itemgetter as itg
