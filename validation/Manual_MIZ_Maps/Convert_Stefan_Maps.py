import sys,os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import rtree.index	as Rindex
from datetime import datetime
import shapely.geometry as SHgeom

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import fns_plotting as Fplt
import geometry_planar as GP
import fns_Stefan_Maps as FSM

if 0:
   fmon        = '201402'
   fix_invalid = True
   check_ends  = True
   use_thresh  = False
   if 0:
      vsn   = 'SWARP_MIZ_maps_def'
   elif 0:
      vsn   = 'SWARP_MIZ_maps_def_v2'
   else:
      vsn   = 'SWARP_MIZ_maps_def_new'
else:
   fmon        = '201309'
   vsn         = 'SWARP_MIZ_maps_def'
   fix_invalid = True
   check_ends  = False
   use_thresh  = True

if 0:
   rootdir  = '/Volumes/sim/tim/Projects/SWARP/Stefan'
else:
   rootdir  = '/Users/timill/Documents/SWARP/Stefan'
tdir     = rootdir+'/'+fmon+'/'+vsn+'/txt'
idir     = rootdir+'/'+fmon+'/'+vsn+'/pdf'
outdir   = rootdir+'/'+fmon+'/'+vsn+'/polygons'

# location outputs
if not os.path.exists(outdir):
   os.mkdir(outdir)
figdir   = outdir+'/png'
if not os.path.exists(figdir):
   os.mkdir(figdir)
txtdir   = outdir+'/txt'
if not os.path.exists(txtdir):
   os.mkdir(txtdir)

#######################################################################
# make standard stereographic basemap
# - for plotting and also for limiting search area
# lat_ts is latitude of true scale.
# (lon_0,lat_0) is central point -> (x,y)=(0,0)
rad   = 10.          # approx radius of image (degrees)
xmax  = rad*111.e3   # half width of image [m]
ymax  = rad*111.e3   # half height of image [m]
cres  = 'i'          # resolution of coast (c=coarse,l=low,i=intermediate,h)
#
if fmon=='201402':
   lon_0    = 36.    # deg E
elif fmon=='201309':
   lon_0    = 5.    # deg E
lat_0    = 80.    # deg N
lat_ts   = lat_0  # deg N
bmap  = Basemap(width=2*xmax,height=2*ymax,\
                resolution=cres,projection='stere',\
                lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)
#######################################################################

if 0:
   day0  = 1
   day1  = 28
   show  = False
else:
   day   = 1
   # day   = 5
   # day   = 8
   # day   = 17
   # day   = 18
   # day   = 25
   # day   = 28
   day0  = day
   day1  = day
   show  = True
   # fix_invalid = False
   # check_ends  = False

refno       = {}
check_endsD = {}
for day in range(1,29):
   refno.update({day:0})
   check_endsD.update({day:check_ends})

if fmon=='201402':
   for day in [1,15]:
      refno[day]  = 1
   for day in [1,28]:
      check_endsD[day]  = False

for iday in range(day0,day1+1):
   cday  = '%2.2d' %(iday)
   cdate = fmon+cday
   fname = tdir+'/'+cdate+'_lonlat.txt'
   if os.path.exists(fname):
      print('\n'+60*'--'+'\nReading '+fname+'\n')
      Lins     = FSM.read_txt_file(fname)
      Linos    = []

      #################################################################
      if len(Lins)>1:
         #################################################################
         # convert to list of line_info objects
         for n,Lin in enumerate(Lins):
            # only ever 1 part, Lin[0]
            if n==refno[iday]:
               # 1st is usually MIZ but not always
               # - define manually
               func_val = 1 #MIZ
               # Lino_ref = FSM.line_info(cdate,Lin[0],bmap,func_val)
               Lino_ref = FSM.line_info(Lin[0],bmap,cdate=cdate,func_val=func_val)
            else:
               func_val = 0 #ice edge
               #Lino     = FSM.line_info(cdate,Lin[0],bmap,func_val)
               Lino     = FSM.line_info(Lin[0],bmap,cdate=cdate,func_val=func_val)
               Linos.append(Lino)
         #################################################################
               
         Polys,Linos2 = FSM.Linos2polys(Lino_ref,Linos,\
               fix_invalid=fix_invalid,check_ends=check_endsD[iday],use_thresh=use_thresh)

         # write text files with polygons' coords
         fout  = txtdir+'/'+cdate+"_polys.txt"
         FSM.write_polys(fout,Polys)
         print('>'+fout)

         if 1:
            # TEST PLOT - with polygons
            #
            fig   = plt.figure()
            ax1   = fig.add_subplot(1,1,1)
            lstil = ['-','--']

            for Poly in Polys:
               x,y   = np.array(Poly.xy_coords).transpose()
               shp   = SHgeom.Polygon(Poly.xy_coords)
               if shp.is_valid:
                  bmap.plot(x,y,'k',linewidth=2,ax=ax1)
               else:
                  bmap.plot(x,y,'c',linewidth=2,ax=ax1)


            # MIZ line
            x,y   = np.array(Lino_ref.xy_coords).transpose()
            bmap.plot(x,y,linestyle=lstil[Lino_ref.func_val],ax=ax1)
            bmap.plot(x[:1],y[:1],'v',ax=ax1) # 1st point

            # ice edges
            for Lino in Linos:
               x,y   = np.array(Lino.xy_coords).transpose()
               bmap.plot(x,y,linestyle=lstil[Lino.func_val],ax=ax1)
               bmap.plot(x[:1],y[:1],'^',ax=ax1) # 1st point


            Fplt.finish_map(bmap,ax=ax1)
            if show:
               plt.show(fig)
            else:
               figname  = figdir+'/'+cdate+"_polys.png"
               print('>'+figname)
               fig.savefig(figname)
            ax1.cla()
            plt.close(fig)
