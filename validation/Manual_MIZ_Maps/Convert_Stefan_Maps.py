import sys,os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.basemap import Basemap, cm
import rtree.index	as Rindex
from datetime import datetime
import shapely.geometry as SHgeom

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import fns_plotting as Fplt
import geometry_planar as GP

#######################################################################
def write_polys(fname,Polys):
   count = -1
   fid   = open(fname,'w')
   fid.write('Poly   lon   lat   MIZ\n')
   for Poly in Polys:
      count = count+1
      llc   = Poly.ll_coords
      for i,ll in enumerate(llc):
         ss = '%d   %f   %f   %d\n' %(count,ll[0],ll[1],Poly.func_vals[i])
         fid.write(ss)
   fid.close()
   return
#######################################################################

#######################################################################
class line_info:
   # make object with helpful info:
   def __init__(self,cdate,ll_coords,bmap,func_val):

      if type(cdate)==type('hi'):
         self.datetime  = datetime.strptime(cdate,'%Y%m%d')
      else:
         self.datetime  = cdate # already a datetime object

      self.func_val  = func_val
      llc            = 1*ll_coords

      if type(bmap)==type([]):
         # already given list of x,y
         coords   = 1*bmap
      else:
         # get x,y from basemap
         lon,lat  = np.array(llc).transpose()
         x,y      = bmap(lon,lat)
         coords   = list(np.array([x,y]).transpose())

      # make sure coords is a list of tuples
      coords   = [tuple(xy) for xy in coords]

      ###################################################
      # check for repeats
      xyc2 = []
      llc2 = []
      for i,cc in enumerate(coords):
         if cc not in xyc2:
            xyc2.append(cc)
            llc2.append(tuple(llc[i]))

      self.ll_coords = llc2
      self.xy_coords = xyc2
      self.length    = len(llc2)
      self.perimeter = GP.calc_perimeter(xyc2)
      ###################################################

      ###################################################
      # make rtree index
      idx   = Rindex.Index()
      for i in range(self.length):
         xp,yp = xyc2[i]
         idx.insert(i,(xp,yp,xp,yp)) # a point is a rectangle of zero side-length

      self.index  = idx
      ###################################################

      return
#######################################################################

#######################################################################
class poly_info:
   # make object with helpful info:
   def __init__(self,dt,ll_coords,xy_coords,func_vals):

      self.datetime  = dt
      self.length    = len(xy_coords)
      #
      fvals = list(func_vals)
      xyc2  = [tuple(xyp) for xyp in xy_coords] # list of tuples (otherwise logical ops are difficult)
      llc2  = [tuple(llp) for llp in ll_coords] # list of tuples (otherwise logical ops are difficult)

      critter  =     (len(ll_coords)!=len(xy_coords))\
                  or (len(ll_coords)!=len(func_vals))\
                  or (len(xy_coords)!=len(func_vals))
      if critter:
         print('lengths (ll,xy,f):')
         print(len(ll_coords))
         print(len(xy_coords))
         print(len(func_vals))
         raise ValueError('Inconsistent lengths of inputs')

      # close if necessary
      if xyc2[0]!=xyc2[-1]:
         llc2.append(llc2[0])
         xyc2.append(xyc2[0])
         fvals.append(fvals[0])

      x,y   = np.array(xyc2).transpose()
      area  = GP.area_polygon_euclidean(x,y)
      if area<0:
         # reverse order (want anti-clockwise ordering)
         fvals.reverse()
         xyc2.reverse()
         llc2.reverse()

      self.area      = abs(area)
      self.func_vals = fvals
      self.ll_coords = llc2
      self.xy_coords = xyc2
      self.perimeter = GP.calc_perimeter(xyc2,closed=True)

      # make rtree index
      idx   = Rindex.Index()
      for i in range(self.length):
         xp,yp = self.xy_coords[i]
         idx.insert(i,(xp,yp,xp,yp)) # a point is a rectangle of zero side-length
      self.index  = idx

      return
#######################################################################

#######################################################################
def read_txt_file(fname):
   # reads in text file from Stefan,
   # outputs list of polygons

   fid   = open(fname,'r')
   lin   = fid.readline()    # 1 line
   fid.close()
   #
   lins     = lin.strip('\n').split('Line')[1:] # line 1 = ice edge, line 2 = MIZ-pack line
   Lins     = []

   # loop over line types
   for ll in lins:
      ll2   = ll.split('Part ')[1:]
      parts = []

      ################################################################################
      # loop over parts
      for part in ll2:
         ll3         = part.split('n')[1:-1] # now a list with members ['lon1 lat1','lon2 lat2',...]

         # loop over coords
         ll_coords   = [] # this will be a list of tuples [(lon1,lat1),(lon2,lat2),...]
         for cc in ll3:
            cc2      = cc.split()
            lonlat   = (float(cc2[0]),float(cc2[1]))
            ll_coords.append(lonlat)

         parts.append(ll_coords)# list of list of tuples
      ################################################################################

      Lins.append(parts)

   return Lins
#######################################################################

#######################################################################
def shapely2poly(shp,Lino,Lino_ref):
   # if shp.buffer(0) has worked (shp.is_valid = True),
   # extract info from shp

   xy_coords   = list(shp.exterior.coords)
   ll_coords   = []
   func_vals   = []
   for cc in xy_coords:
      x,y   = cc
      i0    = list(Lino.index.nearest(cc))[0]
      i1    = list(Lino_ref.index.nearest(cc))[0]
      # print(i0,i1,Lino.length,Lino_ref.length)

      d0    = np.sqrt( pow(x-Lino    .xy_coords[i0][0],2) + pow(y-Lino    .xy_coords[i0][1],2) )
      d1    = np.sqrt( pow(x-Lino_ref.xy_coords[i1][0],2) + pow(y-Lino_ref.xy_coords[i1][1],2) )
      if d0<d1:
         # part of "Lino" line (usually ice edge)
         ll_coords.append(Lino.ll_coords[i0])
         func_vals.append(Lino.func_val)
      else:
         # part of "Lino_ref" line (usually MIZ edge)
         ll_coords.append(Lino_ref.ll_coords[i1])
         func_vals.append(Lino_ref.func_val)

   Poly  = poly_info(Lino.datetime,ll_coords,xy_coords,func_vals)
   return Poly
#######################################################################

#######################################################################
def check_refpoints(Poly,Lino_ref,Lino):
   # check end points of Lino_ref are not too far from their neighbours:
   coords   = Lino_ref.xy_coords
   fv       = np.array(Poly.func_vals)
   nref     = np.arange(len(fv),dtype=int)[fv==Lino_ref.func_val]
   # print(fv)
   # print(nref)
   # print(len(fv),len(Poly.xy_coords))

   endc  = []
   endc2 = []
   nv1   = []
   nv2   = []
   endi  = nref[[0,-1]]
   endi2 = nref[[1,-2]]

   LR = SHgeom.LineString(Poly.xy_coords)
   for n in range(2):
      m  = endi[n]
      cc = tuple(Poly.xy_coords[m])
      nv1 .append(list(Lino_ref.index.nearest(cc))[0])
      endc.append(cc)
      #
      m  = endi2[n]
      cc = tuple(Poly.xy_coords[m])
      nv2  .append(list(Lino_ref.index.nearest(cc))[0])
      endc2.append(cc)

   ##################################################################
   for m,n in enumerate(nv1):
      # print(n,m)
      # loop over endpoints of Lino_ref

      # make line_info object from Poly
      # => want no repeats (not closed curve now)
      lon0,lat0   = np.array(Poly.ll_coords[:-1]).transpose()
      PolyL       = line_info(Poly.datetime,Poly.ll_coords[:-1],Poly.xy_coords[:-1],-1)
      Pfunc_vals  = Poly.func_vals[:-1]

      cc0   = endc[m]                           # own coords
      in0   = list(PolyL.index.nearest(cc0))[0] # own index in polygon
      cc1   = endc2[m]                          # coords of MIZ neighbour in polygon
      in1   = list(PolyL.index.nearest(cc1))[0] # index of MIZ neighbour in polygon
      #
      in2   = list(Lino.index.nearest(cc0))[0]  # index of nearest point on ice edge
      cc2   = Lino.xy_coords[in2]               # coords of nearest point on ice edge
      in2   = list(PolyL.index.nearest(cc2))[0] # index of this point in polygon

      # dx = cc2[0]-cc0[0]
      # dy = cc2[1]-cc0[1]
      # LL = SHgeom.LineString([(cc0[0]+.01*dx,cc0[1]+.01*dy),\
      #                         (cc0[0]+.99*dx,cc0[1]+.99*dy)])
      # print(cc0,cc2)
      # print(LL)
      # if LL.intersects(LR):
      #    print('bad')
      #    P1 = LL.intersection(LR)
      #    
      #    print(P1)

      # print(cc0,in0)
      # print(cc1,in1)
      # print(cc2,in2)

      if in1<in0:
         # increasing index moves away from MIZ line
         if in2>in0:
            # remove points between indices
            # - get points <= in0
            xy_coords2  = PolyL.xy_coords[:in0+1]
            ll_coords2  = PolyL.ll_coords[:in0+1]
            fvals2      = Pfunc_vals[:in0+1]

            # - get points >= in2
            xy_coords2.extend(PolyL.xy_coords[in2:])
            ll_coords2.extend(PolyL.ll_coords[in2:])
            fvals2    .extend(Pfunc_vals[in2:])
         else:
            # remove points >in0 & <in2
            # ie keep points in2<=index<=in0
            xy_coords2  = PolyL.xy_coords[in2:in0+1]
            ll_coords2  = PolyL.ll_coords[in2:in0+1]
            fvals2      = Pfunc_vals[in2:in0+1]
      else:
         # decreasing index moves away from MIZ line
         if in2<in0:
            # remove points between indices
            # - get points <= in2
            xy_coords2  = PolyL.xy_coords[:in2+1]
            ll_coords2  = PolyL.ll_coords[:in2+1]
            fvals2      = Pfunc_vals[:in2+1]

            # - get points >= in0
            xy_coords2.extend(PolyL.xy_coords[in0:])
            ll_coords2.extend(PolyL.ll_coords[in0:])
            fvals2    .extend(Pfunc_vals[in0:])
         else:
            # remove points >in2 & <in0
            # ie keep points in0<=index<=in2
            xy_coords2  = PolyL.xy_coords[in0:in2+1]
            ll_coords2  = PolyL.ll_coords[in0:in2+1]
            fvals2      = Pfunc_vals[in0:in2+1]

      #######################################################################
      # make new polygon
      Poly  = poly_info(Lino.datetime,ll_coords2,xy_coords2,fvals2)
      #######################################################################

   # end loop over endpoints of Lino_ref
   ##########################################################################

   return Poly
#######################################################################

#######################################################################
def Linos2polys(Lino_ref,Linos,fix_invalid=True,check_ends=True,use_thresh=False):

   ####################################################################
   def make_poly(Lino_ref,Lino,use_thresh=False):
      xy_coords   = Lino.xy_coords
      ll_coords   = Lino.ll_coords
      xy_coords1  = Lino_ref.xy_coords
      ll_coords1  = Lino_ref.ll_coords

      ########################################################################
      # first remove points from Lino that are already in Lino_ref
      # - should no longer be duplicate points 
      #   as repeats in indiv lines are removed in line_info
      xyc   = []
      llc   = []
      for i,cc in enumerate(xy_coords):
         if cc not in xy_coords1:
            xyc.append(cc)
            llc.append(ll_coords[i])

      lon,lat  = np.array(llc).transpose()
      Lino     = line_info(Lino.datetime,llc,xyc,Lino.func_val)
      ########################################################################

      cc0   = Lino.xy_coords[0]
      cc1   = Lino.xy_coords[-1]
      N0    = Lino.length
      j0    = list(Lino_ref.index.nearest(cc0))[0]
      j1    = list(Lino_ref.index.nearest(cc1))[0]
      jj    = [j0,j1]
      ccs   = [cc0,cc1]

      # print(j0,j1,N0,len(Lino.xy_coords))
      # print(cc0,Lino_ref.xy_coords[j0])
      # print(cc1,Lino_ref.xy_coords[j1])
      # print('\n')

      #####################################################################################
      if use_thresh:
         # if over a certain threshold,
         # connect to nearest end point
         x0,y0    = Lino_ref.xy_coords[0]
         x1,y1    = Lino_ref.xy_coords[-1]
         thresh   = 150.e3
         # print((x0,y0),(x1,y1))

         # loop over ends of Lino
         for i,cc in enumerate(ccs):
            j     = jj[i]
            x,y   = Lino_ref.xy_coords[j]
            dd    = np.sqrt( pow(x -cc[0],2) + pow(y -cc[1],2) ) # length to nearest point
            d0    = np.sqrt( pow(x0-cc[0],2) + pow(y0-cc[1],2) ) # length to 1st end point
            d1    = np.sqrt( pow(x1-cc[0],2) + pow(y1-cc[1],2) ) # length to last end point
            # print(j,cc)
            # print(dd/1.e3,thresh/1.e3,d0/1.e3,d1/1.e3)
            # print(jj)
            if dd>thresh:
               if d0<d1:
                  jj[i] = 0
               else:
                  jj[i] = len(Lino_ref.xy_coords)-1
            # print(jj)
         j0,j1 = jj
      #####################################################################################

      # print(j0,j1,N0,len(Lino.xy_coords))
      # print(cc0,Lino_ref.xy_coords[j0])
      # print(cc1,Lino_ref.xy_coords[j1])
      # print('\n')
      #
      fvals = list(np.zeros(N0,dtype=int)+Lino.func_val)
      if j1>j0:
         xy_coords2  = Lino_ref.xy_coords[j0:j1+1]
         xy_coords2.reverse()
         xyc.extend(xy_coords2)
         #
         ll_coords2  = Lino_ref.ll_coords[j0:j1+1]
         ll_coords2.reverse()
         llc.extend(ll_coords2)
         #
         fvals2   = list(np.zeros(len(xy_coords2),dtype=int)+Lino_ref.func_val)
         fvals.extend(fvals2)
      else:
         xy_coords2  = Lino_ref.xy_coords[j1:j0+1]
         xyc.extend(xy_coords2)
         #
         ll_coords2  = Lino_ref.ll_coords[j1:j0+1]
         llc.extend(ll_coords2)
         #
         fvals2   = list(np.zeros(len(xy_coords2),dtype=int)+Lino_ref.func_val)
         fvals.extend(fvals2)
      
      Poly  = poly_info(Lino.datetime,llc,xyc,fvals)
      return Poly,Lino
   ####################################################################

   Polys    = []
   Linos2   = []
   for Lino in Linos:
      Poly,Lino   = make_poly(Lino_ref,Lino,use_thresh=use_thresh)
      Linos2.append(Lino)

      # print(Poly.func_vals)
      ##########################################################################
      if fix_invalid:
         shp   = SHgeom.Polygon(Poly.xy_coords)
         if not shp.is_valid:
            # make valid polygon with shp.buffer(0)
            print('invalid - correcting with shapely')
            shp   = shp.buffer(0)
            if hasattr(shp,'geoms'):
               print('> multi-polygon')
               # get largest poly from multipolygon:
               L  = 0
               nn = 0
               for n,shpp in enumerate(shp.geoms):
                  # print(shpp.length)
                  # print(shpp.area)
                  if shpp.length>L:
                     nn = n
                     L  = shpp.length
               Poly  = shapely2poly(shp.geoms[nn],Lino,Lino_ref)
            else:
               print('> one polygon')
               Poly  = shapely2poly(shp,Lino,Lino_ref)
      ##########################################################################

      ##########################################################################
      if check_ends:
         # check end points of Lino_ref are not too far from their neighbours:
         Poly  = check_refpoints(Poly,Lino_ref,Lino)
         if 0:
            shp   = SHgeom.Polygon(Poly.xy_coords)
            if not shp.is_valid:
               # make valid polygon with shp.buffer(0)
               print('invalid - correcting with shapely')
               shp   = shp.buffer(1).buffer(-1)
               if hasattr(shp,'geoms'):
                  print('> multi-polygon')
                  # get largest poly from multipolygon:
                  L  = 0
                  nn = 0
                  for n,shpp in enumerate(shp.geoms):
                     print(shpp.length)
                     print(shpp.area)
                     if shpp.length>L:
                        nn = n
                        L  = shpp.length
                  Poly  = shapely2poly(shp.geoms[nn],Lino,Lino_ref)
               else:
                  print('> one polygon')
                  Poly  = shapely2poly(shp,Lino,Lino_ref)
      ##########################################################################

      Polys.append(Poly)
   return Polys,Linos2
#######################################################################

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

rootdir  = '/Volumes/sim/tim/Projects/SWARP/Stefan'
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

if 1:
   day0  = 1
   day1  = 28
   show  = False
else:
   # day   = 1
   # day   = 5
   # day   = 8
   # day   = 17
   # day   = 18
   # day   = 25
   day   = 28
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
      Lins     = read_txt_file(fname)
      Linos    = []
      Lens     = []

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
               Lino_ref = line_info(cdate,Lin[0],bmap,func_val)
            else:
               func_val = 0 #ice edge
               Lino     = line_info(cdate,Lin[0],bmap,func_val)
               Linos.append(Lino)
         #################################################################
               
         Polys,Linos2 = Linos2polys(Lino_ref,Linos,\
               fix_invalid=fix_invalid,check_ends=check_endsD[iday],use_thresh=use_thresh)

         if 0:
            for n,Lin in enumerate(Lins):
               print('n (0=MIZ edge,1=ice edge)), no of parts')
               print(n,len(Lin))
         elif 0:
            # TEST PLOT
            #
            fig   = plt.figure()
            ax1   = fig.add_subplot(1,1,1)
            lstil = ['-','--','-.']
            for n,Lin in enumerate(Lins):
               for part in Lin:
                  lon,lat  = np.array(part).transpose()
                  x,y      = bmap(lon,lat)
                  bmap.plot(x,y,linestyle=lstil[n],ax=ax1)
                  # bmap.plot(x[:1],y[:1],'^',ax=ax1) # 1st point

            Fplt.finish_map(bmap,ax=ax1)
            plt.show(fig)
            ax1.cla()
            plt.close(fig)
         elif 1:
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

            fout  = txtdir+'/'+cdate+"_polys.txt"
            print('>'+fout)
            write_polys(fout,Polys)

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
