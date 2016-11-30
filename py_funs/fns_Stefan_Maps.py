import sys,os
import numpy as np
import rtree.index	as Rindex
from datetime import datetime
import shapely.geometry as SHgeom
from matplotlib import pyplot as plt

SR = os.getenv('SWARP_ROUTINES')
sys.path.append(SR+'/py_funs')
import fns_plotting as Fplt
import geometry_planar as GP
import geometry_sphere as GS
import mod_reading as MR
import MIZchar as MC

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
   def __init__(self,ll_coords,bmap,cdate=None,func_val=None):

      if type(cdate)==type('string'):
         self.datetime  = datetime.strptime(cdate,'%Y%m%d')
      else:
         self.datetime  = cdate # already a datetime object (or is None)

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
def read_txt_file(fname):
   # reads in text file (polygons from Stefan's data)
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

   Poly  = poly_info(ll_coords,xy_coords,cdate=Lino.datetime,func_vals=func_vals)
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
      # PolyL       = line_info(Poly.datetime,Poly.ll_coords[:-1],Poly.xy_coords[:-1],-1)
      PolyL       = line_info(Poly.ll_coords[:-1],Poly.xy_coords[:-1],cdate=Poly.datetime,func_val=-1)
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
      #Poly  = poly_info(Lino.datetime,ll_coords2,xy_coords2,fvals2)
      Poly  = poly_info(ll_coords2,xy_coords2,cdate=Lino.datetime,func_vals=fvals2)
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
      # Lino     = line_info(Lino.datetime,llc,xyc,Lino.func_val)
      Lino     = line_info(llc,xyc,cdate=Lino.datetime,func_val=Lino.func_val)
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
      fvals = N0*[Lino.func_val]
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
      
      # Poly  = poly_info(Lino.datetime,llc,xyc,fvals)
      Poly  = poly_info(llc,xyc,cdate=Lino.datetime,func_vals=fvals)
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
