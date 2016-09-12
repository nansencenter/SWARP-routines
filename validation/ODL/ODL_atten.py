import mod_ODL as ODL
from matplotlib import pyplot as plt
import fns_plotting as Fplt
import mod_reading as mr
import test_HYCOM_diag as THD

# parameters in WIM mode
#  TODO put in as a range of values ?


################################################################
def atten_nondim(h,T,young,visc_rp):
   # attenuation per floe
   # - for scalars (ie single grid cell)

   # h        = inputs(1)
   # om       = inputs(2)
   # young    = inputs(3)
   # visc_rp  = inputs(4)

   om       = 2*np.pi/T
   inputs   = np.array([h,om,young,visc_rp])

   outs  = Mwim.atten_youngs(inputs) # calls interface to fortran code
   # outputs  = (/damping,kice,kwtr,int_adm,
   #              ac,modT,argR,argT/)

   k_visc   = outs[0]
   kice     = outs[1]
   kwtr     = outs[2]
   alp_scat = outs[4]

   return k_visc,kice,kwtr,alp_scat
################################################################


################################################################
def Dmax2Dmean(Dmax):
   # TODO put model's formula into THD
   fsd   = THD.fsd_info_smooth(Dmax)
   return fsd.Dmean
################################################################


################################################################
def atten_dim(k_visc,alp_scat,c,Dmax):
# TODO change Dmax to Dmin, should be able to calculate from Dmax
# from floe scaling in topaz code ...
   Dmean = Dmax2Dmean(Dmax)
   alp   = 2*c*k_visc+alp_scat*c/Dmean
   return alp
################################################################


################################################################
def get_atten_coeffs(T,fice,hice,Dmax,mask):
   young    = 5.49e9
   visc_rp  = 13.
   nx,ny    = fice.shape
   alp_dim  = np.zeros((nx,ny))

   for i in range(nx):
      for j in range(ny):
         if mask[i,j]:
            c     = fice[i,j]
            h     = hice[i,j]
            dmax  = Dmax[i,j]

            # do calc's
            k_visc,kice,kwtr,alp_scat  = atten_nondim(h,T,young,visc_rp)
            alp_dim[i,j]               = atten_dim(k_visc,alp_scat,c,dmax)

   return alp_dim
################################################################


################################################################
polys = ODL.get_polys()

# loop over polys
# for poly in polys:

# object with:
# - poly and date info
# - mapping to local coords
wavdir      = None # if unknown, use PCA
size_factor = 2.# multiply bounds of polygon by this
poly        = polys[4]
opi         = ODL.ODL_poly_info(poly,wavdir=wavdir)
date        = opi.datetime

if 0:
   # NB: basemap x,y origin is not at centre of basemap
   print(opi.basemap(90,0))

if 0:
   # ===================================================================
   # test plot
   fig   = plt.figure()
   ax    = fig.add_subplot(1,1,1)

   # plot polygon
   bmap  = opi.plotmap
   bmap.plot(opi.lon_poly,opi.lat_poly,latlon=True,ax=ax)

   # plot corners of potential domain
   lon_cnr,lat_cnr   = opi.get_rectangle(latlon=True,factor=size_factor)
   bmap.plot(lon_cnr,lat_cnr,'r',latlon=True,ax=ax)
   Fplt.finish_map(bmap)

   ttl   = date.strftime('%Y%m%d %H:%M:%S')
   ax.set_title(ttl)
   plt.show(fig)
   # ===================================================================


if 1:
   # ===================================================================
   # check waves in domain (are they getting in?) 
   year  = opi.datetime.year
   bmap  = opi.plotmap

   ################################################################
   # wave data from ERA_Interim
   eraidir  = '/work/shared/nersc/ERA-I_0.125deg/'
   eraifile = 'SWH_'+str(year)+'.nc'
   ncfil    =  eraidir+eraifile
   print(ncfil)

   nci      =  mr.nc_getinfo(ncfil)
   ################################################################

   
   ################################################################
   # find nearest date
   nci_nearestdate   = opi.nearestDate(nci.datetimes)
   rec               = nci.datetimes.index(nci_nearestdate)
   print(nci_nearestdate,rec)
   ################################################################
   
   for idx in range(rec-8,rec+8):
      pobj,bmap   = nci.plot_var('swh',time_index=idx,\
            bmap=bmap,show=False)

      # plot polygon
      bmap.plot(opi.lon_poly,opi.lat_poly,'k',latlon=True,ax=pobj.ax)

      # plot corners of potential domain
      lon_cnr,lat_cnr   = opi.get_rectangle(latlon=True,factor=size_factor)
      bmap.plot(lon_cnr,lat_cnr,'r',latlon=True,ax=pobj.ax)

      dts      = date.strftime('%Y%m%d %H:%M:%S')
      pobj.ax.set_title(dts)
      figname  = 'out/swh'+dts[:8]+'_test_%2.2i.png' %(idx)

      print('saving '+figname)
      pobj.fig.savefig(figname)

      pobj.ax.cla()
      plt.close(pobj.fig)
   # ===================================================================

if 1:
   bindir   = '.' # folder with files
   gridpath = '.' # folder with regional.grid.a etc
   flist    = mr.file_list(bindir,'archv_wav','.a',gridpath=gridpath) # 

   ################################################################
   # find nearest date
   nearestdate = opi.nearestDate(flist.datetimes)
   rec         = flist.datetimes.index(nearestdate)
   print(nearestdate,rec)
   ################################################################


   ################################################################
   model_attens         = []
   model_attens_logmean = []
   model_attens_logstd  = []
   for idx in range(rec-8,rec+8):
      fobj  = flist.file_objects[idx]
      Hs    = fobj.get_var('swh')
      fice  = fobj.get_var('fice')
      Dmax  = fobj.get_var('Dmax')
      
      # find if waves-in-ice
      good           = np.array(1-fice.values.mask  ,dtype='bool')   # not land
      ice            = np.logical_and(good,fice>.15)                 # fice>0.15
      waves          = np.logical_and(good,Hs>.05)                   # Hs>5cm
      waves_in_ice   = np.logical_and(waves,ice)

      #see if inside rectangle
      mlon,mlat         = fobj.get_lonlat()
      lon_cnr,lat_cnr   = opi.get_rectangle(latlon=True,factor=size_factor)
      coords_cnr        = [(lon_cnr[i],lat_cnr[i]) for i in range(len(lon_cnr))]
      inside            = GP.maskgrid_outside_polygon(mlon,mlat,coords_cnr)

      # waves-in-ice & inside rectangle
      # - calc atten
      atten_area           = np.logical_and(waves_in_ice,inside)
      in_area              = np.zeros_like(mlon)
      atten_coeff          = np.zeros_like(mlon)
      in_area[atten_area]  = 1.

      # TODO get T from Dr. Fabs wave periods
      # TODO loop over all values of T?
      T  = 10

      # all the atten's in the reduced model waves-in-ice area
      # - each element of model_attens will be 1 model date around the polygon time
      # 1. TODO compare each element with mean & spread of Dr Fab's one
      # 2. TODO reduce area to each sub-polygon from Dr Fab
      alp_dim  = get_atten_coeffs(T,fice,hice,Dmax,in_area)[in_area]
      model_attens.append(1*alp_dim)
      model_attens_logmean.append(np.mean(np.log10(alp_dim)))
      model_attens_logstd.append(np.std(np.log10(alp_dim)))

      # outline of waves-in-ice area
      pobj,bmap   = fobj.plot_var('swh',\
            bmap=bmap,show=False)
      bmap.pcontour(mlon,mlat,in_area,'--k',levels=[.5],ax=pobj.ax)

      # plot polygon
      bmap.plot(opi.lon_poly,opi.lat_poly,'k',latlon=True,ax=pobj.ax)

      # plot corners of potential domain
      lon_cnr,lat_cnr   = opi.get_rectangle(latlon=True,factor=size_factor)
      bmap.plot(lon_cnr,lat_cnr,'r',latlon=True,ax=pobj.ax)

      dts      = date.strftime('%Y%m%d %H:%M:%S')
      pobj.ax.set_title(dts)
      figname  = 'out2/swh'+dts[:8]+'_test_%2.2i.png' %(idx)

      print('saving '+figname)
      pobj.fig.savefig(figname)

      pobj.ax.cla()
      plt.close(pobj.fig)

   # TODO plot mean & error bars vs time
   # TODO same for ODL
   # ===================================================================

