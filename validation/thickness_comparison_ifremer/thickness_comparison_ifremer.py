import os,sys
if 1:
   from   mpl_toolkits.basemap import pyproj
else:
   import pyproj
import mod_reading          as mr
import numpy                as np
import geometry_sphere      as GS
import fns_plotting         as FP
from matplotlib import pyplot as plt
import datetime as dtm


###############################################
def wgs84_ellipsoid_mat():
   a,fi = 6378137,298.257223563
   # f=flattening=(a-b)/a=1/fi
   f    = 1./fi
   b    = (1-f)*a
   ecc  = np.sqrt(1-pow(b/a,2)) # eccentricity
   # print(a,fi,f,b,ecc)
   # print(a/(a-b))
   return a,ecc
###############################################


###############################################

bmap    = FP.start_HYCOM_map('Arctic')
ODLmap  = pyproj.Proj('+init=EPSG:3413')
"""
srs='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45'+\
    ' +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
"""


# directories
compname  = os.getenv('HOSTNAME')
if 'hexagon' in compname:
    #hexagon
    HEX        = True
    odir       = '/work/shared/nersc/msc/cs2_smos_ifremer/' # path to observations
    mdir       = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice_erai2/'
    gridpath   = None
else:
    HEX        = False
    mdir       = '/mnt/sda1/work/Model-outputs/thickness_comp_ifremer/TP4/'
    odir       = '/mnt/sda1/work/Model-outputs/thickness_comp_ifremer/cs2_smos_ifremer'
    gridpath   = mdir+'/topo'


lonlat_file = 'lonlat_grid.nc'
if 0:
    olist       = ['cs2_smos_ice_thickness_20150302_20150308.nc']
    Mdir        = ['2015_060']
    datetimes   = [dtm.datetime(2015,3,2)]
else:
    olist = os.listdir(odir)
    olist.remove(lonlat_file)
    dates   = []
    Olist   = []
    for ofil in olist:
       cdate1   = ofil[-20:-12]
       dt   = dtm.datetime.strptime(cdate1,'%Y%m%d')
       if dt<dtm.datetime(2015,3,2):
           print('Removing '+ofil)
       else:
           Olist.append(ofil)
           dates.append(dt)

    Lst         = sorted([(dt,i) for i,dt in enumerate(dates)])
    olist       = [Olist[tup[1]] for tup in Lst]
    datetimes   = [tup[0] for tup in Lst]
    Mdir        = []
    for i,dt in enumerate(datetimes):
        jday    = '_%3.3i' %(int(dt.strftime('%j'))-1)
        if HEX:
           Mdir.append(dt.strftime('%Y')+'/'+dt.strftime('%Y')+jday+'/binaries/DAILY/')
        else:
           Mdir.append(dt.strftime('%Y')+jday)
        print(dt,olist[i],Mdir[i])
    # sys.exit()



ofil    = olist[-1]
outdir   = 'out/'

# make eg plots
PLOT_EG  = 0
if not os.path.exists(outdir):
   os.mkdir(outdir)

hmax     = 4
po       = FP.plot_object()
shorts   = ['BackThick'              ,'CS2'                   ,'SMOS'             ,'AnThick'              ,'AnThickErr']
clabs    = ['Background Thickness, m','Cryosat-2 thickness, m','SMOS thickness, m','Analysis Thickness, m','Analysis Thickness Error, m']
vlist    = ['background_thickness'   ,'cs2_thickness'         ,'smos_thickness'   ,'analysis_thickness'   ,'analysis_thickness_err']
ts_data  = {'RMSE_both_ice':[],'Bias_both_ice':[],\
            'RMSE_either_ice':[],'Bias_either_ice':[]}

for i,ofil in enumerate(olist):
   # =======================================================================
   # get arrays for later
   print(mdir+'/'+Mdir[i])
   flist = mr.file_list(mdir+'/'+Mdir[i],'DAILY','.a',gridpath=gridpath)
   nci   = mr.nc_getinfo(odir+'/'+ofil,lonlat_file=odir+'/'+lonlat_file)
   # nci.plot_var('analysis_thickness',clim=[0,hmax])
   DOCB  = False
   if i==0:
       olon,olat    = nci.get_lonlat()
       oX,oY        = ODLmap(olon,olat)
       mlon,mlat    = flist.get_lonlat()
       mX,mY        = ODLmap(mlon,mlat)
       DOCB       = True
   hobs = nci.get_var('analysis_thickness')
   # =======================================================================


   # =======================================================================
   if PLOT_EG:
      # plot obs
      if not os.path.exists(outdir):
         os.mkdir(outdir)

      for i,vhobs in enumerate(vlist):
         if 'smos' in vhobs:
             clim   = [0,.5]
         elif 'err' in vhobs:
            clim   = None
         else:
            clim   = [0,hmax]
         po,bmap    = nci.make_png(vhobs,clabel=clabs[i],clim=clim,figdir=outdir,pobj=po)
         po.ax.cla()
         #po.fig.clear()
      sys.exit()
   # =======================================================================

   cdate1   = ofil[-20:-12]
   cdate2   = ofil[-11:-3]

   print('Calculating weekly average')
   print(cdate1+' - '+cdate2+'\n')

   # time average
   hav_m    = flist.time_average('hice',\
                start_date=cdate1+'T120000Z',\
                end_date=cdate2+'T120000Z',\
                inner_points=False)

   # spatial smoothing
   if 0:
      radius   = 20.e3
      hav_m    = mr.smooth(mlon,mlat,hav_m.data,radius,mask=hav_m.mask)
      figname  = 'TP4weekly_'+cdate1+'_smoothed.png'
   else:
      figname  = 'TP4weekly_'+cdate1+'.png'

   # interpolate onto observation grid
   hav_m2   = mr.reproj_mod2obs(mX,mY,hav_m,oX,oY,mask=1*hobs.values.mask)
   hdiff    = hav_m2-hobs.values
   good     = np.logical_not(hdiff.mask)

   # 1st estimate (lower bound)
   # - only consider pixels with ice in both datasets
   both_ice       = np.copy(good)
   both_ice[good] = np.logical_and(hav_m2[good]>0.,hobs[good]>0.)
   RMSEb          = np.sqrt(np.mean(hdiff[both_ice]**2))
   BIASb          = np.mean(hdiff[both_ice])
   ts_data['RMSE_both_ice'].append(RMSEb)
   ts_data['Bias_both_ice'].append(BIASb)

   # 2nd estimate (upper bound)
   # - consider pixels with ice in one of the datasets
   either_ice        = np.copy(good)
   either_ice[good]  = np.logical_or(hav_m2[good]>0.,hobs[good]>0.)
   RMSEe             = np.sqrt(np.mean(hdiff[either_ice]**2))
   BIASe             = np.mean(hdiff[either_ice])
   ts_data['RMSE_either_ice'].append(RMSEe)
   ts_data['Bias_either_ice'].append(BIASe)

   np.savez(outdir+'/thickness_anomaly_'+cdate1+'.npz',lon=olon,lat=olat,\
           anomaly=hdiff.data,mask=hdiff.mask)

   print('\n-------------------')
   print('RMSE (m) = '+str(RMSEe))
   print('Bias (m) = '+str(BIASe))
   print('-------------------\n')

   # mask water for plotting
   hdiff.mask  = np.logical_not(either_ice)
   if 0:
       print('Plotting weekly average')
       PC = bmap.pcolor(mlon,mlat,hav_m,ax=po.ax,vmin=0,vmax=hmax,latlon=True)
       if DOCB:
          po.fig.colorbar(PC)
       FP.finish_map(bmap,ax=po.ax)
       # po.ax.set_title('Mod')
   elif 0:
       print('Plotting weekly average - interpolated')
       PC = bmap.pcolor(olon,olat,hav_m2,ax=po.ax,vmin=0,vmax=hmax,latlon=True)
       if DOCB:
          po.fig.colorbar(PC)
       FP.finish_map(bmap,ax=po.ax)
       figname  = figname.replace('.png','_interp.png')
   else:
       print('Plotting weekly average - anomaly')
       PC = bmap.pcolor(olon,olat,hdiff,ax=po.ax,latlon=True,vmin=-1.5,vmax=1.5)
       #
       if DOCB:
          cbar = po.fig.colorbar(PC)
          cbar.set_label('Thickness anomaly, m',rotation=270,labelpad=20,fontsize=16)
       #
       if 0:
           levels  = [-2,-1.5,-1,-.5,0,.5,1,1.5,2]
           CS   = bmap.contour(olon,olat,hdiff,levels=levels,ax=po.ax,latlon=True)
           po.ax.clabel(CS, fontsize=10, inline=1,fmt = '%1.1f',ticks=levels)
       # po.ax.set_title(cdate1+' - '+cdate2)
       
       FP.finish_map(bmap,ax=po.ax)
       figname  = figname.replace('.png','_anomaly.png')

   print('Saving '+outdir+'/'+figname)
   po.fig.savefig(outdir+'/'+figname)
   po.ax.cla()
   tso  = mr.time_series(datetimes[:i+1],ts_data,\
           filename=outdir+'/thickness_error.txt',overwrite=True)
   # sys.exit()
plt.close(po.fig)
