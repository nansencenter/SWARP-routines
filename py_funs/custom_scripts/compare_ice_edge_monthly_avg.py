import mod_reading as mr
import fns_plotting as Fplt
import mod_HYCOM_utils as MHU
import datetime as DTM
from matplotlib import pyplot as plt
from matplotlib import gridspec as GrdSpc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os,sys
import numpy as np

PLOT_MOD = False
PLOT_OBS = False
year     = 2006

# model data
print('Load model daily files...\n')
ddir  = '/work/timill/RealTime_Models/TP4a0.12/expt_01.4/data'
fli1  = mr.file_list(ddir,'DAILY','.a')

# monthly averaged conc
print('Load OSISAF monthly files...\n')
osidir   = '/work/shared/nersc/OSISAF/'
osipat   = 'osisaf-nh_aggregated_ice_concentration_nh_polstere-100'
fli2     = mr.file_list(osidir,osipat,'.nc')

# basemap
gs    = GrdSpc.GridSpec(1,1)
bmap  = Fplt.start_HYCOM_map('Arctic')
fig   = plt.figure()
ax    = fig.add_subplot(gs[0,0])

dto0        = DTM.datetime(year,1,1)
dto0_,idx0  = fli2.nearestDate(dto0)
# print(idx0)

hgi         = MHU.HYCOM_grid_info(gridpath=ddir)
cell_areas  = hgi.get_areas()

ice_areas   = []
ice_vols    = []
ice_vols2   = []
dates       = []

HAVE_LONLAT = False
HAVE_CBAR   = False
xyann       = (0.05,.925)
for idx in range(idx0,idx0+12):
   start_date  = fli2.datetimes[idx]
   end_date    = fli2.datetimes[idx+1]-DTM.timedelta(1)
   print('Start date:')
   print(start_date)
   print('End date:')
   print(end_date)

   if not HAVE_LONLAT:
      mlon,mlat   = fli1.get_lonlat()
      olon,olat   = fli2.get_lonlat()
      HAVE_LONLAT = True
      print('\nHave lon/lat\n')

   if 0:
      # quick test to check the plotting
      fice_mod = fli1.get_var('fice',time_index=0).values   #masked array
   else:
      print("Getting model's monthly average of conc")
      fice_mod = fli1.time_average('fice',start_date=start_date,end_date=end_date)#masked array
      print("Getting model's monthly average of thickness")
      hice_mod = fli1.time_average('hice',start_date=start_date,end_date=end_date)#masked array
      print("Getting model's monthly average of volume")
      vol_mod  = fli1.time_average('Volume',start_date=start_date,end_date=end_date)#masked array
      print(fice_mod.mean())


   # save time series of ice area
   ice_areas.append(np.sum(cell_areas*fice_mod))
   ice_vols.append(np.sum(cell_areas*vol_mod))
   ice_vols2.append(np.sum(cell_areas*fice_mod*hice_mod))
   dates.append(start_date)

   # convert OSISAF % to fraction
   fice_obs    = .01*fli2.get_var('ice_conc_avg',time_index=idx).values   #masked array


   if PLOT_MOD:
      # ================================================================================
      # plot model monthly av + OSISAF ice edge
      PC = bmap.pcolor(mlon,mlat,fice_mod,ax=ax,latlon=True,vmin=0,vmax=1)
      if not HAVE_CBAR:
         divider     = make_axes_locatable(ax)
         cax         = divider.append_axes("right", size="5%", pad=0.30)
         cbar        = fig.colorbar(PC,cax=cax)
         cbar.set_label("Sea ice concentration",rotation=270,labelpad=20,fontsize=16)
         HAVE_CBAR   = True

      bmap.contour(olon,olat,fice_obs,levels=[.15],colors='g',linewidths=2,ax=ax,latlon=True)
      Fplt.finish_map(bmap,ax=ax)

      # annotation
      tlabel   = start_date.strftime('%b %Y')
      ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18,color='r')

      figname  = 'out/model_monthly_avg_'+start_date.strftime('%Y%m')+'.png'
      print('Saving '+figname)

      fig.savefig(figname)
      ax.cla()
      # ================================================================================


   if PLOT_OBS:
      # ================================================================================
      # plot OSISAF monthly av
      PC = bmap.pcolor(olon,olat,fice_obs,ax=ax,latlon=True,vmin=0,vmax=1)
      Fplt.finish_map(bmap,ax=ax)

      # annotation
      tlabel   = start_date.strftime('%b %Y')
      ax.annotate(tlabel,xy=xyann,xycoords='axes fraction',fontsize=18,color='r')

      figname  = 'out/OSISAF_monthly_avg_'+start_date.strftime('%Y%m')+'.png'
      print('Saving '+figname+'\n')

      fig.savefig(figname)
      ax.cla()
      # ================================================================================


   # sys.exit()

plt.close(fig)

data  = {'ice_area':ice_areas,
         'ice_volume':ice_vols,
         'ice_volume2':ice_vols2}
TS    = mr.time_series(dates,data,filename='out/time_series_'+str(year)+'.txt')

out   = TS.plot('ice_area')
plt.show(out[0].fig)
