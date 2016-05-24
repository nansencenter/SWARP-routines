import os,sys
import mod_reading as MR
import time
import datetime as DTM

# get yesterday's date
fmt   = '%Y%m%d'
tday  = time.strftime(fmt)
dtm0  = DTM.datetime.strptime(tday,fmt)
dtm   = dtm0-DTM.timedelta(1)
yday  = dtm.strftime(fmt)
idx   = 0

if 1:
   # scalar plot of concentration
   vobj  = 'fice'
else:
   # plot ice speed with unit vectors showing direction
   vobj = MR.make_plot_options('uice',\
      vec_opt=2,layer=0,conv_fac=1,wave_mask=False,ice_mask=True,dir_from=True)

FCdir = os.getenv('TP4_REALTIME_RES')+'/ice_only/'+yday
if 0:
   # binary file
   jday  = int(dtm.strftime('%j'))-1  #
   dts   = '%4.4i_%3.3i' %(dtm.year,jday)
   ddir  = FCdir+'/bin/'
   afil  = ddir+'/TP4archv.'+dts+'_120000Z.a'
   print(afil)
   fobj  = MR.HYCOM_binary_info(afil)
elif 0:
   # single-record netcdf
   ddir  = FCdir+'/netcdf/'
   dts   = yday+'_%2.2i' %(dtm.hour)
   ncfil = ddir+'/TP4archv_'+dts+'.nc'
   print(ncfil)
   fobj  = MR.nc_getinfo(ncfil)
   vobj  = 'fice'
elif 1:
   # multi-record netcdf
   ddir  = FCdir+'/final_output/'
   dts   = yday+'T000000Z'
   ncfil = ddir+'/SWARPiceonly_forecast_start'+dts+'.nc'
   print(ncfil)
   fobj  = MR.nc_getinfo(ncfil)
   vobj  = 'fice'
   idx   = 2
else:
   # file_object_list
   pattern  = 'archv'
   if 1:
      # list of binary files
      ddir  = FCdir+'/bin'
      ext   = '.a'
   elif 1:
      # list of netcdf files
      ddir     = FCdir+'/netcdf'
      ext      = '.nc'

   fobj  = MR.file_list(ddir,pattern,ext)
   vobj  = 'fice'

# 2nd variable to plot: ice velocity as a quiver plot
vobj2 = MR.make_plot_options('uice',\
   vec_opt=3,layer=0,conv_fac=1,wave_mask=False,ice_mask=True,dir_from=True)


   # =========================================================================
if 0:
   # use imshow for a fast plot
   fobj.imshow(vobj,time_index=idx,date_label=2)
elif 0:
   # use plot_var for a projected plot
   fobj.plot_var(vobj,time_index=idx,date_label=2)
elif 0:
   # plot pair of variables
   fobj.plot_var_pair(vobj,vobj2,time_index=idx,date_label=2)
elif 0:
   # plot variable & save to png
   fobj.make_png(vobj,figdir='out',date_label=2,time_index=idx)
elif 0:
   # plot pair of variables & save to png
   fobj.make_png_pair(vobj,vobj2,figdir='out2',date_label=2,time_index=idx)
elif 0:
   # plot fice & OSISAF ice edge
   # save to png if figname is given
   fobj.compare_ice_edge_obs(date_label=2,time_index=idx,figname=fobj.basename+'.OSISAF.png')
   # =========================================================================


   # =========================================================================
elif 0:
   # MIZ width diagnostic
   fobj.MIZmap(var_name='fice',do_sort=True,time_index=idx,outdir='out')
elif 0:
   # distance to ice edge diagnostic
   fobj.areas_of_disagreement(outdir='out',time_index=idx)
   # =========================================================================


   # =========================================================================
elif 0:
   # make png's of 1 variable for all time records
   fobj.make_png_all(vobj,date_label=2,figdir='out')

elif 1:
   # make png's of 2 variables for all time records
   fobj.make_png_pair_all(vobj,vobj2,date_label=2,figdir='out')
   # =========================================================================

# compare_ice_edge_obs_all(
