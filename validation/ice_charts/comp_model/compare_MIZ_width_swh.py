import matplotlib
matplotlib.use('Agg')

import MIZchar as mc
import mod_reading as mr
import matplotlib.pyplot as plt
import numpy as np
import geometry_sphere as GS
import os,sys
import datetime
import geometry_planar as GP


# ============================================================================= 
def average_width_model(cdate,fli,outdir='.'):
   # average MIZ width for a given date
   # use binary files since more likely to be there
   # returns None if nothing there


   # loop over 6-h intervals:
   daily_areas    = []
   daily_widths   = []
   for hr in range(0,24,6):
      ss    = cdate+(" %2.2i" %(hr))
      dto   = datetime.datetime.strptime(ss,"%Y%m%d %H")

      #print(dto)
      if dto in fli.datetimes:
         idx   = fli.datetimes.index(dto)
         mil   = fli.get_solutions(time_index=idx)
         ofil0 = fli.file_list[idx].split('/')[-1]
         ofil  = ofil0.replace('.txt','_summary.txt')
         out   = mil.save_summary(outdir+'/'+ofil)

         if out is not None:
            summ,long_names   = out
            daily_areas.append(summ['tot_area'])
            daily_widths.append(summ['int_width_mean'])


   if len(daily_areas)==0:
      return
   else:
      daily_areas    = np.array(daily_areas)
      daily_widths   = np.array(daily_widths)

      # take daily average
      A_av  = np.mean(daily_areas)
      W_av  = np.sum(daily_areas*daily_widths)/np.sum(daily_areas)

      print('\nAvg MIZ width: '+str(W_av/1.e3)+'\n')
      return W_av,A_av
# ============================================================================= 


# ============================================================================= 
def get_width_chart(dto,fli,vertices,outdir='.'):
   # MIZ area for a given date

   idx   = fli.datetimes.index(dto)
   mil   = fli.get_solutions(time_index=idx,vertices=vertices)

   ofil0 = fli.file_list[idx].split('/')[-1]
   ofil  = ofil0.replace('.txt','_summary.txt')
   out   = mil.save_summary(outdir+'/'+ofil)
   if out is None:
      return
   else:
      summ,long_names   = out
      A  = summ['tot_area']
      W  = summ['int_width_mean']
      return W,A
# ============================================================================= 
   
   

# for waves: scan \pm4 records
before         = -4
after          = 4
starting_date  = datetime.datetime(2015,5,7)
cyear          = starting_date.strftime("%Y")

print("starting date:")
print(starting_date)
print("\n")


# =====================================================================================
# set paths
wavefile = '/work/shared/nersc/msc/WAVES_INPUT/ERA-I_1deg/SWH_'+cyear+'.nc'          # yearly file of ERA-I waves
chartdir = "/work/shared/nersc/msc/DMI_icecharts/MIZ_FA/MIZ_polys_classified/"+cyear # path to ice chart files (FA criteria)
FCdir    = "out/model_polygons"                                                      # where forecast results are untarred to
outdir_m = "out/model_polygons_summaries"                                            # where to save polygons from model
outdir_c = "out/chart_polygons_summaries"                                            # where to save polygons from model
outdir   = "out"                                                                     # where to save polygons from model
for odir in [outdir,outdir_m,outdir_c]:
   if not os.path.exists(odir):
      os.mkdir(odir)
# =====================================================================================

TSfil = outdir+'/test_width.txt'

# =====================================================================================
# computations
if 1:
   w     = open(TSfil,'w')

   ADD_WAVES   = 0
   if ADD_WAVES:
      w.write('date\tmodel_width\tchart_width\t max_swh\n')
   else:
      w.write('date\tmodel_width\tchart_width\n')
   w.close()

   vertices    = [(-50,50),(20,50),(20,84),(-50,84)] # region of analysis (gre)
   dates       = []
   flist       = os.listdir(chartdir)
   nci         = mr.nc_getinfo(wavefile)
   wlon,wlat   = nci.get_lonlat() #lon,lat of wave grid
   mask        = GP.maskgrid_outside_polygon(wlon,wlat,vertices) # if True, these are wave points inside vertices

   chart_list  = mr.polygon_file_list(chartdir,prefix='',\
                                      suffix='_Greenland_WA_MIZpolys.txt',\
                                      date_format='%Y%m%d%H%M')
   model_list  = mr.polygon_file_list(FCdir,prefix='TP4archv_wav.',\
                                      suffix='_dmax.txt',
                                      date_format="%Y_%j_%H%M%S",add_day=True)


   for day in chart_list.datetimes:
           if day < starting_date: 
                   continue
           print("\n")
           print day

           # model area
           out = average_width_model(day.strftime("%Y%m%d"),model_list,outdir=outdir_m)
           if out is not None:
              model_width,model_area  = out
           else:
              print('skipping date '+day.strftime("%Y%m%d")+' no model output')
              # sys.exit()
              continue

           # ========================================================================================================
           # chart area
           chart_width,chart_area   = get_width_chart(day,chart_list,vertices,outdir=outdir_c)
           # ========================================================================================================


           if ADD_WAVES:
              # ========================================================================================================
              # max swh
              DT,idx=nci.nearestDate(day)
              swh_max=[]
              #w.write(day.strftime("%Y%m%d")+'\t'+str(tot_area)+'\t')
              
              for i in range(before,after):
                      
                      print swh_max
                      swh=nci.get_var('swh',time_index=idx+i)
                      good=np.logical_not(swh.values.mask) # swh.values is masked array
                      good=np.logical_and(good,mask)       # restricts to area from vertices
                      #swh_max=max(swh.values.data[good].max(),swh_max)	
                      swh_max=swh.values.data[good].max()
                      #w.write(str(swh_max)+'\t')
              #w.write('\n')

              ss  = day.strftime("%Y%m%d")+'\t'+str(model_width)+'\t'+str(chart_width)+'\t'+str(swh_max)+'\n'
              # ========================================================================================================
           else:
              ss  = day.strftime("%Y%m%d")+'\t'+str(model_width)+'\t'+str(chart_width)+'\n'
              
              
           print('date/model/chart')
           print(ss)
           w  = open(TSfil,'a')
           w.write(ss)
           w.close()

           continue
           # break
# sys.exit()
# =====================================================================================


# =====================================================================================
# load time series file & plot
print('Time series in '+TSfil)

def convert_datetime(dto,info):
   if info['time_units']=='days':
      xfac  = 24*3600. # seconds in 1 day
   elif info['time_units']=='hours':
      xfac  = 3600. # seconds in 1h
   elif info['time_units']=='minutes':
      xfac  = 60. # seconds in 1min
   else:
      xfac  = 1. # seconds in 1min

   return (dto-info['refdate']).total_seconds()/xfac

ts       = mr.read_time_series(TSfil)
yscaling = 1.e-3 # m -> km

lines       = []
po,lin,info = ts.plot('model_width',time_units='days',yscaling=yscaling)
lines.append(lin)
po,lin,info = ts.plot('chart_width',time_units='days',yscaling=yscaling,pobj=po,color='r')
lines.append(lin)



# po.ax.set_xlabel('Date')
po.ax.set_ylabel('MIZ width, km')
po.ax.legend(lines,['Model','Chart'])

# xticks
xt = []
XT = []
for mon in range(1,13):
   if mon>starting_date.month:
      dto   = datetime.datetime(int(cyear),mon,1)
      xt.append(convert_datetime(dto,info))
      XT.append(dto.strftime('%D'))

po.ax.set_xticks(xt)#,rotation='vertical')
po.ax.set_xticklabels(XT,rotation='vertical')

TSfig = TSfil.replace('.txt','.png')
print('Saved time series plot to '+TSfig)
po.fig.savefig(TSfig,bbox_inches='tight')
# =====================================================================================
