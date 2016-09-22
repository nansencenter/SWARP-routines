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
def average_area_model(cdate,vertices,fcdir0,outdir='.'):
   # average MIZ area for a given date
   # use binary files since more likely to be there
   # returns -100 if nothing there

   fcdir    = fcdir0+'/'+cdate+'/'
   if not os.path.exists(fcdir):
      return
   lst      = os.listdir(fcdir)
   if cdate in lst:
      # check if need to add another cdate
      fcdir+= cdate+'/'

   lst   = os.listdir(fcdir)
   if 'bin' in lst:
      bindir   = fcdir+'bin'
   else:
      bindir   = fcdir+'binaries'


   # make file_list object from binary files
   # - treat in same way as multi-record netcdf
   fli   = mr.file_list(bindir,'archv_wav','.a')
   if fli.number_of_time_records==0:
      return


   # loop over 6-h intervals:
   daily_areas = []
   for hr in range(0,24,6):
      dto   = datetime.datetime(int(cdate[:4]),int(cdate[4:6]),int(cdate[6:8]),hr)

      #print(dto)
      if dto in fli.datetimes:
         idx   = fli.datetimes.index(dto)
         out   = fli.MIZmap(no_width=True,vertices=vertices,time_index=idx,outdir=outdir)
         # out   = fli.MIZmap(vertices=vertices,time_index=idx,outdir=outdir)
         # sys.exit()
         #
         tfil  = out[out.keys()[0]]
         pil   = mc.single_file(tfil)

         tot_area=0
         for pio in pil.poly_info_objects:
            #pio.area		# approximate area (Euclidean after projection using NP as center)
            #pio.ll_coords #list of coords of boundaries
            lon,lat=np.array(pio.ll_coords).transpose()
            area=GS.area_polygon_ellipsoid(lon,lat)
            # print(area)
            tot_area+=area

         print('\nTot area: '+str(tot_area)+'\n')
         daily_areas.append(tot_area)

   if len(daily_areas)==0:
      return
   else:
      # take daily average
      TotArea  = np.mean(daily_areas)
      print('\nAvg tot area: '+str(TotArea)+'\n')
      return TotArea
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
FCdir    = "/work/timill/old_forecasts/TP4wamnsea/"+cyear                            # where forecast results are untarred to
outdir_m = "out/model_polygons"                                                          # where to save polygons from model
outdir   = "out"                                                                     # where to save polygons from model
for odir in [outdir,outdir_m]:
   if not os.path.exists(odir):
      os.mkdir(odir)
# =====================================================================================

TSfil = outdir+'/test_area.txt'
w     = open(TSfil,'w')
w.write('date\tmodel_area\tchart_area\t max_swh\n')

vertices    = [(-50,50),(20,50),(20,84),(-50,84)] # region of analysis (gre)
dates       = []
flist       = os.listdir(chartdir)
nci         = mr.nc_getinfo(wavefile)
wlon,wlat   = nci.get_lonlat() #lon,lat of wave grid
mask        = GP.maskgrid_outside_polygon(wlon,wlat,vertices) # if True, these are wave points inside vertices


for tf in flist:
	# print tf
#	if (tf != '201505101200_Greenland_WA_MIZpolys.txt'):
#		continue
	cdate=tf[0:8]
	t=datetime.datetime(int(cdate[0:4]),int(cdate[4:6]),int(cdate[6:8]))
	dates+=[t]
dates.sort()

for day in dates:
	if day < starting_date: 
		continue
	print("\n")
	print day

        # model area
	model_area = average_area_model(day.strftime("%Y%m%d"),vertices,FCdir,outdir=outdir_m)
	if model_area is None:
	   print('skipping date '+day.strftime("%Y%m%d")+' no model output')
           # sys.exit()
           continue

        # ========================================================================================================
        # chart area
	pil=mc.single_file(chartdir+'/'+day.strftime("%Y%m%d")+'1200_Greenland_WA_MIZpolys.txt')
	#po=pil.plot_all(show=False,latlon=True)
	#plt.savefig('testfig.png')
	
	#po.fig.show()

	pil_reduced=pil.reduce_area(vertices)
	#- closing polygon
	#po0=pil_reduced.plot_all(show=False,latlon=True)
	#po0.fig.show()
	#po0=pil_reduced.plot_all(show=False,latlon=True,check_flags=True)
	
	lon,lat=np.array(vertices).transpose()
	print('rectangle area\t'+str(GS.area_polygon_ellipsoid(lon,lat)))
	
	chart_area=0
	for pio in pil_reduced.poly_info_objects:
		#pio.area		# approximate area (Euclidean after projection using NP as center)
		#pio.ll_coords #list of coords of boundaries
		lon,lat=np.array(pio.ll_coords).transpose()
		area=GS.area_polygon_ellipsoid(lon,lat)
		print(area)
		chart_area+=area
        # ========================================================================================================


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
        # ========================================================================================================
	
	
        ss  = day.strftime("%Y%m%d")+'\t'+str(model_area)+'\t'+str(chart_area)+'\t'+str(swh_max)+'\n'
        print(ss)
        w.write(ss)

	continue	
	
	
w.close()


# load time series file & plot
print('Saved time series to '+TSfil)

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
yscaling = 1.e-6 # m^2 -> km^2

lines       = []
po,lin,info = ts.plot('model_area',time_units='days',yscaling=yscaling)
lines.append(lin)
po,lin,info = ts.plot('chart_area',time_units='days',yscaling=yscaling,pobj=po)
lines.append(lin)



# po.ax.set_xlabel('Date')
po.ax.set_ylabel('MIZ area, km$^2$')
po.ax.legend(lines,['Model','Chart'])

# xticks
xt = []
XT = []
for mon in range(1,13):
   if mon>starting_date.month:
      dto   = datetime.datetime(int(cyear),mon,1)
      xt.append(convert_datetime(dto,info))
      XT.append(dto.strftime('%d %b'))

po.ax.set_xticks(xt)#,rotation='vertical')
po.ax.set_xticklabels(XT,rotation='vertical')

TSfig = TSfil.replace('.txt','.png')
print('Saved time series to '+TSfig)
po.fig.savefig(TSfig,bbox_inches='tight')
