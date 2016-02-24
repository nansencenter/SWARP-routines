import mod_reading as mr
# echo echo ???
#reg='TP4'
#y1='2015'
#d1='347'
#y2='2015'
#d2='351'
#''+reg+'DAILY_'+y1+'_'+d1+'_'+y2+'_'+d2+'.a'
#
# if show figures: 'True' or 'False'
showfigs='False'

# week year
wy='2015'
# first week day
wd='347'
# plot year
yd='2015'
# plot day
dd='351'
# Area of plot TP4, FR1, BS1, gre, bar, ....?
hycomreg='TP4'

dres='/work/users/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice/analysis/binaries_all/DAILY/'
# resdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/'+wy+'_GOOD/'+wy+'_'+yd+'/binaries/'
outdir='/work/timill/RealTime_Models/results_hindcasts/Figures/'
osdir='/work/shared/nersc/msc/OSI-SAF/'+wy+'_nh_polstere/'

afiles=[ ''+dres+'TP4DAILY_2015_347_2015_351.a' ]
#,\
#         'TP4DAILY_2016_010_2016_011.a',\
#         'TP4DAILY_2016_010_2016_011.a',\
#         'TP4DAILY_2016_010_2016_011.a' ]

for afil in afiles:
   vo=mr.make_plot_options('hice',ice_mask=True)
   vo2=mr.make_plot_options('fice',ice_mask=True)
   # fice vs OSI-SAF
   print "afil: ", afil
   hi=mr.HYCOM_binary_info(afil)
   dtm=hi.datetime
   cdate=dtm.strftime('%Y%m%d')
   figname='comp_OSISAF_'+hycomreg+'_'+cdate+'.png'
   hi.compare_ice_edge_obs(figname=figname,show='False',HYCOMreg=hycomreg)
   # only OSI-SAF fice
   osfil=''+osdir+'ice_conc_nh_polstere-100_multi_'+cdate+'1200.nc'
   print "osfil: ", osfil
   nci=mr.nc_getinfo(osfil)
   figname='OSISAF_'+hycomreg+'_'+cdate+'.png'
   #nci.plot_var(vo,show=showfigs,figname=figname)
   nci.make_png(vo2,clim=[0,1],show='False',HYCOMreg=hycomreg)
      
   
   # hi.plot_var(vo,clim=[0,4])
   




