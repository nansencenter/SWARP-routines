from matplotlib import pyplot as plt
from matplotlib import gridspec as GrdSpc
import mod_reading as mr
import numpy as np

nci   = mr.nc_getinfo('SWARPwavesice_WW3_forecast_start20151216T000000Z.nc')
fig   = plt.figure(figsize=(20,20))
gs    = GrdSpc.GridSpec(1,2)
gs.update(left=0.075, right=0.9, wspace=0.25)

Hs_max   = 6
dH       = 1
Dmax_max = 250
dD       = 50
reg      = 'FR1' # 'gre'
dtl      = 2

for idx in range(nci.number_of_time_records):

   dto   = nci.datetimes[idx]
   cdt1  = dto.strftime('%Y%m%dT%H%M%SZ')
   cdt2  = dto.strftime('%d %b %Y %H:%M')
   if 0:
      mon   = dto.strftime('%b')
      cdt2  = cdt2.replace(mon,mon.upper())

   if dtl==0:
      fig.text(.5,.77,cdt2,\
            horizontalalignment='center',\
            fontsize=18,\
            color='k')

   ax1      = fig.add_subplot(gs[0,0])
   po1      = mr.plot_object(fig=fig,ax=ax1)
   po1,bmap = nci.plot_var('Hs'  ,time_index=idx,pobj=po1,\
               clabel='Significant wave height, m',\
               clim=[0,Hs_max],\
               smoothing=0,\
               date_label=dtl,date_color='r',\
               HYCOMreg=reg,show=False)
   po1.cbar.set_ticks(np.arange(0,Hs_max+dH,dH)) 

   ax2      = fig.add_subplot(gs[0,1])
   po2      = mr.plot_object(fig=fig,ax=ax2)
   po2,bmap = nci.plot_var('Dmax',time_index=idx,pobj=po2,\
               clabel='Maximum floe size, m',\
               clim=[0,Dmax_max],\
               smoothing=0,\
               date_label=dtl,date_color='b',\
               HYCOMreg=reg,show=False)
   po2.cbar.set_ticks(np.arange(0,Dmax_max+dD,dD))

   figname  = 'Hs_Dmax_'+cdt1+'.png'
   print('\nSaving '+figname+'\n')
   fig.savefig(figname,bbox_inches='tight',dpi=100)# ,pad_inches=1,bbox_inches=None)
   # plt.show(fig)

   ax1.cla()
   ax2.cla()
   fig.clear()

