import os,sys
from datetime     import datetime   as dtm
from datetime     import timedelta  as tdl
from matplotlib   import pyplot     as plt
import mod_reading                  as Mr

def parse_src(fil):
   f  = open(fil)
   L  = f.readlines()
   f.close()
   print(f)
   print(L)
   #
   lut   = {}
   for lin in L:
      ss = lin.split('=')
      print(ss)
      lut.update({ss[0]:ss[1]})
      print(lut)
   return lut

class forecast_info:

   def __init__(self,datestr=None,fctype='ice_only'):

      fmt   = '%Y%m%d'
      if datestr is not None:
         self.startdate = dtm.strptime(datestr,fmt)
      else:
         self.startdate = dtm.today()
         datestr        = self.startdate.strftime(fmt)

      SR    = os.getenv('SWARP_ROUTINES')
      tfil  = SR+'/forecast_scripts/'
      if fctype=='ice_only':
         nc0   = 'SWARPiceonly_forecast'
      elif fctype=='wavesice':
         nc0   = 'SWARPwavesice_forecast'
      elif fctype=='wavesice_ww3arctic':
         nc0   = 'SWARPwavesiceWW3_forecast'
      else:
         raise ValueError('Invalid value of fctype: '+fctype)

      ncbase   = nc0+'_start'+datestr+'T000000Z'
      fcdir    = '/work/timill/RealTime_Models/results/TP4a0.12/'+fctype+'/'+datestr+'/final_output'
      ncfil    = fcdir+'/'+ncbase+'.nc'

      print('\nOpening '+ncfil+' ...\n')
      nci      = Mr.nc_getinfo(ncfil)
      #
      self.netcdf_info  = nci
      self.basename     = ncbase
      self.datetimes    = []

      for tval in nci.timevalues:
         dtmo  = nci.timeval_to_datetime(tval)
         self.datetimes.append(dtmo)

      return

   def fc2png(self,dtmo=None,vname='icec',figdir='.',date_label=True):

      dts   = self.datetimes
      nci   = self.netcdf_info

      if dtmo is None:
         # last day
         idx   = len(dts)-1
         dtmo  = dts[idx]
      if dtmo in self.datetimes:
         idx   = dts.index(dtmo)

      fig,ax,bmap = nci.plot_var(vname,HYCOMreg='TP4',time_index=idx,show=False)

      if date_label:
         tlabel   = dtmo.strftime('%d %b %Y %H:%M')
         ax.annotate(tlabel,xy=(0.05,.925),xycoords='axes fraction',fontsize=18)

      datestr  = dtmo.strftime('%Y%m%dT%H%M%SZ')
      outname  = figdir+'/'+self.basename+'_dump'+datestr+'_'+vname'.png'

      print('\nSaving '+outname+' ...\n')
      fig.savefig(outname)
      ax.cla()
      plt.close(fig)

      return

   def fc2png_all(self,figdir='.',**kwargs):
      # plot all snapshots for a given variable in a forecast
      # can make a gif with convert for example

      subdir   = figdir+'/'+self.basename
      if not os.path.exists(subdir):
         os.mkdir(subdir)

      for dtmo in self.datetimes:
         self.fc2png(dtmo=dtmo,figdir=subdir,**kwargs)

      return

