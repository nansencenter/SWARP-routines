import os,sys
import mod_reading as mr
from datetime import datetime as DTM
from datetime import timedelta as TDEL
from matplotlib import pyplot as plt

rootdir  = '/work/timill/DANY_PAUL_tests'

# DMI ice charts
TS = {}
ts = mr.read_time_series(rootdir+'/time_series_DMI_2005_2006.txt') # DMI time-series 2005-2006 merged
TS.update({'DMI':ts})
# print(ts.variables)

SETOPT   = 1
VBLOPT   = 1
if len(sys.argv)>=2:
   SETOPT   = int(sys.argv[1])
if len(sys.argv)>=3:
   VBLOPT   = int(sys.argv[2])

# model results
ddirs = {'M1' :'CONST_REFR_1days_2005_242_730',\
         'M2' :'CONST_REFR_2days_2005_242_730',\
         'M3' :'CONST_REFR_3days_2005_242_730',\
         'M4' :'CONST_REFR_4days_2005_242_730',\
         'M2L':'CONST_REFR_INIT_LARGE_2days_2005_242_730',\
         'M2D':'CONST_REFR_2days_adv_dfloe_2005_242_730/'}

colors = {'DMI':'r',\
          'M1' :'b',\
          'M2' :'k',\
          'M3' :'c',\
          'M4' :'g',\
          'M2L':'c',\
          'M2D':'b'}

lstyles  = {'DMI':'-',\
            'M1' :'-',\
            'M2' :'-',\
            'M3' :'-',\
            'M4' :'-',\
            'M2L':'-',\
            'M2D':'-'}


for ddir in ddirs:
   tfil  = rootdir+'/analysis/'+ddirs[ddir]+'/time_series/time_series_MIZ_dmax_custom.txt'
   ts    = mr.read_time_series(tfil)
   TS.update({ddir:ts})

po       = None
lines    = []
refdate  = TS['DMI'].dates[0]

if VBLOPT==1:
   vbl0     = 'width'
   vbl      = 'Mean_intersecting_width'
   yscaling = 1.e-3 # m -> km
   ylab     = 'MIZ width, km'

elif VBLOPT==2:
   vbl0     = 'area'
   vbl      = 'Total_area'
   yscaling = 1.e-6 # m^2 -> km^2
   ylab     = 'MIZ area, km$^2$'

if SETOPT==1:
   # DMI, 2days + adv dfloe + large init
   labs  = ['DMI','M2','M2L','M2D']
   print(labs)
   TSfig = 'out/DMI_model_methods_'+vbl0+'.png'
elif SETOPT==2:
   # DMI, 2days + adv dfloe + large init
   labs  = ['DMI','M1','M2','M3','M4']
   print(labs)
   TSfig = 'out/DMI_model_time_scale_'+vbl0+'.png'

for lab in labs:
   ts = TS[lab]
   po,lin,info = ts.plot(vbl,\
         refdate=refdate,time_units='days',yscaling=yscaling,\
         pobj=po,color=colors[lab],linestyle=lstyles[lab],linewidth=1.25)
   lines.append(lin)


# po.ax.set_xlabel('Date')
po.ax.set_ylabel(ylab)
po.ax.legend(lines,labs,loc='center left', bbox_to_anchor=(1, 0.5))

# =============================================================
# xticks
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

xt = []
XT = []
yr = 2005
for mon in range(2,25,2):

   if mon<=12:
      dto   = DTM(yr,mon,1)
   else:
      dto   = DTM(yr+1,mon-12,1)

   xt.append(convert_datetime(dto,info))
   XT.append(dto.strftime('%m/%y'))

po.ax.set_xticks(xt)#,rotation='vertical')
po.ax.set_xticklabels(XT,rotation='vertical')
# =============================================================

xmin  = convert_datetime(DTM(2005,1,1),info)
xmax  = convert_datetime(DTM(2007,1,3),info)
po.ax.set_xlim([xmin,xmax])

print('Saved time series plot to '+TSfig)
po.fig.savefig(TSfig,bbox_inches='tight')

# plt.show(po.fig)
