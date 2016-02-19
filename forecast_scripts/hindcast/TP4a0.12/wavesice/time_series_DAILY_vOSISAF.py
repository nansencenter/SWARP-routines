import mod_reading as mr
from getopt import getopt
import os,sys
from datetime import datetime as Dtm
import numpy as np
from matplotlib import pyplot as plt

# ==========================================================================
# options
indir    = None
outdir   = '.'

opts,args   = getopt(sys.argv[1:],"",["indir=","outdir="])
for opt,arg in opts:
   if opt=='--indir':
      indir = arg
   if opt=='--outdir':
      outdir   = arg

if indir is None:
   raise ValueError('specify indir with --indir="..."')
# ==========================================================================

def dblank(keys):
   out   = {}
   for key in keys:
      out.update({key:[]})
   return out

dates = os.listdir(indir) # directories corresponding to different dates
fmt   = '%Y%m%d'
tref  = Dtm(2000,1,1)

tids     = []
regions  = ['gre','bar','lab']
cats     = ['Over','Under']
widths   = {}
perims   = {}
for cat in cats:
   widths.update({cat:dblank(regions)})
   perims.update({cat:dblank(regions)})

print('Reading...')
N  = 0
for date in dates:
   N     = N+1
   dtm   = Dtm.strptime(date,fmt)
   dt    = dtm-tref
   tids.append(dt.total_seconds()/float(24*3600))
   # print(tids)
   # print(N,date)

   ddir  = indir+'/'+date+'/vOSISAF/'
   Lst   = os.listdir(ddir)
   for cat in cats:
      for reg in regions:
         ss = cat+'_'+reg+'_summary.txt'
         P  = 0
         W  = 0
         # print(N,date,cat,reg)
         for fil in Lst:
            if ss in fil:
               mps   = mr.read_MIZpoly_summary(ddir+fil)
               # print(mps.filename)
               P     = mps.Total_perimeter/1.e3
               W     = mps.Mean_intersecting_width/1.e3
               break
         perims[cat][reg].append(P)
         widths[cat][reg].append(W)
         # print(perims[cat][reg])
         # print(widths[cat][reg])

tids        = sorted([(e,i) for i,e in enumerate(tids)])
tids,isort  = np.array(tids).transpose()
isort       = np.array(isort,dtype='int')

print('Plotting...')
yr          = 2015
xticks      = []
xtick_labs  = []
for mon in range(3,14):
   if mon>12:
      dtm   = Dtm(yr+1,mon-12,1)
   else:
      dtm   = Dtm(yr,mon,1)
   dt    = dtm-tref
   xticks.append(dt.total_seconds()/float(24*3600))
   Mon   = dtm.strftime('%B')
   # xtick_labs.append('1 '+Mon[:4])
   xtick_labs.append(Mon[0])


sgn   = {'Over':-1,'Under':1}
for reg in regions:
   fig   = plt.figure()
   ax    = fig.add_subplot(1,1,1)
   Ptot  = 0
   Wavg  = 0
   lins  = []
   legs  = []
   for cat in cats:
      P     = np.array(perims[cat][reg])           [isort]
      W     = sgn[cat]*np.array(widths[cat][reg])  [isort]
      lin  ,= ax.plot(tids,W)
      #
      lins.append(lin)
      legs.append(cat)
      #
      Ptot  = Ptot+P
      Wavg  = Wavg+P*W

   Wavg  = Wavg/Ptot
   lin  ,= ax.plot(tids,Wavg)
   lins.append(lin)
   legs.append('Bias')
   ax.set_ylabel('Distance to ice edge, km')
   ax.legend(lins,legs,loc='upper left')

   if 0:
      ax.set_xlabel('Time, days since '+dates[0])
   else:
      ax.xaxis.set_ticks(xticks)
      ax.xaxis.set_ticklabels(xtick_labs)
      ax.set_xlabel('Time, starting 1 March 2015')

   figname  = outdir+'/time_series_DAILY_vOSISAF_'+reg+'.png'
   print('Saving to '+figname)
   fig.savefig(figname)
   if 0:
      plt.show(fig)
   else:
      ax.cla()
      fig.clf()

plt.close(fig)
