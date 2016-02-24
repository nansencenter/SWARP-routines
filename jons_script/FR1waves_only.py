import matplotlib
matplotlib.use('Agg')
import os
import mod_reading as mr

# ddir     = '/work/timill/RealTime_Models/FR1a0.03/expt_01.4/data/' # dir with files of interest
# ddir = '/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/wavesice/analysis/binaries_all/archv_wav/'
ddir = '/work/timill/RealTime_Models/results_hindcasts/Figures/MIZ/'
outdir   = ddir
# '/work/timill/RealTime_Models/results_hindcasts/Figures/MIZ/'       # where to save
if not os.path.exists(outdir):
   os.mkdir(outdir)

print('Sorting files...')
afiles   = []
for fil in os.listdir(ddir):
   if fil[-2:]=='.a' and 'TP4' in fil:
      afiles.append(ddir+fil)

for afil in afiles:
   print('Opening '+afil)
   hbi   = mr.HYCOM_binary_info(afil)

   #########################################################
   if 1:
      # usual way to get time
      cdt   = hbi.datetime.strftime('%Y%m%dT%H%M%SZ')
   else:
      # waves-only run has incorrect 'model day'
      # -read it manually
      bfil  = open(hbi.bfile)
      lins  = bfil.readlines()
      bfil.close()

      # print(afil)
      # print(lins)
      cdate = lins[13].split()[0]
      ctime = lins[14].split()[0]
      cdt   = cdate+'T'+ctime+'Z'
   #########################################################

   odir  = outdir+cdt
   if not os.path.exists(odir):
      os.mkdir(odir)
      print('Making MIZ map...')
      hbi.MIZmap(do_sort=True,plotting=True,outdir=odir)
