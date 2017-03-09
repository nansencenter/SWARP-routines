import os,sys
import mod_reading as mr
from datetime import datetime as DTM
from datetime import timedelta as TDEL

rootdir  = '/work/timill/DANY_PAUL_tests/data'

ddirs = {'M1' :'CONST_REFR_1days_2005_242_730',\
         'M2' :'CONST_REFR_2days_2005_242_730',\
         'M3' :'CONST_REFR_3days_2005_242_730',\
         'M4' :'CONST_REFR_4days_2005_242_730',\
         'M2L':'CONST_REFR_INIT_LARGE_2days_2005_242_730',\
         'M2D':'CONST_REFR_2days_adv_dfloe_2005_242_730/'}

me    = sys.argv[0]
opts  = sys.argv[1:]
if len(opts)==1:
   ddir  = rootdir+'/'+ddirs[opts[0]]
else:
   print('Usage: '+me+' run_name')
   print('options for run_name')

   for s in ddirs:
      d_dir = rootdir+'/'+ddir[s]
      print(s+' = '+d_dir)

   print('\n')
   raise ValueError('\n\nNot enough inputs to '+me+'\n')

# ddir  = sys.argv[1]
print('sort files in '+ddir+'...')
fli   = mr.file_list(ddir,'TP4archv_wav','.a')
# print(fli.datetimes)

if 0:
   # just do a few times for testing
   start_date  = DTM.strptime('20050831','%Y%m%d')# +TDEL(0.5)
   end_date    = DTM.strptime('20050901','%Y%m%d')
   plotting    = True
   show        = False
else:
   # do all times
   start_date  = None
   end_date    = None
   plotting    = False
   show        = False

# define Greenland Sea according to Dany's definition
verts = [(-44,50),(-44,89),(8,89),(8,50)]
verts.append(verts[0])

outdir   = ddir.replace('data','analysis')
print('\noutput files in '+outdir+'\n')
if not os.path.exists(outdir):
   os.mkdir(outdir)
fli.MIZmap_all(start_date=start_date,end_date=end_date,outdir=outdir,\
      vertices=verts,plotting=plotting,show=show)
