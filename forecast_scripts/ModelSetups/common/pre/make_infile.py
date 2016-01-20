import sys
from getopt import getopt

class getinfo:
   def __init__(self,infile):

      f     = open(infile)
      lines = f.readlines()
      f.close()

      self.days         = None
      self.restart_opts = None
      self.refyear      = None
      self.expt_dir     = None
      self.nest_outer   = "F"
      self.nest_inner   = "F"
      self.final_hour   = "00"

      for lin in lines:
         Lin   = lin.split()
         name  = Lin[0]
         if len(Lin)>2:
            # keep as list
            setattr(self,name,Lin[1:])
         else:
            setattr(self,name,Lin[1])

      if self.restart_opts is None:
         self.restart_opts = len(self.days)*['F']

      return

# ==============================================================
# options
opts,args   = getopt(sys.argv[1:],"",["infile="])
for opt,arg in opts:
   if opt=='--infile':
      infile   = arg

if infile is None:
   raise ValueError('specify "infile" input (--infile=)')
# ==============================================================

info  = getinfo(infile)
for name in vars(info).keys():
   val   = getattr(info,name)
   # print(name,val)
   if val is None:
      raise ValueError(name+' not specified in '+infile)

outfile  = info.expt_dir+'/infile.in'
f        = open(outfile,'w')

f.write('3.1          # version of inputfile (limits.dat)            DO NOT CHANGE\n')
f.write(info.rungen+'          # (rungen)   Version number of run\n')
f.write(info.refyear+'         # Reference year for simulation  (365. for spinup)\n')

# Limits
ss    = 4*' '
ss2   = 6*' '
day   = info.days[0]
f.write(ss[:4-len(day)]+day+' 00      # (nday1)    First day of integration NOTE FORMAT F9.2\n')
lin0  = ss2[:6-len(day)]+day+'.00  F     T        '+info.restart_opts[0]+'\n'

day   = info.days[-1]
f.write(ss[:4-len(day)]+day+' '+info.final_hour+'      # (nday2)    Last day of integration NOTE FORMAT F9.2\n')
lin1  = ss2[:6-len(day)]+day+'.00  F     T        '+info.restart_opts[-1]

# forcing etc
f.write('ecncF era40  # forcing option, month, ecmwf, ncepr, ecmo, ecnc\n')
f.write(' 200.0       # temperature relaxation time scale (F in blkdat.input)\n')
f.write(' 200.0       # salinity    relaxation time scale\n')
f.write('F  2         # laverage n Accumulate monthly averages every n hours\n')
f.write('T  F         # Switch on and accumulate (T) or overwrite (F) daily averages\n')

# nesting
f.write(info.nest_outer+'  6         # lnesto, nestdto - saves nesting bnd cond at nestdto intervals\n')
f.write(info.nest_inner+'  6         # lnesti, nestdti - read and apply nesting bnd cond at nestdto intervals\n')
f.write('F CSR F F    # Tides (true,  CSR/FES, apply currents)\n')
f.write('F            # lgridp     Activate storage of gridpoint information\n')

# diag days
f.write('# Days     Assi  Diagno.  Restart            (1992 for first synoptic experiment)\n')
f.write(lin0)


for n,day in enumerate(info.days[1:-1]):
   ss=6*' '
   f.write(ss[:6-len(day)]+day+'.00  F     T        '+info.restart_opts[n+1]+'\n')

f.write(lin1)
f.close()
