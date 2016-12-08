import MIZchar as mc
import fns_plotting as FP

if 0:
   bmap  = FP.start_HYCOM_map('gre')
else:
   lon0  = -10.
   lat0  = 77.
   lon1  = 2.
   lat1  = 83.
   bbox  = [lon0,lat0,lon1,lat1]
   bmap  = FP.start_map(bbox)

ddir  = '/home/nersc/lucia/MIZvalid'
tfil  = ddir+'/TP4archv_wav.2015_351_180000_dmax.txt'
tfo   = mc.single_file(tfil)

METH     = 5
Psolns   = tfo.get_solutions(METH=METH)
Psolns.plot_solutions(bmap)
