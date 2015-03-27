inputs to TP4 model in this folder
1. infile.mal.outer: edited by make_infile4forecast.sh to become infile.in

2. blkdat.input: set frequency of archive outputs
> frequency (days) is 'dsurfq'
> needs to be a multiple of 'batrop'       (barotropic time step)
> batrop needs to be a divisor of 'baclin' (baroclinic time step)
eg:
dsurfq=9999.0 (no archive),   batrop=800, baclin=25
dsurfq=0.25   (6h),           batrop=800, baclin=25
dsurfq=0.125  (3h),           batrop=400, baclin=12.5
