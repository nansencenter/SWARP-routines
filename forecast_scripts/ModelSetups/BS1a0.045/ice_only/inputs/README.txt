Inputs for TP4 ice-only forecast in this folder
===========================================================================

1. infile.mal.outer: edited by make_infile4forecast.sh to become infile.in

2. blkdat.input:
*set frequency of archive outputs
> frequency (days) is 'dsurfq'
> needs to be a multiple of 'batrop'       (barotropic time step)
> batrop needs to be a divisor of 'baclin' (baroclinic time step)
eg:
dsurfq=9999.0 (no archive),   batrop=800, baclin=25
dsurfq=0.25   (6h),           batrop=800, baclin=25
dsurfq=0.125  (3h),           batrop=400, baclin=12.5

*Also check experiment number is correct

3. archv.extract:
*set variables to dump to TP4archv*.[ab] files

4. preprocess.sh:
* make sure archv.extract is copied to SCRATCH at runtime

5. pbsjob.sh:
* make sure $SWARP_PP post-processing option is correct

6. flags:
* keep these up-to-date in case need to recompile
* restore_build should recover them if Build is cleaned

7. THISFC.src:
* has experiment number, base directories for scripts, results and migrate backup
