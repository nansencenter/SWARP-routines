% read in reastart and daily mean
% calculate principal component analysis PCA of the forecast and compare to restart

m_proj('stereographic','lat',90,'long',-45,'radius',50);

homedir=pwd
sdir=[homedir,'/SCRATCH/']

rsd='/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data/';
rsfile='TP4restart2016_087_00.a';
rsf=[rsd,rsfile]
rsicefile='TP4restart2016_087_00ICE.uf';
rsice=[rsd,rsicefile]
rsdate=rsfile(11:18)

fcd='/work/timill/RealTime_Models/results/TP4a0.12/ice_only/20160330/final_output/'
file='SWARPiceonly_forecast_start20160330T000000Z.nc'
fcf=[fcd,file];
date=file(28:35);
% FC Variables: icec, icetk, usurf, vsurf, sst, ssh, sss, uice, vice

% Read in restart grid
idm=800;
jdm=880;
[lonrs,dumb,dumb,dumb]=loada([rsd,'regional.grid.a'],1,idm,jdm);
[latrs,dumb,dumb,dumb]=loada([rsd,'regional.grid.a'],2,idm,jdm);


% Read restart files
% [fld,lon,lat,depths]=loadrestart(rstfile,varname,layer1,layer2);
%cd(rsd)
frs=read_restart_ice(rsice,'ficem',idm,jdm,1);
hrs=read_restart_ice(rsice,'hicem',idm,jdm,1);
I=find(frs==0); frs(I)=NaN;
I=find(hrs==0); hrs(I)=NaN;

figure(100)
m_pcolor(lonrs,latrs,frs);
caxis([0 1]);
shading flat;
colorbar;
m_gshhs_l('patch',[.2 .2 .2]);
m_grid;
title(['Sea ice fraction',rsdate])

figure(200)
m_pcolor(lonrs,latrs,hrs);
caxis([0 5]);
shading flat;
colorbar;
m_gshhs_l('patch',[.2 .2 .2]);
m_grid;
title(['Sea ice thickness - ',rsdate])


fcrec=1

%cd(homedir)

% Read in forecast variables (ncread need matlab 2012)
% VARDATA = NCREAD(FILENAME,VARNAME)
lonfc=ncread(fcf,'longitude');
latfc=ncread(fcf,'latitude');
ffc=ncread(fcf,'icec');
hfc=ncread(fcf,'icetk');
I=find(ffc==0); ffc(I)=NaN;
I=find(hfc==0); hfc(I)=NaN;


figure(300)
m_pcolor(lonfc,latfc,ffc(:,:,fcrec));
caxis([0 1]);
shading flat;
colorbar;
m_gshhs_l('patch',[.2 .2 .2]);
m_grid;
title(['Sea ice fraction FC - ',date])


figure(400)
m_pcolor(lonfc,latfc,hfc(:,:,fcrec));
caxis([0 5]);
shading flat;
colorbar;
m_gshhs_l('patch',[.2 .2 .2]);
m_grid;
title(['Sea ice thickness FC - ',date])




