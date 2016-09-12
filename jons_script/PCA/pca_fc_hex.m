% read in reastart and daily mean
% calculate principal component analysis PCA of the forecast and compare to restart

m_proj('stereographic','lat',90,'long',-45,'radius',50);
plotm=0;

%rsb='/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/2015/'

homedir=pwd;
fcdir='/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data/';
%pcadir=['/work/timill/PCA/']
rsdir='/work/timill/RealTime_Models/TP4a0.12/expt_01.5/data/';
cd(rsdir)

list=dir('TP4restart*ICE.uf')

%rsd='/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data/';
%rsfile='TP4restart2016_087_00.a';
%rsf=[rsd,rsfile]
%rsicefile='TP4restart2016_087_00ICE.uf';
%rsice=[rsd,rsicefile]
%rsdate=rsfile(11:18)

%fcd='/work/timill/RealTime_Models/results/TP4a0.12/ice_only/20160330/final_output/'
%file='SWARPiceonly_forecast_start20160330T000000Z.nc'
%fcf=[fcd,file];

% FC Variables: icec, icetk, usurf, vsurf, sst, ssh, sss, uice, vice
st=1;
% Read in restart grid
idm=800;
jdm=880;
[lon,dumb,dumb,dumb]=loada([fcdir,'regional.grid.a'],1,idm,jdm);
[lat,dumb,dumb,dumb]=loada([fcdir,'regional.grid.a'],2,idm,jdm);
LON=lon(1:st:end,1:st:end);
LAT=lat(1:st:end,1:st:end);
[IDM,JDM]=size(LON);

% for now only use member 1
mem=1
N = length(list)
FRS=zeros(IDM*JDM,N);
HRS=zeros(IDM*JDM,N);
year=zeros(1,1:N);
tyd=zeros(1,1:N);

for n=1:N
 if n==1
  plotm=1
 else
  plotm=0
 end

 file=list(n).name
 tpdate=file(11:18);
 year(n)=str2num(file(11:14));
 tyd(n)=str2num(file(16:18));
 
 % Read restart files
 % [fld,lon,lat,depths]=loadrestart(rstfile,varname,layer1,layer2);
 %cd(rsd)
 frs=read_restart_ice(file,'ficem',idm,jdm,mem);
 hrs=read_restart_ice(file,'hicem',idm,jdm,mem);
% I=find(frs==0); frs(I)=NaN;
% I=find(hrs==0); hrs(I)=NaN;

 if plotm==1
  figure(100+n); clf;
  m_pcolor(lon,lat,frs);
  caxis([0 1]);
  shading flat;
  colorbar;
  m_gshhs_l('patch',[.2 .2 .2]);
  m_grid;
  title(['Sea ice fraction - ',tpdate])
 
  figure(200+n); clf;
  m_pcolor(lon,lat,hrs);
  caxis([0 5]);
  shading flat;
  colorbar;
  m_gshhs_l('patch',[.2 .2 .2]);
  m_grid;
  title(['Sea ice thickness - ',tpdate])
 end

 frs=frs(1:st:end,1:st:end);
 hrs=hrs(1:st:end,1:st:end);
 FRS(:,n)=frs(:)';
 HRS(:,n)=hrs(:)';

 if plotm==1
 figure(1000+n); clf;
 m_pcolor(LON,LAT,frs);
 caxis([0 1]);
 shading flat;
 colorbar;
 m_gshhs_l('patch',[.2 .2 .2]);
 m_grid;
 title(['Sea ice fraction - ',tpdate])
  
 figure(2000+n); clf;
 m_pcolor(LON,LAT,hrs);
 caxis([0 5]);
 shading flat;
 colorbar;
 m_gshhs_l('patch',[.2 .2 .2]);
 m_grid;
 title(['Sea ice thickness - ',tpdate])
 pause
 end

end

save FRS20152016 FRS LON LAT list year tyd
save HRS20152016 HRS LON LAT list year tyd


% Find all grid points that have ice all period?
% "detrend", remove time mean for each grid point
%I=find(sum(isfinite(mFRS),2));
%If=find(sum(isfinite(FRS),2));
%Ih=find(sum(isfinite(HRS),2));
% If or Ih ?
%If=find(sum(FRS,2))
%FRS=FRS(If,:);
%HRS=HRS(If,:);
%mFRS=mean(FRS,2);
%mHRS=mean(HRS,2);
%
%IN=length(If)
%for n=1:N
% dFRS(:,n)=FRS(:,n)-mFRS;
% dHRS(:,n)=HRS(:,n)-mHRS;
%end
%
%for j=1:IN
%[UF,DF,VF]=svd(dFRS);
%[UH,DH,VH]=svd(dHRS);


cd(homedir)






