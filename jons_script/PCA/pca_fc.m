% read in reastart and daily mean
% calculate principal component analysis PCA of the forecast and compare to restart

m_proj('stereographic','lat',90,'long',-45,'radius',50);
plotm=0;

rsb='/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/2015/'

homedir=pwd;
pcadir=['/work/timill/PCA/']

cd(pcadir)

list=dir('TP4restart2015*ICE.uf')

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

% Read in restart grid
idm=800;
jdm=880;
[lonrs,dumb,dumb,dumb]=loada(['regional.grid.a'],1,idm,jdm);
[latrs,dumb,dumb,dumb]=loada(['regional.grid.a'],2,idm,jdm);

% for now only use member 1
mem=1
N = length(list)
FRS=zeros(idm*jdm,N);
HRS=zeros(idm*jdm,N);

for n=1:N
 file=list(n).name
 tpdate=file(11:18);

 % Read restart files
 % [fld,lon,lat,depths]=loadrestart(rstfile,varname,layer1,layer2);
 %cd(rsd)
 frs=read_restart_ice(file,'ficem',idm,jdm,mem);
 hrs=read_restart_ice(file,'hicem',idm,jdm,mem);
% I=find(frs==0); frs(I)=NaN;
% I=find(hrs==0); hrs(I)=NaN;

 if plotm==1
  figure(100+n); clf;
  m_pcolor(lonrs,latrs,frs);
  caxis([0 1]);
  shading flat;
  colorbar;
  m_gshhs_l('patch',[.2 .2 .2]);
  m_grid;
  title(['Sea ice fraction - ',tpdate])
 
  figure(200+n); clf;
  m_pcolor(lonrs,latrs,hrs);
  caxis([0 5]);
  shading flat;
  colorbar;
  m_gshhs_l('patch',[.2 .2 .2]);
  m_grid;
  title(['Sea ice thickness - ',tpdate])
 end

FRS(:,n)=frs(:)';
HRS(:,n)=hrs(:)';

end

% Find all grid points that have ice all period?
% "detrend", remove time mean for each grid point
%I=find(sum(isfinite(mFRS),2));
%If=find(sum(isfinite(FRS),2));
%Ih=find(sum(isfinite(HRS),2));
% If or Ih ?
If=find(sum(FRS,2))
FRS=FRS(If,:);
HRS=HRS(If,:);
mFRS=mean(FRS,2);
mHRS=mean(HRS,2);

IN=length(If)
for n=1:N
 dFRS(:,n)=FRS(:,n)-mFRS;
 dHRS(:,n)=HRS(:,n)-mHRS;
end

for j=1:IN
[UF,DF,VF]=svd(dFRS);
[UH,DH,VH]=svd(dHRS);


cd(homedir)






