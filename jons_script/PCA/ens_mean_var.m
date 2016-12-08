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
% Read in restart grid
IDM=800;
JDM=880;
[lon,dumb,dumb,dumb]=loada([fcdir,'regional.grid.a'],1,IDM,JDM);
[lat,dumb,dumb,dumb]=loada([fcdir,'regional.grid.a'],2,IDM,JDM);

% in ...ICE.uf files with 5 different fields, size of one
% memeber 
recnr=5*IDM*JDM*8;

% get all ensemble members
mem=100
N = length(list);
MEM=zeros(N,1);
% final results
Fm=zeros(IDM,JDM,N);
Fv=zeros(IDM,JDM,N);
Hm=zeros(IDM,JDM,N);
Hv=zeros(IDM,JDM,N);
CVFv=zeros(IDM,JDM,N);
CVHv=zeros(IDM,JDM,N);

% ensemble field
F=zeros(IDM,JDM,mem);
H=zeros(IDM,JDM,mem);

year=zeros(1,N);
tyd=zeros(1,N);

for n=1:N

 file=list(n).name
 tpdate=file(11:18);
 year(n)=str2num(file(11:14));
 tyd(n)=str2num(file(16:18));

 dfile=dir(file);
 sfile=dfile.bytes;
 mem=sfile/recnr;

 disp([file,' mem:',num2str(mem)])
 MEM(n)=mem;
 if mem==100

 % loop over all memmbers
 for m=1:mem
 % Read restart files
 % [fld,lon,lat,depths]=loadrestart(rstfile,varname,layer1,layer2);
 F(:,:,m)=read_restart_ice(file,'ficem',IDM,JDM,m);
 H(:,:,m)=read_restart_ice(file,'hicem',IDM,JDM,m);
 end % for m=1:mem
 
 % set all values below 0 to NaN
 F(F<=0)=NaN;
 H(H<=0)=NaN;
 % calculate mean and variance over dim 3, the ensemble 
 Fm(:,:,n)=nanmean(F,3);
 Hm(:,:,n)=nanmean(H,3);
 Fv(:,:,n)=nanvar(F,0,3);
 Hv(:,:,n)=nanvar(H,0,3);
 CVFv(:,:,n)=Fv(:,:,n)./Fm(:,:,n);
 CVHv(:,:,n)=Hv(:,:,n)./Hm(:,:,n);


 % new flipped hot map from white to dark
 hot=colormap(hot);
 hot2=flipud(hot);

 if plotm==1
 
 figure(1000); clf;
 m_pcolor(lon,lat,Fm(:,:,n));
 caxis([0 1]);
 shading flat;
 colormap(jet);
 colorbar;
 m_gshhs_l('patch',[.2 .2 .2]);
 m_grid;
 title(['Ensemble mean Sea ice fraction - ',tpdate,'mem=',num2str(mem)])
  
 figure(2000); clf;
 m_pcolor(lon,lat,Hm(:,:,n));
 caxis([0 5]);
 shading flat;
 colormap(jet);
 colorbar;
 m_gshhs_l('patch',[.2 .2 .2]);
 m_grid;
 title(['Ensemble mean Sea ice thickness - ',tpdate,' mem=',num2str(mem)])
 
 figure(3000); clf;
 m_pcolor(lon,lat,CVFv(:,:,n));
 caxis([0 1]);
 shading flat;
 colormap(hot2);
 colorbar;
 m_gshhs_l('patch',[.2 .2 .2]);
 m_grid;
 title(['CV of Sea ice fraction - ',tpdate,' mem=',num2str(mem)])
  
 figure(4000); clf;
 m_pcolor(lon,lat,CVHv(:,:,n));
 caxis([0 1]);
 shading flat;
 colormap(hot2);
 colorbar;
 m_gshhs_l('patch',[.2 .2 .2]);
 m_grid;
 title(['CV of Sea ice thickness - ',tpdate,' mem=',num2str(mem)]) 
  
 %pause
 end 

 clear F H 

 end % if mem==100
end % for n=1:N


% save FRS20152016 FRS LON LAT list year tyd
% save HRS20152016 HRS LON LAT list year tyd


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






