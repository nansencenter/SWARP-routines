% collect results from "final_output" Netcdf files
% Extrect first (same as restart) and last (+1week or +57 rec  hindcast/forecast)

m_proj('stereographic','lat',90,'long',-45,'radius',50);

rsdir='/work/timill/RealTime_Models/TP4a0.12/expt_01.5/data/';

hcdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bug_old/';
name='bug_old';
%hcdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bug/';
%name='bug';
%hcdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bugfix2/';
%name='bugfix2';


plotfig=0

% /work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bug/2015/2015_060/final_output

YEAR={'2015'} %, '2016'}
idm=609
jdm=881
ny=length(YEAR)


% 
n=0

for y=1:ny
ydir=[hcdir,YEAR{y},'/']

files=dir(ydir);
dirfiles=[files.isdir];
fold=files(dirfiles); 
% dont use folder "." and ".." the 2 first
fold=fold(3:end);

D=length(fold)
for d = 1:D

   n=n+1; 
   % EX. SWARP_hindcast_ice_only_start20150302T000000Z.nc
   files1=dir([ydir,fold(d).name,'/final_output/SWARP_hindcast_*Z.nc']);
   %files2=dir([ydir,fold(d+1).name,'/final_output/SWARP_hindcast_*Z.nc']);
   if (length(files1)==1)
      file1=[ydir,fold(d).name,'/final_output/',files1.name];
   %   file2=[ydir,fold(d+1).name,'/final_output/',files2.name];
   else
      display 'PROBLEM - to many files in final_output folder !!'
   end
   
   %date1=file1(30:37)
   %date2=file2(30:37)
   % VARDATA = NCREAD(FILENAME,VARNAME)
   time1=ncread(file1,'time',1,1);
   time2=ncread(file1,'time',57,1);
   ts1=datestr(time1./(24*3600)+datenum(1970,01,01,00,00,00),'yyyymmddHHSS');
   ts2=datestr(time2./(24*3600)+datenum(1970,01,01,00,00,00),'yyyymmddHHSS');
   lat=ncread(file1,'latitude');
   %lat2=ncread(file2,'latitude');
   lon=ncread(file1,'longitude');
   %lon2=ncread(file2,'longitude');
   icec1=ncread(file1,'icec',[1 1 1],[idm jdm 1]);
   icec2=ncread(file1,'icec',[1 1 57],[idm jdm 1]);
   iceh1=ncread(file1,'icetk',[1 1 1],[idm jdm 1]);
   iceh2=ncread(file1,'icetk',[1 1 57],[idm jdm 1]);
   
   
   sst1=ncread(file1,'sst',[1 1 1],[idm jdm 1]);
   sst2=ncread(file1,'sst',[1 1 57],[idm jdm 1]);
   sss1=ncread(file1,'sss',[1 1 1],[idm jdm 1]);
   sss2=ncread(file1,'sss',[1 1 57],[idm jdm 1]);
   
   uice1=ncread(file1,'uice',[1 1 1],[idm jdm 1]);
   uice2=ncread(file1,'uice',[1 1 57],[idm jdm 1]);
   vice1=ncread(file1,'vice',[1 1 1],[idm jdm 1]);
   vice2=ncread(file1,'vice',[1 1 57],[idm jdm 1]);


   display([' Step ',num2str(n),'  Restart:',num2str(ts1),....
            '  Hindcast:',num2str(ts2)])
    
   TIME2(n)=str2num(ts2);
   ICEC1(:,:,n)=icec1; 
   ICEC2(:,:,n)=icec2;
   ICEH1(:,:,n)=iceh1;
   ICEH2(:,:,n)=iceh2;
   
   % TODO moore fields sst ssh sss etc ?
   SST1(:,:,n)=sst1;
   SST2(:,:,n)=sst2;
   SSS1(:,:,n)=sss1;
   SSS2(:,:,n)=sss2;
   
   UICE1(:,:,n)=uice1;
   UICE2(:,:,n)=uice2;
   VICE1(:,:,n)=vice1;
   VICE2(:,:,n)=vice2;


end

end

% Remove ICEC and ICEH values < 0
ICEC1(ICEC1<0)=0;
ICEC2(ICEC2<0)=0;
ICEH1(ICEH1<0)=0;
ICEH2(ICEH2<0)=0;

N=n;
for t=1:N-1
   TIME(t)=TIME2(t);
   DIFFC(:,:,t)   =ICEC2(:,:,t)-ICEC1(:,:,t+1);
   DIFFH(:,:,t)   =ICEH2(:,:,t)-ICEH1(:,:,t+1);
   DIFFSST(:,:,t) =SST2(:,:,t) -SST1(:,:,t+1);
   DIFFSSS(:,:,t) =SSS2(:,:,t) -SSS1(:,:,t+1);
   DIFFUICE(:,:,t)=UICE2(:,:,t)-UICE1(:,:,t+1);
   DIFFVICE(:,:,t)=VICE2(:,:,t)-VICE1(:,:,t+1);

end

save(['DIFFC',name],'DIFFC','ICEC1','ICEC2','TIME','lon','lat') 
save(['DIFFH',name],'DIFFH','ICEH1','ICEH2','TIME','lon','lat')
save(['DIFFSST',name],'DIFFSST','SST1','SST2','TIME','lon','lat')
save(['DIFFSSS',name],'DIFFSSS','SSS1','SSS2','TIME','lon','lat')
save(['DIFFUICE',name],'DIFFUICE','UICE1','UICE2','TIME','lon','lat')
save(['DIFFVICE',name],'DIFFVICE','VICE1','VICE2','TIME','lon','lat')


if plotfig==1
 for t=1:N   
   % Figures
   % t=1
   figure(1000+t); clf;
   m_pcolor(lon,lat,DIFFC(:,:,t));
   shading flat;
   colorbar;
   m_gshhs_l('patch',[.2 .2 .2]);
   m_grid;
   title(['Concentration Hindcast - Ass.Resatart Date: ',num2str(TIME(t))])
   print(gcf,'-dpng','-r300',['Figures/fice_Hind-AssRest_',num2str(TIME(t))]);

   figure(2000+t); clf;
   m_pcolor(lon,lat,DIFFH(:,:,t));
   shading flat;
   colorbar;
   m_gshhs_l('patch',[.2 .2 .2]);
   m_grid;
   title(['Thickness (m) Hindcast - Ass.Resatart Date: ',num2str(TIME(t))])
   print(gcf,'-dpng','-r300',['Figures/hice_Hind-AssRest_',num2str(TIME(t))]);

 end
end



