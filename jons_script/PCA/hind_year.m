% collect results from "final_output" Netcdf files
% Extrect first (same as restart) and last (+1week or +57 rec  hindcast/forecast)

m_proj('stereographic','lat',90,'long',-45,'radius',50);

rsdir='/work/timill/RealTime_Models/TP4a0.12/expt_01.5/data/';

% add names to .mat file
hcdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bug_old/';
name='bug_old';
%hcdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bugfix/';
%name='bugfix';
%hcdir='/work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bugfix2/';
%name='bugfix2';

head='hindall'

getdata=1
plotyy=0
plotfig=0
savemat=1
% /work/timill/RealTime_Models/results_hindcasts/TP4a0.12/ice_only/bug/2015/2015_060/final_output

% Calculate Area in regular grid
DX=12.5e3; % m
DY=12.5e3; % m
A=DX*DY;

% time step
dt=3;       % hour
ddt=24/dt;  % steps per day 

YEAR={'2015', '2016'}
idm=609
jdm=881
ny=length(YEAR)

% if run through years and get data
if getdata==1
n=0;
wn=0;

% create empty result matrix n in N steps
time=zeros(1,ddt*ny*365);
datenr=zeros(1,ddt*ny*365);
VOL=zeros(1,ddt*ny*365);
ICEA=zeros(1,ddt*ny*365);
AWMH=zeros(1,ddt*ny*365);

% only weekly data dw in DW steps
twd=zeros(1,ny*52);
DV=zeros(1,ny*52);
DA=zeros(1,ny*52);
DH=zeros(1,ny*52);

%TIME=zeros(1,ddt*ny*365);
%CICE=zeros(idm,jdm,ddt*ny*365);
%HICE=zeros(idm,jdm,ddt*ny*365);
%SST=zeros(idm,jdm,ny*365);
%SSS=zeros(idm,jdm,ny*365);
%UICE=zeros(idm,jdm,ny*365);
%VICE=zeros(idm,jdm,ny*365);

icec2=zeros(idm,jdm);
iceh2=zeros(idm,jdm);

for y=1:ny
ydir=[hcdir,YEAR{y},'/']

% assume only hindcast day folders (YYYY_DDD)
files=dir(ydir);
dirfiles=[files.isdir];
fold=files(dirfiles); 
% dont use folder "." and ".." the 2 first
fold=fold(3:end);

D=length(fold)
for d = 1:D
   wn=wn+1; 
   % EX. SWARP_hindcast_ice_only_start20150302T000000Z.nc
   files=dir([ydir,fold(d).name,'/final_output/SWARP_hindcast_*Z.nc']);
   %files2=dir([ydir,fold(d+1).name,'/final_output/SWARP_hindcast_*Z.nc']);
   if (length(files)==1)
      file=[ydir,fold(d).name,'/final_output/',files.name];
   %   file2=[ydir,fold(d+1).name,'/final_output/',files2.name];
   else
      display 'PROBLEM - to many files in final_output folder !!'
      exit
   end
   
   %date1=file1(30:37)
   %date2=file2(30:37)
   % VARDATA = NCREAD(FILENAME,VARNAME)
   time=ncread(file,'time');
    
   if length(time)<56
      display('Days missing in hindcast: time<56 rec')
      display(['file: ',file])
      error('ERROR - STOP % CHECK YOUR DATA')
   end
   lat=ncread(file,'latitude');
   lon=ncread(file,'longitude');
   % DONT NEED THIS
   %x=ncread(file,'x');
   %y=ncread(file,'y'); 
   %X=x.*100; % (km) x/y in 100 km
   %Y=y.*100; % (km)
   %A=ones(IDM,JDM).*NaN;
   %for i=1:IDM-1
   %   for j=1:jdm-1
   %      A(i,j)=(x(i+1)-x(i))
   
   % compare last forecast day and initial assimilated field
   icec1=ncread(file,'icec',[1 1 1],[idm jdm 1]);
   iceh1=ncread(file,'icetk',[1 1 1],[idm jdm 1]);
   icec1(icec1<=0)=NaN;
   iceh1(iceh1<=0)=NaN;
   icec2(icec2<=0)=NaN;
   iceh2(iceh2<=0)=NaN;
   %dc(wn)=icec2-icec1;
   %dh(wn)=iceh2-iceh1;
   v1=icec1.*iceh1.*A;
   v2=icec2.*iceh2.*A;
   a1=icec1.*A;
   a2=icec2.*A;
   dv=v2-v1;
   da=a2-a1;
   h1=nansum(v1(:))./nansum(a1(:));
   h2=nansum(v2(:))./nansum(a2(:));
   DV(wn)=nansum(dv(:));
   DA(wn)=nansum(da(:));
   DH(wn)=h2-h1;
   tw=ncread(file,'time',1,1);
   tws=datestr(tw./(24*3600)+datenum(1970,01,01,00,00,00),'yyyy-mm-dd-HH:SS');
   twd(wn)=datenum(tws);
   % get last forecast fields (7 days) for next step
   icec2=ncread(file,'icec',[1 1 57],[idm jdm 1]);
   iceh2=ncread(file,'icetk',[1 1 57],[idm jdm 1]);

   for t=1:56

   d=ncread(file,'time',t,1);
   ts=datestr(d./(24*3600)+datenum(1970,01,01,00,00,00),'yyyy-mm-dd-HH:SS');
   td=datenum(ts);
   icec=ncread(file,'icec',[1 1 t],[idm jdm 1]);
   iceh=ncread(file,'icetk',[1 1 t],[idm jdm 1]);
   %sst=ncread(file,'sst',[1 1 t],[idm jdm t]);
   %sss=ncread(file,'sss',[1 1 t],[idm jdm t]);
   %uice=ncread(file,'uice',[1 1 t],[idm jdm t]);
   %vice=ncread(file,'vice',[1 1 t],[idm jdm t]);
   
   % TODO check that all dates are recorded ? +3h, so datenum(ts(n-1)) +0.125 =
   % datenum(ts(n)) ???
   n=n+1;
   %TIME(n)=str2num(ts);
   datenr(n)=td;
   disp(['Step:',num2str(n),' (',num2str(N),') Date:',datestr(datenr(n),'yyyy-mm-dd-HH:SS')])
   %CICE(:,:,n)=icec; 
   %HICE(:,:,n)=iceh;
   %SST(:,:,n)=sst;
   %SSS(:,:,n)=sss;
   %UICE(:,:,n)=uice;
   %VICE(:,:,n)=vice;

   % calculate average sea ice volume, ice area and area weighted thickness
   % VOL(n)=sum(sum(CICE(:,:,n).*HICE(:,:,n).*A));
   % ICEA(n)=sum(sum(CICE(:,:,n).*A));
   % AWMH(n)=VOL(n)./ICEA(n);
   
   icec(icec<=0)=NaN;
   iceh(iceh<=0)=NaN;
   vol=icec.*iceh.*A;
   VOL(n)  = nansum(vol(:));
   icea=icec.*A;
   ICEA(n) = nansum(icea(:));
   AWMH(n) = VOL(n)./ICEA(n);
   
   clear icec iceh 
   end

end

end

% NUMBER OF TOTAL TIME REC
N=n;
WN=wn;
% resize output to nr of recordings (N)
datenr=datenr(1:N)
VOL=VOL(1:N);
ICEA=ICEA(1:N);
AWMH=AWMH(1:N);
% weekly data, skip initial step no forecast to compare
DV=DV(2:WN);
DA=DA(2:WN);
DH=DH(2:WN);
twd=twd(2:WN)

%TIME=TIME(1:N)
%CICE=CICE(:,:,1:N);
%HICE=HICE(:,:,1:N);
%SST=SST(:,:,1:N);
%SSS=SSS(:,:,1:N);
%UICE=UICE(:,:,1:N);
%VICE=VICE(:,:,1:N);

% Remove ICEC and ICEH values < 0
%ICEC(ICEC<0)=0;
%ICEH(ICEH<0)=0;

if savemat==1
   save(['ICEVOL_20152016_',name],'datenr','VOL','ICEA','AWMH','DV','DA','DH','twd') 

%save(['CICE',head,'_',name],'CICE','ICEC','TIME','lon','lat') 
%save(['HICE',head,'_',name],'HICE','HICE','TIME','lon','lat')
%save(['SST',head,'_',name],'SST','TIME','lon','lat')
%save(['SSS',head,'_',name],'SSS','TIME','lon','lat')
%save(['UICE',head,'_',name],'UICE','TIME','lon','lat')
%save(['VICE',head,'_',name],'VICE','TIME','lon','lat')
end
end % if getdata==1

if plotyy==1
% Plot calculate average sea ice volume etc
km2=10e6;
km3=10e9;

figure(10)
[ax,ha1,ha2]=plotyy(datenr,VOL./km3,twd,DV./km3);
set(ha1,'color','k','linewidth',2);
set(ha2,'color','r','linestyle','--','linewidth',2);
set(ax(2),'Xtick',[],'Xticklabel',[],'ycolor','r');
set(ax(1),'Ycolor','k');
datetick('x','m')
xlabel('Year 2015 ansd 2016')
set(get(ax(1),'Ylabel'),'String','Tot. Sea ice Volume (km^3)')
set(get(ax(2),'Ylabel'),'String','Diff. in Sea ice Volume (km^3)')

figure(20)
[ax,ha1,ha2]=plotyy(datenr,ICEA./km2,twd,DA./km2)
set(ha1,'color','k','linewidth',2);
set(ha2,'color','r','linestyle','--','linewidth',2);
set(ax(2),'Xtick',[],'Xticklabel',[],'ycolor','r');
set(ax(1),'Ycolor','k');
datetick('x','m');
xlabel('Year 2015 ansd 2016')
set(get(ax(1),'Ylabel'),'String','Tot. Sea ice Area (km^2)');
set(get(ax(2),'Ylabel'),'String','Diff. in Sea ice Area (km^2)');

figure(30)
[ax,ha1,ha2]=plotyy(datenr,AWMH,twd,DH);
set(ha1,'color','k','linewidth',2);
set(ha2,'color','r','linestyle','--','linewidth',2);
set(ax(2),'Xtick',[],'Xticklabel',[],'ycolor','r');
set(ax(1),'Ycolor','k');
datetick('x','m');
xlabel('Year 2015 ansd 2016')
set(get(ax(1),'Ylabel'),'String','Area weighted mean height (m)');
set(get(ax(2),'Ylabel'),'String','Diff. in area weighted mean height (m)');
end

if plotfig==1

   % check tontinuity of Time
   figure; plot(1:N,td); ylabel('TIME axis');
   
   for t=1 %:N   
   % Figures
   % t=1
   figure(1000+t); clf;
   m_pcolor(lon,lat,CICE(:,:,t));
   shading flat;
   colorbar;
   m_gshhs_l('patch',[.2 .2 .2]);
   m_grid;
   title(['CICE: ',num2str(TIME(t))])
   print(gcf,'-dpng','-r300',[figdir,'/cice_',num2str(TIME(t))]);

   figure(2000+t); clf;
   m_pcolor(lon,lat,HICE(:,:,t));
   shading flat;
   colorbar;
   m_gshhs_l('patch',[.2 .2 .2]);
   m_grid;
   title(['HICE: ',num2str(TIME(t))])
   print(gcf,'-dpng','-r300',[figdir,'/hice_',num2str(TIME(t))]);

 end
end



