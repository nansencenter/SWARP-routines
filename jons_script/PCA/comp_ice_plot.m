% read in yearly sea ice data from "hind_year" and 
% compare sea ice vol, area, thickness
% compare the "jump" from the freerun 7 days forecast to the ensemble sim. from met.no

figdir='/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/jons_script/PCA/'
savefigs=1

name={'bug_old','bugfix','bugfix2'};
c={'r','b','k'};
l={'-','--','-'};
N=3;

mai1=datenum(2015,05,01);
jun1=datenum(2015,06,01);
jul1=datenum(2015,07,01);
aug1=datenum(2015,08,01);

for i=1:N

% load data datenr, VOL, ICEA, AWMH, DV, DA, DH
load(['ICEVOL_20152016_',name{i},'.mat']);

% loop over experiments
for z=3
   if z==1
      tz1=1;
      tz2=length(datenr);
      dz1=datenr(tz1);
      dz2=datenr(tz2);
      zname='';
elseif z==2
      I=find(datenr>=jul1); tz1=I(1);
      I=find(datenr>=aug1); tz2=I(1);
      tz1=tz1;
      tz2=tz2;
      dz1=datenr(tz1);
      dz2=datenr(tz2);
      zname='zoom';
   elseif z==3
      I=find(datenr>=mai1); tz1=I(1);
      I=find(datenr>=jun1); tz2=I(1);
      tz1=tz1;
      tz2=tz2;
      dz1=datenr(tz1);
      dz2=datenr(tz2);
      zname='zoom2';
   end

% Plot calculated averages of sea ice volume etc
km2=10e6;
km3=10e9;

figure(10); hold on;
h=plot(datenr(tz1:tz2),VOL(tz1:tz2)./(1e3*km3));
set(h,'color',c{i},'linestyle',l{i},'linewidth',2);
datetick('x','m')
xlabel('Year 2015 and 2016','Fontsize',12,'fontweight','bold')
ylabel('Tot. Sea ice Volume (1x10^3 km^3)','Fontsize',12,'fontweight','bold')
set(gca,'Fontsize',12,'fontweight','bold')

figure(11); hold on;
h=plot(twd,DV./(1e3*km3));
set(h,'color',c{i},'linestyle',l{i},'linewidth',2);
datetick('x','m')
xlabel('Year 2015 and 2016','Fontsize',12,'fontweight','bold')
ylabel('Diff. in Sea ice Volume (1x10^3 km^3)','Fontsize',12,'fontweight','bold')
set(gca,'Fontsize',12,'fontweight','bold')

figure(20); hold on;
h=plot(datenr(tz1:tz2),ICEA(tz1:tz2)./(1e6*km2));
set(h,'color',c{i},'linestyle',l{i},'linewidth',2);
datetick('x','m');
xlabel('Year 2015 and 2016','Fontsize',12,'fontweight','bold')
ylabel('Tot. Sea ice Area (1x10^6 km^2)','Fontsize',12,'fontweight','bold');
set(gca,'Fontsize',12,'fontweight','bold')

figure(21); hold on;
h=plot(twd,DA./(1e6*km2));
set(h,'color',c{i},'linestyle',l{i},'linewidth',2);
datetick('x','m');
xlabel('Year 2015 and 2016','Fontsize',12,'fontweight','bold')
ylabel('Diff. in Sea ice Area (1x10^6 km^2)','Fontsize',12,'fontweight','bold');
set(gca,'Fontsize',12,'fontweight','bold')

figure(30); hold on;
h=plot(datenr(tz1:tz2),AWMH(tz1:tz2));
set(h,'color',c{i},'linestyle',l{i},'linewidth',2);
datetick('x','m');
xlabel('Year 2015 and 2016','Fontsize',12,'fontweight','bold')
ylabel('Area weighted mean height (m)','Fontsize',12,'fontweight','bold');
set(gca,'Fontsize',12,'fontweight','bold')

figure(31); hold on;
h=plot(twd,DH);
set(h,'color',c{i},'linestyle',l{i},'linewidth',2);
datetick('x','m');
xlabel('Year 2015 and 2016','Fontsize',12,'fontweight','bold')
ylabel('Diff. in area weighted mean height (m)','Fontsize',12,'fontweight','bold');
set(gca,'Fontsize',12,'fontweight','bold')

end % if z

 if i<N
  clear AWMH DA DH DV ICEA VOL datenr twd
 end

end % for i=1:N

% write legends in emtpy figures
% TODO not 0 !!!!!
t1=datenr(1);
t2=datenr(end);
figure(10); hold on;
plot([dz1 dz2],[0 0],'r--',[dz1 dz2],[0 0],'b--',[dz1 dz2],[0 0],'k--');
h=legend('bug old','bugfix','bugfix2','Location','SouthEast')
set(h,'Fontsize',10,'fontweight','bold');
xlim([dz1 dz2]);
if z==2; datetick('x','dd'); xlabel('July 2015','Fontsize',12,'fontweight','bold'); end;
if z==3; datetick('x','dd'); xlabel('May 2015','Fontsize',12,'fontweight','bold'); end;
figure(11); hold on;
plot([t1 t2],[0 0],'r--',[t1 t2],[0 0],'b--',[t1 t2],[0 0],'k--');
h=legend('bug old','bugfix','bugfix2','Location','SouthEast')
set(h,'Fontsize',10,'fontweight','bold');
figure(20); hold on;
plot([dz1 dz2],[0 0],'r--',[dz1 dz2],[0 0],'b--',[dz1 dz2],[0 0],'k--');
h=legend('bug old','bugfix','bugfix2','Location','SouthEast')
set(h,'Fontsize',10,'fontweight','bold');
xlim([dz1 dz2]);
if z==2; datetick('x','dd'); xlabel('July 2015','Fontsize',12,'fontweight','bold'); end;
if z==3; datetick('x','dd'); xlabel('May 2015','Fontsize',12,'fontweight','bold'); end;
figure(21); hold on;
plot([t1 t2],[0 0],'r--',[t1 t2],[0 0],'b--',[t1 t2],[0 0],'k--');
h=legend('bug old','bugfix','bugfix2','Location','SouthEast')
set(h,'Fontsize',10,'fontweight','bold');
figure(30); hold on;
plot([dz1 dz2],[0 0],'r--',[dz1 dz2],[0 0],'b--',[dz1 dz2],[0 0],'k--');
h=legend('bug old','bugfix','bugfix2','Location','SouthEast')
set(h,'Fontsize',10,'fontweight','bold');
xlim([dz1 dz2]);
if z==2; datetick('x','dd'); xlabel('July 2015','Fontsize',12,'fontweight','bold'); end;
if z==3; datetick('x','dd'); xlabel('May 2015','Fontsize',12,'fontweight','bold'); end;
figure(31); hold on;
plot([t1 t2],[0 0],'r--',[t1 t2],[0 0],'b--',[t1 t2],[0 0],'k--');
h=legend('bug old','bugfix','bugfix2','Location','SouthEast')
set(h,'Fontsize',10,'fontweight','bold');

if savefigs==1
% save figures
figure(10); print(gcf,'-dpng','-r300',[figdir,'/Vol_20152016',zname]);
figure(11); print(gcf,'-dpng','-r300',[figdir,'/DiffVol_20152016']);
figure(20); print(gcf,'-dpng','-r300',[figdir,'/Area_20152016',zname]);
figure(21); print(gcf,'-dpng','-r300',[figdir,'/DiffArea_20152016']);
figure(30); print(gcf,'-dpng','-r300',[figdir,'/Thick_20152016',zname]);
figure(31); print(gcf,'-dpng','-r300',[figdir,'/DiffThick_20152016']);

end



