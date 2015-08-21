yr    = 2015;
vlist = {'D2M','MSL','T2M','TCC','U10M','V10M'};
Nv    = length(vlist);
%%
cyr   = num2str(yr);
fmt   = 'yyyymmdd';

%%matlab date is days since 1 Jan 0 AD
dref  = datenum([cyr   ,'0101'],fmt);%%start of file (start of year)
dref2 = datenum(['1950','0101'],fmt);%%nc reference

for j=1:Nv
   vbl   = vlist{j};
   ncfil = ['ec_atmo_geo_la_',vbl,'_',cyr,'.nc'];
   %%
   disp(' ');
   disp(['Opening ',ncfil, '...']);
   %%
   nc    = netcdf(ncfil);
   time  = nc{'time',1}(:);%%hours since 1950-01-01
   close(nc);
   nrec  = length(time);

   if 1
      %% check even spacing of records
      time2 = time(1)+(0:nrec-1)'*6;
      jbad  = find(time2~=time);
      if  ~isempty(jbad)
         tbad        = time(jbad);
         tbad_days   = dref2+tbad/24;
         disp(' ');
         disp('   Bad time values:');
         disp(tbad);
         disp('   Corresponding dates:');
         disp(datestr(tbad_days,fmt));
         disp(' ');
      end
   end



   d_today  = floor(now);
   j_today  = d_today+1-dref; %% julian day
   Nfc      = 8*4+3;          %% forecast stops at 12:00 on the 9th day from today

   disp(['Variable                                 : ',vbl]);
   disp(['Number of records in file                : ',num2str(nrec         )]);
   disp(['Day number                               : ',num2str(j_today      )]);
   disp(['Number of records up to end of today     : ',num2str(j_today*4    )]);
   disp(['Number of records up to end of forecast  : ',num2str(j_today*4+Nfc)]);
   disp(' ');
end
