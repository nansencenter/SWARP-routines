%%changed calc of dmean in mod_wavesice.F to be more readable
%%-simple script to check it
dmax  = 80
dmin  = 20;
f     = .9;
xi    = 2;
mom   = 1

%%new method
mm    = 0;
r     = dmax/dmin;
while ( r > xi )
   r  = r/xi;
   mm = mm+1;
end

mm

if ( mm > 0 )

   nsum  = 0.0;
   ndsum = 0.0;

   for m = 0:mm
      nm    = (1.0-f)*(f*xi^2)^m;
      dm    = dmax/(xi^m);
      nsum  = nsum +nm;
      ndsum = ndsum+nm*dm^mom;
   end
   dmean = ndsum/nsum;
else
   dmean = dmin^mom;
end

dmean2   = dmean;

%%old method
if ( mm > 0 )

  n     = 0.0;
  nsum  = 0.0;
  nd    = 0.0;
  ndsum = 0.0;

  for m = 0:mm
     n     = (1.0-f)*(f*xi^2)^m;
     nd    = n/(xi^m);
     nsum  = nsum +n;
     ndsum = ndsum+nd;
     dfac  = ndsum/nsum;
  end
  dmean = dfac*dmax

else
  dmean = dmin
end

if mom==2
   Sbot     = dmean2;
   Slat     = 4*dmean;
   alp_lat  = Slat/(Slat+Sbot)
else
   test_method = dmean2/dmean
end
