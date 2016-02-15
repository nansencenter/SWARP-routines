format long

if 0
   a              = 6378273;  % semi-major axis (m)
   ecc            = .1;       % eccentricity of ellipse
   ellipsoid_mat  = [a,ecc];  % ellipsoid formatted ala matlab

   if 1
      lons  = [40,50,50,40];
      lats  = [60,60,70,70];
   end

   A  = areaint(lats,lons,ellipsoid_mat,'degrees')
else
   ecc   = .2;
   phi   = 70
   chi   = convertlat([1,ecc], phi, 'geodetic', 'conformal', 'degrees')
   phi2  = convertlat([1,ecc], chi, 'conformal','geodetic', 'degrees')
end
