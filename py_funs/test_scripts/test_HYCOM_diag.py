import numpy as np

class struct:
   def __init__(self,lut):
      for name in lut.keys():
         setattr(self,name,lut[name])
      return

class diag_info:

   #################################################################
   def __init__(self,tfil):

      f     = open(tfil)
      lines = f.readlines()
      f.close()

      ##########################################
      for i,lin in enumerate(lines):

         if '#' in lin:
            spl   = lin.split('#')
            name  = spl[1].strip()
            val   = spl[0].strip()
            if '*' not in val:
               val   = float(val)
               setattr(self,name,val)
            else:
               setattr(self,name,np.nan)

         elif 'Omega' in lin:
            ifreq = i+2
            # print(ifreq)

         elif 'dirn' in lin:
            idir  = i+2
            # print(idir)
        
         elif 'spectrum' in lin:
            ispec = i+2
            # print(ispec)
        
         elif 'atten' in lin:
            iatt = i+2
            # print(ispec)
      ##########################################

      # integer fields
      for name in ['nfreq','ndir','i','j']:
         val   = getattr(self,name)
         setattr(self,name,int(val))

      # vector
      omega = []
      nfreq = self.nfreq
      for i in range(ifreq,ifreq+nfreq):
         omega.append(float(lines[i]))
      self.omega  = np.array(omega)

      # vector
      atten = []
      for i in range(iatt,iatt+nfreq):
         atten.append(float(lines[i]))
      self.atten  = np.array(atten)

      # vector
      wavdir   = []
      ndir     = self.ndir
      for i in range(idir,idir+ndir):
         wavdir.append(float(lines[i]))
      self.wavdir = np.array(wavdir)

      # array
      Sdir  = []
      for i in range(ispec,ispec+ndir*nfreq):
         Sdir.append(float(lines[i]))

      Sdir        = np.array(Sdir)
      Sdir        = Sdir.reshape(ndir,nfreq,order='fortran')
      self.Sdir   = Sdir

      # wave stress (magnitude)
      self.tau_wav  = np.sqrt(self.taux_wav**2+self.tauy_wav**2)

      self.get_prams()
      self.do_ints()

      return
   #################################################################

   #################################################################
   def get_prams(self):
      
      lut   = {}

      lut.update({'radian':np.pi/180.})
      lut.update({'radinv':180./np.pi})
      lut.update({'rho_wtr':1025})
      lut.update({'gravity':9.81})

      self.prams  = struct(lut)
      return
   #################################################################

   #################################################################
   def do_ints(self):

      wt_omega          = 2+0*self.omega
      wt_omega[[0,-1]]  = 1
      wt_omega[1::2]    = 4
      #
      dom      = abs(self.omega[1]-self.omega[0])
      wt_omega = dom/3.*wt_omega

      radian   = self.prams.radian
      radinv   = 1/radian
      dtheta   = radian*abs(self.wavdir[1]-self.wavdir[0])
      wt_theta = dtheta+0*self.wavdir

      # get rotated version directly
      mth   = (self.wavdir*wt_theta).dot(self.Sdir).dot(wt_omega)

      # tau_wav (mag of wave stress)
      # - these are relative to a lon-lat grid
      # - model outputs are on native grid
      # - model outputs are also for the previous time step
      adv_dir  = -radian*(90.+self.wavdir)
      cos      = np.cos(adv_dir)
      sin      = np.sin(adv_dir)
      cp       = self.prams.gravity/self.omega           # phase vel: om/k=om/(om^2/g)=g/om
      cg       = cp/2.                                   # group vel: dom/dk=1/(2*om/g)=g/2/om=cp/2
      fac      = self.prams.gravity*self.prams.rho_wtr   # m/s^2*kg/m^3=N/m^3=Pa/m
      source   = self.Sdir.dot(np.diag(-fac*self.atten)) # m^2s/m*Pa/m=Pa.s

      taux_wav = -(cos*wt_theta).dot(source).dot(wt_omega*cg/cp) # (W-E)
      tauy_wav = -(sin*wt_theta).dot(source).dot(wt_omega*cg/cp) # (S-N)
      tau_wav  = np.sqrt(taux_wav**2+tauy_wav**2)

      Sfreq = wt_theta.transpose().dot(self.Sdir)
      m0    = Sfreq.dot(wt_omega)
      m2    = Sfreq.dot(wt_omega*self.omega**2)

      lut   = {}

      lut.update({'Hs':4*np.sqrt(m0)})
      lut.update({'Tp':2*np.pi*np.sqrt(m0/m2)})
      lut.update({'mwd':mth/m0})
      lut.update({'m0':m0})
      lut.update({'m2':m2})
      lut.update({'mth':mth})
      lut.update({'taux_wav':taux_wav})
      lut.update({'tauy_wav':tauy_wav})
      lut.update({'tau_wav':tau_wav})

      self.test      = struct(lut)
      self.wt_omega  = wt_omega
      self.wt_theta  = wt_theta

      return
   #################################################################

   #################################################################
   def compare(self):

      print(' ')
      vv = vars(self.test)

      for name in vv.keys():
         v1 = vv[name]
         v0 = getattr(self,name)
         print(name+' (here/HYCOM): '+str(v0)+' | '+str(v1))

      print(' ')
      return
   #################################################################

class therm_HYCOM:

   ################################################################
   # def __init__(self,tml,sml,hml,hice,fice,hsnw):
   def __init__(self,hice,fice,hsnw):

      self.units  = {}

      # variables
      # self.tml    = tml    # deg C     temperature of mixed layer
      # self.sml    = sml    # 1e-3      salinity of mixed layer
      # self.hml    = hml    # m         thickness of mixed layer
      # self.units.update({'tml':'deg C'})
      # self.units.update({'sml':'1e-3'})
      # self.units.update({'hml':'m'})

      self.fice   = fice   # -         ice fraction
      self.hice   = hice   # m         ice thickness
      self.hsnw   = hsnw   # m         snow thickness
      self.fsnw   = fice   # -         snow fraction = ice fraction
      self.units.update({'fice':''})
      self.units.update({'hice':'m'})
      self.units.update({'hsnw':'m'})
      self.units.update({'fsnw':''})

      # other parameters
      self.fusi      = 3.02e8 # J/m^3     heat of fusion of ice
      self.fuss      = 1.10e8 # J/m^3     heat of fusion of snow
      self.cpsw      = 3987.  # J/K/kg    specific heat of seawater
      self.t0deg     = 273.15 # K         zero deg celcius in K
      self.fice_max  = .995   # -         maximum sea ice fraction
      self.units.update({'fusi':'J/m^3'})
      self.units.update({'fuss':'J/m^3'})
      self.units.update({'cpsw':'J/K/kg'})
      self.units.update({'t0deg':'K'})

      # dependent quantities
      # self.tice_f = 273.216-.057*sml-self.t0deg          # deg C     freezing point of sea water
      # self.rhosw  = (1000.+self.sig(tml-self.t0deg,sml)) # kg/m^3    density of sea water
      # self.cpfac  = self.cpsw*self.rhosw*hml             # J/K/m^2   energy to heat 1m^2 of sw 1K  
      # self.units.update({'tice_f':'deg C'})
      # self.units.update({'rhosw':'kg/m^3'})
      # self.units.update({'cpfac':'J/K/m^2'})

      # ice
      self.vol_ice   = fice*hice                # m^3/m^2 = m  Volume density of ice
      self.heat_ice  = self.fusi*self.vol_ice   # J/m^2        Heat density of ice
         # flux=heat/dt; <0 freezing; >0 melting
         # vol=-flux*dt/fusi: >0 freezing, <0 melting
      self.units.update({'vol_ice':'m'})
      self.units.update({'heat_ice':'J/m^2'})

      # snow
      self.vol_snw   = self.fsnw*hsnw           # m^3/m^2 = m  Volume density of snow
      self.heat_snw  = self.fuss*self.vol_snw   # J/m^2        Heat density of snow
      self.units.update({'vol_snw':'m'})
      self.units.update({'heat_snw':'J/m^2'})

      # # mixed layer
      # self.heat_ml   = self.cpfac*(self.tml-self.tice_f)# J/m^2        Heat density of mixed layer
      # self.units.update({'heat_ml':'J/m^2'})
      return               
   ################################################################


   ################################################################
   def sig(self,t,s):
      # (from stmt_fns.h)
      # --- sigma-theta as a function of temp (deg c) and salinity (psu)
      # --- (friedrich-levitus, polynomial fit that is cubic in T and linear in S)

      # coefficients for sigma-0 (based on Brydon & Sun fit)
      c1    = -1.36471E-01    # const. coefficent
      c2    =  4.68181E-02    # T      coefficent
      c3    =  8.07004E-01    # S      coefficent
      c4    = -7.45353E-03    # T^2    coefficent
      c5    = -2.94418E-03    # T*S    coefficent
      c6    =  3.43570E-05    # T^3    coefficent
      rc6   =  1.0/c6
      c7    =  3.48658E-05    # T^2S   coefficent
      pref  =  0.0            # reference pressure

      sig   = (c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
      return sig
   ################################################################


   ################################################################
   def vol_melted(self,heat,matl='ice'):
      # heat=q=Q*dt (J/m^2)
      if matl == 'ice':
         vol   = heat/self.fusi # melted ice; <0 implies frozen
      elif matl == 'snow':
         vol   = max(0,heat/self.fuss) # melted snow; <0 implies frozen (impossible)
      elif matl == 'ice_snow':
         # for lat melt:
         # dv_lat=(hi+hs)*df=q/(hi*Fi+hs*Fs)
         fac   = (self.hice+self.hsnw)/(self.hice*self.fusi+self.hsnw*self.fuss) # m/J*m^2
         vol   = max(0,fac*heat) # m=m^3/m^2
      return vol
   ################################################################


   ################################################################
   def time_step(self,region='TP4'):
      # time step for thermodynamics (baclin from blkdat.input)

      if region=='TP4':
         baclin   = 400.0  # s   baroclinic time step
      elif region=='BS1':
         baclin   = 300.0  # s   baroclinic time step
      elif region=='FR1':
         baclin   = 300.0  # s   baroclinic time step

      return baclin
   ################################################################


   ################################################################
   def evolve_ice_lat(self,h,f,dv_lat):

      # lat melt/freeze
      df = dv_lat/h
      f2 = f+df
      h2 = h

      if f2<0:
         # all ice has melted
         return 0,0
      elif f2>self.fice_max:
         # ice covers grid cell
         dv_vrt   = (f2-self.fice_max)*h
         f2       = self.fice_max
         h2       = h+dv_vrt/f2

      return f2,h2
   ################################################################

   ################################################################
   def evolve_ice_vrt(self,h,f,dv_vrt):

      # vrt melt/freeze
      dh = dv_vrt/f
      h2 = h+dh
      f2 = f

      if h2<0:
         # all ice has melted
         return 0,0

      return f2,h2
   ################################################################
   

   ################################################################
   def freeze_vol_ice(self,dv,alp_lat):

      f        = self.fice
      h        = self.hice
      v        = f*h
      v_final  = v+dv
      if dv<0:
         raise ValueError('dv should be >=0')

      its   = 0
      print('iteration '+str(its)+'; error '+str(dv))
      while abs(dv)>1e-12:
         dv_lat   = alp_lat*dv
         df       = dv_lat/h 
         f        = f+df

         dv_vrt   = dv-dv_lat
         dh       = dv_vrt/f
         h        = h+dh

         v  = f*h
         dv = v_final-v

         its   = its+1
         print('iteration '+str(its)+'; error '+str(dv))

      if f>self.fice_max:
         print('fice too big')
         f        = self.fice_max
         v        = f*h
         dv_vrt   = v_final-v
         dh       = dv_vrt/f
         h        = h+dh

      v  = f*h
      dv = v_final-v # final error
      return f,h,dv
   ################################################################


   ################################################################
   def evolve_ice(self,heat_vrt,heat_lat,alp_lat):

      f  = self.fice
      h  = self.hice
      hs = self.hsnw
      v  = f*h

      heat_rest   = (heat_vrt+heat_lat)-(self.heat_ice+self.heat_snw)
      if heat_rest>0.:
         # everything melts, heat_rest to warm mixed layer
         return 0,0,heat_rest
      elif heat_vrt>self.heat_ice: 
         # only snow left
         # - sweep into ocean (heat from this is negligible?)
         return 0,0,heat_rest
      else:
         heat_rest   = 0.

      ################################################################
      # special cases
      # - no vert, or no lat
      dv_vrt   = -self.vol_melted(self,heat_vrt,matl='ice')
      if heat_lat==0 and heat_vrt==0:
         # no heat
         f  = self.fice
         h  = self.hice
         return f,h,heat_rest

      elif heat_lat==0:
         # only vert freezing
         f  = self.fice
         h  = self.hice+dv_vrt/f
         return f,h,heat_rest

      elif heat_vrt==0:
         # only lat freezing
         h  = self.hice
         if heat_lat<0:
            # freeze ice
            dv_lat   = -self.vol_melted(self,heat_lat,matl='ice')
            f        = min(self.fice_max,self.fice+dv_lat/h)
            dv       = f*h-self.vol_ice
            h        = h+(dv_lat-dv)/f
         else:
            # melt ice + snow
            dv_lat   = -self.vol_melted(self,heat_lat,matl='ice_snow')
            f        = max(0,self.fice+dv_lat/(h+self.hsnw))
            dv       = f*h-self.vol_ice
            h        = h+(dv_lat-dv)/f

         return f,h,heat_rest
      ################################################################


      if heat_lat>0:
         # need to melt both snow and ice
         dv_lat   = -self.vol_melted(self,heat_lat,matl='ice_snow')
      else:
         # need to freeze only ice
         dv_lat   = -self.vol_melted(self,heat_lat,matl='ice')


      ##################################################################
      # do vrt/lat changes in both orders, then take average:
      # -1- vrt,lat:
      h2a,f2a        = self.evolve_ice_vrt(h,f,dv_vrt)
      h3a,f3a        = self.evolve_ice_lat(h2a,f2a,dv_lat)

      # effective vrt vol change
      dv_vrt_eff_a   = dv_vrt+f3a*(h3a-h2a)
      Qvrt_eff_a     = (dv_vrt_eff_a*self.fusi)/dt/f            # m*J/m^3/s=J/s/m^2=W/m^2

      # effective lat volume change and flux
      # for effective fluxes (W/m^2)
      # - divide by surface area
      # - [lateral area] = [bottom area]*bet_lat
      dv_lat_eff_a   = h2a*(f3a-f2a)
      if dv_lat_eff_a<0:
         # melting ice & snow
         Qlat_eff_a     = (h2a*self.fusi+hs*self.fuss)*(f3a-f2a)/dt/(bet_lat*f)  # m*J/m^3/s=J/s/m^2=W/m^2
      else:
         # freezing ice
         Qlat_eff_a     = (h2a*self.fusi)*(f3a-f2a)/dt/(bet_lat*f)  # m*J/m^3/s=J/s/m^2=W/m^2
      ##################################################################

      ##################################################################
      # do vrt/lat changes in both orders, then take average:
      # -2- lat,vrt:
      h2b,f2b        = self.evolve_ice_lat(h,f,dv_lat)
      h3b,f3b        = self.evolve_ice_vrt(h2b,f2b,dv_vrt)

      # effective lat volume change and flux
      # for effective fluxes (W/m^2)
      # - divide by surface area
      # - [lateral area] = [bottom area]*bet_lat
      bet_lat        = alp_lat/(1-alp_lat)
      dv_lat_eff_b   = h*(f2b-f)
      if dv_lat_eff_b<0:
         # melting ice & snow
         Qlat_eff_b     = (h*self.fusi+hs*self.fuss)*(f2b-f)/dt/(bet_lat*f)  # m*J/m^3/s=J/s/m^2=W/m^2
      else:
         # freezing ice
         Qlat_eff_b     = (h*self.fusi)*(f2b-f)/dt/(bet_lat*f)  # m*J/m^3/s=J/s/m^2=W/m^2

      # effective vrt volume change and flux
      dv_vrt_eff_b   = f2b*(h2b-h)+f2b*(h3b-h2b)
      Qvrt_eff_b     = (dv_vrt_eff_b*self.fusi)/dt/f            # m*J/m^3/s=J/s/m^2=W/m^2
      ##################################################################


      ##################################################################
      # average:
      h2          = .5*(h2a+h2b)
      f2          = .5*(f2a+f2b)
      dv_vrt_eff  = .5*(dv_vrt_eff_a+dv_vrt_eff_b)
      dv_lat_eff  = .5*(dv_lat_eff_a+dv_lat_eff_b)
      Qvrt_eff    = .5*(Qvrt_eff_a+Qvrt_eff_b)
      Qlat_eff    = .5*(Qlat_eff_a+Qlat_eff_b)
      ################################################################


      return f2,h2,heat_rest,Qvrt_eff,Qlat_eff
   ###################################################################


   ###################################################################
   def split_flux(self,Qother,Qcool,Qatm,dt,alp_lat=0.):
      # Qother: vertical heat that shouldn't be converted to 
      #         lateral ones (conduction from top, extra heat from snow, penetration of light)
      # Qcool>0: heat into ice from cooling the mixed layer
      # Qatm>0: heat into open water through leads
      # NB Q* are in J/m^2; Q=flux*dt
      # TODO: use Dmax to get alp_lat

      # initial partition:
      q_dist = dt*(  Qcool+(1-self.fice)*Qatm            ) # heat J/m^2
      q_vrt  = dt*(  self.fice*Qother+(1-alp_lat)*Qdist  ) # heat J/m^2
      q_lat  = dt*(  alp_lat*Qdist                       ) # heat J/m^2

      return q_vrt,q_lat
   ###################################################################

###################################################################
class fsd_info_smooth:
   # fsd with continuous PDF (unlike RG method)
   # * PDF function is A*D^{-(1+alpha)} if D\in[Dmin,Dmax] else 0
   # * cumulative probability function is P(d>D) = (D^{-alpha}-D_max^{-alpha})/(D_min^{-alpha}-D_max^{-alpha})
   def __init__(self,Dmax,Dmin=20,fragility=.9,xi=2):
      self.Dmax      = Dmax
      self.Dmin      = Dmin
      self.fragility = fragility
      self.xi        = xi
      self.alpha     = np.log(xi*xi*f)/np.log(xi) #Exponent from Toyota
      A              = self.alpha/(pow(self.Dmin,-self.alpha)-pow(self.Dmax,-self.alpha))

      self.normalisation_constant   = A
      self.Dmean                    = self.moment(1)
      self.Dsq_mean                 = self.moment(2)

      return

   def moment(self,m):
      A     = self.normalisation_constant
      mom   = A/(self.alpha-m)*(pow(self.Dmin,m-self.alpha)-pow(self.Dmax,m-self.alpha))
      return mom

   def surf_area_frac(self,h):
      # for a square: lateral area is Nfloes*h*perim=Nfloes*h*(4*Dmean), bottom area is Nfloes*Dsq_mean
      # beta_lat  = S_lat/S_bot = 4*h*Dmean/Dsq_mean
      beta_lat    = 4*h*self.Dmean/self.Dsq_mean
      alpha_lat   = beta_lat/(1+beta_lat)
      return beta_lat,alpha_lat
###################################################################
