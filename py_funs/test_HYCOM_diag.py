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
      mth   = (self.wavdir*wt_theta).transpose().dot(self.Sdir).dot(wt_omega)

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
         print(name+' : '+str(v0)+' | '+str(v1))

      print(' ')
      return
   #################################################################
