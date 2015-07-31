def smoother(fvals_array):
	import numpy as np
	fvals = np.array(fvals_array,dtype=float)
	n = 0
	N = len(fvals)
	if 1 not in fvals and 2 in fvals:
		fvals[fvals==2] = 1
	if 0 not in fvals and 2 in fvals:
		fvals[fvals==2] = 0
	if len(fvals[fvals==1])!=0 and len(fvals[fvals==0])!=0:
		while n < N:
			en = fvals[n]
			if (en == 0 or en == 1) and en != fvals[n-1] and (en != fvals[n+1] if n < N-1 else en != fvals[0]):
				# this is done to delete singularities
				fvals[n] = fvals[n-1]
				n += 1
			elif en == 2:
				# if it's 2 measure the 2's string
				for m,em in enumerate(fvals[n:]):
					# check the end to decide the sign of the gradient
					if fvals[m] == 1:
						end_check = 1
						break
					elif fvals[m] == 0:
						end_check = 0
						break
					if m == N:
						dbl = 1
						end_check = fvals[dbl]
						while end_check != 0 or end_check != 1:
							dbl += 1
							end_check = fvals[dbl]
							if dbl >= N:
								fvals[fvals==2] = .5
								return(fvals)
				ncont = m + 1
				while fvals[n-1] == 2:
					n -= 1
					ncont += 1
				# calculate the unit for the gradient
				unit = 1/float(ncont)
				if fvals[n-1] == end_check:
					# in this case flatten the unknown to same value
					# generally polynias
					for l in range(ncont):
						if n < 0:
							fvals[N+n] = fvals[N+n-1]
						else:
							fvals[n] = fvals[n-1]
						n += 1
				elif fvals[n-1] == 1:
					# negative gradient
					for l in range(ncont):
						if n < 0:
							fvals[N+n] = fvals[N+n-1] - unit
						else:
							fvals[n] = fvals[n-1] - unit
						if fvals[n] < 0.05:
							fvals[n] = 0
						n += 1
				elif fvals[n-1] == 0:
					# positive gradient
					for l in range(ncont):
						if n < 0:
							fvals[N+n] = fvals[N+n-1] + unit
						else:
							fvals[n] = fvals[n-1] + unit
						if fvals[n] > 0.95:
							fvals[n] = 1
						n += 1
			else:
				# not a singularity nor a 2
				n += 1
	else:
		if len(fvals[fvals==1]) == 0:
			fvals[fvals==2] = 1
		elif len(fvals[fvals==0]) == 0:
			fvals[fvals==2] = 0
		else:
			for n in range(N/2):
				fvals[n] = 1
			fvals[fvals!=1] = 0	
	return(fvals)

################################################################################################
class pca_mapper:

   def __init__(self,xy_coords):
      import numpy as np
      #
      x,y   = np.array(xy_coords).transpose()
      #
      self.x0 = np.mean(x)
      self.y0 = np.mean(y)
      self.x  = x
      self.y  = y
      #
      xy_rel   = np.array([x-self.x0,y-self.y0]) # 2xN matrix
      cov      = xy_rel.dot(xy_rel.transpose()) # covariance (2x2 matrix)
      #
      self.evals,self.evecs   = np.linalg.eig(cov)
      self.X,self.Y           = self.map(x,y,inverse=False)

      return
      
   def map(self,x,y,inverse=False):

      import numpy as np

      if not inverse:
         # in:  basemap coordinates
         # out: coordinates relative to principal components
         xy_rel   = np.array([x-self.x0,y-self.y0]) # 2xN matrix
         X,Y      = self.evecs.transpose().dot(xy_rel)
      else:
         # in:  coordinates relative to principal components
         # out: basemap coordinates
         xy    = np.array([x,y])
         X,Y   = self.evecs.dot(xy)
         X     = X+self.x0
         Y     = Y+self.y0

      return X,Y

   def set_func_vals(self):

      import geometry_planar as GP
      import numpy as np

      Nc       = len(self.X)
      coords   = np.array([self.X,self.Y]).transpose()
      ss       = GP.arc_length(coords,closed=True)
      P        = ss[-1] # perimeter
      ss       = ss[:-1] # drop last val
      #
      nvec     = range(Nc)
      fvals    = 0*ss

      # longest directions
      # - these can be the 2 zeros of a sine function
      Xsort = sorted([(x_,i) for i,x_ in enumerate(self.X)])
      i0    = Xsort[0][1]
      i1    = Xsort[-1][1]

      # orientation doesn't matter
      # - swap indices if inconvenient
      if i1<i0:
         i0,i1 = i1,i0

      # 1st half of polygon
      s_top       = ss[i0:i1]
      ntop        = range(i0,i1)
      L_top       = ss[i1]-ss[i0]
      fvals[ntop] = np.sin((np.pi/L_top)*(s_top-s_top[0]))

      # 2nd half of polygon
      s_bot = list(ss[i1:])
      nbot  = range(i1,Nc)
      s_bot.extend(list(P+ss[:i0]))
      nbot.extend(range(i0))
      L_bot       = ss[i0]+P-ss[i1]
      fvals[nbot] = np.sin(np.pi+(np.pi/L_bot)*(s_bot-s_bot[0]))

      return fvals
################################################################################################
