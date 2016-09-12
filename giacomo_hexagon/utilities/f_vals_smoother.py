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

