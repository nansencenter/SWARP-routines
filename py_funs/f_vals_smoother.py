def smoother(fvals_array):
	import numpy as np
	if len(fvals_array) > 20:
		fvals = np.array(fvals_array, dtype=float)
		n = 0
		N = len(fvals)
		while n < N:
			en = fvals[n]
			if en == 2:
				for m,em in enumerate(fvals[n:]):
					if em == 1:
						end_check = 1
						break
					elif em == 0:
						end_check = 0
						break
				ncont = m + 1
				unit = 1/float(ncont)
				if fvals[n-1] == end_check:
					for l in range(ncont):
						fvals[n] = fvals[n-1]
						n += 1
				elif fvals[n-1] == 1:
					for l in range(ncont):
						fvals[n] = fvals[n-1] - unit
						n += 1
				elif fvals[n-1] == 0:
					for l in range(ncont):
						fvals[n] = fvals[n-1] + unit
						n += 1
			else:
				n += 1
	else:
		fvals=fvals_array
	return(fvals)
