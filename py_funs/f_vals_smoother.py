def smoother(fvals_array):
	import numpy as np
	fvals = np.array(fvals_array, dtype=float)
	n = 0
	N = len(fvals)
	while n < N:
		en = fvals[n]
		if en != 2 and en != fvals[n-1] and en != fvals[n+1]:
			# this is done to delete singularities
			fvals[n] = fvals[n-1]
			n += 1
		elif en == 2:
			# if it's 2 measure the 2's string
			for m,em in enumerate(fvals[n:]):
				# check the end to decide the sign of the gradient
				if em == 1:
					end_check = 1
					break
				elif em == 0:
					end_check = 0
					break
			ncont = m + 1
			# calculate the unit for the gradient
			unit = 1/float(ncont)
			if fvals[n-1] == end_check:
				# in this case flatten the unknown to same value
				# generally polynias
				for l in range(ncont):
					fvals[n] = fvals[n-1]
					n += 1
			elif fvals[n-1] == 1:
				# negative gradient
				for l in range(ncont):
					fvals[n] = fvals[n-1] - unit
					if fvals[n] < 0.1:
						fvals[n] = 0
					n += 1
			elif fvals[n-1] == 0:
				# positive gradient
				for l in range(ncont):
					fvals[n] = fvals[n-1] + unit
					if fvals[n] > 0.9 and fvals[n] < 1:
						fvals[n] = 1
					n += 1
		else:
			# not a singularity nor a 2
			n += 1
	return(fvals)
