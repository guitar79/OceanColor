#!/usr/bin/env python

def myfunction (all_fnames):
	results = []
	for fname in all_fnames:
		if '2003' in fname: 	# if 2003 is found in fname
			results.append(fname)
		return results
		
def myprocedure(n):
	a,b = 0, 1
	for i in range(n):
		print(a)
		a,b = b,a+b
		