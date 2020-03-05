#!/usr/bin/env python3

import msprime, pyslim, sys, getopt

def main():
	ne = 50
	L = 3e7
	m = 1e-7
	r = 1e-8
	outname = "test.vcf"

	try:
		opts, args = getopt.getopt(sys.argv[1:], "hn:L:m:r:o:") # colon indicates argument needed for an option
	except getopt.GetoptError as err:
		print(err)
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('slim_makeburnin.py -n <Ne diploid> -L <length in bp> -m <mutation rate per site> -r <recombination rate per site> -o <outputfile>')
			sys.exit()
		elif opt == '-n':
			ne = int(arg)
		elif opt == '-L':
			L = int(arg)
		elif opt == '-m':
			m = float(arg)
		elif opt == '-r':
			r = float(arg)
		elif opt == '-o':
			outname = arg

	ts = msprime.simulate(sample_size=2*ne, Ne=ne, length=L, mutation_rate=0, recombination_rate=r) 
	ts = msprime.mutate(ts, rate=m, random_seed=1, model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))
	with open(outname
, 'w') as outfile:
		ts.write_vcf(outfile)

if __name__ == "__main__":
	main()

