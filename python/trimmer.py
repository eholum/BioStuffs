"""
Trims primers from reads.
"""
import os

def trim(f, o, trimmed):
	
	infile = open(f)
	outfile = open(o,"w")
	
	line = "."
	
	while line:
		line = infile.readline()
		seq = infile.readline()
		
		outfile.write(line)
		
		seq = seq[19:]
		
		if not trimmed:
			seq = seq[0:-19]
			outfile.write(seq + "\n")
		else:
			outfile.write(seq)
		
	infile.close()
	outfile.close()
	
	del infile
	del outfile
			

def trim_files(inpath, outpath, numfiles):
	
	if not outpath[-1] == "/" :
		outpath = outpath + "/"
		
	fastqs = []
	for dirname, dirnames, filenames in os.walk(inpath):
		for name in filenames:
			if ".fastq" in name:
				fastqs += [os.path.join(dirname, name)]
		
	count = 1
	for f in fastqs:
	
		trimmed = "trimmed" in f
		
		temp = f.split("/")[-1][0:-6] + "_NOPRIMER.fastq"
		outfile = outpath + temp
		
		print "trimming: " + f + " into: " + outfile
		trim(f, outfile, trimmed)
		
		count += 1
		if (count > numfiles):
			return
	
	
if __name__ == "__main__":
	
	import sys
	argv = sys.argv
	
	inpath = argv[1]
	outpath = argv[2]
	
	numfiles = 10000
	
	if len(argv) >= 4:
		numfiles = int(argv[3])
	
	trim_files(inpath, outpath, numfiles)
	
	
