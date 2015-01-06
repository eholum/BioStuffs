"""
Format those damn files according to @dirtygeews specifications.
"""
def clean(f, o, s):
	
	lines = f.readlines()
	
	sample = ''
	
	for line in lines:
		if (line[0] == ">"):
			o.write((sample + "\n").replace(":",''))
			sample = ''
			o.write(line[0] + s + line[1:])
		else:
			sample = sample + line[0:-2]
			
	o.write((sample + "\n").replace("-",''))
	o.close()

if __name__ == "__main__":
	
	import sys
	argv = sys.argv
	
	f = open(argv[1])
	o = open(argv[2],'w')
	s = argv[3]
	
	clean(f,o,s)
