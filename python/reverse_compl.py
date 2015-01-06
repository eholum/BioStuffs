d = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}

"""
Simple script to compute reverse complements for raw tufts data reads.
"""
def complement(string):
	s = ""
	for i in range(len(string)-1, -1, -1):
		s += d[string[i]]
	return s

def reverse_compl(lines, o):
	string = ""
	for i in range(len(lines)):
		if lines[i][0] == '>':
			o.write(complement(string) + "\n")
			string = ""
			o.write(lines[i])
		else:
			string += lines[i][0:-1]
             
	o.write(complement(string) + "\n")
	
if __name__=="__main__":
	
	import sys
		
	infile = sys.argv[1]
	outfile = sys.argv[2]
	
	lines = open(infile).readlines()
	o = open(outfile, 'w')
	reverse_compl(lines, o)
