"""
Reads and filters a fasta file for read lengths
Params:
infile outfile (optional)minlength (optional)maxlength
example:
"python filter.py infile.fasta outfile.fasta 50 70" 
will read in infile.fasta and write all lines between 50 and 70 to outfile.fasta
"""
def filter(infile, outfile, minlength, maxlength):
    
    print "reading from " + infile.name
    print "keeping lines of length in range " + str((minlength,maxlength))

    line1 = infile.readline()

    while line1:
        if line1[0] == ">":
            line2 = infile.readline()
            if minlength <= len(line2) <= maxlength:
                outfile.write(line1)
                outfile.write(line2)

        line1 = infile.readline()

    outfile.close()

    print "results written to " + outfile.name

    return outfile


if __name__ == "__main__":

    import sys

    infile = open(sys.argv[1],"r")
    outfile = open(sys.argv[2],"w")

    # default values if they aren't given
    minlength = 50
    maxlength = 70

    if len(sys.argv) > 3:
        minlength = int(sys.argv[3])
    if len(sys.argv) > 4:
        maxlength = int(sys.argv[4])

    filter(infile, outfile, minlength, maxlength)
