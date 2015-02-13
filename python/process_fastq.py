# Left side barcodes
lbc = {"AACCA":"V6L-0"}

# Right side barcodes
rbc = {"TGGTT":"V6R-0"}

samples = {"V6L-0V6R-0":"TE.D0B1"}

files = {}

import gc
from argparse import ArgumentParser
from Bio import SeqIO
from time import time

# If you don't need it delete it
def expand_maps():
    """
    Just add each primer string with a hamming distance of 1 to the barcode maps rather
    than computing it for each line/string. Faster this way.
    """
    for key in lbc.keys():
        value = lbc[key]
        
        for i in xrange(len(key)):
            lbc[key[0:i] + "N" + key[i+1:]] = value
            lbc[key[0:i] + "A" + key[i+1:]] = value
            lbc[key[0:i] + "C" + key[i+1:]] = value
            lbc[key[0:i] + "T" + key[i+1:]] = value
            lbc[key[0:i] + "G" + key[i+1:]] = value
            
    for key in rbc.keys():
        value = rbc[key]
        
        for i in xrange(len(key)):
            rbc[key[0:i] + "N" + key[i+1:]] = value
            rbc[key[0:i] + "A" + key[i+1:]] = value
            rbc[key[0:i] + "C" + key[i+1:]] = value
            rbc[key[0:i] + "T" + key[i+1:]] = value
            rbc[key[0:i] + "G" + key[i+1:]] = value
            
            
def prepare_files(o):
    """
    Creates a set of files for Demultiplexing.
    """
    for key in samples:
        name = samples[key]
        
        ######
        # RENAME THESE
        ######
        t = o + name + "_SOMETHING.fasta"
        files[name] = (t,[])
        
    files["NO_MATCH"] = (o + "NO_MATCH.fasta",[])
    files["LOW_QUALITY_READS"] = (o + "LOW_QUALITY_READS.fasta",[])
    
    for key in files:
        f = open(files[key][0], 'w')
        f.close()
    
def dump_files():
    """
    Dumps queued content into the appropriate file and garbage collects.
    """
    for key in files:
        t = files[key]
        f = open(t[0],'a')
        
        for line in t[1]:
            f.write(line)

        files[key] = (t[0],[])
        f.close()
    
    print gc.collect()


def get_sample_name(rec, filter_limit):
    """
    Returns the filename in which to dump the record
    """
    d = str(rec.seq)
    leftBC = d[0:5]
    rightBC = d[-5:]
    
    if not ((leftBC in lbc)    and (rightBC in rbc)):
        return "NO_MATCH"
    
    if compute_phred_score(rec) <= filter_limit:
        return "LOW_QUALITY_READS"
        
    ids = lbc[leftBC] + rbc[rightBC]
    
    return samples[ids]
 
 
def compute_phred_score(rec):
    """
    Computes the phred score of the record by determning the lowest score in the barcodes
    """
    phreds = rec.letter_annotations["phred_quality"]
    
    lps = phreds[0:5]
    rps = phreds[-5:]
    
    # Magic number
    low_score = 60
    
    for i in xrange(len(lps)):
        if lps[i] > 2:
            low_score = min(lps[i],low_score)
        if rps[i] > 2:
            low_score = min(rps[i],low_score)

    return low_score
        
        
def trim_adaptors(record, fadaptor, radaptor):
    """
    Removes the adaptors from the front and end of the read, if they exist.
    """
    indexf = record.seq.find(fadaptor)
    indexr = record.seq.find(radaptor)
        
    if indexf >= 0: 
        len_fadaptor = len(fadaptor)
        record = record[indexf+len_fadaptor:]
    if indexr >= 0:
        record = record[:indexr]

    return record
         

def process(infile, outdir, fadaptor, radaptor, filter_limit):
    """
    Process the fastq file. For each record (read), we do the following:
    
    1) Trim the adaptors
    2) Remove the barcodes
    3) Lookup the appropriate file, filter if necessary
    4) Queue the file for writing
    """
    
    print "reading from: " + infile
    print "writing to: " + outdir
    
    # If you only want perfect matches remove this line.
    expand_maps()
    ##########
    prepare_files(outdir)
    
    count = 0
    filtered = 0
    data = SeqIO.parse(infile,"fastq")
    start = time()
    
    for read in data:
        
        # Inefficient
        if count % 25000 == 0 and count > 1:
            print "processed {0} reads in {1} seconds ({2}% have been filtered)".format(\
                str(count), str(time() - start), str(((filtered+0.0)/count)*100)) 
        
        # Trim the adaptor
        read = trim_adaptors(read, fadaptor, radaptor)
        
        # Filter and find file
        f = get_sample_name(read, filter_limit)
        
        # Remove the barcodes?
        # Does something need to go here?
        
        # Inefficient
        if (f == "LOW_QUALITY_READS"):
            filtered += 1

        temp = files[f]
        print(format(read,"fasta"))
        temp[1].append(format(read,"fasta"))
        
        count += 1
        
        # Arbitrary... but seems to prevent memory overloads.
        if count % 500000 == 0:
            print "DUMPING FILES... "
            dump_files()
            
    print "analysis complete, successfully processed {0} reads... writing to files".format(str(count))
    dump_files()
    
    print "files written, filtered {0} reads {1}%".format(str(filtered),str(((filtered+0.0)/count)*100))
 
        
if __name__ == "__main__":
    
    parser = ArgumentParser(description="Process a fastq file according to gws.")
    parser.add_argument("infile", help="Input file. Should be fastq", type=str)
    parser.add_argument("front_adaptor", help="front adaptor for reads", type=str)
    parser.add_argument("rear_adaptor", help="rear adaptor for reads", type=str)
    parser.add_argument("-o", "--outdir", help="Output filename directory", type=str, default="")
    parser.add_argument("-f", "--filter_limit", help="filter limit", type=int, default=20)
    args = parser.parse_args()
    
    infile = args.infile
    outdir = args.outdir if args.outdir[-1] == "/" else args.outdir + "/"
    filter_limit = args.filter_limit
    fadaptor = args.front_adaptor
    radaptor = args.rear_adaptor
    
    process(infile, outdir, fadaptor, radaptor, filter_limit)