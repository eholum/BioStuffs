from Bio import SeqIO
import gc
import os
import time

# Left side barcodes
lbc = {"CCTAAACTACGG":"F-01",
"TGCAGATCCAAC":"F-02",
"CCATCACATAGG":"F-03",
"GTGGTATGGGAG":"F-04",
"ACTTTAAGGGTG":"F-05",
"GAGCAACATCCT":"F-06",
"TGTTGCGTTTCT":"F-07",
"ATGTCCGACCAA":"F-08",
"AGGTACGCAATT":"F-09",
"ACAGCCACCCAT":"F-10",
"TGTCTCGCAAGC":"F-11",
"GAGGAGTAAAGC":"F-12",
"GTTACGTGGTTG":"F-13",
"TACCGCCTCGGA":"F-14",
"CGTAAGATGCCT":"F-15",
"TACCGGCTTGCA":"F-16",
"ATCTAGTGGCAA":"F-17",
"CCAGGGACTTCT":"F-18",
"CACCTTACCTTA":"F-19",
"ATAGTTAGGGCT":"F-20",
"GCACTTCATTTC":"F-21",
"TTAACTGGAAGC":"F-22",
"CGCGGTTACTAA":"F-23",
"GAGACTATATGC":"F-24",}

rbc = {"CCGTAGTTTAGG":"R-01",
"GTTGGATCTGCA":"R-02",
"CCTATGTGATGG":"R-03",
"CTCCCATACCAC":"R-04",
"CACCCTTAAAGT":"R-05",
"AGGATGTTGCTC":"R-06",
"AGAAACGCAACA":"R-07",
"TTGGTCGGACAT":"R-08",
"AATTGCGTACCT":"R-09",
"ATGGGTGGCTGT":"R-10",
"GCTTGCGAGACA":"R-11",
"GCTTTACTCCTC":"R-12",
"CAACCACGTAAC":"R-13",
"TCCGAGGCGGTA":"R-14",
"AGGCATCTTACG":"R-15",
"TGCAAGCCGGTA":"R-16",}

"""
# Right side barcodes
rbc = {"CCTAAACTACGG":"R-01",
"TGCAGATCCAAC":"R-02",
"CCATCACATAGG":"R-03",
"GTGGTATGGGAG":"R-04",
"ACTTTAAGGGTG":"R-05",
"GAGCAACATCCT":"R-06",
"TGTTGCGTTTCT":"R-07",
"ATGTCCGACCAA":"R-08",
"AGGTACGCAATT":"R-09",
"ACAGCCACCCAT":"R-10",
"TGTCTCGCAAGC":"R-11",
"GAGGAGTAAAGC":"R-12",
"GTTACGTGGTTG":"R-13",
"TACCGCCTCGGA":"R-14",
"CGTAAGATGCCT":"R-15",
"TACCGGCTTGCA":"R-16",}"""

samples = {
"F-01R-01":"XA1_T0",
"F-01R-02":"XA2_T0",
"F-01R-03":"XA3_T0",
"F-01R-04":"XA4_T0",
"F-01R-05":"XB1_T0",
"F-01R-06":"XB2_T0",
"F-01R-07":"XB3_T0",
"F-01R-08":"XB4_T0",
"F-02R-01":"XC1_T0",
"F-02R-02":"XC2_T0",
"F-02R-03":"XC3_T0",
"F-02R-04":"XC4_T0",
"F-02R-05":"DA1_T0",
"F-02R-06":"DA2_T0",
"F-02R-07":"DA3_T0",
"F-02R-08":"DA4_T0",
"F-03R-01":"DB1_T0",
"F-03R-02":"DB2_T0",
"F-03R-03":"DB3_T0",
"F-03R-04":"DB4_T0",
"F-03R-05":"DC1_T0",
"F-03R-06":"DC2_T0",
"F-03R-07":"DC3_T0",
"F-03R-08":"DC4_T0",
"F-04R-01":"QA1_T0",
"F-04R-02":"QA2_T0",
"F-04R-03":"QA3_T0",
"F-04R-04":"QA4_T0",
"F-04R-05":"QB1_T0",
"F-04R-06":"QB2_T0",
"F-04R-07":"QB3_T0",
"F-04R-08":"QB4_T0",
"F-05R-01":"QC1_T0",
"F-05R-02":"QC2_T0",
"F-05R-03":"QC3_T0",
"F-05R-04":"QC4_T0",
"F-05R-05":"DQA1_T0",
"F-05R-06":"DQA2_T0",
"F-05R-07":"DQA3_T0",
"F-05R-08":"DQA4_T0",
"F-06R-01":"DQB1_T0",
"F-06R-02":"DBQ2_T0",
"F-06R-03":"DBQ3_T0",
"F-06R-04":"DQB4_T0",
"F-06R-05":"DQC1_T0",
"F-06R-06":"DQC2_T0",
"F-06R-07":"DQC3_T0",
"F-06R-08":"DQC4_T0",
"F-07R-01":"XA1_T1",
"F-07R-02":"XA2_T1",
"F-07R-03":"XA3_T1",
"F-07R-04":"XA4_T1",
"F-07R-05":"XB1_T1",
"F-07R-06":"XB2_T1",
"F-07R-07":"XB3_T1",
"F-07R-08":"XB4_T1",
"F-08R-01":"XC1_T1",
"F-08R-02":"XC2_T1",
"F-08R-03":"XC3_T1",
"F-08R-04":"XC4_T1",
"F-08R-05":"DA1_T1",
"F-08R-06":"DA2_T1",
"F-08R-07":"DA3_T1",
"F-08R-08":"DA4_T1",
"F-09R-01":"DB1_T1",
"F-09R-02":"DB2_T1",
"F-09R-03":"DB3_T1",
"F-09R-04":"DB4_T1",
"F-09R-05":"DC1_T1",
"F-09R-06":"DC2_T1",
"F-09R-07":"DC3_T1",
"F-09R-08":"DC4_T1",
"F-10R-01":"QA1_T1",
"F-10R-02":"QA2_T1",
"F-10R-03":"QA3_T1",
"F-10R-04":"QA4_T1",
"F-10R-05":"QB1_T1",
"F-10R-06":"QB2_T1",
"F-10R-07":"QB3_T1",
"F-10R-08":"QB4_T1",
"F-11R-01":"QC1_T1",
"F-11R-02":"QC2_T1",
"F-11R-03":"QC3_T1",
"F-11R-04":"QC4_T1",
"F-11R-05":"DQA1_T1",
"F-11R-06":"DQA2_T1",
"F-11R-07":"DQA3_T1",
"F-11R-08":"DQA4_T1",
"F-12R-01":"DQB1_T1",
"F-12R-02":"DBQ2_T1",
"F-12R-03":"DBQ3_T1",
"F-12R-04":"DQB4_T1",
"F-12R-05":"DQC1_T1",
"F-12R-06":"DQC2_T1",
"F-12R-07":"DQC3_T1",
"F-12R-08":"DQC4_T1",
"F-13R-01":"XA1_T2",
"F-13R-02":"XA2_T2",
"F-13R-03":"XA3_T2",
"F-13R-04":"XA4_T2",
"F-13R-05":"XB1_T2",
"F-13R-06":"XB2_T2",
"F-13R-07":"XB3_T2",
"F-13R-08":"XB4_T2",
"F-14R-01":"XC1_T2",
"F-14R-02":"XC2_T2",
"F-14R-03":"XC3_T2",
"F-14R-04":"XC4_T2",
"F-14R-05":"DA1_T2",
"F-14R-06":"DA2_T2",
"F-14R-07":"DA3_T2",
"F-14R-08":"DA4_T2",
"F-15R-01":"DB1_T2",
"F-15R-02":"DB2_T2",
"F-15R-03":"DB3_T2",
"F-15R-04":"DB4_T2",
"F-15R-05":"DC1_T2",
"F-15R-06":"DC2_T2",
"F-15R-07":"DC3_T2",
"F-15R-08":"DC4_T2",
"F-16R-01":"QA1_T2",
"F-16R-02":"QA2_T2",
"F-16R-03":"QA3_T2",
"F-16R-04":"QA4_T2",
"F-16R-05":"QB1_T2",
"F-16R-06":"QB2_T2",
"F-16R-07":"QB3_T2",
"F-16R-08":"QB4_T2",
"F-17R-01":"QC1_T2",
"F-17R-02":"QC2_T2",
"F-17R-03":"QC3_T2",
"F-17R-04":"QC4_T2",
"F-17R-05":"DQA1_T2",
"F-17R-06":"DQA2_T2",
"F-17R-07":"DQA3_T2",
"F-17R-08":"DQA4_T2",
"F-18R-01":"DQB1_T2",
"F-18R-02":"DBQ2_T2",
"F-18R-03":"DBQ3_T2",
"F-18R-04":"DQB4_T2",
"F-18R-05":"DQC1_T2",
"F-18R-06":"DQC2_T2",
"F-18R-07":"DQC3_T2",
"F-18R-08":"DQC4_T2",
"F-19R-01":"BacHo_DA",
"F-19R-02":"BacHo_DB",
"F-19R-03":"BacHo_DC",
"F-19R-04":"BacHo_DQA",
"F-19R-05":"BacHo_DQB",
"F-19R-06":"BacHo_DQC",
"F-19R-07":"H2_Metatrans14",
"F-19R-08":"H3_Metatrans14",
"F-20R-01":"H4_Metatrans14",
"F-20R-02":"H5_Metatrans14",
"F-20R-03":"H6_Metatrans14",
"F-20R-04":"H7_Metatrans14",
"F-20R-05":"H8_Metatrans14",
"F-20R-06":"D2_Metatrans14",
"F-20R-07":"D3_Metatrans14",
"F-20R-08":"D4_Metatrans14",
"F-21R-01":"D5_Metatrans14",
"F-21R-02":"D6_Metatrans14",
"F-21R-03":"D7_Metatrans14",
"F-21R-04":"D8_Metatrans14",
"F-21R-05":"PE1_doseH",
"F-21R-06":"PE1_doseD",
"F-21R-07":"PE2_doseH",
"F-21R-08":"PE2_doseD",
"F-22R-01":"619_D_Homo",
"F-22R-02":"619_H_Homo",
"F-22R-03":"R_CK4_G4_F1",
"F-22R-04":"R_CK4_G6_F2",
"F-22R-05":"R_POPA_G2_F1",
"F-22R-06":"R_POPA_G3_F2",
"F-22R-07":"S_CK6_G1_F1",
"F-22R-08":"S_POPA_G6_F2",
"F-23R-01":"S_POPA_G7_F2",
"F-23R-02":"623_D_Homo",
"F-23R-03":"623_H_Homo",
"F-23R-04":"R_CK4_G4_F2",
"F-23R-05":"R_CK4_G6_F3",
"F-23R-06":"R_POPA_G2_F2",
"F-23R-07":"R_POPA_G3_F3",
"F-23R-08":"S_CK6_G1_F2",
"F-24R-01":"S_CK6_G7_F4",
"F-24R-02":"S_POPA_G6_F3",
"F-24R-03":"S_POPA_G7_F5",
"F-24R-04":"R_CK4_G4_F3",
"F-24R-05":"R_CK4_G6_F5",
"F-24R-06":"R_POPA_G3_F4",
"F-24R-07":"S_CK6_G1_F3",
"F-24R-08":"S_CK6_G7_F5",}

# Can also be passed as arguments
fadaptor = "ACGCTCTTCCGATCT"
radaptor = "AGATCGGAAGAGC"

def process(left, right, outdir, filter_limit, filter_length):
    """
    Process the fastq file. For each record (read), we do the following:
    
    1) Trim the adaptors
    2) Remove the barcodes
    3) Lookup the appropriate file, filter if necessary
    4) Queue the file for writing
    """
    
    print("reading from: ")
    print(("left: " + left))
    print("right: " + right)
    print("writing to: " + outdir)
    print("limit: " + str(filter_limit))
    
    
    # Create a set of file names, a list of lines to be written, and line counts
    left_files = {}
    right_files = {}

    prepare_files("left", o, left_files)
    prepare_files("right", o, right_files)
    
    count = 0
    left_filtered = 0
    right_filtered = 0
    ftrimmed = 0
    rtrimmed = 0
    
    right_primer_dimers = 0
    left_primer_dimers = 0
    
    left_data = SeqIO.parse(left,"fastq")
    right_data = SeqIO.parse(right,"fastq")
    
    start = time.time()
    
    left_read = next(left_data);
    right_read = next(right_data);
    while (left_read or right_read):
        
        if (not left_read) or (not right_read):
            raise Exception("This isn't right? More reads in one that the other")
       
        count += 1
        if count % 25000 == 1 and count > 1:
            print("processed {0} reads in {1} seconds ({2}% have been filtered from left, {3}% have been filtered from right)".format(\
                str(count),\
                str(time.time() - start),\
                str(((left_filtered+0.0)/count)*100),\
                str(((right_filtered+0.0)/count)*100)))
            print("trimmed {0} front {1} rear".format(ftrimmed, rtrimmed))
        
        # Trim the adaptor
        left_read, fcheck, rt = trim_adaptors(left_read)
        left_primer_dimers += rt
            
        right_read, ft, rcheck = trim_adaptors(right_read)
        right_primer_dimers += ft
            
        ftrimmed += fcheck
        rtrimmed += rcheck
        
        # Filter and find file
        sample_name = get_sample_name(left_read, right_read)
        left_sample_name = sample_name
        right_sample_name = sample_name
        
        # Remove barcodes
        left_read = left_read[12:]
        right_read = right_read[:-12]
    
        # Also remove the primers
        left_read = left_read[20:]
        right_read = right_read[:-17]
        
        # Trim the reads
        left_read = trim_left_read(left_read, filter_limit)
        right_read = trim_right_read(right_read, filter_limit)
        
        # If the reads are too short throw them away
        if (len(left_read) < filter_length):
            left_sample_name = "LOW_QUALITY_READS"
            left_filtered += 1
        if (len(right_read) < filter_length):
            right_sample_name = "LOW_QUALITY_READS"
            right_filtered += 1

        # format and append
        format_lines(left_sample_name, right_sample_name, left_read, right_read, left_files, right_files)
        
        # Arbitrary... but seems to prevent memory overloads.
        if count % 500000 == 0:
            print("DUMPING FILES... ")
            dump_files(left_files)
            dump_files(right_files)
            
        try:
            left_read = next(left_data)
            right_read = next(right_data)
        except:
            break
            
    print("analysis complete, successfully processed {0} reads... writing to files".format(str(count)))
    dump_files(left_files)
    dump_files(right_files)
        
    print("trimmed {0} front {1} rear".format(ftrimmed, rtrimmed))
    print("primer dimers: {0} left, {1} right".format(left_primer_dimers, right_primer_dimers))
    
    print("files written, {0} reads processed, left filtered {1}, right filtered {2}".format(\
                str(count),str(left_filtered),str(right_filtered)))
    
 
def prepare_files(direction, outputfolder, file_map):
    """
    Creates a set of files for Demultiplexing.
    """
    loc = outputfolder + direction + "/"
    
    # create the folder directory if it doesn't exist
    try:
        os.stat(outputfolder)
    except:
        os.mkdir(outputfolder)
        
    try:
        os.stat(loc)
    except:
        os.mkdir(loc)
        
    for key in samples:
        name = samples[key]
        t = loc + name + ".fasta"
        
        # File name, lines list, lines written
        file_map[name] = [t,[], 1]
        
    file_map["NO_MATCH"] = [loc + "NO_MATCH.fasta", [], 1]
    file_map["LOW_QUALITY_READS"] = [loc + "LOW_QUALITY_READS.fasta", [], 1]
    
    for key in file_map:
        f = open(file_map[key][0], 'w')
        f.close()


def dump_files(file_map):
    """
    Dumps queued content into the appropriate file and garbage collects.
    """
    for key in file_map:
        t = file_map[key]
        f = open(t[0],'a')
        
        for line in t[1]:
            f.write(line)

        file_map[key][1] = []
        f.close()
    
    gc.collect()
    

def get_trim_index(rec):
    """
    Computes the phred score of the record 
    """
    phreds = rec.letter_annotations["phred_quality"]
    
    # Maximum possible score
    low_score = 60
    
    for i in range(len(phreds)):
        low_score = min(phreds[i], low_score)

    return low_score


def get_sample_name(leftread, rightread):
    """
    Returns the filename in which to dump the record
    """
    leftBC = leftread.seq[0:12]
    rightBC = rightread.seq[-12:]
    
    if not ((leftBC in lbc) and (rightBC in rbc)):
        return "NO_MATCH"
        
    ids = lbc[leftBC] + rbc[rightBC]
    
    if not ids in samples:
        return "NO_MATCH"
        
    return samples[ids]


def trim_adaptors(record):
    """
    Removes the adaptors from the front and end of the read, if they exist.
    """
    indexf = record.seq.find(fadaptor)
    indexr = record.seq.find(radaptor)
    
    f = 0
    r = 0
    
    # Trim the rear adaptor first if it's there, obvi
    if indexr >= 0:
        record = record[:indexr]
        r = 1
    if indexf >= 0: 
        len_fadaptor = len(fadaptor)
        record = record[indexf+len_fadaptor:]
        f = 1

    return [record, f, r]


def format_lines(left_sample_name, right_sample_name, left_read, right_read, left_files, right_files):
    """
    Converts a fastq line into a fasta line, removes the adapters, and appends the header name to
    the line.
    """
    # Link to files
    left_tupl = left_files[left_sample_name]
    right_tupl = right_files[right_sample_name]
    
    left_line = format(left_read, "fasta")
    left_line = "{0}{1}_{2}:{3}".format(left_line[0], left_sample_name, left_tupl[2], left_line[1:])
    left_tupl[1].append(left_line)
    left_tupl[2] += 1
    
    right_line = format(right_read, "fasta")
    right_line = "{0}{1}_{2}:{3}".format(right_line[0], right_sample_name, right_tupl[2], right_line[1:])
    right_tupl[1].append(right_line)
    right_tupl[2] += 1
    
    
def trim_left_read(left_read, filter_limit):  
    
    phreds = left_read.letter_annotations["phred_quality"]
    
    for i in range(len(phreds)):
        if (phreds[i] < filter_limit):
            return left_read[:i]

    return left_read


def trim_right_read(right_read, filter_limit):
    
    phreds = right_read.letter_annotations["phred_quality"]
    
    for i in range(len(phreds) - 1, 0, -1):
        if (phreds[i] < filter_limit):
            return right_read[i:]

    return right_read
    
# 
# BECCA SET THESE:
#

# output directory
o = "test_out/"

# Left/right files
left = "R2_test.fastq"
right = "R1_RevComp_test.fastq"

# Filter limit
filter_limit = 20

# Seq length limit
seq_limit = 100

process(left, right, o, filter_limit, seq_limit)