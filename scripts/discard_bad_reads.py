import os.path
import itertools
import pandas
import re
from collections import OrderedDict
import math
import sys, getopt

class FastqRead(object):
    def __init__(self, read):
        self.desc, self.seq, self.qual = read
    
    def __repr__(self):
        return self.desc + "\n" + self.seq + "\n+\n" + self.qual + "\n"

def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return zip(*args)

def parse_fastq(f):
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield desc, seq, qual

def main(argv):

  work_dir = ''

  try:
    opts, args = getopt.getopt(argv,"h:d:",["dir="])
  except getopt.GetoptError:
    print('To run: python discard_bad_reads.py -d <files_directory>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print('python discard_bad_reads.py -d <files_directory>')
      sys.exit()
    elif opt in ("-d", "--dir"):
      work_dir = arg

  R1 = os.path.join(work_dir, "Undetermined_S0_L001_R1_001.fastq")
  I1 = os.path.join(work_dir, "Undetermined_S0_L001_I1_001.fastq")
  R2 = os.path.join(work_dir, "Undetermined_S0_L001_R2_001.fastq")
  R_OUT = os.path.join(work_dir, "Undetermined_S0_L001_R1_001_trimmed.fastq")
  I_OUT = os.path.join(work_dir, "Undetermined_S0_L001_I1_001_trimmed.fastq")
  R2_OUT = os.path.join(work_dir, "Undetermined_S0_L001_R2_001_trimmed.fastq")

  I1_handle = open(I1)
  indices = (FastqRead(x) for x in parse_fastq(I1_handle))

  R1_handle = open(R1)
  fwds = (FastqRead(x) for x in parse_fastq(R1_handle))

  R2_handle = open(R2)
  revs = (FastqRead(x) for x in parse_fastq(R2_handle))

  P2rev = "GCATCGATGAAGAACGCAGCGGCTGACTGACT"
  P2rev = "GCATCGATGAAGAACGCAGC"

  counter = 0
  trimmed_seqs = 0
  no_rev = 0
  with open(R_OUT, "w") as r_out:
    with open(R2_OUT, "w") as r2_out:
      with open(I_OUT, "w") as i_out:
        for fwd, rev, idx in zip(fwds, revs, indices):
          if P2rev in fwd.seq:
            strindex = fwd.seq.index(P2rev)
            if strindex > 150:
              r_out.write("@%s\n%s\n+\n%s\n" % (fwd.desc, fwd.seq, fwd.qual))
              r2_out.write("@%s\n%s\n+\n%s\n" % (rev.desc, rev.seq, rev.qual))
              i_out.write("@%s\n%s\n+\n%s\n" % (idx.desc, idx.seq, idx.qual))
              counter += 1
            else:
              r_out.write("@%s\n%s\n+\n%s\n" % (fwd.desc, fwd.seq, fwd.qual))
              r2_out.write("@%s\n%s\n+\n%s\n" % (rev.desc, rev.seq, rev.qual))
              i_out.write("@%s\n%s\n+\n%s\n" % (idx.desc, idx.seq, idx.qual))
              no_rev += 1
          else:
            trimmed_seqs += 1

  I1_handle.close()
  R1_handle.close()

  print("Sequences with no rvs primers:", no_rev)
  print("Sequences with rvs primer and removed:", trimmed_seqs)
  print("Sequences with rvs primer and kept:", counter)
  print("Percent of reads with kept primers:", counter/(counter+trimmed_seqs+no_rev)*100, "%")

if __name__ == "__main__":
  main(sys.argv[1:])
