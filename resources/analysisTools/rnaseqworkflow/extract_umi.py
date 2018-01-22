#!/usr/bin/env python

# Extract UMI from read 1 and append it to both reads. With PyPy, this script is ~5x faster than UMI-tools.
# author: Jeongbin Park (j.park@dkfz.de)
# input: read1, read2, outread1, outread2, pattern (same as UMI-tools)
# usage: extract_umi.py read1 read2 outread1 outread2 pattern

import sys
from datetime import datetime

def main(argv):
    fn1 = argv[1]
    fn2 = argv[2]
    fon1 = argv[3]
    fon2 = argv[4]
    pattern = argv[5]

    class Seq:
        def read(self, f):
            self.l = str(f.readline())
            if len(self.l) == 0:
                return False
            self.s = str(f.readline())
            self.p = str(f.readline())
            self.q = str(f.readline())
            return True

        def insert_tag(self, tag):
            bpos = self.l.find(' ')
            self.l = self.l[:bpos] + '_' + tag + self.l[bpos:]

        def extract_umi(self, umi):
            umiseq = self.s[umi[0]:umi[1]]
            self.s = self.s[:umi[0]] + self.s[umi[1]:]
            self.q = self.q[:umi[0]] + self.q[umi[1]:]
            return umiseq

        def write(self, fo):
            fo.write(self.l + self.s + self.p + self.q)

    rpos = 0
    nstart = -1
    umis = []
    for c in pattern:
        if c == 'N':
            if nstart == -1:
                nstart = rpos
        elif c == 'X':
            if nstart != -1:
                umis.append( (nstart, rpos+1, ) )
                nstart = -1
        else:
            print("Error!")
            exit()
        rpos += 1

    if nstart != -1:
        umis.append( (nstart, rpos, ) )

    print("Pattern: %s"%pattern)


    with open(fn1) as f1, open(fn2) as f2, open(fon1, "w+") as fo1, open(fon2, "w+") as fo2:
        s1 = Seq()
        s2 = Seq()

        proc_lines = 0
        while s1.read(f1) and s2.read(f2):
            for umi in umis:
                umiseq = s1.extract_umi(umi)
                s1.insert_tag(umiseq)
                s2.insert_tag(umiseq)
            s1.write(fo1)
            s2.write(fo2)

            proc_lines += 1
            if proc_lines%100000 == 0:
                t = datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[:-3]
                print("%s INFO Parsed %d reads"%(t, proc_lines))

if __name__ == "__main__":
    main(sys.argv)