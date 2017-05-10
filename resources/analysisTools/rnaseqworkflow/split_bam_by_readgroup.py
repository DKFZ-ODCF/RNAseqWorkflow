#!/usr/bin/env python

# Split BAM file by read group.
# author: Jeongbin Park (j.park@dkfz.de)
# usage: split_bam_by_readgroup.py BAMFILE [OUTDIR]

import sys, pysam, os

def main(args):
    outdir = args[2] if len(args) >= 3 else "."
    samfile = pysam.Samfile(args[1], "rb")

    rgs = {rg['ID']: rg for rg in samfile.header['RG']}

    outsamfiles = {}
    for read in samfile:
        rg = dict(read.tags)['RG']
        outfile = outsamfiles.get(rg, None)
        if outfile == None:
            outfilename = rgs[rg]['LB'] + "_" + '_'.join(rg.split("_")[-3:]) + "_merged.bam"
            outfile = outsamfiles[rg] = pysam.Samfile(os.path.join(outdir, outfilename), "wb", template=samfile)
            print(outfilename)
        outfile.write(read)

    for outsamfile in outsamfiles.values():
        outsamfile.close()
    samfile.close()

if __name__ == "__main__":
    main(sys.argv)