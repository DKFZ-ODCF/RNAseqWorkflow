#!/usr/bin/env python

# Parse raw barcode file from various single cell sequencer.
# author: Jeongbin Park (j.park@dkfz.de)
# usage: parse_barcodefile.py [-h] [-t [RETURNTYPE]] [-s [SAMPLE]] [-n [MAXCNT]] barcode_file [barcode_file ...]

from __future__ import print_function
from collections import defaultdict
import argparse, csv, sys, os

HEADERS = ["sample", "row", "col", "barcode"] # Order is important

def main(args):
    samples = {}
    rgs = []
    for fn in args.fns:
        with open(fn, 'r') as csvfile:
            is_tsv = '\t' in csvfile.readline()

        with open(fn, 'rb') as csvfile:
            r = csv.reader(csvfile, delimiter='\t' if is_tsv else ',')
            hdr_idx_dic = {}
            dbgcnt = 0
            for cnt, row in enumerate(r):
                if len(hdr_idx_dic) == 0:
                    found_headers = []
                    for idx, hdr in enumerate(row):
                        hdr = hdr.lower()
                        if hdr in HEADERS:
                            hdr_idx_dic[hdr] = idx
                    found_headers = filter(lambda e: e in hdr_idx_dic, HEADERS)
                    if not all(e in found_headers for e in ['sample', 'barcode']):
                        print("Error: not all headers present in the input file!", file=sys.stderr)
                        exit(1)
                else:
                    data = [row[hdr_idx_dic[hdr]] for hdr in found_headers]
                    sample_name = data[0].replace(" ", "_")
                    if not sample_name in samples:
                        samples[sample_name] = [len(samples)+1, 0]
                    samples[sample_name][1] += 1
                    if len(data) == 4:
                        fc = "%s-R%sC%s"%(sample_name, data[1], data[2])
                        barcode = data[3]
                    elif len(data) == 2:
                        fc = "%s-%02d"%(sample_name, cnt)
                        barcode = data[1]
                    if len(args.sample) == 0 or (len(args.sample) > 0 and args.sample == sample_name):
                        if args.returntype == 0:
                            fnhead = "%s-%s-C%d-H%03d_%s_R"%(os.environ["PID"], sample_name, samples[sample_name][0], samples[sample_name][1], barcode)
                            print('\t'.join([fc, barcode, fnhead+"1.fastq.gz", fnhead+"2.fastq.gz"]))
                            dbgcnt += 1
                        elif args.returntype == 1:
                            rgs.append("ID:{0}_C{1}_H{2:03d}_{3} LB:{4}_{5} PL:ILLUMINA SM:sample_{4}_{5} PU:{6}".format(os.environ["RUN_ID"], samples[sample_name][0], samples[sample_name][1], barcode, sample_name, os.environ["PID"], os.environ["RUN_ID"].split("_")[-1]))
                            dbgcnt += 1
                        if args.maxcnt == dbgcnt:
                            break
    if args.returntype == 1:
        print(' , '.join(rgs))
    elif args.returntype == 2:
        print('\n'.join(samples))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse raw barcode files (welllist, etc) from single cell sequencing platforms.")
    parser.add_argument("-t", "--returntype", dest="returntype", nargs="?", type=int, default=0,
                              help='Specify output type (0: barcode file for Je, 1: Read Groups, 2: Samples).')
    parser.add_argument("-s", "--sample", dest="sample", nargs="?", type=str, default="",
                              help='Print list only with specified sample. Neglected with return type 2.')
    parser.add_argument("-n", "--maxcount", dest="maxcnt", nargs="?", type=int, default=-1,
                              help='Limit number of returned entries per file (except for samples). Useful for debug purpose.')
    parser.add_argument('fns', metavar='barcode_file', type=str, nargs='+',
                              help='Raw barcode information file from sequencing machine.')
    args = parser.parse_args()
    main(args)
