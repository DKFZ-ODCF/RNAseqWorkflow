#!/usr/bin/env python

#######################################
# bsnp.py
# v1.1
# 3 nov 2015
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# extract allele frequencies from alignment file for in-silico genotyping (fingerprinting).
# returns base counts for ref and alt and betavalue-like allele frequency.
#
# automatically detects if contigs are named "chr1" or "1" in alignment file.
# in strand-specific WGBS sequencing data only the strand not affected by conversion is considered.

def mod_bsnp(sysargv):
	import sys
	import pysam
	import argparse
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[bsnp %Y-%m-%d %H:%M:%S]")
	
			
	#######################################
	# arguments, filehandles

	parser = argparse.ArgumentParser(prog="methylCtools bsnp", version="1.0", description="extract allele frequencies from alignment file for in-silico genotyping (fingerprinting).")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")
	
	parser.add_argument("-q", "--mapQ", metavar="INT", dest="mapQ", type=int, default=1, action="store", help="skip alignments with mapQ smaller than INT [1]")
	parser.add_argument("-Q", "--baseQ", metavar="INT", dest="baseQ", type=int, default=20, action="store", help="skip bases with baseQ smaller than INT (Phred+33) [20]")
	parser.add_argument("-c", "--cov", metavar="INT", dest="minCOV", type=int, default=0, action="store", help="skip positions with coverage smaller INT [0]")
	
	parser.add_argument("-r", "--splitRG", dest="allRG", action="store_true", help="split by read group")
	parser.add_argument("-t", "--trimPE", dest="trimPE", action="store_true", help="only consider first read in overlapping PE reads")
	parser.add_argument("-b", "--bs", dest="isBS", action="store_true", help="input is WGBS sequencing data")

	parser.add_argument("-n", "--name", metavar="STR", dest="sname", default="input.bam", action="store", help="sample name to be used in output file [input.bam]")	
	
	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inPOS", metavar="positions.bed", action="store", default=False, help="list of SNP positions to be assessed")
	groupinput.add_argument("inBAM", metavar="alignments.bam", action="store", default=False, help="alignment file (sorted, indexed)")
	
	args = parser.parse_args(sysargv)

	try:
		posIN = open(args.inPOS, "r")
		samfileIN = pysam.Samfile(args.inBAM, "rb")
		out = sys.stdout
	except IOError as strerror:
		sys.exit("methylCtools bsnp: error: %s" % strerror)

	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))


	######################################
	# prepare

	if args.sname == "input.bam": sname = args.inBAM.split("/")[-1]
	else: sname = args.sname

	h = samfileIN.header													# remove "chr" to match bam file, undo later
	#print h
	flagCHR = h["SQ"][1]["SN"][:3] == "chr"

	rgs = ["combined"]														# check read groups
	if args.allRG:
		if "RG" in h:
			for rg in h["RG"]: rgs += [rg["ID"]]
		if args.qf: sys.stderr.write("%s status: %i read groups in input file\n" % (nicetime(), len(rgs)-1))

	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

				
	######################################
	# main

	snps = dict()															# dict of all snp positions
	snpsv = []																# vector of all snp names
	c = dict()																# c[readgroup][snp] = [red, green]
	for rg in rgs: c[rg] = dict()
	
	n = 0
	for l in posIN:
		if l[0] == "#": continue											# process SNPs
		ls = l.rstrip().split("\t")
		snp = [ls[0], int(ls[2]), ls[5], ls[6], ls[7], ls[3]]
		if snp[2] == "-":													# convert to + strand (as in pileup)
			snp[2] = "+"
			snp[3] = complement[snp[3]]
			snp[4] = complement[snp[4]]
		if not flagCHR: snp[0] = snp[0][3:]
		
		snps[ls[3]] = snp
		snpsv += [ls[3]]
		for rg in rgs: c[rg][snp[-1]] = [0, 0]
		
		for pileupcolumn in samfileIN.pileup(snp[0], start=snp[1]-5, end=snp[1]+3):
			if pileupcolumn.pos == snp[1]-1:

				if args.trimPE:												# trim overlapping PE reads: check if qname occurs more than once
					once, twice = [], []
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.qname in once: twice += [pileupread.alignment.qname]
						else: once += [pileupread.alignment.qname]	

				for pileupread in pileupcolumn.pileups:
					if not pileupread.is_del:								# check if deletion at pos
						
						if args.allRG:										# get read group
							patags = dict(pileupread.alignment.tags)
							if "RG" in patags: rg = patags["RG"]
							else: rg = False
						else: rg = False
								
						if args.trimPE:										# trim overlapping PE reads
							if pileupread.alignment.qname in twice and pileupread.alignment.is_read2: continue
					
						if args.isBS:										# only use informative strand in BS mode (C could be converted to T)
							if snp[3] == "C" or snp[4] == "C":
								if pileupread.alignment.is_read1 != pileupread.alignment.is_reverse: continue
							elif snp[3] == "G" or snp[4] == "G":
								if pileupread.alignment.is_read1 == pileupread.alignment.is_reverse: continue
			
						if pileupread.alignment.qual[pileupread.qpos] < (args.baseQ+33): continue
						if pileupread.alignment.mapq < args.mapQ: continue	# check for minimum mapQ and baseQ
						
						if pileupread.alignment.seq[pileupread.qpos] == snp[3]:
							if rg: c[rg][snp[-1]][0] += 1					# count
							c["combined"][snp[-1]][0] += 1
						elif pileupread.alignment.seq[pileupread.qpos] == snp[4]:
							if rg: c[rg][snp[-1]][1] += 1
							c["combined"][snp[-1]][1] += 1
						
		n += 1
		if n%100 == 0:
			if args.qf: sys.stderr.write("%s status: %i positions processed\n" % (nicetime(), n))
	
	out.write("# bsnp v1.1 output\n")
	out.write("# call: %s\n" % " ".join(sys.argv))
	
	for s in snpsv:
		if not flagCHR: snps[s][0] = "chr" + snps[s][0]
	
	for rg in rgs:
		n = 0
		if rg == "combined": snamerg = sname
		else: snamerg = sname + ":" + rg
		for s in snpsv:
			snp = snps[s]
			if sum(c[rg][s]) >= args.minCOV:
				if c[rg][s][0]+c[rg][s][1] > 0: out.write("%s\t%s\t%i\t%s\t%s\t%i\t%i\t%.3f\n" % (snamerg, snp[0], snp[1], "*", snp[-1], c[rg][s][0], c[rg][s][1], c[rg][s][0]/float(c[rg][s][0]+c[rg][s][1])))
				else: out.write("%s\t%s\t%i\t%s\t%s\t%i\t%i\tNA\n" % (snamerg, snp[0], snp[1], "*", snp[-1], c[rg][s][0], c[rg][s][1]))
				n += 1
		if args.qf:
			if rg == "combined": sys.stderr.write("%s status: %i SNPs covered in all read groups (coverage >= %i)\n" % (nicetime(), n, args.minCOV))
			else: sys.stderr.write("%s status: %i SNPs covered in read group %s (coverage >= %i)\n" % (nicetime(), n, rg, args.minCOV))
		

if __name__ == "__main__":
    import sys
    mod_bsnp(sys.argv[1:])

