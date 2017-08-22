#!/usr/bin/python

from __future__ import division

from glob import glob
from os import system, path, mkdir
from math import ceil
from multiprocessing import Pool

DEBUG = False

def run_command(cmd):
    print("Running command:\n" + cmd)
    if DEBUG:
        return 0
    else:
        return system(cmd)

def main(options):
    PRIMER = options.adapter # "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG"
    QUAL = options.qual # 25 # Quality filter on 3' side
    NEXTSEQ = options.nextseq # True # well known problem with NextSeq data due to 2 color chemistry
                              # (https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases, http://seqanswers.com/forums/showthread.php?t=45130)
                              # Should be used together with QUAL parameter. This could be used instead of POLYG trimming.

    POLYA = options.polya3 # True # 3'
    POLYT = options.polyt3 # False # 3'
    POLYG = options.polyg3 # False # 3'
    POLYN = options.polyn3 # False # 3'

    POLY_MINBASES = options.polyminbases # 6
    MIN_LEN = options.minlen #20 # minimum read length after trimming

    CUTADAPT_BIN = options.cutadapt #"/home/parkj/.local/bin/cutadapt"
    OUT_DIR = options.OUT_DIR[0]

    if len(options.tmpdir) > 0:
        tmpdir = options.tmpdir
    else:
        tmpdir = path.join(OUT_DIR, "tmp")

    run_command("mkdir -p " + tmpdir)

    with open(options.read1) as f:
        read1 = filter(None, (line.rstrip() for line in f))

    if options.read2:
        with open(options.read2) as f:
            read2 = filter(None, (line.rstrip() for line in f))
        files = zip(read1, read2)
    else:
        files = [[fn] for fn in read1]

    cmd_lines = []
    fifos = []

    first_cutadapt = options.adapter or QUAL > 0
    second_cutadapt = POLYA or POLYT or POLYN or POLYG

    for fns in files:
        runid = fns[0].split('/')[-3]
        prefix = runid + "_" + "'.'.join(path.basename(fns[0]).split(".")[:-2])
        out_file = path.join(OUT_DIR, prefix + ".fastq.gz")
        fifo1 = path.join(tmpdir, prefix + ".fifo.fastq")
        fifos.append(fifo1)

        if options.read2:
            runid = fns[0].split('/')[-3]
            prefix = runid + "_" + '.'.join(path.basename(fns[1]).split(".")[:-2])
            out_file2 = path.join(OUT_DIR, prefix + ".fastq.gz")
            fifo2 = path.join(tmpdir, prefix + ".fifo.fastq")
            fifos.append(fifo2)

        cmds = [ CUTADAPT_BIN, ' '.join(fns) ]
        if first_cutadapt:
            if options.adapter:
                cmds.append("-a %s"%PRIMER)
                if options.read2:
                    cmds.append("-A %s"%PRIMER)
            if QUAL > 0:
                if NEXTSEQ:
                    cmds.append("--nextseq-trim=%d"%QUAL)
                else:
                    cmds.append("-q %d"%QUAL)
            if second_cutadapt:
                cmds.append("-o " + fifo1)
                if options.read2:
                    cmds.append("-p " + fifo2)

        if second_cutadapt:
            if first_cutadapt:
                cmds += ["&", CUTADAPT_BIN, fifo1]
                if options.read2:
                    cmds.append(fifo2)

            cmds.append("--overlap %d"%POLY_MINBASES)

            if POLYA:
                cmds.append("-a A{1000}")
            if POLYT:
                cmds.append("-a T{1000}")
            if POLYN:
                cmds.append("-a N{1000}")
            if POLYG:
                cmds.append("-a G{1000}")

        if MIN_LEN > 0:
            cmds.append("-m %d"%MIN_LEN)

        cmds.append("-o " + out_file)
        if options.read2:
            cmds.append("-p " + out_file2)

        cmd_lines.append(' '.join(cmds))

    if first_cutadapt and second_cutadapt:
        for fifo in fifos:
            run_command("mkfifo " + fifo)

    p = Pool(options.cores)
    status = p.map(run_command, cmd_lines)

    if first_cutadapt and second_cutadapt:
        for fifo in fifos:
            run_command("rm -f " + fifo)

    if any(status):
        print("Error: There was at least one error while running cutadapt.")


    run_command("rm -fr " + tmpdir)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run cutadapt with appropriate options.")

    parser.add_argument('-a', '--adapter', dest='adapter', type=str)
    parser.add_argument('-q', '--qual', dest='qual', type=int, default=25)
    parser.add_argument('-n', '--nextseq', dest='nextseq', action='store_true')
    parser.add_argument('-A3', '--polya3', dest='polya3', action='store_true')
    parser.add_argument('-T3', '--polyt3', dest='polyt3', action='store_true')
    parser.add_argument('-G3', '--polyg3', dest='polyg3', action='store_true')
    parser.add_argument('-N3', '--polyn3', dest='polyn3', action='store_true')
    parser.add_argument('-P', '--polyminbases', dest='polyminbases', type=int, default=6)
    parser.add_argument('-l', '--minlen', dest='minlen', type=int, default=20)
    parser.add_argument('-b', '--cutadapt', dest='cutadapt', type=str, default='cutadapt')
    #parser.add_argument('-j', '--jobs', dest='jobs', type=int, default=1)
    #parser.add_argument('-k', '--chunk', dest='chunk', type=int, default=0)
    parser.add_argument('-c', '--cores', dest='cores', type=int, default=1)
    parser.add_argument('-1', '--read1', dest='read1', type=str, required=True)
    parser.add_argument('-2', '--read2', dest='read2', type=str)
    parser.add_argument('-t', '--tmpdir', dest='tmpdir', type=str, default="")
    parser.add_argument('OUT_DIR', type=str, nargs=1)

    options = parser.parse_args()
    main(options)