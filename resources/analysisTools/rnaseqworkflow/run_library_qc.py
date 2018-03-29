# run_library_qc.py
# Purpose: run library QC for single cell RNA seq data. Produces statistics of 1) total reads, 2) number of genes that expresses at least 1 read, 3) Mitochondrial reads fraction, 4) RNA reads, 5) non-genic reads.
# Usage: run_library_qc.py sample pid featurecount_path outdir
import sys, os, tarfile
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

coding_types = "IG_C_gene;IG_D_gene;IG_J_gene;IG_LV_gene;IG_V_gene;TR_C_gene;" \
               "TR_J_gene;TR_V_gene;TR_D_gene;IG_pseudogene;IG_C_pseudogene;" \
               "IG_J_pseudogene;IG_V_pseudogene;TR_V_pseudogene;TR_J_pseudogene;" \
               "protein_coding;pseudogene;polymorphic_pseudogene"
# rrna_types = "rRNA;Mt_rRNA;rRNA_pseudogene"
rrna_types = "rRNA;rRNA_pseudogene"

model = os.environ.get("GENE_MODELS", "")
coding_types = os.environ.get("CODING_TYPES", coding_types).split(";")
if model == "":
    print("Error: environment variable GENE_MODELS is not set!")

coding_genes = []
mito_genes = []
rrna_genes = []
with open(model) as f:
    for line in f:
        if line[0] == "#": continue
        entries = line.rstrip().split("\t")
        if entries[2] == "gene":
            features = dict([i.strip().replace('"', '').split() for i in entries[-1].rstrip(';').split(";")])
            if features.get("gene_type", features["gene_biotype"]) in coding_types:
                coding_genes.append(features["gene_id"])
            if features.get("gene_type", features["gene_biotype"]) in rrna_types:
                rrna_genes.append(features["gene_id"])
            if entries[0].upper().lstrip("CHR") == "MT" or entries[0].upper().lstrip("CHR") == "M":
                mito_genes.append(features["gene_id"])

sample = sys.argv[1]
pid = sys.argv[2]
count_file = os.path.join(sys.argv[3], '_'.join([sample, pid, "featureCounts.count.tsv"]))
print "Processing count file..."
with open(count_file) as f:
    for line in f:
        if line[0] == "#":
            cell_ids = line.strip().split("\t")[2:]
            total_read_cnt = [0] * len(cell_ids)
            coding_gene_cnt = [0] * len(cell_ids)
            rrna_read_cnt = [0] * len(cell_ids)
            mito_read_cnt = [0] * len(cell_ids)
            non_genic_read_cnt = [0] * len(cell_ids)
            continue
        entries = line.strip().split("\t")
        gene_id = entries[0].split(".")[0]
        cnts = [int(e) for e in entries[2:]]
        for i, c in enumerate(cnts):
            total_read_cnt[i] += c
        if gene_id in mito_genes:
            for i, c in enumerate(cnts):
                mito_read_cnt[i] += c
        if gene_id in rrna_genes:
            for i, c in enumerate(cnts):
                rrna_read_cnt[i] += c
        if gene_id in coding_genes:
            for i, c in enumerate(cnts):
                if c > 0:
                    coding_gene_cnt[i] += 1

if len(sys.argv) > 5:
    ngcidx = 0
    fns = sys.argv[5:]
    for fn in fns:
        with open(fn) as f:
            for line in f:
                entries = line.strip().split()
                if entries[0] == "Unassigned_NoFeatures":
                    ngcnts = map(int, entries[1:])
                    for ngcnt in ngcnts:
                        non_genic_read_cnt[ngcidx] = ngcnt
                        ngcidx += 1
                    break

rrna_read_ratio = [0 if b == 0 else a/float(b)*100 for a, b in zip(rrna_read_cnt, map(sum, zip(total_read_cnt, non_genic_read_cnt)))]
mito_read_ratio = [0 if b == 0 else a/float(b)*100 for a, b in zip(mito_read_cnt, map(sum, zip(total_read_cnt, non_genic_read_cnt)))]
non_genic_read_ratio = [0 if b == 0 else a/float(b)*100 for a, b in zip(non_genic_read_cnt, map(sum, zip(total_read_cnt, non_genic_read_cnt)))]

with open(os.path.join(sys.argv[4], "qc_%s_%s.txt"%(sample, pid)), "w") as fo:
    fo.write('qc_metric\t' + '\t'.join(cell_ids) + '\n')
    fo.write('total_read_cnt\t' + '\t'.join(map(str, total_read_cnt)) + '\n')
    fo.write('coding_gene_cnt\t' + '\t'.join(map(str, coding_gene_cnt)) + '\n')
    fo.write('rrna_read_ratio\t' + '\t'.join(map(str, rrna_read_ratio)) + '\n')
    fo.write('mito_read_ratio\t' + '\t'.join(map(str, mito_read_ratio)) + '\n')
    fo.write('non_genic_read_ratio\t' + '\t'.join(map(str, non_genic_read_ratio)) + '\n')

def plot_counts(entries, cell_ids, title, axis, isratio):
    sorted_cell_ids, sorted_cnts = zip(*sorted(zip(cell_ids, entries), key=lambda e: e[1]))
    y_pos = np.arange(len(sorted_cell_ids))
    axis.bar(y_pos, sorted_cnts, align='center')
    axis.set_xticks([])
    axis.set_ylabel(title, fontsize=10.0)
    axis.set_xlim([-1, len(sorted_cnts)])
    axis.set_ylim([0, sorted_cnts[-1] + sorted_cnts[-1]*0.1])

fig, axes = plt.subplots(5, 1, figsize=(10, 15))
fig.suptitle("QC plot for " + sample + " of " + pid)
for idx, (title, cnts) in enumerate(zip(["Total reads", "Coding genes"], [total_read_cnt, coding_gene_cnt])):
    plot_counts(cnts, cell_ids, title, axes[idx], False)
for idx, (title, cnts) in enumerate(zip(["rRNA reads %", "Mito. reads %", "Non-genic reads %"], [rrna_read_ratio, mito_read_ratio, non_genic_read_ratio])):
    plot_counts(cnts, cell_ids, title, axes[idx+2], True)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(sys.argv[4], "qc_%s_%s.pdf"%(sample, pid)))
