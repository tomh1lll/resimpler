#! /usr/bin/env python
# -*- coding: utf-8 -*-

"run as python estimate_SFS_stats.py $inputfile"
"Appends to the given output file"

import dendropy
from dendropy.calculate import popgenstat
import os
import sys

seqs = dendropy.DnaCharacterMatrix.get(path=sys.argv[1],schema="fasta")

out = open("popgen_stats.txt","a")

pop = sys.argv[1].split("/")[1].split("_")[2].split(".")[0]
gene = sys.argv[1].split("/")[1].split("_")[0]
fbtr = sys.argv[1].split("/")[1].split("_")[1]
td = popgenstat.tajimas_d(seqs)
tw = popgenstat.wattersons_theta(seqs)
tp = popgenstat.average_number_of_pairwise_differences(seqs)
ss = popgenstat.num_segregating_sites(seqs)

out.write("\t".join([str(pop),str(gene),str(fbtr),str(td),str(tw),str(tp),str(ss)]) + "\n")
