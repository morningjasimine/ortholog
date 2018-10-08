## Introduction
This project is develop for single copy orthologs  detection based on colinear of reference genome.

## pareration
1) dowload gtf, pep file of interesting species from Ensembl database
2) extract longest genes on each locus as representative genes of each species
3) and rename gene name with a unique species-specific prefix, like HOMO_ENSP00000220809
4) cat all fasta format protein file into one file, like all.pep.fa 
5) generate a config file as follow format:
```
	#prefix type    identity	coverage
	HOMO	ref	    0.5			0.5
	CITE	nonref	0.5			0.5
	RAT     nonref	0.5			0.5
	#--END--
```

## Run program step by step
#  blast
0) blastp all.pep.fa to all.pep.fa, and generate m8 format result.

#  filter blast result according to conditions in config file, and result file will be writen into 1_blastp directory
1) perl Step01.pl <m8 file> <pep file> <config>

#  detect single copy orthologs, and result file will be writen into 2_ortholog directory
2) perl Step02.pl <config> <gff dir>

test
