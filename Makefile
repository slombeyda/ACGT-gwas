ACGTtoPGM: ACGTtoPGM.cpp
	g++ -O3 -o ACGTtoPGM ACGTtoPGM.cpp

ACGT.sample.tsv:
	cat snp_cric_maf005_chr1.txt | awk -F '\t' '{s=""; for (i=1; i<=22; i++) { s=s$i"\t"; } for (i=NF-16; i<NF; i++) { s=s$i"\t"; } s=s$NF; print s;}' > ! ACGT.sample.tsv

clean:
	\rm -f ACGTtoPGM ACGT.sample.tsv
