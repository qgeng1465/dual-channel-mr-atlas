NR==FNR{keep[$1]=1; next}
FNR==1{print "SNP\tSNPChr\tSNPPos\tZscore\tAssessedAllele\tOtherAllele\tGene\tGeneChr\tGenePos\tNrSamples"; next}
($8 in keep || $9 in keep){print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$13}
