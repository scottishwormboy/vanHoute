LIMS8810
Contains rough code for the anaysis of DMS3vir experiment. 
Outline as follows:
Download 2x 250bp MiSeq reads. 
Create high quality reads with flash
map with bwa mem to both host and phage
Find SNPs in DMS3vir
Pass to R
Select SNPs in 10bp region (8bp of 'seed' spacer + 2 PAM)
sum allele frequencies for each spacer in each population