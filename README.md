# Mutation_Calling
Nextflow pipeline for mutation Calling 

### Creating a reference db 

MutationCall_DNA_MakingReferenceSNPdb.nf 

A script that creates a reference database for "known" snps that you can use for the SNP calling. 

important: The currect version the script stops after hard filtering and sorting the snps. The processing of extracting snps shared in all samples for now needs to be made seperately. 

run it as follows 

```

~/Programs/Nextflow/nextflow run /home/xabras/Scripts/Mutation_Calling/MutationCall_DNA_MakingReferenceSNPdb.nf --bam "../*.bam" --ref /jumbo/db/Staphylococcus_epidermidis_ATCC12228/Fasta/GCF_000007645.1_ASM764v1_genomic.fna --outdir /jumbo/WorkingDir/B19-057/Data/Meta/Alignment/Epidermidis/NextflowPipe

```

### SNP calling pipeline DNA

MutationCall_DNA_pipe.nf

Script that process the alignment files and calls the snps 

```
~/Programs/Nextflow/nextflow run /home/xabras/Scripts/Mutation_Calling/MutationCall_DNA_pipe.nf --bam "../*.bam" --ref /jumbo/db/Staphylococcus_aureus_NCTC_8325/Fasta/GCF_000013425.1_ASM1342v1_genomic.fna --outdir /jumbo/WorkingDir/B19-057/Data/Meta/Alignment/Aureus/NextflowPipe_SNPCalling --known /jumbo/WorkingDir/B19-057/Data/Meta/Alignment/Aureus/SNPDB/Merged.Sorted.sharedSNPS.Creatingdb_Aureus.vcf

```


### SNP calling pipeline RNA 

MutationCall_RNA_pipe.nf

This is how you run it for example GRCH38 reference genome. The only thing that differs this script compared to the DNA is that you need to splitCigarN as the mutationcalling step cannot handle this yet. 

```

~/Programs/Nextflow/nextflow run /home/xabras/Scripts/Mutation_Calling/MutationCall_RNA_pipe.nf --bam "/jumbo/WorkingDir/B19-053/Data/Meta/SNP_calling/Two-passAlignment/*.bam" --ref /jumbo/db/Homo_sapiens/Ensembl/GRCh38.90/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.toplevel.canonical.fa --outdir /jumbo/WorkingDir/B19-053/Data/Meta/SNP_calling/NextFlow --known /jumbo/db/Homo_sapiens/Ensembl/GRCh38.90/Homo_sapiens_assembly38.dbsnp_withoutchr.vcf

```

### Important for Roy-Users

You need to specify a core that can handle qsub, for the moment not everynode can do this which will make the master worker crash (For a fact node 3 cannot qsub). Otherwise you can run it on the starting node but in that case you need to screen (and pray that the start node does not crash...). 

Make sure you have the newest nextflow installed, atleast version 19.07.0, there is a bug where the scripts might crash by random for the master worker in the older versions. 

### Additional Scripts 

CreateCladogramBasedOnVCF.py

This script creates a cladogram picture based on the vcf files that you have. You just specify which range you want to screen and which VCFs you want to use. One important note is that the script it fairly slow. I have tested it at 2M but it takes a couple of ours. 
Output is the cladogram figure but also the R code so you can check and plot the data manually. 

run it as follows

```
python CreateCladogramBasedOnVCF.py NC_007795.1:138161-140000 1.vcf 2.vcf 3.vcf 

```
