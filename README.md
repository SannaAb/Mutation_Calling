# Mutation_Calling
Nextflow pipeline for mutation Calling 

### Creating a reference db 

MutationCall_DNA_MakingReferenceSNPdb.nf 

A script that creates a reference database for "known" snps that you can use for the SNP calling. 

important: The currect version the script stops after hard filtering and sorting the snps. The processing of extracting snps shared in all samples for now needs to be made seperately. 

run it as follows 

```

/home/xabras/.conda/envs/Nextflow/bin/nextflow run /home/xabras/Scripts/Mutation_Calling/MutationCall_DNA_MakingReferenceSNPdb.nf --bam "../*.bam" --ref /jumbo/db/Staphylococcus_epidermidis_ATCC12228/Fasta/GCF_000007645.1_ASM764v1_genomic.fna --outdir /jumbo/WorkingDir/B19-057/Data/Meta/Alignment/Epidermidis/NextflowPipe

```

### SNP calling pipeline

MutationCall_DNA_pipe.nf

Script that process the alignment files and calls the snps 

```



```