#!/usr/bin/env nextflow

params.bam = ''
params.ref = ''
params.known = '' 

sequences = Channel.fromPath( params.bam ).map { file -> tuple(file.baseName, file) }
sequencesintomutationalcalling = Channel.fromPath( params.bam ).map { file -> tuple(file.baseName, file) } // You cannot use the same input twice in a channel, therefore we are copying it in twice

ref_index = file(params.ref)
ref_known = file(params.known)

params.outdir = ''

//Something is of with the transfer of the files

//----------------------Samtools-------------------------------

process run_samtools {
        publishDir params.outdir, mode: 'copy', overwrite: true
        //errorStrategy 'ignore'

        clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node4-bfx.medair.lcl'

        module 'samtools/1.9'

        input:
	set file_ID, file(bamfile) from sequences
	       
	output:
	set file_ID, "${file_ID}.*" into samtools_out
 
        script:
        """
	samtools index ${bamfile}
	samtools flagstat ${bamfile} > ${file_ID}.flagstat
	samtools idxstats ${bamfile} > ${file_ID}.idxstat
        """
}

//----------------------AddReadsGroups-------------------------------

process run_addingReadsGroups {
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
	executor 'sge'
	queue 'bfxcore.q@node3-bfx.medair.lcl'
	queue 'bfxcore.q@node2-bfx.medair.lcl'

	input:
	set file_ID, file(bamfile) from sequencesintomutationalcalling

	output:
        set file_ID, "${file_ID}.With_RG.bam", "${file_ID}.With_RG.bai" into rg_out
	
	
	script:
	"""
	java -jar /apps/bio/apps/picard/2.1.0/picard.jar AddOrReplaceReadGroups I=${bamfile} O=${file_ID}.With_RG.bam SORT_ORDER=coordinate RGID=${file_ID} RGLB=${file_ID} RGPL=illumina RGPU=${file_ID} RGSM=${file_ID} CREATE_INDEX=True
	"""
}


//----------------------RemovePCRduplicates-------------------------------

process runRemovePCRdup { 
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
	executor 'sge'
	queue 'bfxcore.q@node3-bfx.medair.lcl'
	queue 'bfxcore.q@node2-bfx.medair.lcl'

	input:
        set file_ID, file(rgbam), file(rgbai) from rg_out

	output: 
	set file_ID, "${file_ID}.rmdup.bam","${file_ID}.rmdup.bai" into rmdup_out

	script: 
	"""
	java -Xmx4g -jar /apps/bio/apps/picard/2.1.0/picard.jar MarkDuplicates I=${rgbam} O=${file_ID}.rmdup.bam CREATE_INDEX=True M=${file_ID}.marked_dup_metrics.txt
	"""

}


//----------------------SplitNCigar-------------------------------

// SplitNcigar is important when having RNA seq data 
// Important to Note is that this is extremely memory consuming. Run it with many cores on your largest node! 

process runSplitNCigar { 
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 10'
	executor 'sge'
	queue 'bfxcore.q@node7-bfx.medair.lcl'
	queue 'bfxcore.q@node6-bfx.medair.lcl'

	input:
        set file_ID, file(rdupbam), file(rdupbai) from rmdup_out

	output: 
	set file_ID, "${file_ID}.Split.bam" into split_bam1, split_bam2

	script: 
	"""
	java -Xmx50g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${ref_index} -I ${rdupbam} -o ${file_ID}.Split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
	"""

}


//----------------------RealignTargetCreator-------------------------------

process run_realignTargetCreator {
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	
	clusterOptions='-pe mpi 4'
	executor 'sge'
	queue 'bfxcore.q@node7-bfx.medair.lcl'
	queue 'bfxcore.q@node6-bfx.medair.lcl'

	input:
	set file_ID, file(splitbam) from split_bam1
	
	output: 
	set file_ID, "${file_ID}.intervals" into intervals

	script:
	"""
	/apps/bio/apps/samtools/1.6/samtools index ${splitbam}

	java -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${splitbam} -o ${file_ID}.intervals -R ${ref_index} --known ${ref_known}
	"""
}


//----------------------IndelRealignment-------------------------------

process run_indelrealigner{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 4'
	executor 'sge'
	queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'

	input:
        set file_ID, file(interval) from intervals
	set file_ID, file(splitbam) from split_bam1

	output:
	set file_ID, "${file_ID}.realigned.bam" into realigned1, realigned2

	script:
	"""
	/apps/bio/apps/samtools/1.6/samtools index ${splitbam}

	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref_index} --targetIntervals ${interval} -known ${ref_known} -I ${splitbam} -o ${file_ID}.realigned.bam

	"""
}

//----------------------BaseRecalibration-------------------------------

process run_baserecalibrat{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'

	input:
	set file_ID, file(realbam) from realigned1

	output:
	set file_ID, "${file_ID}.grp" into grps

	script:
	"""
	
	/apps/bio/apps/samtools/1.6/samtools index ${realbam}

	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref_index} -knownSites ${ref_known} -I ${realbam} -o ${file_ID}.grp
	"""
}



//----------------------PrintReads-------------------------------

process run_PrintReads{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'

	input:
	set file_ID, file(group) from grps
	set file_ID, file(realbam) from realigned2

	output:
	set file_ID, "${file_ID}.recal.bam" into recal

	script:
	"""
	/apps/bio/apps/samtools/1.6/samtools index ${realbam}

	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T PrintReads -R ${ref_index} -I ${realbam} -BQSR ${group} -o ${file_ID}.recal.bam
	"""
}


//----------------------SNPCalling-------------------------------

process run_SNPcalling{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'

	input:
	set file_ID, file(recalbam) from recal

	output:
	set file_ID, "${file_ID}.vcf" into rawsnps

	script:
	"""
	/apps/bio/apps/samtools/1.6/samtools index ${recalbam}
	
	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref_index} -I ${recalbam} -o ${file_ID}.vcf
	"""
}

 
//----------------------HardFilter-------------------------------


process run_HardFilter{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'

	input:
	set file_ID, file(unfilteredvcf) from rawsnps

	output:
	set file_ID, "${file_ID}.flag.vcf" into hardfilteredsnps

	script:
	"""
	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref_index} --variant ${unfilteredvcf} -o ${file_ID}.flag.vcf --filterExpression "DP < 50" --filterName "LowDP" --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "MQ < 40.0" --filterName "MQ"
	"""
}

//----------------------PassedHardFilter-------------------------------

process run_greppassed{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'

	input:
	set file_ID, file(filteredvcf) from hardfilteredsnps

	output:
	set file_ID, "${file_ID}.Passed.vcf" into passedfilteredsnps

	script:
	"""
	grep -E '^#|PASS' ${filteredvcf} > ${file_ID}.Passed.vcf
	"""
}

//----------------------SortingVCF-------------------------------

process run_sortVCF{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'
	queue 'bfxcore.q@node7-bfx.medair.lcl'
	
	input:
	set file_ID, file(passedvcf) from passedfilteredsnps

	output:
	set file_ID, "${file_ID}.Sorted.vcf" into sortedvcf
	
	script:
	"""
	java -jar /apps/bio/apps/picard/2.1.0/picard.jar SortVcf I=${passedvcf} O=${file_ID}.Sorted.vcf
	"""
}

