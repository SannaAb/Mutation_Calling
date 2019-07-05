#!/usr/bin/env nextflow

params.bam = ''
params.ref = ''

sequences = Channel.fromPath( params.bam ).map { file -> tuple(file.baseName, file) }
sequencesintomutationalcalling = Channel.fromPath( params.bam ).map { file -> tuple(file.baseName, file) } // You cannot use the same input twice in a channel, therefore we are copying it in twice

ref_index = file(params.ref)

//params.workingDir = '.'
//params.outdir = '/jumbo/WorkingDir/B19-057/Data/Meta/Alignment/Epidermidis/testNextflow'

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

	input:
        set file_ID, file(rgbam), file(rgbai) from rg_out

	output: 
	set file_ID, "${file_ID}.rmdup.bam","${file_ID}.rmdup.bai" into rmdup_out1, rmdup_out2

	script: 
	"""
	java -Xmx4g -jar /apps/bio/apps/picard/2.1.0/picard.jar MarkDuplicates I=${rgbam} O=${file_ID}.rmdup.bam CREATE_INDEX=True M=${file_ID}.marked_dup_metrics.txt
	"""

}

//----------------------RealignTargetCreator-------------------------------

process run_realignTargetCreator {
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	
	clusterOptions='-pe mpi 4'
	executor 'sge'
	queue 'bfxcore.q@node6-bfx.medair.lcl'

	input:
	set file_ID, file(rdupbam), file(rdupbai) from rmdup_out1
	
	output: 
	set file_ID, "${file_ID}.intervals" into intervals

	script:
	"""
	java -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${rdupbam} -o ${file_ID}.intervals -R ${ref_index}
	"""
}


//----------------------IndelRealignment-------------------------------

process run_indelrealigner{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 4'
	executor 'sge'
	queue 'bfxcore.q@node6-bfx.medair.lcl'
	

	input:
        set file_ID, file(interval) from intervals
	set file_ID, file(rdupbam), file(rdupbai) from rmdup_out2

	output:
	set file_ID, "${file_ID}.realigned.bam" into realigned

	script:
	"""
	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref_index} --targetIntervals ${interval} -I ${rdupbam} -o ${file_ID}.realigned.bam
	"""
}

//----------------------SNPCalling-------------------------------

process run_SNPcalling{
	publishDir params.outdir, mode: 'copy', overwrite: true
	//errorStrategy 'ignore'
	clusterOptions='-pe mpi 1'
        executor 'sge'
        queue 'bfxcore.q@node6-bfx.medair.lcl'

	input:
	set file_ID, file(realbam) from realigned

	output:
	set file_ID, "${file_ID}.vcf" into rawsnps

	script:
	"""
	/apps/bio/apps/samtools/1.6/samtools index ${realbam}
	
	java -Xmx4g -jar /apps/bio/apps/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref_index} -I ${realbam} -o ${file_ID}.vcf
	"""
}



 


