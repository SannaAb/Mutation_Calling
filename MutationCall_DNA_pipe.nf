#!/usr/bin/env nextflow

params.fastq = ''
params.ref = ''

//sequences = Channel.fromPath(params.fastq)

sequences = Channel
                .fromPath(params.fastq)
                .map { file -> tuple(file.baseName, file) }

ref_index = file(params.ref)
params.workingDir = '/home/xabras/NextflowPractice'



//----------------------Bowtie2-------------------------------

process run_Bowtie2 {
        publishDir params.workingDir, mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        clusterOptions='-pe mpi 4'
        executor 'sge'
        queue 'bfxcore.q@node4-bfx.medair.lcl'

        module 'bowtie2/2.2.9'

        input:
		set file_ID, file(fq) from sequences

        output:
        set file_ID, file("${file_ID}_bowtie2_alignment.sam") into bowtie_out

        script:

        """
		bowtie2 --sensitive -x ${ref_index} -U ${fq} -S ${file_ID}_bowtie2_alignment.sam
        """
}


//----------------------Samtools-------------------------------

process run_Samtools {
        publishDir params.workingDir, mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	executor 'sge'
	queue 'bfxcore.q@node4-bfx.medair.lcl'
	cpus '1'

        module 'samtools/1.9'

        input:
		set file_ID, sam from bowtie_out

        output:
        set file_ID, "${file_ID}_*" into samtools_out

        script:

        """
		samtools view -hb ${sam} > ${file_ID}_bowtie2_alignment.bam
        samtools sort ${file_ID}_bowtie2_alignment.bam -o ${file_ID}_bowtie2_sorted.bam
        samtools index ${file_ID}_bowtie2_sorted.bam
        """
}


//----------------------Samtools-------------------------------

process sort_files {
        publishDir params.workingDir, mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		set file_ID, bam from samtools_out

        script:

        """
        cd ${params.workingDir}
        mkdir -p ${file_ID}
        mv ${file_ID}*.* ${file_ID}/
        """
}






