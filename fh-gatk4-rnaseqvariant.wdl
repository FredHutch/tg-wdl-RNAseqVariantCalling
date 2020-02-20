####### RNA Seq from paired end sequencing - variant calling
# Modeled after https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels
# requires a bed file containing the regions you want to call variants on rather than doing whole transcriptome
workflow rnaSeqVariantCalling {

  # sample specific parameters
  String base_file_name
  Array[File] R1Fastq
  Array[File] R2Fastq

  #Software modules for gizmo
  String GATKModule
  String samtoolsModule
  String starModule
  String perlModule
  
  # chemistry specific parameters
  File bedFile

  # genome specific parameters
  String ref_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  # Annovar Inputs
  # Note:  For Annovar, please reference: Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010
  String annovarDIR 
  String annovar_protocols
  String annovar_operation

  # Inputs for STAR
  String STARgenomeDir  # with default 100bp read length, for local
  String referenceGenome 

  # Prepare bed file and check sorting
  call SortBed {
    input:
      unsorted_bed = bedFile,
      ref_dict = ref_dict,
      modules = GATKModule
  }

  call concatenateFastqs {
    input:
    base_file_name = base_file_name,
    fastqR1String = R1Fastq,
    fastqR2String = R2Fastq
  }

  # Align reads using STAR for the entire sample
	call StarAlign { 
		input: 
      genomeDir = STARgenomeDir,
			fastq1 = concatenateFastqs.R1fastq,
      fastq2 = concatenateFastqs.R2fastq,
			base_file_name = base_file_name,
      referenceGenome = referenceGenome,
      modules = starModule + " " + samtoolsModule,
      threads = 4
	}

  # Split, split reads
  call SplitNCigarReads {
    input:
      input_bam = StarAlign.alignedBam,
      input_bam_index = StarAlign.alignedBamBai,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      interval_list = SortBed.intervals,
      modules = GATKModule + " " + samtoolsModule
  }
  # Generate the recalibration model by interval and apply it
  call ApplyBaseRecalibrator {
    input:
      input_bam = SplitNCigarReads.splitBam,
      input_bam_index = SplitNCigarReads.splitBamIndex,
      base_file_name = base_file_name,
      intervals = SortBed.intervals,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      modules = GATKModule + " " + samtoolsModule
  }

    # Generate haplotype caller vcf
    call HaplotypeCaller {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        intervals = SortBed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        modules = GATKModule
    }

	call VariantFiltration {
		input:
			input_vcf = HaplotypeCaller.output_vcf,
			input_vcf_index = HaplotypeCaller.output_index,
			base_file_name = base_file_name + ".variant_filtered.vcf.gz",
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_dict = ref_dict,
			modules = GATKModule
	}

    # Annotate variants
    call annovar {
      input:
        input_vcf = VariantFiltration.output_vcf,
        ref_name = ref_name,
        annovarDIR = annovarDIR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        modules = perlModule
    }

# Outputs that will be retained when execution is complete
output {
    File analysisReadyBam = SplitNCigarReads.splitBam
    File analysisReadyBamIndex = SplitNCigarReads.splitBamIndex
    File GATK_vcf = VariantFiltration.output_vcf
    File GATK_annotated_vcf = annovar.output_annotated_vcf
    File GATK_annotated = annovar.output_annotated_table

  }
}# End workflow

#### TASK DEFINITIONS


# Prepare bed file and check sorting
task SortBed {
  File unsorted_bed
  File ref_dict
  String modules
  command {
    set -eo pipefail

    echo "Sort bed file"
    sort -k1,1V -k2,2n -k3,3n ${unsorted_bed} > sorted.bed

    echo "Transform bed file to intervals list with Picard----------------------------------------"
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      BedToIntervalList \
      -I=sorted.bed \
      -O=sorted.interval_list \
      -SD=${ref_dict}
  }
  runtime {
    modules: "${modules}"
    docker: "broadinstitute/gatk:4.1.0.0"
  }
  output {
    File intervals = "sorted.interval_list"
    File sorted_bed = "sorted.bed"
  }
}


# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  File input_bam
  File intervals 
  File input_bam_index
  String base_file_name
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String modules
  command {
    set -e

  gatk --java-options "-Xms4g" \
    BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${base_file_name}.recal_data.csv \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      --intervals ${intervals} \
      --interval-padding 0 

  gatk --java-options "-Xms4g" \
    ApplyBQSR \
      -bqsr ${base_file_name}.recal_data.csv \
      -I ${input_bam} \
      -O ${base_file_name}.recal.bam \
      -R ${ref_fasta} \
      --intervals ${intervals} \
      --interval-padding 0 
  

    #finds the current sort order of this bam file
  samtools view -H ${base_file_name}.recal.bam|grep @SQ|sed 's/@SQ\tSN:\|LN://g' > ${base_file_name}.sortOrder.txt

  }
  runtime {
    memory: 8000
    cpu: 4
    partition: "campus"
    walltime: "01:00:00"
    modules: "${modules}"
  }
  output {
    File recalibrated_bam = "${base_file_name}.recal.bam"
    File recalibrated_bai = "${base_file_name}.recal.bai"
    File sortOrder = "${base_file_name}.sortOrder.txt"
  }
}


# HaplotypeCaller per-sample
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String base_file_name
  File intervals
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf
  File dbSNP_vcf_index
  String modules 
  command {
    set -e

    gatk --java-options "-Xms2g" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${base_file_name}.vcf.gz \
      -dont-use-soft-clipped-bases \
      --intervals ${intervals} \
      --interval-padding 0 \
      --dbsnp ${dbSNP_vcf}

  }
  runtime {
    memory: 8000
    walltime: "00:45:00"
    modules: "${modules}"
  }
  output {
    File output_vcf = "${base_file_name}.vcf.gz"
    File output_index = "${base_file_name}.vcf.gz.csi"
  }
}




# Align to the genome using STAR
task StarAlign {
  String genomeDir
	File fastq1
  File fastq2
	String base_file_name
  String referenceGenome
  String modules
  Int threads

	command {
		set -e

		STAR \
		--genomeDir ${genomeDir} \
		--runThreadN ${threads} \
		--readFilesIn ${fastq1} ${fastq2} \
		--readFilesCommand "gunzip -c" \
		--sjdbOverhang 100 \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
    --outFileNamePrefix ${base_file_name}.${referenceGenome}. \
    --limitBAMsortRAM 5000000000

    samtools index ${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam

	}
	runtime {
		memory: 8000
		cpu: "${threads}"
		walltime: "02:00:00"
    modules: "${modules}"
	}
  output {
		File alignedBam = "${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam"
    File alignedBamBai = "${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam.bai"
		File output_log = "${base_file_name}.${referenceGenome}.Log.out"
		File output_log_progress = "${base_file_name}.${referenceGenome}.Log.progress.out"
	}
}

# split reads when they align across a junction
task SplitNCigarReads {
  File input_bam
  File input_bam_index
  String base_file_name
  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String modules
  command {
    set -eo pipefail

    gatk --java-options "-Xms4g" \
       SplitNCigarReads \
                -R ${ref_fasta} \
                -I ${input_bam} \
                -O ${base_file_name}.initial.bam 

    gatk --java-options "-Xms4g" \
       AddOrReplaceReadGroups \
       -I=${base_file_name}.initial.bam \
       -O=${base_file_name}.bam  \
       --RGPU=10x \
       --RGLB=10x5prime \
       --RGPL=illumina \
       --RGSM=${base_file_name}

    samtools index ${base_file_name}.bam 
    }
  runtime {
    memory: 16000
    partition: "campus"
    walltime: "02:00:00"
    modules: "${modules}"
  }
  output {
    File splitBam = "${base_file_name}.bam"
    File splitBamIndex = "${base_file_name}.bam.bai"
 }
}


task VariantFiltration {

	File input_vcf
	File input_vcf_index
	String base_file_name

 	File ref_dict
 	File ref_fasta
 	File ref_fasta_index

	String modules

	command {
  set -eo pipefail
		 gatk \
		    VariantFiltration \
			--R ${ref_fasta} \
			--V ${input_vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O ${base_file_name}
  }
  runtime {
    modules: "${modules}"
  }
	output {
    	File output_vcf = "${base_file_name}"
    	File output_vcf_index = "${base_file_name}.tbi"
	}
}




# annotate with annovar
task annovar {
  File input_vcf
  String ref_name
  String annovar_protocols
  String annovar_operation
  String annovarDIR
  String base_vcf_name = basename(input_vcf, ".vcf")
  String modules
  command {
  set -eo pipefail
  
  perl ${annovarDIR}/annovar/table_annovar.pl ${input_vcf} ${annovarDIR}/annovar/humandb/ \
    -buildver ${ref_name} \
    -outfile ${base_vcf_name} \
    -remove \
    -protocol ${annovar_protocols} \
    -operation ${annovar_operation} \
    -nastring . -vcfinput
  }
  runtime {
    docker: "perl:5.28.0"
    modules: "${modules}"
  }
  output {
    File output_annotated_vcf = "${base_vcf_name}.${ref_name}_multianno.vcf"
    File output_annotated_table = "${base_vcf_name}.${ref_name}_multianno.txt"
  }
}

task concatenateFastqs {
  Array[File] fastqR1String
  Array[File] fastqR2String
  String base_file_name

  command {
    set -eo pipefail

    cat "${sep=" " fastqR1String}[@]" > ${base_file_name}_R1.fastq.gz
    cat "${sep=" " fastqR2String}[@]" > ${base_file_name}_R2.fastq.gz
    }
  runtime {
    cpu: 2
    memory: 4000
    partition: "campus"
  }
  output {
    File R1fastq = "${base_file_name}_R1.fastq.gz"
    File R2fastq = "${base_file_name}_R2.fastq.gz"
  }
}