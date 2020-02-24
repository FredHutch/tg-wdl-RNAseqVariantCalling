# tg-wdl-RNAseqVariantCalling
WDL workflow for variant calling from bulk RNA seq for a region of interest using GATK tools.  

Note:  Instructions for running [this workflow](https://github.com/FredHutch/tg-wdl-RNAseqVariantCalling/blob/master/fh-gatk4-rnaseqvariant.wdl) from your local machine are in [this R Markdown](https://github.com/FredHutch/tg-wdl-RNAseqVariantCalling/blob/master/manifestPrep.Rmd)

## Workflow inputs:
This workflow is intended for paired end RNA sequencing data, of 100pb read length or less, from human RNA.  

### Parameter inputs
The [parameter inputs json](https://github.com/FredHutch/tg-wdl-RNAseqVariantCalling/blob/master/rna-variant-gizmo-parameters.json) contains: 
- a set of reference data files for hg38, stored in our local file system, 
- a pre-compiled STAR referebce data set (for STAR 2.7.1, with hg38 and gencode v31) and,
- a pre-cpmared Annovar database using Annovar released 4/16/2018 and rolled up on 1/5/2009 with the specified annotation sources (see: "rnaSeqVariantCalling.annovar_protocols"), and finally
- the Easybuild modules specifying the exact version of the software used to run each of the tools in this workflow.  

This input file will only change if you want to update software or genome reference versions.

### Sample inputs
The [sample inputs json](https://github.com/FredHutch/tg-wdl-RNAseqVariantCalling/blob/master/rna-variant-gizmo-sample.json), will change run-to-run for each of the samples you want to run. This workflow can also be adjusted to run from a batch file for a group of samples at a time (just ask @vortexing).  This inputs file contains:

- the base file name, or sample name you use to identify this sample (and will be the beginning of all output file names),
- a list of the paths in the file system to all the R1/forward reads for this sample in a format like this (order doesn't matter):
```
["/path/R1.001.fastq.gz", "/path/R1.002.fastq.gz",  "/path/R1.003.fastq.gz", "/path/R1.004.fastq.gz"]
```
- a corresponding list of the R2/reverse reads for this sample (order doesn't matter)
- a bed formatted file specifying the regions in the genome over which you'd like to identify variants.  

This one sample, targeted variant calling approach was chosen as for the people wanting to use this workflow, that is likely to be the approach they will take.  If you'd like to edit the workflow to take as an input a batch file so you can scatter over mulitple samples, OR you'd like to send one sample at a time but scatter over each chromosome to do whole genome variant calling, these edits can be made.  Currently (2/2020), doing both at the same time is not possible as scatters within scatters are not yet fully supported in WDL. 


## Outputs
Currently (2/2020) the outputs of this workflow are defined in the output chunk as:
```
output {
    File analysisReadyBam = SplitNCigarReads.splitBam
    File analysisReadyBamIndex = SplitNCigarReads.splitBamIndex
    File GATK_vcf = VariantFiltration.output_vcf
    File GATK_annotated_vcf = annovar.output_annotated_vcf
    File GATK_annotated = annovar.output_annotated_table
}
```
These are the analysis ready bam file that was used for the variant calling, it's index, the vcf created and then filtered, the annotated vcf (via Annovar) and the tabular version of the annotated vcf that is simpler to view in Excel than a vcf but contains identical information.

To add outputs simply add task output references to this workflow chunk. 
