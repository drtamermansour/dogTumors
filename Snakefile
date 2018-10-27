import glob
import os

#fastq_pattern = "data/fastq/*/*_001.fastq.gz"
fastq_pattern = "data/fastq/*/[GMO]N*_001.fastq.gz"
fastq_files = glob.glob(fastq_pattern)
read_ids = [fq[11:-9] for fq in fastq_files]
fragment_ids = list(set([fq[11:-16] for fq in fastq_files]))

normFastq_pattern= "data/fastq/*/[GMO]N*_001.fastq.gz"
normFastq_files = glob.glob(normFastq_pattern)
normFragment_ids = list(set([fq[11:-16] for fq in normFastq_files]))

tumFastq_pattern= "data/fastq/*/[GMO]T*_001.fastq.gz"
#tumFastq_pattern= "data/fastq/*/GT11*_001.fastq.gz"
tumFastq_files = glob.glob(tumFastq_pattern)
tumFragment_ids = list(set([fq[11:-16] for fq in tumFastq_files]))

def get_normal(wildcards):
    # code that returns the normal control of a tumor samle using its tumFragment_ids
    dirName=os.path.dirname(wildcards.tumFragment)
    tumSM=os.path.basename(wildcards.tumFragment).split('_')[0]
    # accomodate samples OT33M OT33P
    if not tumSM[-1].isdigit():
        tumSM=tumSM[:-1]
    normSM=tumSM.replace('T', 'N')
    #print('data/recalib/' + dirName + '/' + tumSM + '_*.bam')
    normAr=glob.glob('data/recalib/' + dirName + '/' + normSM + '_*.bam')
    # accomodate tumor samples without matching control (GT1, GT14, GT15)
    if len(normAr) == 0:
        return 'data/recalib/' + wildcards.tumFragment + '.bam'
    else:
        return normAr[0]

def get_normalTable(wildcards):
    # code that returns the normal control of a tumor samle using its tumFragment_ids
    dirName=os.path.dirname(wildcards.tumFragment)
    tumSM=os.path.basename(wildcards.tumFragment).split('_')[0]
    # accomodate samples OT33M OT33P
    if not tumSM[-1].isdigit():
        tumSM=tumSM[:-1]
    normSM=tumSM.replace('T', 'N')
    #print('data/pileupsum/' + dirName + '/' + tumSM + '_*.table')
    normAr=glob.glob('data/pileupsum/' + dirName + '/' + normSM + '_*.table')
    # accomodate tumor samples without matching control (GT1, GT14, GT15)
    if len(normAr) == 0:
        return 'data/pileupsum/' + wildcards.tumFragment + '.table'
    else:
        return normAr[0]


#class Object(object):
#    pass

#sample = Object()

#for id in tumFragment_ids:
#    sample.tumFragment=id
#    get_normal(sample)


rule all:
    input:
        expand("qc/fastqc/fastq/{read}_fastqc.{ext}", read=read_ids, ext=['html', 'zip']),
        "qc/multiqc/fastq/multiqc_report.html",
        expand("data/uBAM/{fragment}.bam", fragment=fragment_ids),
        expand("data/adap_uBAM/{fragment}.adap.bam", fragment=fragment_ids),
        expand("data/adap_uBAM/{fragment}.adap_metrics.txt", fragment=fragment_ids),
        "refGenome/canFam3_chr.fa",
        "refGenome/BwaIndex/genome.fa", expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        "refGenome/gatkIndex/genome.fa", "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        expand("data/mapped_reads/{fragment}.bam", fragment=fragment_ids),
        expand("data/dedup/{fragment}.bam", fragment=fragment_ids), expand("data/dedup/{fragment}.metrics.txt", fragment=fragment_ids),
        "knowVar/canis_familiaris_SNPs.vcf", "knowVar/canis_familiaris_indels.vcf",
        expand("data/recalib/{fragment}.txt", fragment=fragment_ids),
        expand("data/recalib/{fragment}.bam", fragment=fragment_ids),
        expand("data/PON/{normFragment}.vcf.gz", normFragment=normFragment_ids),
        expand("data/PON/{normFragment}.vcf.gz.tbi", normFragment=normFragment_ids),
        "data/PON.vcf.gz",
        expand("data/rawVC/{tumFragment}.vcf.gz", tumFragment=tumFragment_ids),
        expand("data/realign/{tumFragment}.bam", tumFragment=tumFragment_ids),
        "data/population/commonAlt.vcf",
        expand("data/pileupsum/{fragment}.table", fragment=fragment_ids),
        expand("data/cont/{tumFragment}.table", tumFragment=tumFragment_ids),
        expand("data/seg/{tumFragment}.table", tumFragment=tumFragment_ids),
        expand("data/filteredVC_one/{tumFragment}.vcf.gz", tumFragment=tumFragment_ids),
        expand("data/filteredVC_one/{tumFragment}.stat", tumFragment=tumFragment_ids),
        expand("data/preAdap/{tumFragment}.pre_adapter_detail_metrics", tumFragment=tumFragment_ids),
        expand("data/filteredVC_2nd/{tumFragment}.vcf.gz", tumFragment=tumFragment_ids),
#        expand("data/filteredVC_2nd/{tumFragment}.vcf.gz", tumFragment="Project_TBDY_L1_DTDB05_H1220P_York_1/GT7_S11_L001"),
#        expand("data/cont/{tumFragment}.table", tumFragment="Project_TBDY_L1_DTDB05_H1220P_York_1/GT7_S11_L001"),
#        expand("data/seg/{tumFragment}.table", tumFragment="Project_TBDY_L1_DTDB05_H1220P_York_1/GT7_S11_L001"),
#        expand("data/preAdap/{tumFragment}", tumFragment="Project_TBDY_L1_DTDB05_H1220P_York_1/GT7_S11_L001"),

rule fastqc_pre:
    input:
        "data/fastq/{read}.fastq.gz"
    output:
        html="qc/fastqc/fastq/{read}_fastqc.html",
        zip="qc/fastqc/fastq/{read}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        #module load FastQC/0.11.5 ## version FastQC v0.11.7 was installed in the snakemake env
        source activate snakemake
        fastqc -t {threads} -f fastq -noextract -o qc/fastqc/fastq/$(dirname {wildcards.read}) {input}
        '''

rule multiQC_pre:
    input:
        html=expand("qc/fastqc/fastq/{read}_fastqc.html", read=read_ids)
    output:
        "qc/multiqc/fastq/multiqc_report.html",
        "qc/multiqc/fastq/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 16000000000
    shell:
        "/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o qc/multiqc/fastq qc/fastqc/fastq"

### A) pre-processing
### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
### https://software.broadinstitute.org/gatk/documentation/article?id=6483#step1
### All docs of the picard tools in GATK are still using the Java command line not the gatk. I updated the usage of these tools (see FastqToSam)

## the doc is used from https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_FastqToSam.php
## but the syntex is adjusted as described in https://software.broadinstitute.org/gatk/documentation/quickstart
rule FastqToSam:
    input:
        r1="data/fastq/{fragment}_R1_001.fastq.gz",
        r2="data/fastq/{fragment}_R2_001.fastq.gz"
    output:
        "data/uBAM/{fragment}.bam"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        name=$(basename {input.r1})
        SM=$(echo $name | cut -d "_" -f1)
        LB=$(echo $name | cut -d"_" -f1,2)  ## We use <Index.Sequence> in the Illumina file name as an index to the library
        #batch=$(basename "$(dirname {input.r1})")
        #if [ "$batch" != "trimmed" ];then LB=$batch.$LB;fi
        PL="Illumina"
        ##read Fastq 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>"
        header=$(head -n1 <(zcat {input.r1}) | grep ':*:*:*:*:*:*')
        if [ "$header" != "" ]; then
            RGID=$(echo "$header" | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        else # "make unique ID and PU using checksum"
            checksum=$(shasum {input.r1} | awk '{{ print $1 }}')
            RGID="UnChrPU_"$checksum
        fi
        PU=$RGID.$LB

        module load Java/jdk1.8.0
        source activate gatk
        #java -jar picard.jar FastqToSam
        gatk --java-options "-Xmx6G" FastqToSam \
        -F1={input.r1} \
        -F2={input.r2} \
        -O={output} \
        -SM=$SM \
        -LB=$LB \
        -PL=$PL \
        -RG=$RGID \
        -PU=$PU \
        --TMP_DIR="tmp/{wildcards.fragment}"
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_illumina_MarkIlluminaAdapters.php
rule MarkIlluminaAdapters:
    input:
        "data/uBAM/{fragment}.bam"
    output:
        bam="data/adap_uBAM/{fragment}.adap.bam",
        metrics="data/adap_uBAM/{fragment}.adap_metrics.txt"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        gatk --java-options "-Xmx6G" MarkIlluminaAdapters \
        -I={input} \
        -O={output.bam} \
        -M={output.metrics} \
        --TMP_DIR="tmp2/{wildcards.fragment}"
        '''

rule download_ref:
    output:
        "refGenome/canFam3_chr.fa"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        '''
        mkdir refGenome && cd refGenome
        wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz' -O canFam3.fa.gz
        gunzip canFam3.fa.gz
        cat canFam3.fa | awk '{if($1 ~ ">chrUn_"){f=0;}else if($1 ~ ">chr"){print $0;f=1;}else if(f){print $0;}}' > canFam3_chr.fa
        '''

rule bwa_index:
    input:
        "refGenome/canFam3_chr.fa"
    output:
        "refGenome/BwaIndex/genome.fa",
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    shell:
        '''
        if [ ! -f refGenome/BwaIndex/genome.fa ];then ln -s ../canFam3_chr.fa refGenome/BwaIndex/genome.fa;fi
        module load bwa/0.7.7.r441
        bwa index -p refGenome/BwaIndex/genome -a bwtsw {input}
        '''

rule GATK_index:
    input:
        "refGenome/canFam3_chr.fa"
    output:
        ref="refGenome/gatkIndex/genome.fa",
        index="refGenome/gatkIndex/genome.fa.fai",
        dict="refGenome/gatkIndex/genome.dict",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    shell:
        '''
        if [ ! -f refGenome/gatkIndex/genome.fa ];then ln -s ../canFam3_chr.fa refGenome/gatkIndex/genome.fa;fi
        module load SAMTools/1.5
        module load picardTools/1.89
        samtools faidx "refGenome/gatkIndex/genome.fa"
        java -Xmx4g -jar $PICARD/CreateSequenceDictionary.jar R= {input} O= {output.dict}
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.2.0/picard_sam_SamToFastq.php
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_MergeBamAlignment.php
# With the original MarkDuplicates aproach, MergeBamAlignment runs with SORT_ORDER 'coordinate'. If you run with SORT_ORDER other than 'coordinate' (e.g. in the WDL pipeline implementation), you will need to use SortSam to sort by queryname and SetNmAndUqTags (DEPRECATED: Use SetNmMdAndUqTags instead) to fix these tags
rule align:
    input:
        bam="data/adap_uBAM/{fragment}.adap.bam",
        bwa_ref="refGenome/BwaIndex/genome.fa",
        gatk_ref="refGenome/gatkIndex/genome.fa",
    output:
        bam="data/mapped_reads/{fragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 32
    threads:4
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        module load bwa/0.7.7.r441
        gatk --java-options "-Xmx30G" SamToFastq \
        -I={input.bam} \
        --FASTQ=/dev/stdout \
        --CLIPPING_ATTRIBUTE=XT --CLIPPING_ACTION=2 --INTERLEAVE=true -NON_PF=true \
        --TMP_DIR="tmp3/{wildcards.fragment}" | \
        bwa mem -M -t {threads} -p refGenome/BwaIndex/genome /dev/stdin | \
        gatk --java-options "-Xmx30G" MergeBamAlignment \
        --ALIGNED_BAM=/dev/stdin \
        --UNMAPPED_BAM={input.bam} \
        --OUTPUT={output.bam} \
        -R={input.gatk_ref} --CREATE_INDEX=true --ADD_MATE_CIGAR=true \
        --CLIP_ADAPTERS=false --CLIP_OVERLAPPING_READS=true \
        --INCLUDE_SECONDARY_ALIGNMENTS=true --MAX_INSERTIONS_OR_DELETIONS=-1 \
        --PRIMARY_ALIGNMENT_STRATEGY=MostDistant --ATTRIBUTES_TO_RETAIN=XS \
        --TMP_DIR="tmp4/{wildcards.fragment}"
        '''

# The apprach for MarkDuplicates coded here uses coordinate-sorted and indexed bam files. Recently, the tools added a feature to accept queryname-sorted inputs that in turn by default activates additional features that will give DIFFERENT duplicate flagging results than outlined in this tutorial. Namely, if you provide MarkDuplicates a queryname-sorted BAM, then if a primary alignment is marked as duplicate, then the tool will also flag its (i) unmapped mate, (ii) secondary and/or (iii) supplementary alignment record(s) as duplicate. You can get queryname sorted BAM using SortSam tool.
# You can simultaneously MarkDups and merge RG BAMs by providing the files altogether to MarkDuplicates by specifying each file with I= (picard) or -I (gatk). This works with coordinate-sorted alignments and should also work with queryname-sorted alignment files.

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php
# https://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
# https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/
# https://www.biostars.org/p/229842/
rule MarkDuplicates:
    input:
        bam="data/mapped_reads/{fragment}.bam",
    output:
        bam="data/dedup/{fragment}.bam",
        metrics="data/dedup/{fragment}.metrics.txt",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 24
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        gatk --java-options "-Xmx20G" MarkDuplicates \
        --INPUT={input.bam} \
        --OUTPUT={output.bam} \
        --METRICS_FILE={output.metrics} \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        --CREATE_INDEX=true \
        --TMP_DIR="tmp5/{wildcards.fragment}"
        '''

rule download_knowVar:
    output:
        knownSNPs="knowVar/canis_familiaris_SNPs.vcf",
        knownIndels="knowVar/canis_familiaris_indels.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        '''
        mkdir knowVar && cd knowVar
        wget 'ftp://ftp.ensembl.org/pub/release-94/variation/vcf/canis_familiaris/canis_familiaris.vcf.gz' -O canis_familiaris.vcf.gz
        gunzip canis_familiaris.vcf.gz
        ## change the name of the chromosomes to match the UCSC genome (bet keep the file co-ordinates 1-based)
        grep -v "^#" canis_familiaris.vcf | awk '{{ print "chr"$0 }}' > canis_familiaris_fixedChrNames.vcf
        cd ../
        perl scripts/sortByRef.pl knowVar/canis_familiaris_fixedChrNames.vcf refGenome/gatkIndex/genome.fa.fai > knowVar/canis_familiaris_fixedChrNames_sorted.vcf
        grep "^#" knowVar/canis_familiaris.vcf > knowVar/canis_familiaris_SNPs.vcf
        grep "TSA=SNV" knowVar/canis_familiaris_fixedChrNames_sorted.vcf >> knowVar/canis_familiaris_SNPs.vcf
        grep "^#" knowVar/canis_familiaris.vcf > knowVar/canis_familiaris_indels.vcf
        grep -v "TSA=SNV" knowVar/canis_familiaris_fixedChrNames_sorted.vcf >> knowVar/canis_familiaris_indels.vcf
        module load Java/jdk1.8.0
        source activate gatk
        gatk IndexFeatureFile -F knowVar/canis_familiaris_SNPs.vcf
        gatk IndexFeatureFile -F knowVar/canis_familiaris_indels.vcf
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
# https://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals
# https://software.broadinstitute.org/gatk/documentation/article.php?id=1319
# https://software.broadinstitute.org/gatk/documentation/article.php?id=44
# https://gatkforums.broadinstitute.org/gatk/discussion/11081/base-quality-score-recalibration-bqsr
rule BaseRecalib:
    input:
        bam="data/dedup/{fragment}.bam",
        ref="refGenome/gatkIndex/genome.fa",
        knownSNPs="knowVar/canis_familiaris_SNPs.vcf",
        knownIndels="knowVar/canis_familiaris_indels.vcf",
    output:
        report="data/recalib/{fragment}.txt",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        gatk --java-options "-Xmx15G" BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --use-original-qualities \
        -O {output.report} \
        --known-sites {input.knownSNPs} \
        --known-sites {input.knownIndels} \
        -L refGenome/intervals.bed \
        -ip 100
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php
rule ApplyBQSR:
    input:
        bam="data/dedup/{fragment}.bam",
        ref="refGenome/gatkIndex/genome.fa",
        report="data/recalib/{fragment}.txt",
    output:
        bam="data/recalib/{fragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        gatk --java-options "-Xmx15G" ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.bam} \
        -bqsr {input.report} \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        --add-output-sam-program-record \
        --use-original-qualities \
        -L refGenome/intervals.bed \
        -ip 100
        '''

# https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
rule Mutect2_for_PON:
    input:
        bam="data/recalib/{normFragment}.bam",
        ref="refGenome/gatkIndex/genome.fa",
    output:
        var="data/PON/{normFragment}.vcf.gz",
        tbi="data/PON/{normFragment}.vcf.gz.tbi",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 24,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    threads:1
    shell:
        '''
        name=$(basename {input.bam})
        SM=$(echo $name | cut -d "_" -f1)
        module load Java/jdk1.8.0
        source activate gatk
        gatk --java-options "-Xmx15G" Mutect2 \
        -R {input.ref} \
        -I {input.bam} \
        -tumor $SM \
        -L refGenome/intervals.bed \
        -O {output.var}
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.9.0/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
rule CreateSomaticPON:
    input:
        expand("data/PON/{normFragment}.vcf.gz", normFragment=normFragment_ids),
    output:
        "data/PON.vcf.gz",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        pon_var=(data/PON/*/*.vcf.gz)
        gatk --java-options "-Xmx15G"  CreateSomaticPanelOfNormals \
        $(printf " -vcfs %s" "${{pon_var[@]}}") \
        -O {output}
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php  ## manily for the section about "How GATK4 Mutect2 differs from GATK3 MuTect2"
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.9.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
# remove "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter". This is needed for alt-aware mapping
# Add "--genotype-pon-sites": The panel of normals represents common germline variant sites and also commonly noisy sites in sequencing data. By default, the tool does not reassemble nor emit variant sites that match identically to a PoN variant. To enable genotyping of PoN sites, use the --genotype-pon-sites option. If the match is not exact, e.g. there is an allele-mismatch, the tool reassembles the region, emits the calls and annotates matches in the INFO field with IN_PON.
#

rule somaticVarCalling:
    input:
        tumBam="data/recalib/{tumFragment}.bam",
        normBam=get_normal,
        ref="refGenome/gatkIndex/genome.fa",
        pon="data/PON.vcf.gz",
    output:
        rawVar="data/rawVC/{tumFragment}.vcf.gz",
        realignedBam="data/realign/{tumFragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 24,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        # We need to create these files regardless, even if they stay empty
        #touch {output.realignedBam}
        mkdir -p $(dirname tmp6/{wildcards.tumFragment}.normal_name.txt)
        echo "" > tmp6/{wildcards.tumFragment}.normal_name.txt
     
        # retrieve samples names  
        gatk --java-options "-Xmx7G" GetSampleName -R {input.ref} -I {input.tumBam} \
        -O tmp6/{wildcards.tumFragment}.tumor_name.txt -encode 
        tumor_command_line="-I {input.tumBam} -tumor `cat tmp6/{wildcards.tumFragment}.tumor_name.txt`"
        if [[ "{input.normBam}" != "{input.tumBam}" ]]; then
            gatk --java-options "-Xmx7G" GetSampleName -R {input.ref} -I {input.normBam} \
            -O tmp6/{wildcards.tumFragment}.normal_name.txt -encode 
            normal_command_line="-I {input.normBam} -normal `cat tmp6/{wildcards.tumFragment}.normal_name.txt`"
        else normal_command_line=""
        fi

        gatk --java-options "-Xmx7g" Mutect2 \
        -R {input.ref} \
        $tumor_command_line \
        $normal_command_line \
        -pon {input.pon} \
        --genotype-pon-sites \
        --af-of-alleles-not-in-resource 0.001 \
        -L refGenome/intervals.bed \
        -O {output.rawVar} \
        -bamout {output.realignedBam}
        '''

#            ${"--orientation-bias-artifact-priors " + artifact_prior_table} \

## Create population germline resource containing only common biallelic variants
rule popCommon:
    input:
        ref="refGenome/gatkIndex/genome.fa",
        dogSeq="/mnt/home/mansourt/Tamer/dogSeq/varResults/GenotypeGVCFs_output_max50.pass.combinedFiltered.vcf",
    output:
        "data/population/commonAlt.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        gatk SelectVariants \
        -R {input.ref} \
        -V {input.dogSeq} \
        --restrict-alleles-to BIALLELIC \
        -select 'AN > 50 && AC > 1' \
        -L refGenome/intervals.bed \
        --sites-only-vcf-output=true \
        -O {output}
        '''

rule GetPileupSummaries:
    input:
        bam="data/recalib/{fragment}.bam",
        commonAlt="data/population/commonAlt.vcf",
    output:
        pileupsum="data/pileupsum/{fragment}.table",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 12
    threads:1
    shell:
        '''
        gatk GetPileupSummaries \
        -I {input.bam} \
        -V {input.commonAlt} \
        -L refGenome/intervals.bed \
        -O {output.pileupsum}
        '''

rule CalculateContamination:
    input:
        tumTable="data/pileupsum/{tumFragment}.table",
        normTable=get_normalTable,
    output:
        contTable="data/cont/{tumFragment}.table",
        segTable="data/seg/{tumFragment}.table"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        if [[ "{input.normTable}" != "{input.tumTable}" ]]; then
            NORMAL_CMD="-matched {input.normTable}"
        else NORMAL_CMD=""
        fi

        gatk CalculateContamination \
        -I {input.tumTable} \
        $NORMAL_CMD \
        -O {output.contTable} \
        --tumor-segmentation {output.segTable}
        '''

rule FilterMutectCalls:
    input:
        rawVar="data/rawVC/{tumFragment}.vcf.gz",
        contTable="data/cont/{tumFragment}.table",
        segTable="data/seg/{tumFragment}.table"
    output:
        filterVar="data/filteredVC_one/{tumFragment}.vcf.gz",
        filterStat="data/filteredVC_one/{tumFragment}.stat"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    threads:1
    shell:
        '''
        gatk --java-options "-Xmx12G" FilterMutectCalls \
        -V {input.rawVar} \
        --contamination-table {input.contTable} \
        --tumor-segmentation {input.segTable} \
        --stats {output.filterStat} \
        -O {output.filterVar}
        '''


rule CollectSequencingArtifactMetrics:
    input:
        tumBam="data/recalib/{tumFragment}.bam",
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "data/preAdap/{tumFragment}.pre_adapter_detail_metrics",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 32
    threads:1
    shell:
        '''
        gatk --java-options "-Xmx30G" CollectSequencingArtifactMetrics \
        -R {input.ref} \
        -I {input.tumBam} \
        -O data/preAdap/{wildcards.tumFragment} \
        -VALIDATION_STRINGENCY LENIENT \
        -TMP_DIR="tmp7/{wildcards.tumFragment}"
        '''

rule FilterByOrientationBias:
    input:
        filterVar="data/filteredVC_one/{tumFragment}.vcf.gz",
        preAdap="data/preAdap/{tumFragment}.pre_adapter_detail_metrics",
    output:
        filter2Var="data/filteredVC_2nd/{tumFragment}.vcf.gz",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        gatk --java-options "-Xmx6G" FilterByOrientationBias \
        -V {input.filterVar} \
        --artifact-modes 'G/T' \
        --artifact-modes 'C/T' \
        -P {input.preAdap} \
        -O {output.filter2Var} \
        --tmp-dir="tmp8"
        '''

