mkdir dogTumors && cd dogTumors
work_dir=$(pwd)

mkdir -p $work_dir/{data,scripts}

cd $work_dir/data
wget -r --no-parent -A '*.gz' http://slimsdata.genomecenter.ucdavis.edu/Data/8rhzkpx/Un/Project_TBDY_L1_DTDB05_H1220P_York_1/
wget -r --no-parent -A '*.gz' http://slimsdata.genomecenter.ucdavis.edu/Data/fhnqivtrq/Un/Project_TBDY_L1_DTDB06_H1220P_York_2/
wget -r --no-parent -A '*.gz' http://slimsdata.genomecenter.ucdavis.edu/Data/yoe0kg52dk/Un/Project_TBDY_L1_DTDB07_H1220P_York_3/
wget -r --no-parent -A '*.gz' http://slimsdata.genomecenter.ucdavis.edu/Data/mtun2fifr2/Un/Project_TBDY_L1_DTDB08_H1220P_York_4/

mkdir fastq
for dir in slimsdata.genomecenter.ucdavis.edu/Data/*/Un/*;do
  mv dir fastq/.;
  rmdir -p $(dirname $dir);
done

## remove the undetermined samples from the analysis
for f in data/fastq/*/Undetermined*;do f2=$f.doNotUse;mv $f $f2;done 

## prepare the exome bed file
cat scripts/OID45779_CanFam31_06Mar2017_capture_targets.bed | awk '{if($1 ~ "chrUn_"){next;}else{print $0;}}' > refGenome/intervals.bed

## cp -R $Bovine_seq/hpcc $work_dir/.
export PATH=$HOME/miniconda3/bin:$(pwd)/hpcc:$PATH
source activate snakemake
conda install -c bioconda fastqc  ## to get FastQC v0.11.7 that can work with Novaseq
. hpcc/submit.sh 


## move successful log files into log dir
mkdir -p logs/FastqToSam
for x in snakejob.FastqToSam.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/FastqToSam/.;fi;done 

mkdir -p logs/MarkIlluminaAdapters
for x in snakejob.MarkIlluminaAdapters.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/MarkIlluminaAdapters/.;fi;done 

mkdir -p logs/align
for x in snakejob.align.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/align/.;fi;done

mkdir -p logs/MarkDuplicates
for x in snakejob.MarkDuplicates.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/MarkDuplicates/.;fi;done

mkdir -p logs/BaseRecalib
for x in snakejob.BaseRecalib.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/BaseRecalib/.;fi;done

mkdir -p logs/ApplyBQSR
for x in snakejob.ApplyBQSR.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/ApplyBQSR/.;fi;done

mkdir -p logs/Mutect2_for_PON
for x in snakejob.Mutect2_for_PON.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/Mutect2_for_PON/.;fi;done

mkdir -p logs/somaticVarCalling
for x in snakejob.somaticVarCalling.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/somaticVarCalling/.;fi;done

mkdir -p logs/GetPileupSummaries
for x in snakejob.GetPileupSummaries.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/GetPileupSummaries/.;fi;done

mkdir -p logs/CalculateContamination
for x in snakejob.CalculateContamination.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/CalculateContamination/.;fi;done

mkdir -p logs/CollectSequencingArtifactMetrics
for x in snakejob.CollectSequencingArtifactMetrics.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/CollectSequencingArtifactMetrics/.;fi;done

mkdir -p logs/FilterMutectCalls
for x in snakejob.FilterMutectCalls.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/FilterMutectCalls/.;fi;done

mkdir -p logs/FilterByOrientationBias
for x in snakejob.FilterByOrientationBias.*.e*;do n=$(grep "1 of 1 steps (100%) done" $x | wc -l);if [ $n -eq 1 ];then echo $n;x2=$(echo $x | sed 's/sh\.e/sh\.o/');echo $x $x2;mv $x $x2 logs/FilterByOrientationBias/.;fi;done


## direct run of the gatk engine without snakemake
module load Java/jdk1.8.0
source activate gatk
gatk --java-options "-Xmx15G" AnalyzeCovariates -bqsr "data/recalib/Project_TBDY_L1_DTDB06_H1220P_York_2/GN11_S59_L002.txt" -plots BQSR.pdf

## compare the 100 genome VCF to the dbvar database list you used for recalibration
module load vcftools/0.1.14
dogSeq="/mnt/home/mansourt/Tamer/dogSeq/"
dogSeqVarPass=$dogSeq/varResults/GenotypeGVCFs_output_max50.pass.combinedFiltered.vcf
vcftools --vcf $dogSeqVarPass --diff knowVar/canis_familiaris_fixedChrNames_sorted.vcf --diff-site --out pass_dogSeqVsdbSNP
dogSeqVarRaw=$dogSeq/varResults/GenotypeGVCFs_output_max50.vcf
vcftools --vcf $dogSeqVarRaw --diff knowVar/canis_familiaris_fixedChrNames_sorted.vcf --diff-site --out raw_dogSeqVsdbSNP
#MPS="/mnt/home/mansourt/Tamer/MPS"
#dogSeqVarRawGZ=$MPS/knownVar/control_dogSeq.vcf.gz
#vcftools --gzvcf $dogSeqVarRawGZ --diff knowVar/canis_familiaris_fixedChrNames_sorted.vcf --diff-site --out raw_dogSeqVsdbSNP

## compare common population alleles to ENSEMBL var
vcftools --vcf data/population/commonAlt.vcf --diff knowVar/canis_familiaris_fixedChrNames_sorted.vcf --diff-site --out data/population/commonAltVsdbSNP 
cat data/population/commonAltVsdbSNP.diff.sites_in_files | awk '{if($4=="O")print $0}' > data/population/commonAltVsdbSNP.diff.sites_in_files.bad   

## analysis of the variants file
conda install -c bioconda bcftools
bcftools merge --merge all -Ov -o data/filteredVC_one_all.vcf --force-samples data/filteredVC_one/*/*.vcf.gz
#bcftools merge --merge all -Ov -o data/filteredVC_2nd_all.vcf --force-samples data/filteredVC_2nd/*/*.vcf.gz
#module load vcftools/0.1.15-Feb2018
module load VCFtools/0.1.15-Perl-5.26.1
cat data/filteredVC_one_all.vcf | fill-an-ac > data/filteredVC_one_all_AC.vcf

for type in {GT,MT,OT};do echo $type;#done
  ls -tral data/filteredVC_one/*/${type}*.vcf.gz | wc -l
  label=data/filteredVC_one_${type}
  bcftools merge --merge all -Ov -o ${label}.vcf --force-samples data/filteredVC_one/*/${type}*.vcf.gz ## 39
  cat ${label}.vcf | fill-an-ac > ${label}_AC.vcf
  grep -v "^#" ${label}_AC.vcf | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$4,$5,$8}' | awk -F";" '{print $1}' > ${label}_AC.reduced
  cat ${label}_AC.reduced | awk '{print $5}' | sort | uniq -c | sort -k1,1nr > ${label}_AC.summary
  grep -v "^#" ${label}_AC.vcf | grep "PASS" | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$4,$5,$8}' | awk -F";" '{print $1}' > ${label}_AC.reduced_pass
  cat ${label}_AC.reduced_pass | awk '{print $5}' | sort | uniq -c | sort -k1,1nr > ${label}_AC.summary_pass
done

#grep "PASS" data/filteredVC_one_GT_AC.vcf | grep "AC=3,1,1,1;" > hits/GT.AC.pass.6
#grep "PASS" data/filteredVC_one_GT_AC.vcf | grep "^chr13" | grep "\b467" > hits/GT.AC.pass.6.around
#cat hits/GT.AC.pass.6.around | awk '{print $2,$4,$5,$8}' | awk -F";" '{print $1}'

#grep "PASS" data/filteredVC_one_MT_AC.vcf | grep "AC=5;" > hits/MT.AC.pass.5
#grep "PASS" data/filteredVC_one_MT_AC.vcf | grep "^chr6" | grep "\b890" > hits/MT.AC.pass.5.around

#grep "PASS" data/filteredVC_one_OT_AC.vcf | grep "AC=4;" > hits/OT.AC.pass.4
#grep "PASS" data/filteredVC_one_OT_AC.vcf | grep "^chr11" | grep "\b5234" > hits/OT.AC.pass.4.around1
#grep "PASS" data/filteredVC_one_OT_AC.vcf | grep "^chr5" | grep "\b3256" > hits/OT.AC.pass.4.around2
#cat hits/OT.AC.pass.4.around2 | awk '{print $2,$4,$5,$8}' | awk -F";" '{print $1}'



# install vcfanno for VCF annotation
# vcfanno
# conda install -c bioconda vcfanno 
conda create -n vcfanno vcfanno==0.3.0
source activate vcfanno
conda install -c biobuilds tabix ## 1.6.0

## generate new dog GTF (follow the code from dogAnn)
mkdir -p scripts/UCSC_kent_commands
cp ~/Tamer/dogAnn/scripts/UCSC_kent_commands/gff3ToGenePred scripts/UCSC_kent_commands/.
cp ~/Tamer/dogAnn/scripts/mapGenome.sh scripts/.
cp ~/Tamer/dogAnn/scripts/NCBItoUCSCmapTopsl scripts/.
cp ~/Tamer/dogAnn/scripts/UCSC_kent_commands/pslToChain scripts/UCSC_kent_commands/.
cp ~/Tamer/dogAnn/scripts/UCSC_kent_commands/chainSort scripts/UCSC_kent_commands/.
cp ~/Tamer/dogAnn/scripts/UCSC_kent_commands/liftOver scripts/UCSC_kent_commands/.
cp ~/Tamer/dogAnn/scripts/UCSC_kent_commands/genePredToGtf scripts/UCSC_kent_commands/.
cp ~/Tamer/dogAnn/scripts/genePredToBed scripts/.
cd refGenome
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz ## last modified 9/5/17 (significant change from 9/18/15 in dogAnn)
gunzip ref_CanFam3.1_top_level.gff3.gz
$work_dir/scripts/UCSC_kent_commands/gff3ToGenePred -useName ref_CanFam3.1_top_level.gff3 ref_CanFam3.1_top_level.gpred
mkdir ncbi && cd ncbi
bash $work_dir/scripts/mapGenome.sh $work_dir/refGenome          ## ends by creating ncbi/NCBItoUCSC_map.sorted.chain
cd ../
$work_dir/scripts/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level.gpred $work_dir/refGenome/ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped.gpred unMapped -genePred
#$work_dir/scripts/UCSC_kent_commands/genePredToGtf file ref_CanFam3.1_top_level_mapped.gpred ref_CanFam3.1_top_level_mapped.gtf
cat ref_CanFam3.1_top_level_mapped.gpred | $work_dir/scripts/genePredToBed > ref_CanFam3.1_top_level_mapped.bed

## Extend NCBI annotation (follow the code from dogSeq) 
grep "ID=rna" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{print $9}' | awk -F "[,;]" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; Name = vars["Name"]; gene = vars["gene"]; product = vars["product"]; } { print Dbxref,Name,gene,product }' > NCBI_TransInfo.temp.trans;
grep -v "^#" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{if($3=="gene")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; gene_biotype = vars["gene_biotype"]; } { print Dbxref,gene_biotype }' > NCBI_TransInfo.temp.gene;
module load R/3.5.1-X11-20180131
Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1],header=F);data2=read.delim(args[2],header=F);'\
'dataMerge=merge(data1,data2,by="V1",all.x=F,all.y=T); colnames(dataMerge)=c("Gene_ID","Gene_biotype","Transcript_ID","Gene_Name","Description");'\
'write.table(dataMerge[,c(3,4,2,5)],"NCBI_TransInfo.txt", sep="\t", quote=F, row.names=F, col.names=T);' NCBI_TransInfo.temp.gene NCBI_TransInfo.temp.trans
rm NCBI_TransInfo.temp.*
sed -i 's/%2C/,/g' NCBI_TransInfo.txt; sed -i 's/%3B/;/g' NCBI_TransInfo.txt;
awk 'BEGIN{FS=OFS="\t";}{if(FNR==NR){a[$1]=$2;next}{if(a[$4]!="")print $1,$2,$3,a[$4],$5,$6,$7,$8,$9,$10,$11,$12}}' NCBI_TransInfo.txt ref_CanFam3.1_top_level_mapped.bed > ref_CanFam3.1_top_level_mapped_genes.bed

## annotation
cat ref_CanFam3.1_top_level_mapped.bed | sort -k1,1 -k2,2n | bgzip > ref_CanFam3.1_top_level_mapped.bed.gz
tabix -p bed ref_CanFam3.1_top_level_mapped.bed.gz
cat ref_CanFam3.1_top_level_mapped_genes.bed | sort -k1,1 -k2,2n | bgzip > ref_CanFam3.1_top_level_mapped_genes.bed.gz
tabix -p bed ref_CanFam3.1_top_level_mapped_genes.bed.gz

cd ../
for type in {GT,MT,OT};do echo $type;#done # type=OT
  label=data/filteredVC_one_${type}_AC
  vcfanno scripts/Bed_config.toml $label.vcf > $label.BedAnn.vcf
  grep -v "^#" $label.BedAnn.vcf | awk -F "\t" '{print $8}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } AC = vars["AC"]; transName = vars["transName"]; geneName = vars["geneName"]; } { print AC,transName,geneName }' > $label.BedAnn.Ann; ## AC,transName,geneName
  cat $label.BedAnn.Ann |awk -F "\t" '{print $1}' |awk -F "," '{vars=0;for(i = 1; i <= NF; ++i){vars += $i}} {print vars}' > $label.BedAnn.AC; ## AC_sum
  paste $label.BedAnn.AC $label.BedAnn.Ann > $label.BedAnn.AC_Ann ## AC_sum,AC,transName,geneName
  paste <(grep -v "^#" $label.BedAnn.vcf | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$4,$5,$7,$8}') $label.BedAnn.AC_Ann > $label.BedAnn.AC_Ann.extended 
  cat $label.BedAnn.AC_Ann.extended | awk '{if($10!=""){a[$10]+=$7;b[$10]+=1;}}END{for(i in a)print i,a[i],b[i];}' | sort -k2,2nr > $label.BedAnn.counts ## geneName, AC_sum_inAllPositions, noOfPositions
  cat $label.BedAnn.AC_Ann.extended | grep "PASS" | awk '{if($10!=""){a[$10]+=$7;b[$10]+=1;}}END{for(i in a)print i,a[i],b[i];}' | sort -k2,2nr > $label.BedAnn.pass_counts ## geneName, AC_sum_inAllPositions, noOfPositions (only in PASS records)
  awk '{print $1}' $label.BedAnn.pass_counts | grep -Fwf - $label.BedAnn.counts | sort -k2,2nr > $label.BedAnn.counts2 ## geneName, AC_sum_inAllPositions, noOfPositions (all records for genes with any PASS record in the merged file)
  awk '{if(FNR==NR){a[$1]=$2" "$3;next}{print $1,$2,$3,a[$1];b[$1]+=1;}}END{for(i in a){if(b[i]=="")print i,"- -",a[i]}}' $label.BedAnn.counts2 $label.BedAnn.pass_counts > $label.BedAnn.pass_count_extended
done

## specific gene records
for type in {GT,MT,OT};do echo $type;#done # type=GT
  mkdir -p ${type}
  label=data/filteredVC_one_${type}_AC
  cat $label.BedAnn.pass_counts | awk '{print $1}' | while read gene;do #=PDGFRA #LOC102152242
    grep -w $gene $label.BedAnn.AC_Ann.extended > ${type}/$gene.AC_Ann.extended  ## annotated version of the merged vcf
    chr=$(head -n1 ${type}/$gene.AC_Ann.extended | awk '{print $1}')
    > ${type}/$gene.filteredAfterMereg.inVCF
    > ${type}/$gene.passAfterMereg.inVCF
    for f in data/filteredVC_one/*/${type}*.vcf.gz;do 
      sample=$(echo $f | awk -F"/" '{print $4}' | awk -F"." '{print $1}')
      zcat $f |grep -Fwf <(grep -wv PASS ${type}/$gene.AC_Ann.extended | awk '{print $2}') |grep -w $chr |\
awk -v s=$sample 'BEGIN{FS=OFS="\t";}{print s,$1,$2,$4,$5,$7,$8}' >> ${type}/$gene.filteredAfterMereg.inVCF
      zcat $f |grep -Fwf <(grep -w PASS ${type}/$gene.AC_Ann.extended | awk '{print $2}') |grep -w $chr |\
awk -v s=$sample 'BEGIN{FS=OFS="\t";}{print s,$1,$2,$4,$5,$7,$8}' >> ${type}/$gene.passAfterMereg.inVCF
    done
  done
done


for type in {GT,MT,OT};do echo $type;#done # type=OT
  label=data/filteredVC_one_${type}_AC
  cat $label.BedAnn.pass_counts | awk '{print $1}' | while read gene;do
    grep -w PASS ${type}/$gene.filteredAfterMereg.inVCF | awk '{print $3}' | sort | uniq -c | awk -v g=$gene '{a+=$1}END{print g,a,NR}'
  done > ${type}.FailAfterMerge
  awk '{if(FNR==NR){a[$1]=$2" "$3;next}{print $1,$2,$3,$4,$5,a[$1];}}' ${type}.FailAfterMerge $label.BedAnn.pass_count_extended > $label.BedAnn.pass_count_extended2
done
