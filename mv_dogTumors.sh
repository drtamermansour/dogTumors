mkdir -p Tamer2/dogTumors/data
cp -R Tamer/dogTumors/data/fastq/ Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/uBAM Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/adap_uBAM Tamer2/dogTumors/data/.
mkdir Tamer2/dogTumors/qc
cp -R Tamer/dogTumors/qc/fastqc Tamer2/dogTumors/qc/.
cp -R Tamer/dogTumors/qc/multiqc Tamer2/dogTumors/qc/.
mkdir Tamer2/dogTumors/refGenome
cp Tamer/dogTumors/refGenome/canFam3*.fa Tamer2/dogTumors/refGenome/.
cp -R Tamer/dogTumors/refGenome/*Index Tamer2/dogTumors/refGenome/.
cp -R Tamer/dogTumors/data/mapped_reads Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/dedup Tamer2/dogTumors/data/.
mkdir Tamer2/dogTumors/knowVar
cp -R Tamer/dogTumors/knowVar/canis_familiaris*.vcf Tamer2/dogTumors/knowVar/.
cp -R Tamer/dogTumors/knowVar/canis_familiaris*.vcf.idx Tamer2/dogTumors/knowVar/.
for dir in Tamer/dogTumors/data/recalib/*;do 
 newdir=$(echo $dir | sed 's/Tamer/Tamer2/');
 mkdir -p $newdir;
 cp -R $dir/*.txt $newdir/.
 cp -R $dir/*.bam $newdir/.
 cp -R $dir/*.bai $newdir/.
done
cp -R Tamer/dogTumors/data/PON Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/PON.vcf.gz Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/rawVC Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/realign Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/population Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/pileupsum Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/cont Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/seg Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/filteredVC_one Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/preAdap Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/data/filteredVC_2nd Tamer2/dogTumors/data/.
cp -R Tamer/dogTumors/logs Tamer2/dogTumors/.
cp -R Tamer/dogTumors/failed_logs Tamer2/dogTumors/.
cp Tamer/dogTumors/data/filteredVC_one_all*.vcf Tamer2/dogTumors/data/.
for type in {GT,MT,OT};do echo $type;#done
  label=data/filteredVC_one_${type}
  cp Tamer/dogTumors/${label}.vcf Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/${label}_AC.vcf Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/${label}_AC.reduced Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/${label}_AC.summary Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/${label}_AC.reduced_pass Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/${label}_AC.summary_pass Tamer2/dogTumors/data/.
  label=data/filteredVC_one_${type}_AC
  cp Tamer/dogTumors/$label.BedAnn.vcf Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.Ann Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.AC Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.AC_Ann Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.AC_Ann.extended Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.counts Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.pass_counts Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.counts2 Tamer2/dogTumors/data/.
  cp Tamer/dogTumors/$label.BedAnn.pass_count_extended Tamer2/dogTumors/data/.
  cp -R Tamer/dogTumors/${type} Tamer2/dogTumors/.
  cp -R Tamer/dogTumors/${type}_10genes Tamer2/dogTumors/.
  cp Tamer/dogTumors/${type}.FailAfterMerge Tamer2/dogTumors/.
  cp Tamer/dogTumors/$label.BedAnn.pass_count_extended2 Tamer2/dogTumors/data/.
done
cp -R Tamer/dogTumors/scripts Tamer2/dogTumors/.
cp Tamer/dogTumors/refGenome/ref_CanFam3.1_top_level* Tamer2/dogTumors/refGenome/.
cp Tamer/dogTumors/refGenome/{NCBI_TransInfo.txt,intervals.bed,unMapped} Tamer2/dogTumors/refGenome/.
cp -R Tamer/dogTumors/refGenome/ncbi Tamer2/dogTumors/refGenome/.  
cp Tamer/dogTumors/{raw_dogSeqVsdbSNP.*,pass_dogSeqVsdbSNP.*,BQSR.pdf,contamination_table,all.vcf} Tamer2/dogTumors/.
cp Tamer/dogTumors/{Snakefile,main.sh,main_GATK4.sh} Tamer2/dogTumors/.  
cp -R Tamer/dogTumors/hpcc Tamer2/dogTumors/.  
