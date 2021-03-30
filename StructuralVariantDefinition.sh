#!/bin/bash -x

#----------------------------------------------------------------------------------------------------------------------------
# Variables in input
#----------------------------------------------------------------------------------------------------------------------------

while [[ $# > 0 ]]
do
key="$1"

case $key in
	#-----------------------------
    # samples info
    #-----------------------------
	-c|--configsample)
    CONFIG="$2"
    shift # past argument
    ;;
    
    #-----------------------------
    # file path of project files 
    # i.e /AbsPathTo/StructuralVariantsDefinition/
    # output directory
    #-----------------------------
    -i|--filepath) 
    FILEPATH="$2"
    shift # past argument
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    ;;
    
	#-----------------------------
	# path of the share references
	#-----------------------------
	-b|--bedfile)
    BED="$2"
    shift # past argument
    ;;
    -b_pl|--bedfile_platypus)
    BED_PL="$2"
    shift # past argument
    ;;
    -f|--fasta)
    FASTA="$2"
    shift # past argument
    ;;
    
    #-----------------------------
	# path of the tools
	#-----------------------------
    -fbpath|--freebayespath) 
    FBpath="$2"
    shift # past argument
    ;;
    -gkpath|--gatkpath) 
    GKpath="$2"
    shift # past argument
    ;;
    -vdpath|--vardictpath) 
    VDpath="$2"
    shift # past argument
    ;;
    --fileR)
    FILER="$2"
    shift # past argument
    ;;
    --filePl)
    FILEPL="$2"
    shift # past argument
    ;;
    -plpath|--platypuspath) 
    PLpath="$2"
    shift # past argument
    ;;
    -btpath|--bcftoolspath) 
    BTpath="$2"
    shift # past argument
    ;;
    
    #-----------------------------
    # parameters of the variant callers
    #-----------------------------
    -th|--threads)
    TH="$2"
    shift # past argument
    ;;
    -R|--ram)
    RAM="$2"
    shift # past argument
    ;;
    --MIN_MQ)
    MIN_MQ="$2"
    shift # past argument
    ;;
    --MIN_BQ)
    MIN_BQ="$2"
    shift # past argument
    ;;
    --MIN_AF)
    MIN_AF="$2"
    shift # past argument
    ;;
    --MIN_AO)
    MIN_AO="$2"
    shift # past argument
    ;;
    --MIN_COV)
    MIN_COV="$2"
    shift # past argument
    ;;
    --MAX_DEPTH)
    MAX_DEPTH="$2"
    shift # past argument
    ;;
    
    #-----------------------------
    # settings for features extraction
    #-----------------------------
    --DP_expected_mean)
    EMDP="$2"
    shift # past argument
    ;;
    
    #-----------------------------
	# settings for filter file
	#-----------------------------
    -af|--minallfreq)
    AF="$2"
    shift # past argument
    ;;
    -dp|--mindepth)
    DP="$2"
    shift # past argument
    ;;
    
    #-----------------------------
    # settings for AllData filter
    #-----------------------------
    --AFallData)
    AFallData="$2"
    shift # past argument
    ;;
    --MaxStrBias)
    MaxStrBias="$2"
    shift # past argument
    ;;
    --MinLengthINDELallData)
    MinLengthINDELallData="$2"
    shift # past argument
    ;;
    --NumCallersallData)
    NumCallersallData="$2"
    shift # past argument
    ;;
    --minReadsAllDef)
    MinReadsAllDef="$2"
    shift # past argument
    ;;
    --minLengthAllele)
    MinLengthAllele="$2"
    shift # past argument
    ;;
    --thresholdSimilarity)
    ThrSimilarity="$2"
    shift # past argument
    ;;
    
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done


#----------------------------------------------------------------------------------------------------------------------------
# Creation of the output directories
#----------------------------------------------------------------------------------------------------------------------------

echo -e "\n\n!!!!!!!!PIPELINE STARTS!!!!!!!!\n"

cd $OUTPUT
datestart=$(date)

var=0
while [ -d TEST$var ]; do
    ((var++))
done 

mkdir TEST$var
cd $OUTPUT/TEST$var/
clear

mkdir DepthOfCoverage
mkdir OutputPipeline

#----------------------------------------------------------------------------------------------------------------------------
# Read configuration sample file to create the samplelist.list file and the match.txt file
#----------------------------------------------------------------------------------------------------------------------------

while IFS='' read -r line || [[ -n "$line" ]]; do

#sample_id=$(echo $line | cut -f1 -d:);
#bam=$(echo $line | cut -f2 -d:);
#bam_file_name=${bam##*/};
#sample_name_in_bamfile=$(echo $bam_file_name | cut -f1 -d.);
#echo $sample_name_in_bamfile':'$sample_id >> $OUTPUT/TEST$var/match.txt
#echo $bam >> $OUTPUT/TEST$var/sampleList.list

bam=$(echo $line );
sample_id=$(samtools view -H $bam | grep "^@RG" | cut -f2 -d$'\t' | cut -f2 -d':' );
#echo $sample_id >> $OUTPUT/TEST$var/sampleIDlist.txt
echo $sample_id':'$sample_id >> $OUTPUT/TEST$var/match.txt

done < $CONFIG

#----------------------------------------------------------------------------------------------------------------------------
# Pipeline step (a) Coverage evaluation. For all samples in analysis evaluation of the depth of coverage using GATK
#----------------------------------------------------------------------------------------------------------------------------

echo -e "\n!!!!!!!!!!!!!!!!!! GATK DepthOfCoverage evaluation!!!!!!!!!!!!!!!\n";

java7 -jar $GKpath \
   -T DepthOfCoverage \
   -R $FASTA \
   -o $OUTPUT/TEST$var/DepthOfCoverage/Depth \
   -I $CONFIG \
   -L $BED \
   -mmq $MIN_MQ
   
echo -e "GATK DepthOfCoverage --> DONE\n\n"

#----------------------------------------------------------------------------------------------------------------------------
# Creation of the single sample directories for step (b) and (c)
#----------------------------------------------------------------------------------------------------------------------------

declare -a samples
samples=()

while IFS='' read -r line || [[ -n "$line" ]]; do

cd $OUTPUT/TEST$var/

START=$(date +%s);

bam=$(echo $line );
sample_id=$(samtools view -H $bam | grep "^@RG" | cut -f2 -d$'\t' | cut -f2 -d':' );

#sample_id=$(echo $line | cut -f1 -d:);
#bam=$(echo $line | cut -f2 -d:);

samples=("${samples[@]}" "$sample_id");

mkdir $sample_id;

cd $OUTPUT/TEST$var/$sample_id;	

mkdir variantCalling
mkdir fixedVariantCalling
mkdir filteredVcf

cd $OUTPUT/TEST$var/$sample_id/variantCalling;

#----------------------------------------------------------------------------------------------------------------------------
# Pipeline step (b) Variant Calling for each single sample.
#----------------------------------------------------------------------------------------------------------------------------

echo -e "\n!!!!!!!!!!!!!!!!!! VARIANT CALLING $sample_id!!!!!!!!!!!!!!!\n";

bash $FILEPATH/variant_calling.sh \
-fbpath $FBpath -gkpath $GKpath -plpath $PLpath -vdpath $VDpath --fileR $FILER --filePl $FILEPL \
-s $sample_id -d $OUTPUT/TEST$var/$sample_id/ \
-n $bam -b $BED -b_pl $BED_PL -f $FASTA \
-th $TH --ram $RAM \
--MIN_MQ $MIN_MQ --MIN_BQ $MIN_BQ --MIN_AF $MIN_AF --MIN_AO $MIN_AO --MIN_COV $MIN_COV --MAX_DEPTH $MAX_DEPTH

echo -e "VARIANT CALLING $sample_id --> DONE\n\n"

#----------------------------------------------------------------------------------------------------------------------------
# Pipeline step (c) Post Variant Calling Processing for each single sample.
#----------------------------------------------------------------------------------------------------------------------------

echo -e "\n!!!!!!!!!!!!!!!!!! POST VARIANT CALLING PROCESSING $sample_id!!!!!!!!!!!!!!!\n";

echo -e "\nFixing of the VCF files...\n"
python $FILEPATH/header_fix_freebayes.py -f $OUTPUT/TEST$var/$sample_id/variantCalling/freebayes/freebayes.vcf > $OUTPUT/TEST$var/$sample_id/variantCalling/freebayes/headerFixed.freebayes.vcf
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/freebayes/headerFixed.freebayes.vcf -m -any > $OUTPUT/TEST$var/$sample_id/variantCalling/freebayes/bcf_split.freebayes.vcf
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/freebayes/bcf_split.freebayes.vcf  > $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/freebayes.vcf;

$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/vardict/vardict.vcf -m -any > $OUTPUT/TEST$var/$sample_id/variantCalling/vardict/bcf_split.vardict.vcf 
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/vardict/bcf_split.vardict.vcf  > $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/vardict.vcf 

python $FILEPATH/header_fix_platypus.py -f $OUTPUT/TEST$var/$sample_id/variantCalling/platypus/platypus.vcf > $OUTPUT/TEST$var/$sample_id/variantCalling/platypus/headerFixed.platypus.vcf
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/platypus/headerFixed.platypus.vcf -m -any > $OUTPUT/TEST$var/$sample_id/variantCalling/platypus/bcf_split.platypus.vcf 
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/platypus/bcf_split.platypus.vcf  > $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/platypus.vcf 

python $FILEPATH/header_fix_gatk.py -f $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/gatk_indel.vcf > $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/headerFixed.gatk.indel.vcf
python $FILEPATH/header_fix_gatk.py -f $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/gatk_snp.vcf > $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/headerFixed.gatk.snp.vcf
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/headerFixed.gatk.indel.vcf -m -any > $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/bcf_split.gatk.indel.vcf 
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/bcf_split.gatk.indel.vcf  > $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/gatk.indel.vcf 
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/headerFixed.gatk.snp.vcf -m -any > $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/bcf_split.gatk.snp.vcf 
$BTpath norm -f $FASTA $OUTPUT/TEST$var/$sample_id/variantCalling/gatk/bcf_split.gatk.snp.vcf  > $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/gatk.snp.vcf 

#----------------------------------------------------------------------------------------------------------------------------

cd $OUTPUT/TEST$var/$sample_id/fixedVariantCalling;
mkdir variants

echo -e "\nFeatures Extraction...\n"

python $FILEPATH/feature_extraction_snp.py \
-fb $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/freebayes.vcf \
-vd $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/vardict.vcf \
-pl $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/platypus.vcf \
-gk $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/gatk.snp.vcf \
-s $sample_id -mDP $EMDP -mqVDThreshold $MIN_MQ \
> $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/variants/snp_$sample_id.txt

python $FILEPATH/feature_extraction_indel.py \
-fb $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/freebayes.vcf \
-vd $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/vardict.vcf \
-pl $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/platypus.vcf \
-gk $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/gatk.indel.vcf \
-s $sample_id -mDP $EMDP -mqVDThreshold $MIN_MQ \
> $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/variants/indel_$sample_id.txt

#--------------------------------------------------------------------------------------------------------------------------------------------------------

echo -e "\nVariants filtering...\n"

python $FILEPATH/filters.py -i $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/variants/snp_$sample_id.txt -f $AF -d $DP\
> $OUTPUT/TEST$var/$sample_id/filteredVcf/snp_finalFiltered_$sample_id.txt
python $FILEPATH/filters.py -i $OUTPUT/TEST$var/$sample_id/fixedVariantCalling/variants/indel_$sample_id.txt -f $AF -d $DP\
> $OUTPUT/TEST$var/$sample_id/filteredVcf/indel_finalFiltered_$sample_id.txt

#--------------------------------------------------------------------------------------------------------------------------------------------------------

echo -e "\nVariants filtering for the output results...\n"

python $FILEPATH/filters_forAllData.py -i $OUTPUT/TEST$var/$sample_id/filteredVcf/indel_finalFiltered_$sample_id.txt \
-f $AFallData -nt $MinLengthINDELallData -calls $NumCallersallData -sb $MaxStrBias -vt INDEL\
> $OUTPUT/TEST$var/$sample_id/filteredVcf/indel_allData_$sample_id.txt
python $FILEPATH/filters_forAllData.py -i $OUTPUT/TEST$var/$sample_id/filteredVcf/snp_finalFiltered_$sample_id.txt \
-f $AFallData -calls $NumCallersallData -sb $MaxStrBias -vt SNP\
> $OUTPUT/TEST$var/$sample_id/filteredVcf/snp_allData_$sample_id.txt


echo $OUTPUT/TEST$var/$sample_id/filteredVcf/indel_allData_$sample_id.txt >> $OUTPUT/TEST$var/IndelVariantsPathList_AllData.txt
echo $OUTPUT/TEST$var/$sample_id/filteredVcf/snp_allData_$sample_id.txt >> $OUTPUT/TEST$var/SnpVariantsPathList_AllData.txt

#--------------------------------------------------------------------------------------------------------------------------------------------------------

echo -e "POST VARIANT CALLING PROCESSING $sample_id --> DONE\n\n"
done < $CONFIG

cat $OUTPUT/TEST$var/IndelVariantsPathList_AllData.txt $OUTPUT/TEST$var/SnpVariantsPathList_AllData.txt >> $OUTPUT/TEST$var/VariantPathList_AllData.txt
#----------------------------------------------------------------------------------------------------------------------------
# Pipeline step (d) Allele Determination for all the samples.
#----------------------------------------------------------------------------------------------------------------------------

echo -e "\nAllele Determination...\n"

python $FILEPATH/AlleleDetermination.py \
-GATKdp $OUTPUT/TEST$var/DepthOfCoverage/Depth \
-RegionsInfo $BED \
-MatchID $OUTPUT/TEST$var/match.txt \
-VariantPathList $OUTPUT/TEST$var/VariantPathList_AllData.txt \
-minReadsAllDef $MinReadsAllDef \
-minLengthAllele $MinLengthAllele \
-thresholdSimilarity $ThrSimilarity \
-outputPath $OUTPUT/TEST$var/OutputPipeline/

echo -e "\nVariants statistics...\n"

python $FILEPATH/VariantsStatistics.py \
-pathList $OUTPUT/TEST$var/IndelVariantsPathList_AllData.txt \
-outputPath $OUTPUT/TEST$var/OutputPipeline/

python $FILEPATH/VariantsStatistics.py \
-pathList $OUTPUT/TEST$var/SnpVariantsPathList_AllData.txt \
-outputPath $OUTPUT/TEST$var/OutputPipeline/

python $FILEPATH/AdLocusToVariants.py \
-RegionsInfo $BED \
-IndelVariants $OUTPUT/TEST$var/OutputPipeline/indelStatistics.txt \
-SnpVariants $OUTPUT/TEST$var/OutputPipeline/snpStatistics.txt \
-outputPath $OUTPUT/TEST$var/OutputPipeline/

cat $OUTPUT/TEST$var/OutputPipeline/indelAssignedToLocusStatistic.txt $OUTPUT/TEST$var/OutputPipeline/snpAssignedToLocusStatistic.txt \
>> $OUTPUT/TEST$var/OutputPipeline/All_polimorphisms_AssignToLocusStatistic.txt

echo -e "\nGrouping alleles in categories...\n"

python $FILEPATH/Grouping.py \
-RegionsInfo $BED \
-Variants $OUTPUT/TEST$var/OutputPipeline/All_polimorphisms_AssignToLocusStatistic.txt \
-AlleleFile $OUTPUT/TEST$var/OutputPipeline/Alleles.txt \
-outputPath $OUTPUT/TEST$var/OutputPipeline/

#----------------------------------------------------------------------------------------------------------------------------
# End of the pipeline
#----------------------------------------------------------------------------------------------------------------------------

END=$(date +%s);
echo $((END-START)) | awk '{print int($1/60)":"int($1%60)}'

echo -e "\n\n!!!!!!!!PIPELINE ENDS!!!!!!!!\n"
