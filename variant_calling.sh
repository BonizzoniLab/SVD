#----------------------------------------------------------------------------------------------------------------------------
# Variables in input
#----------------------------------------------------------------------------------------------------------------------------

while [[ $# > 0 ]]
do
key="$1"

case $key in
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
    -plpath|--platypuspath) 
    PLpath="$2"
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
    
  	#-----------------------------
	# path of the working directory
	#-----------------------------
    -d|--workdir)
    WORKDIR="$2"
    shift # past argument
    ;;
    
    #-----------------------------
	# sample and reference files
	#-----------------------------
    -s|--sample_name)
    SAMPLE_NAME="$2"
    shift # past argument
    ;;
    -n|--bamfile)
    BAM="$2"
    shift # past argument
    ;;
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
    # parameters of the variant callers
    #-----------------------------
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
   -th|--threads)
    TH="$2"
    shift # past argument
    ;;
    -R|--ram)
    RAM="$2"
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

mkdir freebayes;
mkdir vardict;
mkdir platypus;
mkdir gatk;

#----------------------------------------------------------------------------------------------------------------------------

echo -e "Run Freebayes...\n"

$FBpath -f $FASTA \
--min-mapping-quality $MIN_MQ --min-base-quality $MIN_BQ \
--min-alternate-fraction $MIN_AF --min-alternate-count $MIN_AO --min-coverage $MIN_COV --genotype-qualities \
--report-genotype-likelihood-max \
 -t $BED $BAM \
> $WORKDIR/variantCalling/freebayes/freebayes.vcf;

#----------------------------------------------------------------------------------------------------------------------------

echo -e "Run Platypus...\n"

export LD_LIBRARY_PATH=/usr/local/lib
export C_INCLUDE_PATH=/usr/local/include

python $PLpath callVariants --bamFiles=$BAM --refFile=$FASTA --output=$WORKDIR/variantCalling/platypus/platypus.vcf \
--minReads=$MIN_AO --minMapQual=$MIN_MQ --minBaseQual=$MIN_BQ --nCPU $TH --regions=$BED_PL;

#----------------------------------------------------------------------------------------------------------------------------

echo -e "Run Gatk indel...\n"

java7 -Xmx$RAM -jar $GKpath -T UnifiedGenotyper -I $BAM \
-R $FASTA --downsample_to_coverage $MAX_DEPTH --min_base_quality_score $MIN_BQ -L $BED --num_threads $TH \
--genotype_likelihoods_model INDEL -o $WORKDIR/variantCalling/gatk/gatk_indel.vcf

echo -e "Run Gatk snp...\n" 

java7 -Xmx$RAM -jar $GKpath -T UnifiedGenotyper -I $BAM \
-R $FASTA --downsample_to_coverage $MAX_DEPTH --min_base_quality_score $MIN_BQ -L $BED --num_threads $TH \
--genotype_likelihoods_model SNP -o $WORKDIR/variantCalling/gatk/gatk_snp.vcf

#----------------------------------------------------------------------------------------------------------------------------

echo -e "Run Vardict..."
 
$VDpath -G $FASTA \
 -f $MIN_AF -N $SAMPLE_NAME -b $BAM -c 1 -S 2 -E 3 -g 4 \
 $BED | $FILER | $FILEPL \
 -N $SAMPLE_NAME -E -f $MIN_AF > $WORKDIR/variantCalling/vardict/vardict.vcf 
 
#----------------------------------------------------------------------------------------------------------------------------
