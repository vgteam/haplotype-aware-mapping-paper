#!/usr/bin/env bash
# test-population-mapping.sh: evaluate the performance impact of population-aware mapping using toil-vg mapeval on AWS

set -ex
# We need extglob support for this script
shopt -s extglob

# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/adamnovak/toil-vg.git@13fd01c03f520c622f729205452ba21ac4b49c07#egg=toil-vg"

# What Docker registry can the corresponding dashboard containers (Grafana, etc.) be obtained from?
TOIL_DOCKER_REGISTRY="quay.io/ucsc_cgl"
# What Toil appliance should we use? Ought to match the locally installed Toil,
# but can't quite if the locally installed Toil is locally modified or
# installed from a non-release Git commit.
TOIL_APPLIANCE_SELF="${TOIL_DOCKER_REGISTRY}/toil:3.18.0a1.dev2629-14af7326ed3531d8fa2a38fe152256e8a208ed52"


# What version of awscli do we use? This has to be compatible with the
# botocore/boto3 that toil-vg and toil can agree on, and each awscli version
# seems to require exactly one botocore version, and pip can't do something
# sensible like find the one that goes with the botocore we need to have when
# we just ask for "awscli". So we work out the right version manually and live
# in hope that <https://github.com/pypa/pip/issues/988> will eventually stop
# plaguing our lab.
# awscli 1.15.85 should pull in botocore 1.10.84 
AWSCLI_PACKAGE="awscli==1.15.85"

# What vg should we use?
# Update this tag to change the Docker that will be used by a restart.
# Just editing the script won't do it; the tag name lives in the Toil job store.
# To update, do something like:
# docker pull quay.io/vgteam/vg:dev-v1.8.0-142-g758c92ec-t190-run
# docker tag quay.io/vgteam/vg:dev-v1.8.0-142-g758c92ec-t190-run quay.io/adamnovak/vg:wholegenome
# docker push quay.io/adamnovak/vg:wholegenome
VG_DOCKER_OPTS=("--vg_docker" "quay.io/adamnovak/vg:wholegenome")

# What node types should we use?
# Comma-separated, with :bid-in-dollars after the name for spot nodes
# We need non-preemptable i3.4xlarge at least to get ~3.8TB storage available so the GCSA indexing jobs will have somewhere to run.
# And we also need more memory (?) than that so some of the later jobs will run.
# Suggested: i3.8xlarge which is worth ~0.70-0.80 (2.something on-demand), and r4.8xlarge which is worth ~0.60
# If using r4 we also need to make sure --nodeStorage, which sets EBS volume size for EBS nodes, is sufficient.
# We can also use r3.8xlarge which are cheaper (~0.50) and come with faster disks anyway.
# But you need a Toil which can tell them apart
NODE_TYPES="i3.8xlarge,r3.8xlarge:0.60"
# How many nodes should we use at most per type?
# Also comma-separated.
MAX_NODES="10,50"
# And at least per type? (Should probably be 0)
# Also comma-separated.
MIN_NODES="0,0"

# What's our unique run ID? Should be lower-case and start with a letter for maximum compatibility.
# See <https://gist.github.com/earthgecko/3089509>
RUN_ID="run$(cat /dev/urandom | LC_CTYPE=C tr -dc 'a-z0-9' | fold -w 32 | head -n 1)"

# What cluster should we use?
CLUSTER_NAME="${RUN_ID}"
# Is our cluster just for this run, or persistent for multiple runs?
PERSISTENT_CLUSTER=0
# Set this to 1 to delete the job store on exit
CLEAN_UP_JOB_TREE=0

# What stage should we restart from, if any?
RESTART_STAGE=""

# What named graph region should we be operating on?
# We can do "21", "whole-genome", "MHC", etc.
INPUT_DATA_MODE="21"

# Set a FASTQ to model reads after
TRAINING_FASTQ="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz"
# And a read simulation seed
READ_SEED="90"

# Define the sample to use for synthesizing reads
SAMPLE_NAME="NA12878"
# And options to filter them out of graphs, along with people related to them
FILTER_OPTS=("--filter_ceph" "--filter_samples" "${SAMPLE_NAME}")

# What min allele frequency limits do we use?
# These come out as minaf, minaf1, minaf2, minaf3 (number = # of zeroes)
MIN_AFS=("0.0335570469" "0.1" "0.01" "0.001")

# What GBWT penalties do we use?
# These come out as numbers in the condition tags.
GBWT_RECOMBINATION_PENALTIES=(5.0 10.0 19.0 20.7 22.0 40.0)

# Put this in front of commands to do or not do them, depending on if we are doing a dry run or not
PREFIX=""

function url_to_store() {
    # Convert an s3:// URL to an output store specifier like aws:us-west-2:bucket/path/to/thing 
    local URL="${1}"
    shift
    echo "${URL}" | sed 's|^s3://|aws:us-west-2:|'
}

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] KEYPAIR_NAME GRAPHS_URL ALIGNMENTS_URL CALLS_URL \n"
    printf "\tKEYPAIR_NAME\tSSH keypair in Amazon us-west-2 region to use for cluster authentication\n"
    printf "\tGRAPHS_URL\tS3 URL where generated graphs, indexes, and reads should be cached\n"
    printf "\tALIGNMENTS_URL\tS3 URL where mapped reads and statistics should be cached\n"
    printf "\tCALLS_URL\tS3 URL where variant calls and statistics should be deposited\n"
    printf "Options:\n\n"
    printf "\t-d\tDo a dry run.\n"
    printf "\t-S\tClean up job store on exit, even if the run failed.\n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg.\n"
    printf "\t-D REGISTRY\tUse the given Docker registry for finding containers (default: ${TOIL_DOCKER_REGISTRY}).\n"
    printf "\t-t CONTAINER\tUse the given Toil container in the cluster (default: ${TOIL_APPLIANCE_SELF}).\n"
    printf "\t-c CLUSTER\tUse the given persistent Toil cluster, which will be created if not present.\n"
    printf "\t-v DOCKER\tUse the given Docker image specifier for vg.\n"
    printf "\t-R RUN_ID\tUse or restart the given run ID.\n"
    printf "\t-s STAGE\tRestart the given stage (construct, construct-sample, sim, bwa-index, map-sim, map-real, call-sim, call-real).\n"
    printf "\t-r REGION\tRun on the given named region (21, MHC, BRCA1).\n"
    exit 1
}

while getopts "hdSp:D:t:c:v:R:s:r:" o; do
    case "${o}" in
        d)
            PREFIX="echo"
            ;;
        S)
            CLEAN_UP_JOB_TREE=1
            ;;
        p)
            TOIL_VG_PACKAGE="${OPTARG}"
            ;;
        D)
            TOIL_DOCKER_REGISTRY="${OPTARG}"
            # TODO: Rewrite TOIL_APPLIANCE_SELF if it isn't going to be overridden.
            ;;
        t)
            TOIL_APPLIANCE_SELF="${OPTARG}"
            ;;
        c)
            CLUSTER_NAME="${OPTARG}"
            PERSISTENT_CLUSTER=1
            ;;
        v)
            VG_DOCKER_OPTS=("--vg_docker" "${OPTARG}")
            ;;
        R)
            # This doesn't change the cluster name, which will still be the old run ID if not manually set.
            # That's probably fine.
            RUN_ID="${OPTARG}"
            ;;
        s)
            RESTART_STAGE="${OPTARG}"
            ;;
        r)
            INPUT_DATA_MODE="${OPTARG}"
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "4" ]]; then
    # Too few arguments
    usage
fi

KEYPAIR_NAME="${1}"
shift
GRAPHS_URL="${1}"
shift
ALIGNMENTS_URL="${1}"
shift
CALLS_URL="${1}"
shift

if [[ "$#" -gt "0" ]]; then
    # Options or extra arguments after positional arguments, which getopts can't handle
    echo 1>&2 "Error: trailing option or unrecognized positional argument: ${1}"
    exit 1
fi

# Also work out the environment-variable-setting strings we need.
TOIL_ENV=(env TOIL_APPLIANCE_SELF="${TOIL_APPLIANCE_SELF}" TOIL_DOCKER_REGISTRY="${TOIL_DOCKER_REGISTRY}")

# Define some cluster management functions

# Return success if the cluster with the given name exists, and failure otherwise
function cluster_exists() {
    local CLUSTER_NAME="${1}"
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" true
}

# Start up the cluster with the given name, using the given keypair, and install everything necessary on it
function start_cluster() {
    local CLUSTER_NAME="${1}"
    shift
    local KEYPAIR_NAME="${1}"
    
    echo "Creating cluster ${CLUSTER_NAME}"
    
    # Start the cluster
    "${TOIL_ENV[@]}" $PREFIX toil launch-cluster "${CLUSTER_NAME}" --leaderNodeType=t2.medium -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"
    
    # We need to manually install git to make pip + git work...
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt update
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt install git -y

    # Ignore the old virtualenv if re-using a cluster

    # For hot deployment to work, toil-vg needs to be in a virtualenv that can see the system Toil
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" virtualenv --system-site-packages venv
}

# Tear down the cluster with the given name
function destroy_cluster() {
    local CLUSTER_NAME="${1}"
    $PREFIX toil destroy-cluster "${CLUSTER_NAME}" -z us-west-2a
}

# Fill in some variables characterizing our input based on the input data mode
case "${INPUT_DATA_MODE}" in
    21)
        # Do a lot of reads
        READ_COUNT="10000000"
        # Evaluate on some of them
        READ_DOWNSAMPLE_PORTION="0.5"
        # Simulate in several chunks
        READ_CHUNKS="32"
        # Mark reads with the feature names from these BEDs, if set
        READ_TAG_BEDS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags.bed")
        # Drop reads whose truth contigs match any of these regexes
        READ_DECOY_REGEXES=()
        # Define a region name to process. This sets the name that the graphs and
        # indexes will be saved/looked for under.
        REGION_NAME="CHR21"
        # Define the contigs we are using
        GRAPH_CONTIGS=("21")
        # And the offsets for calling on each of the contigs (0-based)
        GRAPH_CONTIG_OFFSETS=("0")
        # Define the region to build the graph on, as contig[:start-end], 1-based
        # If a non-whole-contig region is specified for any contig, we will need to convert this to
        # BED, so a region must be specified for all contigs.
        # If this is empty, contig names are sniffed from the input FASTA
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}")
        # Define the VCF and FASTA basenames. We assume each VCF has a TBI.
        # If multiple files are used they are paired up with the contigs in order until the VCFs run out.
        # We now use GRCh38 for compatibility with Platinum Genomes
        CONSTRUCT_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/1kg/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr21_GRCh38.genotypes.20170504.vcf.gz")
        # TODO: This FASTA qhould have the centromeres, but they are hard-masked.
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/GRCh38.only21.fa.gz")
        # What FASTA should we use for BWA mapping and Freebayes calling? It needs to have just the selected regions cut out.
        MAPPING_CALLING_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        # What VCF should we use for the truth? Must be a single VCF.
        # Used for evaluation and also for sample graph construction and read simulation.
        EVALUATION_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/platinum-genomes/2017-1.0/hg38/hybrid/nochr/hg38.hybrid.21.vcf.gz"
        # And a single FASTA
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        # And what high confidence regions should we use there?
        # This can't be specified if we are using regions on GRAPH_REGIONS
        EVALUATION_BED_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/platinum-genomes/2017-1.0/hg38/hybrid/nochr/hg38.hybrid.21.bed.gz"
        
        # We will process an interleaved fastq or a BAM of real reads and evaluate its variant calls too.
        # This can be empty, but these are mutually exclusive
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.21.bam"
        ;;
    MHC)
        # Actually do a smaller test
        READ_COUNT="100000"
        READ_DOWNSAMPLE_PORTION="1.0"
        READ_CHUNKS="2"
        READ_TAG_BEDS=()
        READ_DECOY_REGEXES=()
        REGION_NAME="MHC"
        GRAPH_CONTIGS=("6")
        GRAPH_CONTIG_OFFSETS=("28510119")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:28510120-33480577")
        CONSTRUCT_VCF_URLS=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr6_GRCh38_sites.20170504.vcf.gz")
        # We had "s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-MHC.vcf.gz" but it is malformed
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr6.fa.gz")
        MAPPING_CALLING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/MHC.fa"
        EVALUATION_VCF_URL="${CONSTRUCT_VCF_URLS[0]}"
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/platinum_NA12878_MHC.fq.gz"
        REAL_REALIGN_BAM_URL=""
        ;;
    BRCA1)
        # Do just BRCA1 and a very few reads
        READ_COUNT="20000"
        READ_DOWNSAMPLE_PORTION="1.0"
        READ_CHUNKS="2"
        READ_TAG_BEDS=()
        READ_DECOY_REGEXES=()
        REGION_NAME="BRCA1"
        GRAPH_CONTIGS=("17")
        GRAPH_CONTIG_OFFSETS=("43044293")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:43044294-43125484")
        CONSTRUCT_VCF_URLS=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr17_GRCh38.genotypes.20170504.vcf.gz")
        # We had "s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-BRCA1.vcf.gz" but it is malformed
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr17.fa.gz")
        MAPPING_CALLING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/BRCA1.fa"
        EVALUATION_VCF_URL="${CONSTRUCT_VCF_URLS[0]}"
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/platinum_NA12878_BRCA1.fq.gz"
        REAL_REALIGN_BAM_URL=""
        ;;
    WG38)
        # 10m pairs is enough to test mapping.
        # We don't really want to do sim-based calling anymore.
        READ_COUNT="10000000"
        # Use all of it
        READ_DOWNSAMPLE_PORTION="1.0"
        READ_CHUNKS="32"
        # This BED comes from Karen Miga's centromere and telomere BEDs (see ../real-data in this repo)
        # and Trevor Pesout's LINE and SINE BEDs (from s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/LINE.bed and SINE.bed)
        READ_TAG_BEDS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags.bed")
        # This BED was made by thresholding http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw
        # at 1.0 using bwtool, and applies a "umapAll100" tag at positions where all overlapping 100-mers map.
        # So 100-bp reads touching one of these annotations really should have a mappable 100-mer in them, before variants/errors.
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_umapAll100.bed")
        # This BED was made directly from http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Unique.Mappability.bb 
        # and applies a "umapAny100" tag at positions where at least one overlapping 100-mer maps.
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_umapAny100.bed")
        # These two BEDs were generated by downloading Gencode v24 whole genes and exons from the UCSC table browser,
        # removing the "chr", and then stripping down to contig, start, end and adding tag names.
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_gencodeGene.bed")
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_gencodeExon.bed")
        # These two BEDs were made from https://raw.githubusercontent.com/ENCODE-DCC/encValData/master/GRCh38/GRCh38.chrom.sizes
        # by removing "chr" and turning each record of a "_random" contig into a full-length unlocalized feature, and each
        # "chrUn" contig into a full-length unplaced feature
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_unlocalizedContig.bed")
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_unplacedContig.bed")
        # This pseudoautosomal region BED was manually produced from the PAR table at https://www.ncbi.nlm.nih.gov/grc/human
        # by subtracting 1 from the start posittion and leaving the end position unchanged.
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_PAR.bed")
        # This blacklisted region BED was derived from ENCODE's blacklist BED, converted to our naming scheme.
        # Users should apparently cite:
        # ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome.
        # Nature. 2012 Sep 6;489(7414):57-74. doi: 10.1038/nature11247.
        # Original BED: http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
        # See: https://sites.google.com/site/anshulkundaje/projects/blacklists
        READ_TAG_BEDS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags_blacklist.bed")
        # Drop reads on Y (because NA12878 is XX AFAIK) and also on all the _decoy sequences (because they aren't "real" parts of the reference)
        # TODO: Have a better, graph-internal way to delineate decoys.
        READ_DECOY_REGEXES=("^Y$" "_decoy$")
        REGION_NAME="WG38"
        # We do calling on all the autosomes and also X because NA12878 is XX
        GRAPH_CONTIGS=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
        GRAPH_CONTIG_OFFSETS=("0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0")
        # We don't pass any contig list to graph construction, so we use all the contigs in the FASTA
        GRAPH_REGIONS=()
        CONSTRUCT_VCF_URLS=()
        for CONTIG in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
            # Use the per-chromosome VCFs from 1KG, which we have mirrored in S3
            # Make sure to do all the chromosomes, even Y, for our reference, and in the order they appear in the FASTA.
            CONSTRUCT_VCF_URLS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/1kg/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${CONTIG}_GRCh38.genotypes.20170504.vcf.gz")
        done
        # I built this by dropping the alts from
        # http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz,
        # ordering 1-22, X, Y, M, other, mixing in hs38DH-extra.fa with the
        # decoys that BWA uses (but dropping the HLAs), and dropping the "chr"
        # from all the contig names. It has the centromeres in, but contains
        # lost of soft-masked (lower-case) sequence. I haven't yet checked
        # whether it has the right ref bases for all the 1kg variants.
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/hg38.centromeres.decoys.noAlts.fa.gz")
        MAPPING_CALLING_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        # Evaluate against Platinum Genomes/GIAB hybrid
        EVALUATION_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/platinum-genomes/2017-1.0/hg38/hybrid/nochr/hg38.hybrid.vcf.gz"
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        EVALUATION_BED_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/platinum-genomes/2017-1.0/hg38/hybrid/nochr/hg38.hybrid.bed.gz"
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam"
        ;;
    test)
        # Do some test regions of fake chromosomes
        # Exercise both multi-chromosome and offset capabilities, as well as no-VCF contig support
        READ_COUNT="1000"
        READ_DOWNSAMPLE_PORTION="0.90"
        READ_CHUNKS="2"
        READ_TAG_BEDS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref-tags.bed")
        READ_DECOY_REGEXES=()
        REGION_NAME="test"
        GRAPH_CONTIGS=("ref" "x" "nonvariable")
        GRAPH_CONTIG_OFFSETS=("5" "5" "0")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:6-1133" "${GRAPH_CONTIGS[1]}:6-1001" "${GRAPH_CONTIGS[2]}")
        CONSTRUCT_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref.vcf.gz" "s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/x.vcf.gz")
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.fa")
        MAPPING_CALLING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined-minus-5-bases.fa"
        EVALUATION_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.vcf.gz"
        EVALUATION_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.fa"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL=""
        # Override global sample name with one present in these VCFs, for testing
        SAMPLE_NAME="1"
        ;;
    testFastaRegions)
        # Do some test chromosomes
        # Exercise --fasta_regions contig name inferrence
        READ_COUNT="1000"
        READ_DOWNSAMPLE_PORTION="0.90"
        READ_CHUNKS="2"
        READ_TAG_BEDS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref-tags.bed")
        READ_DECOY_REGEXES=("^x$")
        REGION_NAME="testFastaRegions"
        GRAPH_CONTIGS=()
        GRAPH_CONTIG_OFFSETS=()
        # If empty we pass --fasta_regions to construct
        GRAPH_REGIONS=()
        CONSTRUCT_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref.vcf.gz" "s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/x.vcf.gz")
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.fa")
        MAPPING_CALLING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.fa"
        EVALUATION_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.vcf.gz"
        EVALUATION_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.fa"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL=""
        # Override global sample name with one present in these VCFs, for testing
        SAMPLE_NAME="1"
        ;;
    *)
        echo 1>&2 "Unknown input data set ${INPUT_DATA_MODE}"
        exit 1
esac

# Where do we store our job tree?
# We only need one job tree because only one Toil run runs at a time.
JOB_TREE="aws:us-west-2:${RUN_ID}"

echo "Running run ${RUN_ID} as ${KEYPAIR_NAME} on ${GRAPH_REGIONS[@]} for ${READ_COUNT} read pairs into ${ALIGNMENTS_URL}"

# Decide on sim and real alignment URLS
SIM_ALIGNMENTS_URL="${ALIGNMENTS_URL}/sim"
REAL_ALIGNMENTS_URL="${ALIGNMENTS_URL}/real"

# And for calls
SIM_CALLS_URL="${CALLS_URL}/sim"
REAL_CALLS_URL="${CALLS_URL}/real"

# Make sure we don't leave the cluster running or data laying around on exit.
function clean_up() {
    set +e
    
    if [[ "${PERSISTENT_CLUSTER}" == "0" ]]; then
        # Destroy the cluster
        destroy_cluster "${CLUSTER_NAME}"
    else
        echo "Leaving cluster ${CLUSTER_NAME} running!"
        echo "Destroy with: toil destroy-cluster '${CLUSTER_NAME}' -z us-west-2a"
    fi
    
    if [[ "${CLEAN_UP_JOB_TREE}" == "1" ]] ; then
        echo "Cleaning up job tree"
        toil clean "${JOB_TREE}"
    else 
        echo "Clean with:"
        echo toil clean "${JOB_TREE}"
    fi
}
trap clean_up EXIT

echo "Using cluster ${CLUSTER_NAME}"
echo "Destroy with: toil destroy-cluster '${CLUSTER_NAME}' -z us-west-2a"

if [ "${PERSISTENT_CLUSTER}" == "0" ] || ! cluster_exists "${CLUSTER_NAME}" ; then
    # We need to start up a cluster, because we are using our own or we are using a nonexistent one.
    start_cluster "${CLUSTER_NAME}" "${KEYPAIR_NAME}"
fi

# Install the right versions of all the Python stuff we need at once and hope
# pip is smart enough to figure out mutually compatible dependencies
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install --upgrade \
    pyyaml \
    "${AWSCLI_PACKAGE}" \
    numpy \
    scipy \
    scikit-learn \
    "${TOIL_VG_PACKAGE}"

# We need the master's IP to make Mesos go
MASTER_IP="$($PREFIX toil ssh-cluster --insecure --zone=us-west-2a --logOff "${CLUSTER_NAME}" hostname -i)"

# Strip some garbage from MASTER_IP
MASTER_IP="${MASTER_IP//[$'\t\r\n ']}"

# Work out the Toil cluster options that we pass to each Toil/toil-vg run to
# point it at the cluster. Having --defaultPreemptable makes jobs accept
# preemptable nodes by default. We also need --nodeStorage to specify disk for EBS nodes.
TOIL_CLUSTER_OPTS=(--realTimeLogging --logInfo \
    --batchSystem mesos --provisioner=aws "--mesosMaster=${MASTER_IP}:5050" \
    "--nodeTypes=${NODE_TYPES}" "--maxNodes=${MAX_NODES}" "--minNodes=${MIN_NODES}" \
    --nodeStorage 200 --defaultPreemptable \
    --metrics --retryCount 3 --stats)

########################################################################################################

# Now we are set up to do the actual experiment.

# We express our experimental conditions through a set of functions, keyed on condition name.
# We have one space of condition names for graphs being constructed, and later we will have another set for graph and mapper.

# We use a "haplotypes" condition to make the haplotype xg graphs for simualting from, although we can't really map/call against it.
# We also use a "bwa" condition to represent the BWA indexes, although it's not really a graph.

# Define all the graph conditions we can use
POSSIBLE_GRAPH_CONDITIONS=("snp1kg" "primary" "neg-control" "pos-control" "haplotypes" "bwa")
MIN_AF_NUM=0
for MIN_AF in "${MIN_AFS[@]}" ; do
    # Define all the minaf conditions
    if [[ "${MIN_AF_NUM}" == "0" ]] ; then
        POSSIBLE_GRAPH_CONDITIONS+=("snp1kg-minaf")
    else
        POSSIBLE_GRAPH_CONDITIONS+=("snp1kg-minaf${MIN_AF_NUM}")      
    fi
    MIN_AF_NUM=$((MIN_AF_NUM+1))
done

# Define the conditions with functions

# Produce the input VCF set used to generate the condition (and consequently which toil-vg construct run it needs to be part of).
# The result will be "construct" (for the VCFs used to construct ordinary graphs), or "evaluation" (for the positive control)
function get_graph_condition_step() {
    local CONDITION="${1}"
    if [[ "${CONDITION}" == "pos-control" || "${CONDITION}" == "haplotypes" ]] ; then
        echo "evaluation"
    else
        echo "construct"
    fi
}

# Produce the base URL for aligning against the given graph condition (i.e. .xg without the extension)
# Also works for BWA, in which case it produces the .fa URL without the index extension(s)
function get_graph_condition_base_url() {
    local CONDITION="${1}"
    
    # Work out the base URL for this condition
    case "${CONDITION}" in
        snp1kg)
            local GRAPH_BASE_URL="${GRAPHS_URL}/snp1kg-${REGION_NAME}_filter"
            ;;
        snp1kg-minaf*)
            # Parse out the MIN AF number
            local MIN_AF_NUM="${CONDITION#"snp1kg-minaf"}"
            if [[ -z "${MIN_AF_NUM}" ]] ; then
                MIN_AF_NUM=0
            fi
            # Get the actual threshold for that
            local MIN_AF_THRESHOLD="${MIN_AFS["${MIN_AF_NUM}"]}"
            local GRAPH_BASE_URL=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minaf_${MIN_AF_THRESHOLD}")
            ;;
        primary)
            local GRAPH_BASE_URL="${GRAPHS_URL}/primary"
            ;;
        neg-control)
            local GRAPH_BASE_URL="${GRAPHS_URL}/snp1kg-${REGION_NAME}_minus_${SAMPLE_NAME}"
            ;;
        pos-control)
            # This one gets generated from the evaluation VCF
            local GRAPH_BASE_URL="${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample"
            ;;
        haplotypes)
            # This one is generated from the evaluation VCF and also kind of fake
            local GRAPH_BASE_URL="${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo"
            ;;
        bwa)
            # This one isn't a real graph, it's wherever we should be putting our FASTA indexes
            local GRAPH_BASE_URL="${GRAPHS_URL}/bwa.fa"
            ;;
        *)
            echo 1>&2 "Error: unimplemented condition ${CONDITION}"
            # TODO: we can't necessarily kill the whole shell here, just make the function call fail
            exit 1
    esac
    
    echo "${GRAPH_BASE_URL}"
    
}

# Produce the base URL with optional comma-separated override evaluation index URL, used for toil-vg mapeval
# Doesn't work for BWA
function get_graph_index_base_with_override() {
    local CONDITION="${1}"
    
    if [[ "${CONDITION}" == "bwa" ]] ; then
        echo "Error: The bwa condition can't be fed into mapeval with a graph index base option." 1>&2
        exit 1
    fi
    
    # Look up the real base URL
    local GRAPH_BASE_URL="$(get_graph_condition_base_url "${CONDITION}")"
    
    if [[ "${CONDITION}" == "pos-control" ]] ; then
        # We have an override
        echo "${GRAPH_BASE_URL},${GRAPH_BASE_URL}_withref"
    else
        # No override, pass it along
        echo "${GRAPH_BASE_URL}"
    fi
}

# Return success if the output files for the given condition are ready in S3, and false if not and the condition has to be rerun.
# Works for the bwa condition too, in which case it checks for all the BWA indexes.
function is_graph_condition_done() {
    local CONDITION="${1}"
    
    # Work out the base URL for this condition
    local GRAPH_BASE_URL="$(get_graph_condition_base_url "${CONDITION}")"
    
    local GRAPHS_READY=1
    
    if [[ "${CONDITION}" == "haplotypes" ]] ; then
        # For the haplotypes we need to check the main xg and also the individual haplotype xgs
        for SIM_INDEX_URL in "${GRAPH_BASE_URL}_thread_0.xg" "${GRAPH_BASE_URL}_thread_1.xg" "${GRAPH_BASE_URL}.xg" ; do
            if ! aws s3 ls >/dev/null "${SIM_INDEX_URL}" ; then
                echo "Need to generate haplotype file ${SIM_INDEX_URL}"
                GRAPHS_READY=0
                break
            fi
        done
    else
        if [[ "${CONDITION}" == "bwa" ]] ; then
            # FASTA files get these indexes
            local SUFFIXES=(".amb" ".ann" ".bwt" ".pac" ".sa") 
        else
            # Normal graphs get these indexes
            local SUFFIXES=(".xg" ".gcsa" ".gcsa.lcp")
        fi
        
        
        if [[ "${CONDITION}" == "snp1kg"* || "${CONDITION}" == "neg-control" ]] ; then
            # These all have GBWT and snarl indexes
            SUFFIXES+=(".gbwt" ".snarls") 
        fi
        
        if [[ "${CONDITION}" == "primary" ]] ; then
            # Primary has an (empty) snarl index
            SUFFIXES+=(".snarls")
        fi  
        
        if [[ "${CONDITION}" == "pos-control" ]] ; then
            # The positive control has an associated withref xg that it uses for an override
            SUFFIXES+=("_withref.xg")
        fi
        
        for SUFFIX in "${SUFFIXES[@]}" ; do
            # For each file we want
            if ! aws s3 ls >/dev/null "${GRAPH_BASE_URL}${SUFFIX}" ; then
                # The graphs are not ready yet because this file is missing
                echo "Need to generate index file ${GRAPH_BASE_URL}${SUFFIX}"
                GRAPHS_READY=0
                break
            fi
        done
    fi
    
    if [[ "${GRAPHS_READY}" == "0" ]] ; then
        # Return false
        return 1
    fi
    # Return true
    return 0
}
    
# Add construct options for the given conditions to the global CONDITION_OPTS array.
# You must pass all the conditions at once, because some conditions are generated with the same arguments.
function add_graph_conditions_options() {
    # We only have to deal with the same options for multiple conditions for minaf
    # So we deduplicate the minaf conditions in the input down to just "snp1kg-minaf"
    local DEDUPLICATED_CONDITIONS=()
    local HAVE_MINAF=0
    # Note that $@ mass-expands automatically
    for CONDITION in "$@" ; do
        if [[ "${CONDITION}" == "snp1kg-minaf"* ]] ; then
            if [[ "${HAVE_MINAF}" == "0" ]] ; then
                DEDUPLICATED_CONDITIONS+=("snp1kg-minaf")
                local HAVE_MINAF=1
            fi
        else
            DEDUPLICATED_CONDITIONS+=("${CONDITION}")
        fi
    done

    for CONDITION in "${DEDUPLICATED_CONDITIONS[@]}" ; do
        # Now go through all the deduplicated condition names and produce their options
        case "${CONDITION}" in
            snp1kg)
                CONDITION_OPTS+=("--pangenome")
                ;;
            snp1kg-minaf)
                CONDITION_OPTS+=("--min_af" "${MIN_AFS[@]}")
                ;;
            primary)
                CONDITION_OPTS+=("--primary")
                ;;
            neg-control)
                CONDITION_OPTS+=("--neg_control" "${SAMPLE_NAME}")
                ;;
            pos-control)
                # This one gets generated from the evaluation VCF
                CONDITION_OPTS+=("--sample_graph" "${SAMPLE_NAME}")
                ;;
            haplotypes)
                # This one gets generated from the evaluation VCF and is also kind of fake
                CONDITION_OPTS+=("--haplo_sample" "${SAMPLE_NAME}")
                ;;
            bwa)
                # This one gets generated from a separate FASTA file
                CONDITION_OPTS+=("--bwa_reference" "${MAPPING_CALLING_FASTA_URL}")
                ;;
            *)
                echo 1>&2 "Error: unimplemented condition ${CONDITION}"
                # TODO: we can't necessarily kill the whole shell here, just make the function call fail
                exit 1
        esac
    done
}



# Define, of those, which we will run. This could be all of them, or just one or a few.
# The haplotypes condition always needs to be run if we want any simulated reads.
# And BWA always needs to be run to compare agaisnt bwa
RUN_GRAPH_CONDITIONS=("primary" "haplotypes" "snp1kg" "snp1kg-minaf" "bwa")

# Pass along either all the regions or --fasta_regions to infer them, for construction   
CONSTRUCT_REGION_OPTS=()
if [[ "${#GRAPH_REGIONS[@]}" -eq "0" ]] ; then
    # No regions known.
    CONSTRUCT_REGION_OPTS+=("--fasta_regions")
else
    # Use specified regions
    CONSTRUCT_REGION_OPTS+=("--regions" "${GRAPH_REGIONS[@]}")
fi

for CONSTRUCT_STEP in "construct" "evaluation" ; do
    # For each of the runs we want to do
    
    # Work out what graphs we need
    STEP_GRAPH_CONDITIONS=()
    
    for CONDITION in "${RUN_GRAPH_CONDITIONS[@]}" ; do
        # For each condition that might be involved
        if [[ "$(get_graph_condition_step "${CONDITION}")" == "${CONSTRUCT_STEP}" ]] ; then
            # Do this condition in this step
            STEP_GRAPH_CONDITIONS+=("${CONDITION}")
        fi
    done
    
    # Work out if the graphs are ready
    GRAPHS_READY=1
    # And if not which need to be done
    UNREADY_GRAPH_CONDITIONS=()
    # And whether we will need the VCFs to do them
    NEED_VCFS=0
    for CONDITION in "${STEP_GRAPH_CONDITIONS[@]}" ; do
        # Check if each condition for the step is done
        if ! is_graph_condition_done "${CONDITION}" ; then
            # If any are not, we need to do the step and also actually do those conditions
            GRAPHS_READY=0
            UNREADY_GRAPH_CONDITIONS+=("${CONDITION}")
            
            if [[ "${CONDITION}" != "primary" && "${CONDITION}" != "bwa" ]] ; then
                # If we have to run any condition other than primary or bwa, we need the VCFs
                NEED_VCFS=1
            fi
        fi
    done
    
    if [[ "${GRAPHS_READY}" == "1" ]] ; then
        # Graphs are all ready. Try the next construct step, if any.
        continue
    fi
    
    # Collect together the options we need to make the unmade graphs
    CONDITION_OPTS=()
    add_graph_conditions_options "${UNREADY_GRAPH_CONDITIONS[@]}"
    
    # Determine if we are restarting
    RESTART_OPTS=()
    if [[ ( "${RESTART_STAGE}" == "construct" && "${CONSTRUCT_STEP}" == "construct" ) || \
          ( "${RESTART_STAGE}" == "construct-sample" && "${CONSTRUCT_STEP}" == "evaluation" ) ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi
    
    # Which options vary from one graph construction step to the other
    STEP_OPTS=()
    if [[ "${CONSTRUCT_STEP}" == "construct" ]] ; then
        if [[ "${NEED_VCFS}" == "1" ]] ; then
            # We need input VCFs
            STEP_OPTS+=("--vcf" "${CONSTRUCT_VCF_URLS[@]}")
            # Pass along the filtering options to filter down the input VCFs
            STEP_OPTS+=("${FILTER_OPTS[@]}")
        fi
        STEP_OPTS+=("--fasta" "${CONSTRUCT_FASTA_URLS[@]}")
        STEP_OPTS+=("--out_name" "snp1kg-${REGION_NAME}")
    elif [[ "${CONSTRUCT_STEP}" == "evaluation" ]] ; then
        if [[ "${NEED_VCFS}" == "1" ]] ; then
            # We need input VCFs
            STEP_OPTS+=("--vcf" "${EVALUATION_VCF_URL}")
        fi
        STEP_OPTS+=("--fasta" "${EVALUATION_FASTA_URL}")
        STEP_OPTS+=("--out_name" "platinum-${REGION_NAME}")
    fi
    
    # Construct the graphs
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg construct \
        "${JOB_TREE}" \
        "$(url_to_store "${GRAPHS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --alt_paths \
        "${CONDITION_OPTS[@]}" \
        --handle_unphased arbitrary \
        "${STEP_OPTS[@]}" \
        "${CONSTRUCT_REGION_OPTS[@]}" \
        --merge_graphs \
        --gcsa_index \
        --xg_index \
        --gbwt_index \
        --snarls_index \
        "${RESTART_OPTS[@]}" \
        "${TOIL_CLUSTER_OPTS[@]}"
        
    # Report stats to standard out
    echo "---BEGIN RUN STATS---"
    toil stats "${JOB_TREE}"
    echo "---END RUN STATS---"
    # Clean up on success
    toil clean "${JOB_TREE}"
    
done

# Now work out where in there these simulated reads belong
READS_URL="${GRAPHS_URL}/sim-${SAMPLE_NAME}-${READ_SEED}-${READ_COUNT}-${READ_CHUNKS}"

# Work out if all the reads files are done
READS_READY=1
for READS_FILE_URL in "${READS_URL}/true.pos" "${READS_URL}/sim.fq.gz" ; do
    if ! aws s3 ls >/dev/null "${READS_FILE_URL}" ; then 
        READS_READY=0
        break
    fi
done

if [[ "${READS_READY}" != "1" ]] ; then
    # Now we need to simulate reads from the two haplotypes
    
    RESTART_OPTS=()
    if [[ "${RESTART_STAGE}" == "sim" ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi
    
    BED_TAG_OPTS=()
    for READ_TAG_BED in "${READ_TAG_BEDS[@]}"; do
        BED_TAG_OPTS+=("--tag_bed" "${READ_TAG_BED}")
    done
    
    DECOY_REGEX_OPTS=()
    for DECOY_REGEX in "${READ_DECOY_REGEXES[@]}"; do
        DECOY_REGEX_OPTS+=("--drop_contigs_matching" "${DECOY_REGEX}")
    done
    
    # This will make a "sim.gam".
    # We provide custom sim options with no substitutions over those specified by the FASTQ.
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg sim \
        "${JOB_TREE}" \
        "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_0.xg" \
        "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_1.xg" \
        "${READ_COUNT}" \
        "$(url_to_store "${READS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        "${BED_TAG_OPTS[@]}" \
        "${DECOY_REGEX_OPTS[@]}" \
        --annotate_xg "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
        --gam \
        --fastq_out \
        --seed "${READ_SEED}" \
        --sim_chunks "${READ_CHUNKS}" \
        --fastq "${TRAINING_FASTQ}" \
        --sim_opts "-p 570 -v 165 -i 0.002" \
        "${RESTART_OPTS[@]}" \
        "${TOIL_CLUSTER_OPTS[@]}"
        
    # Report stats to standard out
    echo "---BEGIN RUN STATS---"
    toil stats "${JOB_TREE}"
    echo "---END RUN STATS---"
    # Clean up on success
    toil clean "${JOB_TREE}"
    
fi

# Decide what graphs to do mapeval on
# Note that the sample graph positive control is added separately.

# These are the names of the graph conditions, which will be used as GAM names in mapeval
GAM_NAMES=()
# These are the graph base URLs with optional comma-separated evaluation override URLs
INDEX_BASES=()

for CONDITION in "${RUN_GRAPH_CONDITIONS[@]}" ; do
    if [[ "${CONDITION}" == "haplotypes" ]] ; then
        # Skip this condition for mapping; it is only for haplotype graph generation
        continue
    fi
    
    if [[ "${CONDITION}" == "bwa" ]] ; then
        # Skip this condition because it gets its own option instead.
        continue
    fi

    # For each condition we will run, spit out its GAM name and index_bases entry
    GAM_NAMES+=("${CONDITION}")
    INDEX_BASES+=("$(get_graph_index_base_with_override "${CONDITION}")")
done

################################################################################

# Now do the mapping. We have a system of map conditions, just like the graph conditions.
# Each map condition is a graph condition, then with -gbwt, -gbwt<number>, -mp, and/or -pe tags.
# Map conditions are also defined using functions.


# Produce the graph condition name from a map condition name.
# The graph condition can be used to get the base URL, from which the xg URL can be derived.
function get_map_condition_graph_condition() {
    local CONDITION="${1}"
    # Drop all the suffixes, in order
    local CONDITION="${CONDITION%-pe}"
    local CONDITION="${CONDITION%-mp}"
    # Greedily lop off the biggest gbwt specifier we can get with extglob syntax
    # See <http://wiki.bash-hackers.org/syntax/pattern>
    local CONDITION="${CONDITION%%-gbwt*([0-9.])}"
    echo "${CONDITION}"
}

# Produce the GBWT recombination penalty from a map condition name, or "" if not present
function get_map_condition_gbwt_penalty() {
    local CONDITION="${1}"
    # Grab the number with grep, printing just the match and catching -gbwt with a lookbehind
    # See https://serverfault.com/a/660878/49409
    # Make sure to allow grep to fail
    echo "${CONDITION}" | grep -oP -- '(?<=-gbwt)([0-9.]*)' || true
}

# Work out what map conditions to run
# Doesn't include BWA which is always run.
RUN_MAP_CONDITIONS=("primary-mp-pe" "snp1kg-minaf-mp-pe" "snp1kg-mp-pe")

for PENALTY in "${GBWT_RECOMBINATION_PENALTIES[@]}" ; do
    # Add all the GBWT conditions with the specified penalties
    RUN_MAP_CONDITIONS+=("snp1kg-gbwt${PENALTY}-mp-pe")
done
if [[ "${#GBWT_RECOMBINATION_PENALTIES[@]}" -eq "0" ]] ; then
    # If there are no penalties specified, just add the default-penalty gbwt in
    RUN_MAP_CONDITIONS+=("snp1kg-gbwt-mp-pe")
fi


# Work out the XG URLs that are used for each map condition
XG_URLS=()
for CONDITION in "${RUN_MAP_CONDITIONS[@]}" ; do
    GRAPH_CONDITION="$(get_map_condition_graph_condition "${CONDITION}")"
    XG_URLS+=("$(get_graph_condition_base_url "${GRAPH_CONDITION}").xg")
done

# This will hold the final GAM URLs for simulated reads
SIM_GAM_URLS=()
for CONDITION_NAME in "${RUN_MAP_CONDITIONS[@]}" ; do
    # We generate them from the condition names
    SIM_GAM_URLS+=("${SIM_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}_default.gam") 
done

# Also BWA
SIM_BAM_URLS=("${SIM_ALIGNMENTS_URL}/bwa-mem-pe.bam")
BAM_NAMES=("bwa-pe")
# Don't surject BAMs form simulated reads

# This holds the map conditions we will actually run because they aren't done yet
# TODO: Work out which map conditions are not done and do just them.
UNREADY_MAP_CONDITIONS=()
SIM_ALIGNMENTS_READY=1

# Check if all the expected output alignments exist and only run if not.
for ALIGNMENT_URL in "${SIM_GAM_URLS[@]}" "${SIM_BAM_URLS[@]}" "${SIM_ALIGNMENTS_URL}/position.results.tsv" "${SIM_ALIGNMENTS_URL}/plots/plot-roc.svg"; do
    if ! aws s3 ls >/dev/null "${ALIGNMENT_URL}" ; then
        # The alignments are not ready yet because this file is missing
        echo "Need to generate alignment file ${ALIGNMENT_URL}"
        UNREADY_MAP_CONDITIONS=("${RUN_MAP_CONDITIONS[@]}")
        SIM_ALIGNMENTS_READY=0
        break
    fi
done

if [[ "${SIM_ALIGNMENTS_READY}" != "1" ]] ; then
    
    RESTART_OPTS=()
    if [[ "${RESTART_STAGE}" == "map-sim" ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi
    
    GBWT_PENALTY_OPTS=()
    for CONDITION in "${UNREADY_MAP_CONDITIONS[@]}" ; do
        PENALTY="$(get_map_condition_gbwt_penalty "${CONDITION}")"
        if [[ ! -z "${PENALTY}" ]] ; then
            # This condition we are running uses a GBWT penalty
            GBWT_PENALTY_OPTS+=("${PENALTY}")
        fi
    done
    if [[ "${#GBWT_PENALTY_OPTS[@]}" -ne 0 ]] ; then
        # There are penalties to send, so prefix them
        GBWT_PENALTY_OPTS=("--gbwt-penalties" "${GBWT_PENALTY_OPTS[@]}")
    fi
    
    # Run one big mapeval run that considers all conditions we are interested in
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg mapeval \
        "${JOB_TREE}" \
        "$(url_to_store "${SIM_ALIGNMENTS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --index-bases "${INDEX_BASES[@]}" \
        --gam-names "${GAM_NAMES[@]}" \
        --downsample "${READ_DOWNSAMPLE_PORTION}" \
        --multipath \
        --use-gbwt \
        "${GBWT_PENALTY_OPTS[@]}" \
        --strip-gbwt \
        --use-snarls \
        --bwa --fasta "$(get_graph_condition_base_url bwa)" \
        --fastq "${READS_URL}/sim.fq.gz" \
        --truth "${READS_URL}/true.pos" \
        --plot-sets \
        "Primary vs. BWA:primary-mp-pe,bwa-mem-pe" \
        "GBWT vs. Not:snp1kg-mp-pe,snp1kg-gbwt5.0-mp-pe,snp1kg-minaf-mp-pe" \
        "GBWT Penalty:snp1kg-gbwt5.0-mp-pe,snp1kg-gbwt10.0-mp-pe,snp1kg-gbwt20.7-mp-pe,snp1kg-gbwt40.0-mp-pe" \
        "GBWT Penalty Fine:snp1kg-gbwt19.0-mp-pe,snp1kg-gbwt20.7-mp-pe,snp1kg-gbwt22.0-mp-pe" \
        "${RESTART_OPTS[@]}" \
        "${TOIL_CLUSTER_OPTS[@]}"
        
    # Report stats to standard out
    echo "---BEGIN RUN STATS---"
    toil stats "${JOB_TREE}"
    echo "---END RUN STATS---"
    # Clean up on success
    toil clean "${JOB_TREE}"
fi

# Exit after sim read alignments, so we can compare results so far against other runs.
exit

if [[ ! -z "${REAL_FASTQ_URL}" || ! -z "${REAL_REALIGN_BAM_URL}" ]] ; then
    # We can also do alignments of real data

    # This will hold the final GAM URLs for real reads
    REAL_GAM_URLS=()
    for CONDITION_NAME in "${RUN_MAP_CONDITIONS[@]}" ; do
        # We generate them from the condition names
        REAL_GAM_URLS+=("${REAL_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}_default.gam") 
    done

    # And the BAM URLs
    REAL_BAM_URLS=("${REAL_ALIGNMENTS_URL}/bwa-mem-pe.bam")
    for CONDITION_NAME in "${RUN_MAP_CONDITIONS[@]}" ; do
        # We do surject the real reads.
        # We generate BAM names from the condition names
        REAL_BAM_URLS+=("${REAL_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}-surject.bam")
    done

    # Make sure they exist
    # TODO: If they don't exist, work out specifically which conditions to run.
    # For now we just do all of them.
    # TODO: combine with sim read code
    UNREADY_MAP_CONDITIONS=()
    REAL_ALIGNMENTS_READY=1
    for ALIGNMENT_URL in "${REAL_GAM_URLS[@]}" "${REAL_BAM_URLS[@]}" "${REAL_ALIGNMENTS_URL}/position.results.tsv" "${REAL_ALIGNMENTS_URL}/plots/plot-roc.svg" ; do
        if ! aws s3 ls >/dev/null "${ALIGNMENT_URL}" ; then
            # The alignments are not ready yet because this file is missing
            echo "Need to generate alignment file ${ALIGNMENT_URL}"
            UNREADY_MAP_CONDITIONS=("${RUN_MAP_CONDITIONS[@]}")
            REAL_ALIGNMENTS_READY=0
            break
        fi
    done
    
    # Work out how to send the input reads
    if [[ ! -z "${REAL_FASTQ_URL}" ]] ; then
        # Data is coming in as FASTQ
        DATA_OPTS=(--fastq "${REAL_FASTQ_URL}")
    else
        # Data is coming in as BAM
        DATA_OPTS=(--bam_input_reads "${REAL_REALIGN_BAM_URL}")
    fi

    if [[ "${REAL_ALIGNMENTS_READY}" != "1" ]] ; then
    
        RESTART_OPTS=()
        if [[ "${RESTART_STAGE}" == "map-real" ]] ; then
            # Restart from this stage
            RESTART_OPTS=("--restart")
        fi
        
        GBWT_PENALTY_OPTS=()
        for CONDITION in "${UNREADY_MAP_CONDITIONS[@]}" ; do
            PENALTY="$(get_map_condition_gbwt_penalty "${CONDITION}")"
            if [[ ! -z "${PENALTY}" ]] ; then
                # This condition we are running uses a GBWT penalty
                GBWT_PENALTY_OPTS+=("${PENALTY}")
            fi
        done
        if [[ "${#GBWT_PENALTY_OPTS[@]}" -ne 0 ]] ; then
            # There are penalties to send, so prefix them
            GBWT_PENALTY_OPTS=("--gbwt-penalties" "${GBWT_PENALTY_OPTS[@]}")
        fi
    
        # Run a mapeval run just to map the real reads, under all conditions
        $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg mapeval \
            "${JOB_TREE}" \
            "$(url_to_store "${REAL_ALIGNMENTS_URL}")" \
            --whole_genome_config \
            "${VG_DOCKER_OPTS[@]}" \
            --index-bases "${INDEX_BASES[@]}" \
            --gam-names "${GAM_NAMES[@]}" \
            --downsample "${READ_DOWNSAMPLE_PORTION}" \
            --multipath \
            --use-gbwt \
            "${GBWT_PENALTY_OPTS[@]}" \
            --strip-gbwt \
            --use-snarls \
            --surject \
            --bwa --fasta "${MAPPING_CALLING_FASTA_URL}" \
            "${DATA_OPTS[@]}" \
            --skip-eval \
            "${RESTART_OPTS[@]}" \
            "${TOIL_CLUSTER_OPTS[@]}"
            
            
        # Report stats to standard out
        echo "---BEGIN RUN STATS---"
        toil stats "${JOB_TREE}"
        echo "---END RUN STATS---"
        # Clean up on success
        toil clean "${JOB_TREE}"
    fi
fi

################################################################################

# Now do the calling.

if [[ -z "${EVALUATION_VCF_URL}" ]] ; then
    # Don't calleval without the truth VCF
    echo "No truth VCF; skipping calleval."
    exit 0
fi

if [[ -z "${EVALUATION_BED_URL}" ]] ; then
    # No high confidence bed specified. So we can use a generated region bed to restrict to e.g. BRCA1.

    BED_LINES=()
    for GRAPH_REGION in "${GRAPH_REGIONS[@]}" ; do
        # For each region
        if echo "${GRAPH_REGION}" | grep ":" ; then
            # If it isn't a whole contig, we need a BED line
            # TODO: Validate that no other regions are whole contigs
            
            # Cut it up
            REGION_CONTIG="$(echo ${GRAPH_REGION} | cut -f1 -d':')"
            REGION_RANGE="$(echo ${GRAPH_REGION} | cut -f2 -d':')"
            REGION_START="$(echo ${REGION_RANGE} | cut -f1 -d'-')"
            REGION_END="$(echo ${REGION_RANGE} | cut -f2 -d'-')"
            
            # Convert to 0-based
            ((REGION_START--))
            ((REGION_END--))
            
            # Make a BED line
            BED_LINES+=("$(printf "${REGION_CONTIG}\t${REGION_START}\t${REGION_END}")")
        fi
    done

    # Collect the options we need to use to specify this bed
    BED_OPTS=()
    if [ ${#BED_LINES[@]} -ne 0 ] ; then
        # Upload the BED to the server in the silliest way possible
        TEMP_BED="./temp.bed"
        $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" truncate ${TEMP_BED} --size 0
        for BED_LINE in "${BED_LINES[@]}" ; do
            $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" bash -c "echo \"${BED_LINE}\" >${TEMP_BED}"
        done

        # Remember to use it
        BED_OPTS+=(--vcfeval_bed_regions "${TEMP_BED}")
    fi
else
    # We have a bed we need to use for the high-confidence regions
    BED_OPTS+=(--vcfeval_bed_regions "${EVALUATION_BED_URL}")
fi

# Don't do calls based on simulated data

# What plot sets do we use for calling?
# Each will get a normal (clipped) and an unclipped version.
CALL_PLOT_SETS=("primary-mp-pe-surject-fb,bwa-pe-fb")

if [ ! -z "${REAL_FASTQ_URL}" ] ; then
    # Now the real calls if applicable

    # Now the sim calls
    REAL_CALLS_READY=1
    if ! aws s3 ls >/dev/null "${REAL_CALLS_URL}/plots/roc-weighted.svg" ; then
        # Use just this one file as a marker of calleval done-ness.
        # TODO: use other files/VCFs.
        REAL_CALLS_READY=0
    fi

    # If we ran call we would add:
    #--gams "${REAL_GAM_URLS[@]}" \
    #--gam_names "${RUN_MAP_CONDITIONS[@]}" \
    # --xg_paths "${XG_URLS[@]}" \
    #--call \

    if [[ "${REAL_CALLS_READY}" != "1" ]] ; then
    
        RESTART_OPTS=()
        if [[ "${RESTART_STAGE}" == "call-real" ]] ; then
            # Restart from this stage
            RESTART_OPTS=("--restart")
        fi
    
        $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg calleval \
            "${JOB_TREE}" \
            "$(url_to_store "${REAL_CALLS_URL}")" \
            --whole_genome_config \
            "${VG_DOCKER_OPTS[@]}" \
            --bams "${REAL_BAM_URLS[@]}" \
            --bam_names "${BAM_NAMES[@]}" \
            --chroms "${GRAPH_CONTIGS[@]}" \
            --vcf_offsets "${GRAPH_CONTIG_OFFSETS[@]}" \
            --vcfeval_fasta "${EVALUATION_FASTA_URL}" \
            --vcfeval_baseline "${EVALUATION_VCF_URL}" \
            --caller_fasta "${MAPPING_CALLING_FASTA_URL}" \
            --freebayes \
            "${BED_OPTS[@]}" \
            --sample_name "${SAMPLE_NAME}" \
            --plot_sets "${CALL_PLOT_SETS[@]}" \
            "${RESTART_OPTS[@]}" \
            "${TOIL_CLUSTER_OPTS[@]}"
            
        # Report stats to standard out
        echo "---BEGIN RUN STATS---"
        toil stats "${JOB_TREE}"
        echo "---END RUN STATS---"
        # Clean up on success
        toil clean "${JOB_TREE}"
    fi

fi

echo "Mapping Test Script Successful"


# Cluster (if desired) and trees will get cleaned up by the exit trap

