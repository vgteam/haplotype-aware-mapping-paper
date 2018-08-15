#!/usr/bin/env bash
# test-population-mapping.sh: evaluate the performance impact of population-aware mapping using toil-vg mapeval on AWS

set -ex

# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/adamnovak/toil-vg.git@8c3724bb4226bba261809770ad6befadb4ca7baf#egg=toil-vg"

# What Docker registry can the corresponding dashboard containers (Grafana, etc.) be obtained from?
TOIL_DOCKER_REGISTRY="quay.io/adamnovak"
# What Toil appliance should we use? Ought to match the locally installed Toil,
# but can't quite if the locally installed Toil is locally modified or
# installed from a non-release Git commit.
TOIL_APPLIANCE_SELF="${TOIL_DOCKER_REGISTRY}/toil:3.17.0a1-79c241c0eb273a3af6952202e38bd8d0c4bfe05c"


# What version of awscli do we use? This has to be compatible with the
# botocore/boto3 that toil-vg and toil can agree on, and each awscli version
# seems to require exactly one botocore version, and pip can't do something
# sensible like find the one that goes with the botocore we need to have when
# we just ask for "awscli". So we work out the right version manually and live
# in hope that <https://github.com/pypa/pip/issues/988> will eventually stop
# plaguing our lab.
AWSCLI_PACKAGE="awscli==1.14.70"

# What vg should we use?
# Update this tag to change the Docker that will be used by a restart.
# Just editing the script won't do it; the tag name lives in the Toil job store.
# To update, do something like:
# docker pull quay.io/vgteam/vg:dev-v1.8.0-142-g758c92ec-t190-run
# docker tag quay.io/vgteam/vg:dev-v1.8.0-142-g758c92ec-t190-run quay.io/adamnovak/vg:wholegenome
# docker push quay.io/adamnovak/vg:wholegenome
VG_DOCKER_OPTS=("--vg_docker" "quay.io/vgteam/vg:dev-v1.9.0-123-gd4adc2fd-t217-run")

# What node types should we use?
# Comma-separated, with :bid-in-dollars after the name for spot nodes
# We need non-preemptable i3.4xlarge at least to get ~3.8TB storage available so the GCSA indexing jobs will have somewhere to run.
# And we also need more memory (?) than that so some of the later jobs will run.
# Suggested: i3.8xlarge which is worth ~0.70-0.80 (2.something on-demand), and r4.8xlarge which is worth ~0.60
# If using r4 we also need to make sure --nodeStorage, which sets EBS volume size for EBS nodes, is sufficient.
# We can also use r3.8xlarge which are cheaper (~0.50) and come with faster disks anyway.
# But you need a Toil which can tell them apart
NODE_TYPES="i3.8xlarge,r3.8xlarge:0.50"
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
    printf "\t-s STAGE\tRestart the given stage (construct, construct-sample, sim, map-sim, map-real, call-sim, call-real).\n"
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
        # Mark reads with the feature names from this BED, if set
        READ_TAG_BED="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags.bed"
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
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}")
        # Define the VCF and FASTA basenames. We assume each VCF has a TBI.
        # If multiple files are used they must match up with the contigs.
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
        READ_TAG_BED=""
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
        READ_TAG_BED=""
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
        # 10m pairs would be enough to test mapping, and we could use the real data to test calling
        # But it would be good to have sim calling results too
        # So simulate 500m read pairs, which even at 100 bp/read is ~30x coverage
        READ_COUNT="500000000"
        # Only look at 2% of that (10,000,000 pairs)
        READ_DOWNSAMPLE_PORTION="0.02"
        READ_CHUNKS="32"
        READ_TAG_BED="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/bed/GRCh38_tags.bed"
        REGION_NAME="WG38"
        # We do all the chroms except Y because NA12878 is XX AFAIK
        GRAPH_CONTIGS=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
        GRAPH_CONTIG_OFFSETS=("0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[@]}")
        CONSTRUCT_VCF_URLS=()
        for CONTIG in "${GRAPH_CONTIGS[@]}" ; do
            # Use the per-chromosome VCFs from 1KG, which we have mirrored in S3
            CONSTRUCT_VCF_URLS+=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/1kg/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${CONTIG}_GRCh38.genotypes.20170504.vcf.gz")
            # Use the FASTAs from 
        done
        # I built this by dropping the alts from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz, mixing in
        # hs38DH-extra.fa with the decoys that BWA uses, and dropping the "chr" from all the contig names.
        # It has the centromeres in, but contains lost of soft-masked (lower-case) sequence.
        # I haven't yet checked whether it has the right ref bases for all the 1kg variants.
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
        # Do a couple test regions of fake chromosomes
        # Exercise both multi-chromosome and offset capabilities
        READ_COUNT="1000"
        READ_DOWNSAMPLE_PORTION="0.90"
        READ_CHUNKS="2"
        READ_TAG_BED="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref-tags.bed"
        REGION_NAME="test"
        GRAPH_CONTIGS=("ref" "x")
        GRAPH_CONTIG_OFFSETS=("5" "5")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:6-1133" "${GRAPH_CONTIGS[1]}:6-1001")
        CONSTRUCT_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref.vcf.gz" "s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/x.vcf.gz")
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/ref.fa" "s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/x.fa")
        MAPPING_CALLING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined-minus-5-bases.fa"
        EVALUATION_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.vcf.gz"
        EVALUATION_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/test-data/combined.fa"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL=""
        # Override sample name with one present in these VCFs
        SAMPLE_NAME="1"
        ;;
    *)
        echo 1>&2 "Unknown input data set ${INPUT_DATA_MODE}"
        exit 1
esac

# Where do we store our job tree?
# We only need one job tree because only one Toil run runs at a time.
JOB_TREE="aws:us-west-2:${RUN_ID}"

echo "Running run ${RUN_ID} as ${KEYPAIR_NAME} on ${GRAPH_REGIONS[@]} for ${READ_COUNT} reads into ${ALIGNMENTS_URL}"

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

# Decide what graphs to run
# Note that the sample graph positive control is added separately.
GRAPH_URLS=()
GAM_NAMES=()

# We want the actual graph with its indexes (minus the sample under test)
GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_filter")
GAM_NAMES+=("snp1kg")

MIN_AF_NUM=0
for MIN_AF in "${MIN_AFS[@]}" ; do
    # Add to GRAPH_URLS and GAM_NAMES for all the minaf values
    if [[ "${MIN_AF_NUM}" == "0" ]] ; then
        GAM_NAMES+=("snp1kg-minaf")
    else
        GAM_NAMES+=("snp1kg-minaf${MIN_AF_NUM}")      
    fi
    GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}")
    MIN_AF_NUM=$((MIN_AF_NUM+1))
done

# We want a primary control
GRAPH_URLS+=("${GRAPHS_URL}/primary")
GAM_NAMES+=("primary")

# We want a negative control with no right variants
GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minus_${SAMPLE_NAME}")
GAM_NAMES+=("neg-control")

# Compose actual index bases by zipping
INDEX_BASES=()
GRAPH_NUMBER=0
while [[ "${GRAPH_NUMBER}" -lt "${#GRAPH_URLS[@]}" ]] ; do
    # For each graph, get its graph and override URLs
    GRAPH_URL="${GRAPH_URLS[${GRAPH_NUMBER}]}"
    INDEX_BASES+=("${GRAPH_URL}")
    
    ((GRAPH_NUMBER=${GRAPH_NUMBER}+1))
done

# Check if all the expected output graphs exist and only run if not.
GRAPHS_READY=1
for GRAPH_BASE_URL in "${GRAPH_URLS[@]}" ; do
    # For each graph we want to run
    for SUFFIX in ".xg" ".gcsa" ".gcsa.lcp" ; do
        # For each file absolutely required for the graph
        
        if ! aws s3 ls >/dev/null "${GRAPH_BASE_URL}${SUFFIX}" ; then
            # The graphs are not ready yet because this file is missing
            echo "Need to generate graph file ${GRAPH_BASE_URL}${SUFFIX}"
            GRAPHS_READY=0
            break
        fi
    done
    
    if [[ "${GRAPHS_READY}" == "0" ]] ; then
        break
    fi
done

if [[ "${GRAPHS_READY}" != "1" ]] ; then
    # Graphs need to be generated
    
    RESTART_OPTS=()
    if [[ "${RESTART_STAGE}" == "construct" ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi
    
    # Construct the graphs
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg construct \
        "${JOB_TREE}" \
        "$(url_to_store "${GRAPHS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --vcf "${CONSTRUCT_VCF_URLS[@]}" \
        --fasta "${CONSTRUCT_FASTA_URLS[@]}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --pangenome \
        --min_af "${MIN_AFS[@]}" \
        --primary \
        --pos_control "${SAMPLE_NAME}" \
        --neg_control "${SAMPLE_NAME}" \
        --handle_unphased arbitrary \
        "${FILTER_OPTS[@]}" \
        --regions "${GRAPH_REGIONS[@]}" \
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
        
fi



# Now we run the sample graph positive control, and make the single-haplotype XGs we will use for read simulation.
SAMPLE_GRAPHS_READY=1
# For each graph we want to run
for SUFFIX in ".xg" ".gcsa" ".gcsa.lcp" ; do
    # For each file absolutely required for the graph
    
    if ! aws s3 ls >/dev/null "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample${SUFFIX}" ; then
        # The sample graph is not ready yet because this file is missing
        echo "Need to generate sample graph file ${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample${SUFFIX}"
        SAMPLE_GRAPHS_READY=0
        break
    fi
done

if ! aws s3 ls >/dev/null "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample_withref.xg" ; then
    # The sample graph withref xg is missing
    echo "Need to generate sample graph file ${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample_withref.xg"
    SAMPLE_GRAPHS_READY=0
fi

for SIM_INDEX_URL in "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_0.xg" "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_1.xg" "${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" ; do
    # We also need these indexes for sim
    if ! aws s3 ls >/dev/null "${SIM_INDEX_URL}" ; then
        echo "Need to generate sample graph file ${SIM_INDEX_URL}"
        SAMPLE_GRAPHS_READY=0
        break
    fi
done

if [[ "${SAMPLE_GRAPHS_READY}" != "1" ]] ; then
    # Graphs need to be generated
    
    RESTART_OPTS=()
    if [[ "${RESTART_STAGE}" == "construct-sample" ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi
    
    # Construct the graphs
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg construct \
        "${JOB_TREE}" \
        "$(url_to_store "${GRAPHS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --vcf "${EVALUATION_VCF_URL}" \
        --fasta "${EVALUATION_FASTA_URL}" \
        --out_name "platinum-${REGION_NAME}" \
        --alt_paths \
        --haplo_sample "${SAMPLE_NAME}" \
        --sample_graph "${SAMPLE_NAME}" \
        --handle_unphased arbitrary \
        --regions "${GRAPH_REGIONS[@]}" \
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
        
fi

# Actually use the sample graph positive control, with evaluation XG override
INDEX_BASES+=("${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample,${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample_withref")
GAM_NAMES+=("pos-control")

# Now work out where in there these simulated reads belong
READS_URL="${GRAPHS_URL}/sim-${SAMPLE_NAME}-${READ_SEED}-${READ_COUNT}-${READ_CHUNKS}"

if ! aws s3 ls >/dev/null "${READS_URL}/true.pos" ; then 
    # Now we need to simulate reads from the two haplotypes
    
    RESTART_OPTS=()
    if [[ "${RESTART_STAGE}" == "sim" ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi
    
    BED_TAG_OPTS=()
    if [[ ! -z "${READ_TAG_BED}" ]] ; then
        # We want to tag the reads with overlap against features in this BED.
        BED_TAG_OPTS=("--tag_bed" "${READ_TAG_BED}")
    fi
    
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

# Work out what alignment condition GAM names we hope to generate, with tags
CONDITION_NAMES=()
# And what XGs go with them
XG_URLS=()

# TODO: This is sort of duplicative with GRAPH_URLS and GAM_NAMES above.
# But it is more specific/restrictive for just calling (i.e. we ignore single-ended).

CONDITION_NAMES+=("primary-mp-pe")
XG_URLS+=("${GRAPHS_URL}/primary.xg")

CONDITION_NAMES+=("snp1kg-pe")
XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_filter.xg")

CONDITION_NAMES+=("snp1kg-mp-pe")
XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_filter.xg")

CONDITION_NAMES+=("snp1kg-gbwt-mp-pe")
XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_filter.xg")

MIN_AF_NUM=0
for MIN_AF in "${MIN_AFS[@]}" ; do
    # Make condition names for all the minaf values
    # Make sure to handle normal and gbwt versions
    if [[ "${MIN_AF_NUM}" == "0" ]] ; then
        CONDITION_NAMES+=("snp1kg-minaf-mp-pe")
        CONDITION_NAMES+=("snp1kg-minaf-gbwt-mp-pe")
    else
        CONDITION_NAMES+=("snp1kg-minaf${MIN_AF_NUM}-mp-pe")
        CONDITION_NAMES+=("snp1kg-minaf${MIN_AF_NUM}-gbwt-mp-pe")
    fi
    XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}.xg")
    XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}.xg")
    MIN_AF_NUM=$((MIN_AF_NUM+1))
done

# For the positive control, use the sample + reference xg for calling.
CONDITION_NAMES+=("pos-control-mp-pe")
XG_URLS+=("${GRAPHS_URL}/platinum-${REGION_NAME}_${SAMPLE_NAME}_sample_withref.xg")

CONDITION_NAMES+=("neg-control-mp-pe")
XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minus_${SAMPLE_NAME}.xg")


# This will hold the final GAM URLs for simulated reads
SIM_GAM_URLS=()
for CONDITION_NAME in "${CONDITION_NAMES[@]}" ; do
    # We generate them from the condition names
    SIM_GAM_URLS+=("${SIM_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}_default.gam") 
done

# Also BWA
SIM_BAM_URLS=("${SIM_ALIGNMENTS_URL}/bwa-mem-pe.bam")
BAM_NAMES=("bwa-pe")
# And surjected BAMs
for CONDITION_NAME in "${CONDITION_NAMES[@]}" ; do
    # We generate them from the condition names
    SIM_BAM_URLS+=("${SIM_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}-surject.bam")
    BAM_NAMES+=("${CONDITION_NAME}-surject")
done

# Check if all the expected output alignments exist and only run if not.
SIM_ALIGNMENTS_READY=1
for ALIGNMENT_URL in "${SIM_GAM_URLS[@]}" "${SIM_BAM_URLS[@]}" "${SIM_ALIGNMENTS_URL}/position.results.tsv" "${SIM_ALIGNMENTS_URL}/plots/plot-roc.svg"; do
    if ! aws s3 ls >/dev/null "${ALIGNMENT_URL}" ; then
        # The alignments are not ready yet because this file is missing
        echo "Need to generate alignment file ${ALIGNMENT_URL}"
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
        --strip-gbwt \
        --use-snarls \
        --surject \
        --bwa --fasta "${MAPPING_CALLING_FASTA_URL}" \
        --fastq "${READS_URL}/sim.fq.gz" \
        --truth "${READS_URL}/true.pos" \
        --plot-sets \
        "primary-mp-pe,snp1kg-mp-pe,snp1kg-gbwt-mp-pe,snp1kg-minaf-mp-pe,snp1kg-minaf-gbwt-mp-pe,pos-control-mp-pe,neg-control-mp-pe" \
        "primary-mp,snp1kg-mp,snp1kg-gbwt-mp,snp1kg-minaf-mp,snp1kg-minaf-gbwt-mp,pos-control-mp,neg-control-mp" \
        "bwa-mem-pe,snp1kg-gbwt-mp-pe,snp1kg-pe" \
        "bwa-mem,snp1kg-gbwt-mp,snp1kg" \
        "snp1kg-minaf-mp-pe,snp1kg-minaf1-mp-pe,snp1kg-minaf2-mp-pe,snp1kg-minaf3-mp-pe" \
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
    for CONDITION_NAME in "${CONDITION_NAMES[@]}" ; do
        # We generate them from the condition names
        REAL_GAM_URLS+=("${REAL_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}_default.gam") 
    done

    # And the BAM URLs
    REAL_BAM_URLS=("${REAL_ALIGNMENTS_URL}/bwa-mem-pe.bam")
    for CONDITION_NAME in "${CONDITION_NAMES[@]}" ; do
        # We generate them from the condition names
        REAL_BAM_URLS+=("${REAL_ALIGNMENTS_URL}/aligned-${CONDITION_NAME}-surject.bam")
    done

    # Make sure they exist
    REAL_ALIGNMENTS_READY=1
    for ALIGNMENT_URL in "${REAL_GAM_URLS[@]}" "${REAL_BAM_URLS[@]}" "${REAL_ALIGNMENTS_URL}/position.results.tsv" "${REAL_ALIGNMENTS_URL}/plots/plot-roc.svg" ; do
        if ! aws s3 ls >/dev/null "${ALIGNMENT_URL}" ; then
            # The alignments are not ready yet because this file is missing
            echo "Need to generate alignment file ${ALIGNMENT_URL}"
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

# Now the sim calls
SIM_CALLS_READY=1
if ! aws s3 ls >/dev/null "${SIM_CALLS_URL}/plots/roc-weighted.svg" ; then
    # Use just this one file as a marker of calleval done-ness.
    # TODO: use other files/VCFs.
    SIM_CALLS_READY=0
fi

# It would be nice if we could run genotype, but it is extremely slow (~2.5 hours per chunk on chr21 sim data)
# If we ran call we would add
#--gams "${SIM_GAM_URLS[@]}" \
#--gam_names "${CONDITION_NAMES[@]}" \
# --xg_paths "${XG_URLS[@]}" \
#--call \

# What plot sets do we use for calling?
# Each will get a normal (clipped) and an unclipped version.
CALL_PLOT_SETS=("primary-mp-pe-surject-fb,snp1kg-mp-pe-surject-fb,snp1kg-gbwt-mp-pe-surject-fb,snp1kg-minaf-mp-pe-surject-fb,snp1kg-minaf-gbwt-mp-pe-surject-fb,pos-control-mp-pe-surject-fb,neg-control-mp-pe-surject-fb" \
    "bwa-pe-fb,snp1kg-gbwt-mp-pe-surject-fb,snp1kg-pe-surject-fb,primary-mp-pe-surject-fb" \
    "snp1kg-minaf-mp-pe-surject-fb,snp1kg-minaf1-mp-pe-surject-fb,snp1kg-minaf2-mp-pe-surject-fb,snp1kg-minaf3-mp-pe-surject-fb")

if [[ "${SIM_CALLS_READY}" != "1" ]] ; then

    RESTART_OPTS=()
    if [[ "${RESTART_STAGE}" == "call-sim" ]] ; then
        # Restart from this stage
        RESTART_OPTS=("--restart")
    fi

    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" "${TOIL_ENV[@]}" venv/bin/toil-vg calleval \
        "${JOB_TREE}" \
        "$(url_to_store "${SIM_CALLS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --bams "${SIM_BAM_URLS[@]}" \
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
    #--gam_names "${CONDITION_NAMES[@]}" \
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

