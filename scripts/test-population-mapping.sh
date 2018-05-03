#!/usr/bin/env bash
# test-population-mapping.sh: evaluate the performance impact of population-aware mapping using toil-vg mapeval on AWS

set -ex

# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/adamnovak/toil-vg.git@20d15d0a18ef46ec5ef49a13311362f2da4e56d9#egg=toil-vg"

# What Toil appliance should we use? Ought to match the locally installed Toil,
# but can't quite if the locally installed Toil is locally modified or
# installed from a non-release Git commit.
TOIL_APPLIANCE_SELF="quay.io/ucsc_cgl/toil:3.16.0a1.dev2290-c6d3a2a1677ba3928ad5a9ebb6d862b02dd97998"

# What version of awscli do we use? This has to be compatible with the
# botocore/boto3 that toil-vg and toil can agree on, and each awscli version
# seems to require exactly one botocore version, and pip can't do something
# sensible like find the one that goes with the botocore we need to have when
# we just ask for "awscli". So we work out the right version manually and live
# in hope that <https://github.com/pypa/pip/issues/988> will eventually stop
# plaguing our lab.
AWSCLI_PACKAGE="awscli==1.14.70"

# What vg should we use?
VG_DOCKER_OPTS=("--vg_docker" "quay.io/vgteam/vg:v1.5.0-3156-gfec77632-t159-run")

# What node types should we use?
# Comma-separated, with :bid-in-dollars after the name for spot nodes
# We need non-preemptable i3.4xlarge to get ~3.8TB storage available so the GCSA indexing jobs will have somewhere to run.
NODE_TYPES="r3.8xlarge:0.85,i3.4xlarge"
# How many nodes should we use at most per type?
# Also comma-separated.
# TODO: These don't sort right pending https://github.com/BD2KGenomics/toil/issues/2195
MAX_NODES="2,2"
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

# Should we delete the job store when we exit?
# We do by default, and if the run finishes successfully.
# We don't if we're asked to keep it and Toil errors out.
REMOVE_JOBSTORE=1

# What named graph region should we be operating on?
# We can do "21", "whole-genome", "MHC", etc.
INPUT_DATA_MODE="21"

# TODO: Allow using real input data for speed measurement.

# Set a FASTQ to model reads after
TRAINING_FASTQ="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz"
# And a read simulation seed
READ_SEED="90"

# Define the sample to use for synthesizing reads
SAMPLE_NAME="NA12878"
# And options to filter them out of graphs, along with people related to them
FILTER_OPTS=("--filter_ceph" "--filter_samples" "${SAMPLE_NAME}")

# What min allele frequency limit do we use?
MIN_AF="0.0335570469"

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
    printf "Usage: $0 [Options] KEYPAIR_NAME GRAPHS_URL OUTPUT_URL \n"
    printf "\tKEYPAIR_NAME\tSSH keypair in Amazon us-west-2 region to use for cluster authentication\n"
    printf "\tGRAPHS_URL\tS3 URL where generated graphs, indexes, and reads should be cached\n"
    printf "\tOUTPUT_URL\tS3 URL where output data (mapped reads and statistics) should be deposited\n"
    printf "Options:\n\n"
    printf "\t-d\tDo a dry run\n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg.\n"
    printf "\t-t CONTAINER\tUse the given Toil container in the cluster (default: ${TOIL_APPLIANCE_SELF}).\n"
    printf "\t-c CLUSTER\tUse the given persistent Toil cluster, which will be created if not present.\n"
    printf "\t-v DOCKER\tUse the given Docker image specifier for vg.\n"
    printf "\t-R RUN_ID\tUse or restart the given run ID.\n"
    printf "\t-r REGION\tRun on the given named region (21, MHC, BRCA1).\n"
    exit 1
}

while getopts "hdp:t:c:v:R:r:" o; do
    case "${o}" in
        d)
            PREFIX="echo"
            ;;
        p)
            TOIL_VG_PACKAGE="${OPTARG}"
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
        r)
            INPUT_DATA_MODE="${OPTARG}"
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "3" ]]; then
    # Too few arguments
    usage
fi

KEYPAIR_NAME="${1}"
shift
GRAPHS_URL="${1}"
shift
OUTPUT_URL="${1}"
shift

if [[ "$#" -gt "0" ]]; then
    # Options or extra arguments after positional arguments, which getopts can't handle
    echo 1>&2 "Error: trailing option or unrecognized positional argument: ${1}"
    exit 1
fi

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
    TOIL_APPLIANCE_SELF="${TOIL_APPLIANCE_SELF}" $PREFIX toil launch-cluster "${CLUSTER_NAME}" --leaderNodeType=t2.medium -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"
    
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
        # In several chunks
        READ_CHUNKS="32"
        # Define a region name to process. This sets the name that the graphs and
        # indexes will be saved/looked for under.
        REGION_NAME="CHR21"
        # Define the contigs we are using
        GRAPH_CONTIGS=("21")
        # Define the region to build the graph on, as contig[:start-end]
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}")
        # Define the VCF and FASTA basenames. We assume the VCF has a TBI.
        GRAPH_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg19-CHR21.vcf.gz")
        GRAPH_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/CHR21.fa"
        ;;
    MHC)
        # Actually do a smaller test
        READ_COUNT="100000"
        READ_CHUNKS="2"
        REGION_NAME="MHC"
        GRAPH_CONTIGS=("6")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:28510119-33480577")
        GRAPH_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-MHC.vcf.gz")
        GRAPH_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr6.fa.gz"
        ;;
    BRCA1)
        # Do just BRCA1 and a very few reads
        READ_COUNT="1000"
        READ_CHUNKS="1"
        REGION_NAME="BRCA1"
        GRAPH_CONTIGS=("17")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:43044293-43125483")
        GRAPH_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-BRCA1.vcf.gz")
        GRAPH_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr17.fa.gz"
        ;;
    *)
        echo 1>&2 "Unknown input data set ${INPUT_DATA_MODE}"
        exit 1
esac

# Where do we store our job trees under?
JOB_TREE_BASE="aws:us-west-2:${RUN_ID}"

# What trees should we use for each step? Work it out now so cleanup can clean them later.
# Note that job tree names can't have / so we use -
JOB_TREE_CONSTRUCT="${JOB_TREE_BASE}-construct"
JOB_TREE_SIM="${JOB_TREE_BASE}-sim"
JOB_TREE_MAPEVAL="${JOB_TREE_BASE}-mapeval"

echo "Running run ${RUN_ID} as ${KEYPAIR_NAME} on ${GRAPH_REGIONS[@]} for ${READ_COUNT} reads into ${OUTPUT_URL}"

# Make sure we don't leave the cluster running or data laying around on exit.
function clean_up() {
    set +e
    
    # Delete the Toil intermediates we could have used to restart jobs, since
    # we have a lot of Toil runs and no good way to restart just one
    $PREFIX toil clean "${JOB_TREE_CONSTRUCT}"
    $PREFIX toil clean "${JOB_TREE_SIM}"
    $PREFIX toil clean "${JOB_TREE_MAPEVAL}"
    
    if [[ "${PERSISTENT_CLUSTER}" == "0" ]]; then
        # Destroy the cluster
        destroy_cluster "${CLUSTER_NAME}"
    else
        echo "Leaving cluster ${CLUSTER_NAME} running!"
        echo "Destroy with: toil destroy-cluster '${CLUSTER_NAME}' -z us-west-2a"
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
# preemptable nodes by default. But some jobs still demand non-preemptable
# nodes. So we have some r3.8xlarge:0.85 (with a bit) and some r3.8xlarge (on
# demand).
TOIL_CLUSTER_OPTS=(--realTimeLogging --logInfo \
    --batchSystem mesos --provisioner=aws "--mesosMaster=${MASTER_IP}:5050" \
    "--nodeTypes=${NODE_TYPES}" --defaultPreemptable "--maxNodes=${MAX_NODES}" "--minNodes=${MIN_NODES}" \
    --alphaPacking 2.0)


########################################################################################################

# Now we are set up to do the actual experiment.

# Decide what graphs to run
GRAPH_URLS=()
GAM_NAMES=()

# We want the actual graph with its indexes (minus the sample under test)
GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_filter")
GAM_NAMES+=("snp1kg")

GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}")
GAM_NAMES+=("snp1kg-minaf")

# We want a primary negative control
GRAPH_URLS+=("${GRAPHS_URL}/primary")
GAM_NAMES+=("primary")

GRAPHS_READY=1
for GRAPH_BASE_URL in "${GRAPH_URLS[@]}" ; do
    # For each graph we want to run
    for SUFFIX in ".vg" ".xg" ".gcsa" ".gcsa.lcp" ; do
        # For each file requiredf for the graph
        
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

# TODO: Check if all the expected output graphs exist and only run if not.

if [[ "${GRAPHS_READY}" != "1" ]] ; then
    # Graphs need to be generated
    
    # Construct the graphs
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg construct \
        "${JOB_TREE_CONSTRUCT}" \
        "$(url_to_store "${GRAPHS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --vcf "${GRAPH_VCF_URLS[@]}" \
        --fasta "${GRAPH_FASTA_URL}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --control_sample "${SAMPLE_NAME}" \
        --haplo_sample "${SAMPLE_NAME}" \
        "${FILTER_OPTS[@]}" \
        --regions "${GRAPH_REGIONS[@]}" \
        --min_af "${MIN_AF}" \
        --primary \
        --gcsa_index \
        --xg_index \
        --gbwt_index \
        --snarls_index \
        "${TOIL_CLUSTER_OPTS[@]}"
        
fi

# Now work out where in there these simulated reads belong
READS_URL="${GRAPHS_URL}/sim-${READ_SEED}-${READ_COUNT}-${READ_CHUNKS}"

if ! aws s3 ls >/dev/null "${READS_URL}/true.pos" ; then 
    # Now we need to simulate reads from the two haplotypes
    # This will make a "sim.gam"
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg sim \
        "${JOB_TREE_SIM}" \
        "${GRAPHS_URL}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_0.xg" \
        "${GRAPHS_URL}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_1.xg" \
        "${READ_COUNT}" \
        "$(url_to_store "${READS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --annotate_xg "${GRAPHS_URL}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
        --gam \
        --fastq_out \
        --seed "${READ_SEED}" \
        --sim_chunks "${READ_CHUNKS}" \
        --fastq "${TRAINING_FASTQ}" \
        "${TOIL_CLUSTER_OPTS[@]}"
    
fi

# Run one big mapeval run that considers all conditions we are interested in
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg mapeval \
    "${JOB_TREE_MAPEVAL}" \
    "$(url_to_store "${OUTPUT_URL}")" \
    --whole_genome_config \
    "${VG_DOCKER_OPTS[@]}" \
    --index-bases "${GRAPH_URLS[@]}" \
    --gam-names "${GAM_NAMES[@]}" \
    --multipath \
    --use-gbwt \
    --strip-gbwt \
    --use-snarls \
    --fastq "${READS_URL}/sim.fq.gz" \
    --truth "${READS_URL}/true.pos" \
    --plot-sets primary-mp-pe,primary-mp,snp1kg-mp-pe,snp1kg-mp,snp1kg-gbwt-mp-pe,snp1kg-gbwt-mp,snp1kg-minaf-mp-pe,snp1kg-minaf-mp \
    "${TOIL_CLUSTER_OPTS[@]}"

# Cluster (if desired) and trees will get cleaned up by the exit trap

