#!/usr/bin/env bash
# test-population-mapping.sh: evaluate the performance impact of population-aware mapping using toil-vg mapeval on AWS

set -e

# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/vgteam/toil-vg.git@1d439b2d678e163a35edd5143bc4d1e9ce31990e"

# What Toil should we use?
TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:3.11.0

# What vg should we use?
VG_DOCKER_OPTS=()

# How many nodes should we use at most?
MAX_NODES=6

# What's our unique run ID? Should be lower-case and start with a letter for maximum compatibility.
# See <https://gist.github.com/earthgecko/3089509>
RUN_ID="run$(cat /dev/urandom | LC_CTYPE=C tr -dc 'a-z0-9' | fold -w 32 | head -n 1)"

# What cluster should we use?
CLUSTER_NAME="${RUN_ID}"
MANAGE_CLUSTER=1

# Should we delete the job store when we exit?
# We do by default, and if the run finishes successfully.
# We don't if we're asked to keep it and Toil errors out.
REMOVE_JOBSTORE=1

# Should we delete the outstore at the end of the script, if we're deleting the
# jobstore?
REMOVE_OUTSTORE=1

# What named graph region should we be operating on?
# We can do "21", "whole-genome", "MHC", etc.
INPUT_DATA_MODE="21"

# TODO: Allow using real input data for speed measurement.

# Set a FASTQ to model reads after
TRAINING_FASTQ="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz"
# And a read simulation seed
READ_SEED="90"
# Chunks to simulate in (which affects results)
READ_CHUNKS="32"

# Define the sample to use for synthesizing reads
SAMPLE_NAME="HG00096"

# What min allele frequency limit do we use?
MIN_AF="0.0335570469"

# Put this in front of commands to do or not do them, depending on if we are doing a dry run or not
PREFIX=""

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] KEYPAIR_NAME GRAPHS_PATH OUTPUT_PATH \n"
    printf "\tKEYPAIR_NAME\tSSH keypair in Amazon us-west-2 region to use for cluster authentication\n"
    printf "\tGRAPHS_PATH\tLocal file system path where generated graphs, indexes, and reads should be cached\n"
    printf "\tOUTPUT_PATH\tLocal file syetem path where output data (mapped reads and statistics) should be downloaded\n"
    printf "Options:\n\n"
    printf "\t-d\tDo a dry run\n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg.\n"
    printf "\t-t CONTAINER\tUse the given Toil container in the cluster (default: ${TOIL_APPLIANCE_SELF}).\n"
    printf "\t-c CLUSTER\tUse the given existing Toil cluster.\n"
    printf "\t-v DOCKER\tUse the given Docker image specifier for vg.\n"
    printf "\t-R RUN_ID\tUse or restart the given run ID.\n"
    printf "\t-r REGION\tRun on the given named region (21, MHC, BRCA1).\n"
    printf "\t-k \tKeep the out store and job store in case of error.\n"
    exit 1
}

while getopts "hdp:t:c:v:R:r:k" o; do
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
            MANAGE_CLUSTER=0
            ;;
        v)
            VG_DOCKER_OPTS="--vg_docker ${OPTARG}"
            ;;
        R)
            # This doesn't change the cluster name, which will still be the old run ID if not manually set.
            # That's probably fine.
            RUN_ID="${OPTARG}"
            ;;
        r)
            INPUT_DATA_MODE="${OPTARG}"
            ;;
        k)
            REMOVE_JOBSTORE=0
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
GRAPHS_PATH="${1}"
shift
OUTPUT_PATH="${1}"
shift

if [[ "$#" -gt "0" ]]; then
    # Options or extra arguments after positional arguments, which getopts can't handle
    echo 1>&2 "Error: trailing option or unrecognized positional argument: ${1}"
    exit 1
fi

# Fill in some variables characterizing our input based on the input data mode
case "${INPUT_DATA_MODE}" in
    21)
        # Do a lot of reads
        READ_COUNT="10000000"
        # Define a region name to process. This sets the name that the graphs and
        # indexes will be saved/looked for under.
        REGION_NAME="CHR21"
        # Define the contig we are using
        GRAPH_CONTIG="21"
        # Define the region to build the graph on, as contig[:start-end]
        GRAPH_REGION="${GRAPH_CONTIG}"
        # Define the VCF and FASTA basenames. We assume the VCF has a TBI.
        GRAPH_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg19-CHR21.vcf.gz"
        GRAPH_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/CHR21.fa"
        ;;
    MHC)
        # Actually do a smaller test
        READ_COUNT="100000"
        REGION_NAME="MHC"
        GRAPH_CONTIG="6"
        GRAPH_REGION="${GRAPH_CONTIG}:28510119-33480577"
        GRAPH_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-MHC.vcf.gz"
        GRAPH_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr6.fa.gz"
        ;;
    BRCA1)
        # Do just BRCA1 and a very few reads
        READ_COUNT="1000"
        REGION_NAME="BRCA1"
        GRAPH_CONTIG="17"
        GRAPH_REGION="${GRAPH_CONTIG}:43044293-43125483"
        GRAPH_VCF_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-BRCA1.vcf.gz"
        GRAPH_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr17.fa.gz"
        ;;
    *)
        echo 1>&2 "Unknown input data set ${INPUT_DATA_MODE}"
        exit 1
esac

# Compute where to save our results from the various jobs
OUTPUT_STORE_BASE="aws:us-west-2:cgl-pipeline-inputs/vg_cgl/pop-map-script/runs/${RUN_ID}"
OUTPUT_STORE_URL_BASE="s3://cgl-pipeline-inputs/vg_cgl/pop-map-script/runs/${RUN_ID}"

# Where do we store our jobs under?
JOB_TREE_BASE="aws:us-west-2:${RUN_ID}"

# What trees/outstores should we use for each step? Work it out now so cleanup can fix it
JOB_TREE_CONSTRUCT="${JOB_TREE_BASE}/construct"
OUTPUT_STORE_CONSTRUCT="${OUTPUT_STORE_BASE}/construct"
OUTPUT_STORE_URL_CONSTRUCT="${OUTPUT_STORE_URL_BASE}/construct"
JOB_TREE_SIM="${JOB_TREE_BASE}/sim"
OUTPUT_STORE_SIM="${OUTPUT_STORE_BASE}/sim"
OUTPUT_STORE_URL_SIM="${OUTPUT_STORE_URL_BASE}/sim"
JOB_TREE_MAPEVAL="${JOB_TREE_BASE}/mapeval"
OUTPUT_STORE_MAPEVAL="${OUTPUT_STORE_BASE}/mapeval"
OUTPUT_STORE_URL_MAPEVAL="${OUTPUT_STORE_URL_BASE}/mapeval"

echo "Running run ${RUN_ID} as ${KEYPAIR_NAME} on ${GRAPH_REGION} for ${READ_COUNT} reads into ${OUTPUT_PATH}"

# Make sure we don't leave the cluster running or data laying around on exit.
function clean_up() {
    set +e
    if [[ "${REMOVE_JOBSTORE}" == "1" ]]; then
        # Delete the Toil intermediates we could have used to restart the job
        $PREFIX toil clean "${JOB_TREE_CONSTRUCT}"
        $PREFIX toil clean "${JOB_TREE_SIM}"
        $PREFIX toil clean "${JOB_TREE_MAPEVAL}"
        
        if [[ "${REMOVE_OUTSTORE}" == "1" ]]; then
            # Toil is happily done and we downloaded things OK.
            # Delete the outputs.
            # It is not a problem that we delete simulated reads, etc. because we copied/saved them.
            $PREFIX aws s3 rm --recursive "${OUTPUT_STORE_URL_CONSTRUCT}"
            $PREFIX aws s3 rm --recursive "${OUTPUT_STORE_URL_SIM}"
            $PREFIX aws s3 rm --recursive "${OUTPUT_STORE_URL_MAPEVAL}"
        fi
    fi
    if [[ "${MANAGE_CLUSTER}" == "1" ]]; then
        # Destroy the cluster
        $PREFIX toil destroy-cluster "${CLUSTER_NAME}" -z us-west-2a
    fi
}
trap clean_up EXIT

if [[ "${MANAGE_CLUSTER}" == "1" ]]; then
    # Start up a cluster
    TOIL_APPLIANCE_SELF="${TOIL_APPLIANCE_SELF}" $PREFIX toil launch-cluster "${CLUSTER_NAME}" --leaderNodeType=t2.medium -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"
fi

# We need to manually install git to make pip + git work...
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt update
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt install git -y

# Ignore the old virtualenv if re-using a cluster

# For hot deployment to work, toil-vg needs to be in a virtualenv that can see the system Toil
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" virtualenv --system-site-packages venv

$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install pyyaml
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install aws
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install numpy
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install scipy
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install scikit-learn
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install "${TOIL_VG_PACKAGE}"

# We need the master's IP to make Mesos go
MASTER_IP="$($PREFIX toil ssh-cluster --insecure --zone=us-west-2a --logOff "${CLUSTER_NAME}" hostname -i)"

# Strip some garbage from MASTER_IP
MASTER_IP="${MASTER_IP//[$'\t\r\n ']}"


########################################################################################################

# Now we are set up to do the actual experiment.

if [[ ! -d "${GRAPHS_PATH}" ]]; then
    # Graphs need to be gotten
    
    if [[ -e "${GRAPHS_PATH}" ]]; then
        # It needs to not exist at all
        echo "ERROR: Graph path ${GRAPHS_PATH} is not a directory" 1>&2
        exit 1
    fi
    
    # Construct the graphs
    # Hardcoded constants here go with the snp1kg URLs above.
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg construct \
        "${JOB_TREE_CONSTRUCT}" \
        "${OUTPUT_STORE_CONSTRUCT}" \
        --whole_genome_config \
        ${VG_DOCKER_OPTS} \
        --vcf "${GRAPH_VCF_URL}" \
        --fasta "${GRAPH_FASTA_URL}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --realTimeLogging \
        --control_sample "${SAMPLE_NAME}" \
        --haplo_sample "${SAMPLE_NAME}" \
        --filter_samples "${SAMPLE_NAME}" \
        --regions "${GRAPH_REGION}" \
        --min_af "${MIN_AF}" \
        --primary \
        --gcsa_index \
        --xg_index \
        --gbwt_index \
        --snarls_index
        
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/aws s3 sync --acl public-read "${OUTPUT_STORE_URL_CONSTRUCT}" "${OUTPUT_STORE_URL_CONSTRUCT}"
        
    # Download them
    $PREFIX mkdir -p "${GRAPHS_PATH}"
    $PREFIX aws s3 sync "${OUTPUT_STORE_URL_CONSTRUCT}" "${GRAPHS_PATH}"
fi

# Now work out where in there these simulated reads belong
READS_PATH="${GRAPHS_PATH}/sim-${READ_SEED}-${READ_COUNT}-${READ_CHUNKS}"

if [[ ! -e "${READS_PATH}" ]]; then 
    # Now we need to simulate reads from the two haplotypes
    # This will make a "sim.gam"
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg sim \
        "${JOB_TREE_SIM}" \
        "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_0.xg" \
        "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_1.xg" \
        "${READ_COUNT}" \
        "${OUTPUT_STORE_SIM}" \
        --whole_genome_config \
        ${VG_DOCKER_OPTS} \
        --annotate_xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
        --gam \
        --fastq_out \
        --seed "${READ_SEED}" \
        --sim_chunks "${READ_CHUNKS}" \
        --fastq "${TRAINING_FASTQ}"
        
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/aws s3 sync --acl public-read "${OUTPUT_STORE_URL_SIM}" "${OUTPUT_STORE_URL_SIM}"
        
    # Download them
    $PREFIX mkdir -p "${READS_PATH}"
    $PREFIX aws s3 sync "${OUTPUT_STORE_URL_SIM}" "${GRAPHS_PATH}"
fi

# Now we have both the graphs and reads locally, but built in the cloud.

# Decide what graphs to run
GRAPH_URLS=()
GAM_NAMES=()

# We want the actual graph with its indexes (minus the sample under test)
GRAPH_URLS+=("${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter")
GAM_NAMES+=("snp1kg")

GRAPH_URLS+=("${GRAPHS_PATH}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}")
GAM_NAMES+=("snp1kg-minaf")

# We want a primary negative control
GRAPH_URLS+=("${GRAPHS_PATH}/primary")
GAM_NAMES+=("primary")

# Run one big mapeval run that considers all conditions we are interested in
# TODO: Controls for no haplotype-aware-ness?
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg mapeval \
    "${JOB_TREE_MAPEVAL}" \
    "${OUTPUT_STORE_MAPEVAL}" \
    --whole_genome_config \
    ${VG_DOCKER_OPTS} \
    --index-bases "${GRAPH_URLS[@]}" \
    --gam-names "${GRAPH_NAMES[@]}" \
    --multipath \
    --mpmap_opts "--max-paths 10" \
    --use-gbwt \
    --use-snarls \
    --fastq "${READS_PATH}/sim.fq.gz" \
    --truth "${READS_PATH}/true.pos" \
    --realTimeLogging --logInfo \
    --batchSystem mesos --provisioner=aws "--mesosMaster=${MASTER_IP}:5050" \
    --nodeTypes=r3.8xlarge:0.85 --defaultPreemptable --maxNodes=${MAX_NODES}\
    --alphaPacking 2.0
TOIL_ERROR="$?"
    
# Make sure the output is public
$PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/aws s3 sync --acl public-read "${OUTPUT_STORE_URL_MAPEVAL}" "${OUTPUT_STORE_URL_MAPEVAL}"
    
$PREFIX mkdir -p "${OUTPUT_PATH}"
$PREFIX aws s3 sync "${OUTPUT_STORE_URL_MAPEVAL}" "${OUTPUT_PATH}"
DOWNLOAD_ERROR="$?"

if [[ "${TOIL_ERROR}" == "0" ]]; then
    # Toil completed successfully.
    # We will delete the job store
    REMOVE_JOBSTORE=1
fi

if [[ ! "${DOWNLOAD_ERROR}" == "0" ]]; then
    # Download failed
    # We will keep the out store
    # (we also keep if if Toil fails and we're keeping the jobstore)
    REMOVE_OUTSTORE=0
fi

# Cluster, tree, and output will get cleaned up by the exit trap

