#!/usr/bin/env bash
# test-population-mapping.sh: evaluate the performance impact of population-aware mapping using toil-vg mapeval on AWS

set -ex

# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/adamnovak/toil-vg.git@1f2d99243b1de3d8b6b9ddbbb48ab855213ebdd3#egg=toil-vg"

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
VG_DOCKER_OPTS=("--vg_docker" "quay.io/vgteam/vg:dev-v1.7.0-76-g4f04d04b-t169-run")

# What node types should we use?
# Comma-separated, with :bid-in-dollars after the name for spot nodes
# We need non-preemptable i3.4xlarge at least to get ~3.8TB storage available so the GCSA indexing jobs will have somewhere to run.
# And we also need more memory (?) than that so some of the later jobs will run.
NODE_TYPES="i3.8xlarge,i3.8xlarge:0.90"
# How many nodes should we use at most per type?
# Also comma-separated.
# TODO: These don't sort right pending https://github.com/BD2KGenomics/toil/issues/2195
# We can only get the limits right for preemptable vs. nonpreemptable for the same thing
MAX_NODES="8,8"
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
        # And the offsets for calling on each of the contigs (0-based)
        GRAPH_CONTIG_OFFSETS=("0")
        # Define the region to build the graph on, as contig[:start-end], 1-based
        # If a non-whole-contig region is specified for any contig, we will need to convert this to
        # BED, so a region must be specified for all contigs.
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}")
        # Define the VCF and FASTA basenames. We assume each VCF has a TBI.
        # If multiple files are used they must match up with the contigs.
        # Note that this is hg19 and not GRCh38
        CONSTRUCT_VCF_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg19-CHR21.vcf.gz")
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/CHR21.fa")
        # What FASTA should we use for BWA mapping? It needs to have just the selected regions cut out.
        MAPPING_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        # What VCF should we use for the truth? Must be a single VCF.
        EVALUATION_VCF_URL="${CONSTRUCT_VCF_URLS[0]}"
        # And a single FASTA
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        # And what high confidence regions should we use there?
        # This can't be specified if we are using regions on GRAPH_REGIONS
        EVALUATION_BED_URL=""
        
        # We will process an interleaved fastq or a BAM of real reads and evaluate its variant calls too.
        # This can be empty, but these are mutually exclusive
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL=""
        ;;
    MHC)
        # Actually do a smaller test
        READ_COUNT="100000"
        READ_CHUNKS="2"
        REGION_NAME="MHC"
        GRAPH_CONTIGS=("6")
        GRAPH_CONTIG_OFFSETS=("28510118")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:28510119-33480577")
        CONSTRUCT_VCF_URLS=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr6_GRCh38_sites.20170504.vcf.gz")
        # We had "s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-MHC.vcf.gz" but it is malformed
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr6.fa.gz")
        MAPPING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/MHC.fa"
        EVALUATION_VCF_URL="${CONSTRUCT_VCF_URLS[0]}"
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/platinum_NA12878_MHC.fq.gz"
        REAL_REALIGN_BAM_URL=""
        ;;
    BRCA1)
        # Do just BRCA1 and a very few reads
        READ_COUNT="20000"
        READ_CHUNKS="2"
        REGION_NAME="BRCA1"
        GRAPH_CONTIGS=("17")
        GRAPH_CONTIG_OFFSETS=("43044292")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[0]}:43044293-43125483")
        CONSTRUCT_VCF_URLS=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr17_GRCh38.genotypes.20170504.vcf.gz")
        # We had "s3://cgl-pipeline-inputs/vg_cgl/bakeoff/1kg_hg38-BRCA1.vcf.gz" but it is malformed
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/bakeoff/chr17.fa.gz")
        MAPPING_FASTA_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/BRCA1.fa"
        EVALUATION_VCF_URL="${CONSTRUCT_VCF_URLS[0]}"
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        EVALUATION_BED_URL=""
        REAL_FASTQ_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff/platinum_NA12878_BRCA1.fq.gz"
        REAL_REALIGN_BAM_URL=""
        ;;
    WG38)
        # Do 10m pairs on the whole genome (GRCh38)
        # TODO: Compose a whole genome 38 VCF
        READ_COUNT="10000000"
        READ_CHUNKS="32"
        REGION_NAME="WG38"
        # We do all the chroms except Y because NA12878 is XX AFAIK
        GRAPH_CONTIGS=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
        GRAPH_CONTIG_OFFSETS=("0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0")
        GRAPH_REGIONS=("${GRAPH_CONTIGS[@]}")
        CONSTRUCT_VCF_URLS=()
        for CONTIG in "${GRAPH_CONTIGS[@]}" ; do
            # Use the per-chromosome VCFs from 1KG
            CONSTRUCT_VCF_URLS+=("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${CONTIG}_GRCh38.genotypes.20170504.vcf.gz")
            # Use the FASTAs from 
        done
        # I built this by stripping the "chr" off of ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
        # TODO: Won't that cause trouble with the presence/absence of decoys in the graph being confounded with mapper?
        CONSTRUCT_FASTA_URLS=("s3://cgl-pipeline-inputs/vg_cgl/pop-map/input/GRCh38.fa.gz")
        MAPPING_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        # Evaluate against Platinum Genomes/GIAB hybrid
        EVALUATION_VCF_URL="ftp://platgene_ro:@ussd-ftp.illumina.com/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz"
        EVALUATION_FASTA_URL="${CONSTRUCT_FASTA_URLS[0]}"
        EVALUATION_BED_URL="ftp://platgene_ro:@ussd-ftp.illumina.com/2017-1.0/hg38/hybrid/hg38.hybrid.bed.gz"
        REAL_FASTQ_URL=""
        REAL_REALIGN_BAM_URL="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam"
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
JOB_TREE_CALLEVAL="${JOB_TREE_BASE}-calleval"

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
    
    # Delete the Toil intermediates we could have used to restart jobs, since
    # we have a lot of Toil runs and no good way to restart just one
    $PREFIX toil clean "${JOB_TREE_CONSTRUCT}"
    $PREFIX toil clean "${JOB_TREE_SIM}"
    $PREFIX toil clean "${JOB_TREE_MAPEVAL}"
    $PREFIX toil clean "${JOB_TREE_CALLEVAL}"
    
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
    --alphaPacking 2.0 --metrics)


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

# We want a primary control
GRAPH_URLS+=("${GRAPHS_URL}/primary")
GAM_NAMES+=("primary")

# We want a positive control with just the right variants
GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}")
GAM_NAMES+=("pos-control")


# We want a negative control with no right variants
GRAPH_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minus_${SAMPLE_NAME}")
GAM_NAMES+=("neg-control")

# Check if all the expected output graphs exist and only run if not.
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

if [[ "${GRAPHS_READY}" != "1" ]] ; then
    # Graphs need to be generated
    
    # Construct the graphs
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg construct \
        "${JOB_TREE_CONSTRUCT}" \
        "$(url_to_store "${GRAPHS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --vcf "${CONSTRUCT_VCF_URLS[@]}" \
        --fasta "${CONSTRUCT_FASTA_URLS[@]}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --control_sample "${SAMPLE_NAME}" \
        --haplo_sample "${SAMPLE_NAME}" \
        "${FILTER_OPTS[@]}" \
        --regions "${GRAPH_REGIONS[@]}" \
        --min_af "${MIN_AFS[@]}" \
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
    if [[ "${MIN_AF}" == "0" ]] ; then
        CONDITION_NAMES+=("snp1kg-minaf-mp-pe")
    else
        CONDITION_NAMES+=("snp1kg-minaf${MIN_AF_NUM}-mp-pe")       
    fi
    XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}.xg")
    ((MIN_AF_NUM++_))
done

CONDITION_NAMES+=("pos-control-mp-pe")
XG_URLS+=("${GRAPHS_URL}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}.xg")

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
    SIM_BAM_URLS+=("${SIM_ALIGNMENTS_URL}/${CONDITION_NAME}-surject.bam")
    BAM_NAMES+=("${CONDITION_NAME}-surject")
done

# Check if all the expected output alignments exist and only run if not.
SIM_ALIGNMENTS_READY=1
for ALIGNMENT_URL in "${SIM_GAM_URLS[@]}" "${SIM_BAM_URLS[@]}" ; do
    if ! aws s3 ls >/dev/null "${ALIGNMENT_URL}" ; then
        # The alignments are not ready yet because this file is missing
        echo "Need to generate alignment file ${ALIGNMENT_URL}"
        SIM_ALIGNMENTS_READY=0
        break
    fi
done

if [[ "${SIM_ALIGNMENTS_READY}" != "1" ]] ; then
    
    # Run one big mapeval run that considers all conditions we are interested in
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg mapeval \
        "${JOB_TREE_MAPEVAL}" \
        "$(url_to_store "${SIM_ALIGNMENTS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --index-bases "${GRAPH_URLS[@]}" \
        --gam-names "${GAM_NAMES[@]}" \
        --multipath \
        --use-gbwt \
        --strip-gbwt \
        --use-snarls \
        --surject \
        --bwa --fasta "${MAPPING_FASTA_URL}" \
        --fastq "${READS_URL}/sim.fq.gz" \
        --truth "${READS_URL}/true.pos" \
        --plot-sets \
        "primary-mp-pe,snp1kg-mp-pe,snp1kg-gbwt-mp-pe,snp1kg-minaf-mp-pe,snp1kg-minaf-gbwt-mp-pe,pos-control-mp-pe,neg-control-mp-pe" \
        "primary-mp,snp1kg-mp,snp1kg-gbwt-mp,snp1kg-minaf-mp,snp1kg-minaf-gbwt-mp,pos-control-mp,neg-control-mp" \
        "bwa-mem-pe,snp1kg-gbwt-mp-pe,snp1kg-pe" \
        "bwa-mem,snp1kg-gbwt-mp,snp1kg" \
        "${TOIL_CLUSTER_OPTS[@]}"
fi

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
        REAL_BAM_URLS+=("${REAL_ALIGNMENTS_URL}/${CONDITION_NAME}-surject.bam")
    done

    # Make sure they exist
    REAL_ALIGNMENTS_READY=1
    for ALIGNMENT_URL in "${REAL_GAM_URLS[@]}" "${REAL_BAM_URLS[@]}" ; do
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
    
        # Run a mapeval run just to map the real reads, under all conditions
        $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg mapeval \
            "${JOB_TREE_MAPEVAL}" \
            "$(url_to_store "${REAL_ALIGNMENTS_URL}")" \
            --whole_genome_config \
            "${VG_DOCKER_OPTS[@]}" \
            --index-bases "${GRAPH_URLS[@]}" \
            --gam-names "${GAM_NAMES[@]}" \
            --multipath \
            --use-gbwt \
            --strip-gbwt \
            --use-snarls \
            --surject \
            --bwa --fasta "${MAPPING_FASTA_URL}" \
            "${DATA_OPTS[@]}" \
            --skip-eval \
            "${TOIL_CLUSTER_OPTS[@]}"
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
        BED_OPTS+=(--vcfeval_bed_regions "${TEMP_BED}" --clip_only)
    fi
else
    # We have a bed we need to use for the high-confidence regions
    BED_OPTS+=(--vcfeval_bed_regions "${EVALUATION_BED_URL}" --clip_only)
fi

# Now the sim calls
SIM_CALLS_READY=1
if ! aws s3 ls >/dev/null "${SIM_CALLS_URL}/plots/roc-weighted.svg" ; then
    # Use just this one file as a marker of calleval done-ness.
    # TODO: use other files/VCFs.
    SIM_CALLS_READY=0
fi

# It would be nice if we could run genotype, but it is extremely slow (~2.5 hours per chunk on chr21 sim data)

if [[ "${SIM_CALLS_READY}" != "1" ]] ; then
    $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg calleval \
        "${JOB_TREE_CALLEVAL}" \
        "$(url_to_store "${SIM_CALLS_URL}")" \
        --whole_genome_config \
        "${VG_DOCKER_OPTS[@]}" \
        --gams "${SIM_GAM_URLS[@]}" \
        --gam_names "${CONDITION_NAMES[@]}" \
        --bams "${SIM_BAM_URLS[@]}" \
        --bam_names "${BAM_NAMES[@]}" \
        --xg_paths "${XG_URLS[@]}" \
        --chroms "${GRAPH_CONTIGS[@]}" \
        --vcf_offsets "${GRAPH_CONTIG_OFFSETS[@]}" \
        --vcfeval_fasta "${CONSTRUCT_FASTA_URL}" \
        --vcfeval_baseline "${EVALUATION_VCF_URL}" \
        --call \
        "${BED_OPTS[@]}" \
        --sample_name "${SAMPLE_NAME}" \
        --plot_sets \
            "primary-mp-pe-call,snp1kg-mp-pe-call,snp1kg-gbwt-mp-pe-call,snp1kg-minaf-mp-pe-call,snp1kg-minaf-gbwt-mp-pe-call,pos-control-mp-pe-call,neg-control-mp-pe-call" \
            "bwa-pe-fb,snp1kg-gbwt-mp-pe-call,snp1kg-pe-call" \
        "${TOIL_CLUSTER_OPTS[@]}"
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

    if [[ "${REAL_CALLS_READY}" != "1" ]] ; then
        $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg calleval \
            "${JOB_TREE_CALLEVAL}" \
            "$(url_to_store "${REAL_CALLS_URL}")" \
            --whole_genome_config \
            "${VG_DOCKER_OPTS[@]}" \
            --gams "${REAL_GAM_URLS[@]}" \
            --gam_names "${CONDITION_NAMES[@]}" \
            --bams "${REAL_BAM_URLS[@]}" \
            --bam_names "${BAM_NAMES[@]}" \
            --xg_paths "${XG_URLS[@]}" \
            --chroms "${GRAPH_CONTIGS[@]}" \
            --vcf_offsets "${GRAPH_CONTIG_OFFSETS[@]}" \
            --vcfeval_fasta "${EVALUATION_FASTA_URL}" \
            --vcfeval_baseline "${EVALUATION_VCF_URL}" \
            --call \
            "${BED_OPTS[@]}" \
            --sample_name "${SAMPLE_NAME}" \
            --plot_sets \
                "primary-mp-pe-call,snp1kg-mp-pe-call,snp1kg-gbwt-mp-pe-call,snp1kg-minaf-mp-pe-call,pos-control-mp-pe-call,neg-control-mp-pe-call" \
                "bwa-pe-fb,snp1kg-gbwt-mp-pe-call,snp1kg-pe-call" \
            "${TOIL_CLUSTER_OPTS[@]}"
    fi

fi


# Cluster (if desired) and trees will get cleaned up by the exit trap

