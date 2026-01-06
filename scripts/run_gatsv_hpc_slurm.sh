#!/bin/bash

#SBATCH --job-name=run_gatsv_example
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=log_%J.out
#SBATCH --error=log_%J.err

# =============================================================================
# 1) Configure these paths/values - EDIT THESE FOR YOUR ENVIRONMENT
# =============================================================================

DOCKER_IMG="docker://wchukwu/gatsv_docker:latest"
IMAGE_DIR="/path/to/containers"              # Directory where .sif file is or will be stored
IMAGE="${IMAGE_DIR}/gaTSV.sif"               # Full path to .sif file

INDIR="/path/to/input_data"                  # Must contain example.sv.vcf and example_metadata.txt
OUTDIR="/path/to/results/gaTSV/outputs"      # Outputs will be written here

# Analysis parameters
SAMPLE="example"
REF="hg19"
THREADS="${SLURM_CPUS_PER_TASK:-2}"

# =============================================================================
# 2) Load Singularity module
# =============================================================================

module load singularity

# =============================================================================
# 3) Create necessary directories
# =============================================================================

mkdir -p "${IMAGE_DIR}"
mkdir -p "${OUTDIR}"

# =============================================================================
# 4) Pull Docker image and convert to .sif (only if not already present)
# =============================================================================

if [[ ! -f "${IMAGE}" ]]; then
    echo "Pulling Docker image and converting to SIF..."
    singularity pull --force "${IMAGE}" "${DOCKER_IMG}"
else
    echo "SIF image already exists: ${IMAGE}"
fi

# =============================================================================
# 5) Run the analysis inside the container
# =============================================================================

echo "Starting gaTSV analysis..."
echo "Input directory: ${INDIR}"
echo "Output directory: ${OUTDIR}"
echo "Sample: ${SAMPLE}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}"

singularity exec \
    --cleanenv \
    --bind "${INDIR}:/home/input_data:ro" \
    --bind "${OUTDIR}:/home/output:rw" \
    --bind "${OUTDIR}:/out:rw" \
    "${IMAGE}" \
    /scripts/gaTSV_run.sh \
        /home/input_data/example.sv.vcf \
        /home/input_data/example_metadata.txt \
        "${SAMPLE}" \
        "${REF}" \
        "${THREADS}"

# =============================================================================
# 6) Completion message
# =============================================================================

echo "================================================"
echo "Done. Outputs are in: ${OUTDIR}"
echo "================================================"