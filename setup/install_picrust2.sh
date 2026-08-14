#!/usr/bin/env bash
################################################################################
# Script:  install_picrust2.sh
# Purpose:
#   - Create a dedicated conda environment named "picrust2" and install
#     PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction
#     of Unobserved States) into it from the bioconda channel, following the
#     official installation instructions published by the Huttenhower Lab:
#     https://huttenhower.sph.harvard.edu/picrust/
#   - No specific PICRUSt2 version is pinned; conda will resolve and install
#     the latest version available on bioconda at the time this script is run.
#   - Linux only. The script exits before making changes on other operating
#     systems.
#
# Reference:
#   - Huttenhower Lab PICRUSt2 page: https://huttenhower.sph.harvard.edu/picrust/
#   - PICRUSt2 GitHub wiki (user manual): https://github.com/picrust/picrust2/wiki
#   - Citation: Douglas GM, Maffei VJ, Zaneveld JR, et al. PICRUSt2 for
#     prediction of metagenome functions. Nat Biotechnol 38, 685-688 (2020).
#     https://doi.org/10.1038/s41587-020-0548-6
#
# Prerequisite:
#   - A Linux system with a working conda installation (Miniconda or Anaconda)
#     must already be
#     present on this machine and available on the PATH. This script does
#     NOT install conda itself.
#   - If conda is not yet installed, follow the official installation guide
#     before running this script:
#     https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html
#
# Usage (run from the repository root):
#   chmod +x setup/install_picrust2.sh
#   ./setup/install_picrust2.sh
#
################################################################################


# ==============================================================================
# Shell safety options
# ==============================================================================
# -e : Exit immediately if any command exits with a non-zero status.
# -u : Treat unset variables as an error and exit immediately.
# -o pipefail : Return the exit status of the last command in a pipeline that
#               failed, rather than the exit status of the final command only.
# Together these prevent the script from silently continuing after a failure
# (e.g., if the conda environment creation step fails partway through).
set -euo pipefail


# ==============================================================================
# Configuration parameters
# ==============================================================================
# Name of the conda environment that will be created for PICRUSt2.
# Kept as a variable so it can be changed in one place if needed (e.g. to add
# a version suffix such as "picrust2_2026").
readonly CONDA_ENV_NAME="picrust2"

# Conda channels required to resolve the PICRUSt2 package and its
# dependencies. Order matters: bioconda packages typically depend on
# conda-forge packages, so both channels must be available during solving.
readonly CONDA_CHANNEL_BIOCONDA="bioconda"
readonly CONDA_CHANNEL_CONDAFORGE="conda-forge"

# Name of the package to install. No version is pinned here on purpose, so
# conda will install whatever the latest compatible release on bioconda is
# at the time this script is executed.
readonly PICRUSt2_PACKAGE="picrust2"

# URL shown to the user if conda is not found on PATH, pointing to the
# official installation instructions (Miniconda/Anaconda).
readonly CONDA_INSTALL_HELP_URL="https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html"


# ==============================================================================
# Step 1: Verify the operating system and conda installation
# ==============================================================================
# This installer intentionally supports Linux only. Stop before checking conda
# or changing any environment when invoked on another operating system.
echo "Step 1/4: Checking the operating system and conda installation..."

readonly OPERATING_SYSTEM="$(uname -s)"
if [[ "${OPERATING_SYSTEM}" != "Linux" ]]; then
  echo "ERROR: This PICRUSt2 installer supports Linux only (detected: ${OPERATING_SYSTEM})." >&2
  exit 1
fi

# PICRUSt2 is distributed via bioconda, so a working conda installation is a
# hard prerequisite. `command -v conda` returns the path to the conda
# executable if found, or nothing (with a non-zero exit code) if it is not.

if ! command -v conda >/dev/null 2>&1; then
  # conda was not found on PATH: stop here with an informative error message
  # rather than letting the subsequent `conda create` command fail with a
  # less helpful "command not found" error.
  echo "ERROR: conda was not found on your PATH." >&2
  echo "PICRUSt2 is installed via the bioconda channel, so conda (Miniconda or Anaconda) must be installed first." >&2
  echo "Installation instructions: ${CONDA_INSTALL_HELP_URL}" >&2
  exit 1
fi

# Report the detected conda version for the log/record, useful for
# troubleshooting and reproducibility.
echo "Conda detected: $(conda --version)"


# ==============================================================================
# Step 2: Check whether the target environment already exists
# ==============================================================================
# `conda env list` prints all existing environments. Capture its output before
# matching so grep -q cannot close a pipe early under `set -o pipefail`.
echo "Step 2/4: Checking whether the '${CONDA_ENV_NAME}' environment already exists..."

conda_environment_list="$(conda env list)"
if grep -qE "^${CONDA_ENV_NAME}[[:space:]]" <<< "${conda_environment_list}"; then
  echo "A conda environment named '${CONDA_ENV_NAME}' already exists."
  echo "Step 3/4: Skipping environment creation and verifying the existing installation."
else
  # ==============================================================================
  # Step 3: Create the conda environment and install PICRUSt2
  # ==============================================================================
  # This mirrors the official installation command from the Huttenhower Lab
  # PICRUSt2 page, without pinning a specific version:
  #
  #   conda create -n picrust2 -c bioconda -c conda-forge picrust2
  #
  # -n     : name of the new environment
  # -c     : channel(s) to search for packages, listed in priority order
  # --yes  : automatically confirm the installation plan (non-interactive run)
  echo "Step 3/4: Creating conda environment '${CONDA_ENV_NAME}' and installing PICRUSt2 (this may take several minutes)..."

  conda create \
    --name "${CONDA_ENV_NAME}" \
    --channel "${CONDA_CHANNEL_BIOCONDA}" \
    --channel "${CONDA_CHANNEL_CONDAFORGE}" \
    "${PICRUSt2_PACKAGE}" \
    --yes
fi


# ==============================================================================
# Step 4: Confirm installation and print activation instructions
# ==============================================================================
# The environment cannot be activated from within a non-interactive script
# in a way that persists after the script exits (conda activate modifies the
# current shell session). Instead, print clear instructions for the user to
# activate the environment themselves, and verify installation succeeded by
# checking for the picrust2_pipeline.py entry point inside the new
# environment without needing to activate it first.
echo "Step 4/4: Verifying that PICRUSt2 was installed successfully..."

if conda run -n "${CONDA_ENV_NAME}" picrust2_pipeline.py --version >/dev/null 2>&1 &&
   conda run -n "${CONDA_ENV_NAME}" place_seqs.py --help >/dev/null 2>&1 &&
   conda run -n "${CONDA_ENV_NAME}" hsp.py --help >/dev/null 2>&1; then
  echo "PICRUSt2 installed successfully in the '${CONDA_ENV_NAME}' conda environment."
else
  echo "ERROR: Could not verify the PICRUSt2 commands required by this workflow" >&2
  echo "(picrust2_pipeline.py, place_seqs.py, and hsp.py) in '${CONDA_ENV_NAME}'." >&2
  echo "To reinstall from scratch, run:" >&2
  echo "  conda remove -n ${CONDA_ENV_NAME} --all" >&2
  echo "and then run this installer again." >&2
  exit 1
fi

echo ""
echo "================================================================================"
echo "Installation complete."
echo "Activate the environment in your terminal before running PICRUSt2 with:"
echo "  conda activate ${CONDA_ENV_NAME}"
echo ""
echo "If, once activated, you get an error about the default reference files being"
echo "missing, reinstall PICRUSt2 from source instead of bioconda (see the"
echo "'Install From Source' section at https://huttenhower.sph.harvard.edu/picrust/)."
echo "================================================================================"
