FROM condaforge/miniforge3:latest

# ============================================================================
# Build for CPU (default) or GPU by overriding ENV_FILE at build time.
#
#   CPU:  docker compose build
#   GPU:  docker compose -f docker-compose.yml -f docker-compose.gpu.yml build
#
# DiffSBDD pins its dependencies in a conda environment file. Unlike a simple
# pip torch swap, pytorch-scatter must match the exact torch + CUDA build, so
# we select the whole environment file rather than a single index URL:
#
#   ENV_FILE=/tmp/environment.cpu.yaml  -> CPU-only (default)
#   ENV_FILE=/tmp/environment.gpu.yaml  -> upstream CUDA 11.8 build
#
# To use the GPU at runtime the HOST needs a recent NVIDIA driver and the
# NVIDIA Container Toolkit; the compose GPU override reserves the device.
# ============================================================================
ARG ENV_FILE=/tmp/environment.cpu.yaml

# Layer 1: System packages (openbabel/rdkit come from conda; keep build tools
#          for any source builds and wget/git for fetching checkpoints).
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Layer 2: Bring in both environment files so ENV_FILE can select either.
# Both are vendored under docker/ so the build is self-contained (it does not
# require the upstream environment.yaml to be present in the build context).
#   docker/environment.gpu.yaml (upstream CUDA 11.8) -> /tmp/environment.gpu.yaml
#   docker/environment.cpu.yaml (CPU-only)           -> /tmp/environment.cpu.yaml
COPY docker/environment.gpu.yaml /tmp/environment.gpu.yaml
COPY docker/environment.cpu.yaml /tmp/environment.cpu.yaml

# Layer 3: Create the `diffsbdd` conda env from the selected file and report
#          which torch build landed (CPU-only vs CUDA).
RUN mamba env create -f ${ENV_FILE} && \
    mamba clean -afy && \
    conda run -n diffsbdd python -c "import torch; print('torch', torch.__version__, '| CUDA build:', torch.version.cuda or 'CPU-only')" && \
    conda run -n diffsbdd python -c "import torch_scatter; print('torch_scatter OK')"

# Layer 4: Workspace (DiffSBDD source is bind-mounted here at run time).
WORKDIR /workspace
ENV PYTHONUNBUFFERED=1
RUN echo "source activate diffsbdd" > ~/.bashrc
CMD ["/bin/bash"]
