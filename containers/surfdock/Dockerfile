# SurfDock server container for easydock
#
# Build from the containers/surfdock/ directory:
#   docker build -t surfdock-server .
#
# SurfDock is cloned from the latest master commit at build time.

FROM rapidsai/miniforge-cuda:25.08-cuda12.8.0-base-ubuntu24.04-py3.10@sha256:871560f1d2bc8ba85d1cfe43c97bce786f34bc685b2a9222fe2abd289afd1831

RUN apt-get update && \
    apt-get install --no-install-recommends --yes \
    gcc build-essential git git-lfs wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Clone SurfDock from the latest master commit; pull LFS objects (model weights, tools)
RUN git clone https://github.com/CAODH/SurfDock /app/SurfDock/SurfDock && \
    cd /app/SurfDock/SurfDock && git lfs pull

WORKDIR /app/SurfDock

COPY environment.yaml /app/SurfDock/environment.yaml
RUN conda env create -f /app/SurfDock/environment.yaml

RUN echo "conda activate SurfDock" >> ~/.bashrc
RUN ls SurfDock/comp_surface/tools

COPY msms_i86_64Linux2_2.6.1.tar.gz SurfDock/comp_surface/tools/msms_i86_64Linux2_2.6.1.tar.gz
# Extract bundled surface computation tools (tarballs committed to the repo)
RUN tar -xf /app/SurfDock/SurfDock/comp_surface/tools/APBS_PDB2PQR.tar.gz \
        -C /app/SurfDock/SurfDock/comp_surface/tools/ && \
    tar -xf /app/SurfDock/SurfDock/comp_surface/tools/msms_i86_64Linux2_2.6.1.tar.gz \
        -C /app/SurfDock/SurfDock/comp_surface/tools/

WORKDIR /app/SurfDock/SurfDock
RUN git clone https://github.com/facebookresearch/esm

ENV precomputed_arrays=/app/SurfDock/precomputed/precomputed_arrays
ENV KMP_AFFINITY=disabled
ENV TORCH_HOME=/root/.cache/torch


RUN conda run --no-capture-output -n SurfDock python -c "from utils import so3, torus"
RUN conda run --no-capture-output -n SurfDock python -c "from esm import pretrained; pretrained.load_model_and_alphabet('esm2_t33_650M_UR50D')"

# Copy server files from the build context (containers/surfdock/)
COPY surfdock_server.py     /app/SurfDock/SurfDock/surfdock_server.py
COPY score_dataset_patched.py /app/SurfDock/SurfDock/score_in_place_dataset/score_dataset.py
COPY check_esm_embedding.py /app/SurfDock/SurfDock/bash_scripts/test_scripts/check_esm_embedding.py
COPY global_vars.py         /app/SurfDock/SurfDock/comp_surface/prepare_target/default_config/global_vars.py
COPY init_server.sh         /app/SurfDock/SurfDock/bash_scripts/test_scripts/init_server.sh
COPY dock_server.sh         /app/SurfDock/SurfDock/bash_scripts/test_scripts/dock_server.sh
COPY dock_server_cpu.sh     /app/SurfDock/SurfDock/bash_scripts/test_scripts/dock_server_cpu.sh
RUN chmod +x /app/SurfDock/SurfDock/bash_scripts/test_scripts/init_server.sh \
             /app/SurfDock/SurfDock/bash_scripts/test_scripts/dock_server.sh \
             /app/SurfDock/SurfDock/bash_scripts/test_scripts/dock_server_cpu.sh

WORKDIR /app/SurfDock

ENTRYPOINT ["/bin/sh", "-c", "\
case \"$1\" in \
  server) shift; exec conda run --no-capture-output -n SurfDock python /app/SurfDock/SurfDock/surfdock_server.py \"$@\" ;; \
  help)   echo \"Use: docker run <image> server\" ;; \
  *)      exec \"$@\" ;; \
esac", "--"]
