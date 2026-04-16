FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && apt-get install -y --no-install-recommends \
        wget \
        build-essential \
        python3 \
        python3-pip \
        ocl-icd-opencl-dev \
        opencl-headers \
        libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*

# Build Boost 1.84.0 from source (only the libs needed by the Vina-GPU binary)
RUN cd /tmp \
    && wget -L --tries=3 -O boost_1_84_0.tar.gz \
        "https://sourceforge.net/projects/boost/files/boost/1.84.0/boost_1_84_0.tar.gz/download" \
    && tar xzf boost_1_84_0.tar.gz \
    && cd boost_1_84_0 \
    && ./bootstrap.sh --with-libraries=program_options,system,filesystem --prefix=/usr/local \
    && ./b2 install -j$(nproc) \
    && ldconfig \
    && cd /tmp && rm -rf boost_1_84_0 boost_1_84_0.tar.gz

# Register NVIDIA OpenCL platform so the ICD loader finds it at runtime (--gpus / nvidia-docker)
RUN mkdir -p /etc/OpenCL/vendors \
    && echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# Download AutoDock-Vina-GPU-2.1 (vina-gpu)
ENV BASE=https://raw.githubusercontent.com/DeltaGroupNJUPT/Vina-GPU-2.1/main
RUN mkdir -p /opt/vinagpu/vina-gpu \
    && wget -q -O /opt/vinagpu/vina-gpu/AutoDock-Vina-GPU-2-1  $BASE/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1 \
    && wget -q -O /opt/vinagpu/vina-gpu/Kernel1_Opt.bin         $BASE/AutoDock-Vina-GPU-2.1/Kernel1_Opt.bin \
    && wget -q -O /opt/vinagpu/vina-gpu/Kernel2_Opt.bin         $BASE/AutoDock-Vina-GPU-2.1/Kernel2_Opt.bin \
    && chmod +x /opt/vinagpu/vina-gpu/AutoDock-Vina-GPU-2-1

# Download QuickVina2-GPU-2.1 (qvina-gpu)
RUN mkdir -p /opt/vinagpu/qvina-gpu \
    && wget -q -O /opt/vinagpu/qvina-gpu/QuickVina2-GPU-2-1    $BASE/QuickVina2-GPU-2.1/QuickVina2-GPU-2-1 \
    && wget -q -O /opt/vinagpu/qvina-gpu/Kernel1_Opt.bin        $BASE/QuickVina2-GPU-2.1/Kernel1_Opt.bin \
    && wget -q -O /opt/vinagpu/qvina-gpu/Kernel2_Opt.bin        $BASE/QuickVina2-GPU-2.1/Kernel2_Opt.bin \
    && chmod +x /opt/vinagpu/qvina-gpu/QuickVina2-GPU-2-1

# Download QuickVina-W-GPU-2.1 (qvinaw-gpu)
RUN mkdir -p /opt/vinagpu/qvinaw-gpu \
    && wget -q -O /opt/vinagpu/qvinaw-gpu/QuickVina-W-GPU-2-1  $BASE/QuickVina-W-GPU-2.1/QuickVina-W-GPU-2-1 \
    && wget -q -O /opt/vinagpu/qvinaw-gpu/Kernel1_Opt.bin       $BASE/QuickVina-W-GPU-2.1/Kernel1_Opt.bin \
    && wget -q -O /opt/vinagpu/qvinaw-gpu/Kernel2_Opt.bin       $BASE/QuickVina-W-GPU-2.1/Kernel2_Opt.bin \
    && chmod +x /opt/vinagpu/qvinaw-gpu/QuickVina-W-GPU-2-1

# Download CPU binaries from QVina repo (https://github.com/QVina/qvina)
ENV QVINA_BASE=https://raw.githubusercontent.com/QVina/qvina/master/bin
RUN mkdir -p /opt/vinagpu/qvina \
    && wget -q -O /opt/vinagpu/qvina/qvina2.1  $QVINA_BASE/qvina2.1 \
    && chmod +x /opt/vinagpu/qvina/qvina2.1

RUN mkdir -p /opt/vinagpu/qvinaw \
    && wget -q -O /opt/vinagpu/qvinaw/qvina-w  $QVINA_BASE/qvina-w \
    && chmod +x /opt/vinagpu/qvinaw/qvina-w

RUN mkdir -p /opt/vinagpu/vina \
    && wget -q -O /opt/vinagpu/vina/vina  $QVINA_BASE/vina \
    && chmod +x /opt/vinagpu/vina/vina

COPY vinagpu_server.py /opt/vinagpu/vinagpu_server.py

ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

ENTRYPOINT ["/bin/bash", "-c", "\
case \"$1\" in \
    server) shift; exec python3 /opt/vinagpu/vinagpu_server.py \"$@\" ;; \
    help) echo 'Usage: docker run <image> server'; \
          echo 'GPU programs: vina-gpu, qvina-gpu, qvinaw-gpu'; \
          echo 'CPU programs: vina, qvina, qvinaw'; \
          echo 'Set the program via '\''program'\'' key in the server config (init_server block or top-level).' ;; \
    *) exec \"$@\" ;; \
esac", "--"]
