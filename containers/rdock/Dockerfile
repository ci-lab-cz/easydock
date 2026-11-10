# rDock docking server for EasyDock (--program server).
# rDock (CBDD/rDock, LGPL-3.0) is built from source at a pinned commit; upstream's `make test` gates the build.
# Runtime: Ubuntu 24.04 with the distribution's Python, RDKit and Open Babel (receptor preparation).
#
#   docker build -t rdock .
#   docker run --rm -i rdock help

FROM ubuntu:24.04 AS build
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        make g++ libpopt-dev libpopt0 git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ARG RDOCK_REPO=https://github.com/CBDD/rDock
ARG RDOCK_SHA=ee5226da575306d2b321473a94ff4491e591b347
WORKDIR /src/rDock
RUN git init -q && git fetch -q --depth 1 "$RDOCK_REPO" "$RDOCK_SHA" && git checkout -q FETCH_HEAD \
    && test "$(git rev-parse HEAD)" = "$RDOCK_SHA"

# build_lib first: upstream's `build` target runs build_lib and build_bin as parallel sub-makes, so with -j a binary
# can link against a half-built libRbt.so. RUNPATH $ORIGIN/../lib lets the binaries find the library from any prefix.
RUN set -- RBT_VERSION="v25.02.$(echo "$RDOCK_SHA" | cut -c1-7)" CXX_EXTRA_FLAGS='-Wl,-rpath,\$$ORIGIN/../lib' \
    && make -j"$(nproc)" "$@" build_lib && make -j"$(nproc)" "$@"
RUN make test
RUN make install PREFIX=/opt/rdock RBT_ROOT=/opt/rdock \
    && printf 'repo %s\ncommit %s\n' "$RDOCK_REPO" "$RDOCK_SHA" > /opt/rdock/SOURCE

FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-rdkit python3-numpy openbabel libpopt0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/rdock /opt/rdock
COPY rdock_server.py rdock_prepare_receptor.py /opt/rdock/

ENV RBT_ROOT=/opt/rdock \
    PATH=/opt/rdock/bin:$PATH \
    LD_LIBRARY_PATH=/opt/rdock/lib

ENTRYPOINT ["/bin/bash", "-c", "\
case \"$1\" in \
    server)  shift; exec python3 /opt/rdock/rdock_server.py \"$@\" ;; \
    prepare) shift; exec python3 /opt/rdock/rdock_prepare_receptor.py \"$@\" ;; \
    help) echo 'Usage: docker run --rm -i [-v /data:/data] rdock <command>'; \
          echo '  server            JSON Lines docking server for easydock --program server'; \
          echo '  prepare ...       prepare a receptor: -p protein_H.pdb (-l ref.sdf | --center X Y Z) -o DIR'; \
          echo '  rbdock|rbcavity   the rDock programs; any other command is executed as given' ;; \
    *) exec \"$@\" ;; \
esac", "--"]
