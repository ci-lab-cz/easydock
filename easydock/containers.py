import logging
import os
import platform
import subprocess
import sys
import traceback
from typing import List, Tuple, Union

from easydock.auxiliary import expand_path


def docker_available(container=None):
    try:
        if not container:
            subprocess.run(
                ["docker", "version"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            return True
        else:
            subprocess.run(
                ["docker", "inspect", container],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            return True
    except (subprocess.CalledProcessError, FileNotFoundError, PermissionError) as e:
        # sys.stderr.write(e)
        # sys.stderr.write('Docker is not installed or properly configured. It should be accessible from the command '
        #                  'line without elevated privileges\n')
        return False


def apptainer_available():
    try:
        subprocess.run(
            ["apptainer", "version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, PermissionError) as e:
        # sys.stderr.write(e)
        # sys.stderr.write('Apptainer is not installed or properly configured. It should be accessible from the command '
        #                  'line without elevated privileges\n')
        return False


def singularity_available():
    try:
        subprocess.run(
            ["singularity", "version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, PermissionError) as e:
        # sys.stderr.write(e)
        # sys.stderr.write('Singularity is not installed or properly configured. It should be accessible from the command '
        #                  'line without elevated privileges\n')
        return False


def apptainer_exec(sif_path: str,
                   inner_cmd: List,
                   binds: List[Union[str, Tuple[str, str]]] = None):
    """
    Example of a call.
    apptainer_exec('~/imgs/unipka.sif',
                   ['protonate', '-i', '/path/to/dir/in/input.smi', '-o', '/path/to/dir/out/output.smi'],
                   ['/path/to/dir/in', '/path/to/dir/out'])

    :param sif_path:
    :param inner_cmd:
    :param binds:
    :return:
    """

    system = platform.system()
    sif_path = expand_path(sif_path)
    tmp = []
    for b in binds:
        if isinstance(b, str):
            tmp.append(expand_path(b))
        elif isinstance(b, tuple) and len(b) == 2:
            tmp.append(":".join(map(expand_path, b)))
        else:
            sys.stderr.write(f'Wrong argument type in apptainer_exec function: {b}. '
                             f'It can be only str or a tuple of two strings.\n')
    binds = tmp

    # Normalize bind paths (host_path:container_path)
    apptainer_bind_args = []
    for b in binds:
        if isinstance(b, str):
            apptainer_bind_args.append(b)
        elif isinstance(b, tuple) and len(b) == 2:
            apptainer_bind_args.append(":".join(b))
    if apptainer_bind_args:
        apptainer_bind_args = ['--bind', ','.join(apptainer_bind_args)]

    if system == "Linux":
        # Native Apptainer available
        if apptainer_available():
            cmd = ["apptainer", "run"] + apptainer_bind_args + [sif_path] + inner_cmd
        elif singularity_available():
            cmd = ["singularity", "run"] + apptainer_bind_args + [sif_path] + inner_cmd
        else:
            logging.error('Neither apptainer nor singularity are available.')
            sys.exit(1)

    elif system == "Darwin":
        # macOS: run Apptainer inside Docker
        # Before a Docker image that contains Apptainer named "apptainer:latest" must be built
        #
        # FROM ubuntu:22.04
        #
        # # Install Apptainer
        # RUN apt-get update && apt-get install -y \
        #     wget build-essential squashfs-tools uidmap fuse3 git \
        #     && rm -rf /var/lib/apt/lists/*
        #
        # RUN wget https://github.com/apptainer/apptainer/releases/download/v1.3.4/apptainer_1.3.4_amd64.deb \
        #     && dpkg -i apptainer_1.3.4_amd64.deb
        #
        # ENTRYPOINT ["apptainer"]
        docker_image = "apptainer:latest"

        if not docker_available(docker_image):
            logging.error('Either docker or an image apptainer:latest are not available.')
            sys.exit(1)

        # Convert host paths to absolute so Docker can bind them
        docker_binds = []
        for b in binds:
            if isinstance(b, str):
                host = expand_path(b)
            elif isinstance(b, tuple) and len(b) == 2:
                host = expand_path(b[0])
            if 'host' in locals() and host:
                docker_binds += ["-v", f"{host}:{host}"]  # same path inside docker

        # add dir of sif-file to docker bind
        docker_binds += ["-v", f"{os.path.dirname(sif_path)}:{os.path.dirname(sif_path)}"]

        # Command inside Docker:
        inner = ["run", "--no-mount=bind-paths"] + apptainer_bind_args + [sif_path] + inner_cmd

        cmd = [
                  "docker", "run", "--platform", "linux/amd64", "--privileged",
                  # "-v", f"{os.getcwd()}:{os.getcwd()}",
                  # "-w", os.getcwd()
              ] + docker_binds + [
                  docker_image
              ] + inner

    else:
        logging.error(f'Unsupported system ({system}) to run containerized tools')
        sys.exit(1)

    cmd = ' '.join(cmd)
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

    except subprocess.CalledProcessError as e:

        logging.warning(f'(containers) Error running the command {cmd}\n'
                        f'{traceback.format_exc()}\n'
                        f'STDERR:\n'
                        f'{e.stderr}\n'
                        f'STDOUT:\n'
                        f'{e.stdout}\n')


# apptainer_exec('~/imgs/unipka.sif',
#                ['protonate', '-i', '/home/pavel/python/easydock/test/apptainer/molgpka_test.smi', '-o', '/home/pavel/python/easydock/test/apptainer/molgpka_test_prot.smi'],
#                ['/home/pavel/python/easydock/test/apptainer'])

# apptainer_exec('~/imgs/unipka.sif',
#                # ['protonate', '-i', '/tmp/2.smi', '-o', '/tmp/2out.smi'],
#                ['ls', '-l', '/tmp/2out.smi', '>', '/tmp/text'],
#                ['/tmp'])
