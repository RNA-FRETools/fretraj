FROM ubuntu:21.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends pymol pip && \
    python3 -m pip install -U pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip install fretraj && \
    cp /usr/local/lib/python3.9/dist-packages/fretraj/fretraj_gui.py /usr/lib/python3/dist-packages/pmg_tk/startup
CMD pymol
