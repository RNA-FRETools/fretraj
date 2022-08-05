FROM jupyter/base-notebook
USER root
COPY . /home/jovyan/fretraj
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gcc g++ && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /home/jovyan/fretraj
RUN pip install -e .
WORKDIR /home/jovyan
USER 1000
