FROM ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive
ENV package_dir /usr/local/lib/python3.10/dist-packages/fretraj/
RUN apt-get update && \
    apt-get install -y --no-install-recommends python3-dev pymol pip luakit gcc g++ && \
    python3 -m pip install -U pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip install fretraj && \
    cp "$package_dir"/fretraj_gui.py /usr/lib/python3/dist-packages/pmg_tk/startup && \
    echo {\"root_path\": \"/root\", \"browser\": null, \"local_docs\": null} > "$package_dir"/.fretraj_settings.json
CMD pymol
