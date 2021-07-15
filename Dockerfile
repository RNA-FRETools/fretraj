FROM ubuntu:21.04
ENV DEBIAN_FRONTEND noninteractive
ENV package_dir /usr/local/lib/python3.9/dist-packages/fretraj/
RUN apt-get update && \
    apt-get install -y --no-install-recommends pymol pip midori && \
    python3 -m pip install -U pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip install fretraj && \
    cp "$package_dir"/fretraj_gui.py /usr/lib/python3/dist-packages/pmg_tk/startup && \
    echo {\"root_path\": \"/root\", \"browser\": null, \"local_docs\": null} > "$package_dir"/.fretraj_settings.json
CMD pymol
