# lightweight python base image
# slim contains bash 3.0 - needed for nextflow wrapper
FROM python:3.13-slim-trixie

# show python logs immediately
ENV PYTHONUNBUFFERED=1

# utilities required for nextflow, clean up indexes to keep docker image small
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

# python package dependencies
COPY ../py-requirements.txt /tmp/py-requirements.txt
RUN pip install \
    --no-cache-dir \
    #--require-hashes \
    -r tmp/py-requirements.txt