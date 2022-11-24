FROM        ubuntu:18.04
MAINTAINER  Joon Klaps <joon.klaps@kuleuven.be>

# update and install dependencies
RUN         apt-get update \
                && apt-get install -y \
                    wget \
                    python \
                && apt-get update \
                && wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz \
                && tar -xzf SPAdes-3.15.5-Linux.tar.gz \
                && cd SPAdes-3.15.5-Linux/bin/

