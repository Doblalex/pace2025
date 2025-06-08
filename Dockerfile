FROM debian:latest AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    unzip \
    cmake \
    libboost-all-dev \
    python3 \
    python3-pip \
    python3-venv \
    time \
    ccache \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/lib/ccache:$PATH"


FROM base AS build
WORKDIR /root
RUN git clone https://github.com/Doblalex/pace2025.git
WORKDIR /root/pace2025
RUN git rm instances/official/ || echo "official instances not present"; \
   git submodule update --init --recursive \
RUN mkdir build-release && \
   cd build-release  && \
   cmake .. -DCMAKE_BUILD_TYPE=Release  && \
   make -j $(nproc)
ENV PATH="/root/pace2025/build-release:$PATH"
