FROM phusion/baseimage

# apt-get update to fetch available updates
RUN apt-get update -y 
# apt-get dist-upgrade to apply all packages (with smart handling of obsoletes)
# RUN apt-get dist-upgrade -y

RUN apt-get install build-essential cmake zlib1g-dev git -y  --allow-unauthenticated

RUN git clone https://github.com/pmelsted/bifrost.git  && \
    cd bifrost && mkdir build && cd build  && \
    cmake .. && \
    make && \
    make install

RUN export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/ && \
    export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/ && \
    export PATH=$PATH:/usr/local/lib/
