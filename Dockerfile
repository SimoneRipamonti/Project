FROM ubuntu:latest

LABEL maintainer="simone2.ripamonti@mail.polimi.it"
LABEL version="0.1"
LABEL description="pacs_project"

# Updating Ubuntu packages and install some basic
RUN apt-get update && yes|apt-get upgrade && apt-get install -y wget bzip2 git build-essential vim

# Add a new user named "user"
RUN apt-get -y install sudo

# Add user ubuntu with no password, add to sudo group
RUN adduser --disabled-password --gecos '' user
RUN adduser user sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Add some project related packages
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Rome
RUN apt-get install -y cmake gnuplot libboost-all-dev libeigen3-dev

# Switch to user "user"
USER user
WORKDIR /home/user
RUN chmod a+rwx /home/user/

# install the repo of the project
WORKDIR /home/user
RUN git clone https://github.com/SimoneRipamonti/Project.git project

# darcy case
WORKDIR /home/user/project/case0_example/build
RUN cmake ..; make; ./exe.sh

WORKDIR /home/user/project/case1_darcy/build
RUN cmake ..; make; ./exe.sh

WORKDIR /home/user/project/case2_darcy/build
RUN cmake ..; make; ./exe.sh

WORKDIR /home/user/project/case3_transport/build
RUN cmake ..; make; ./exe.sh

WORKDIR /home/user/project/case4_linear_decay/build
RUN cmake ..; make; ./exe.sh

WORKDIR /home/user/project/case5_2reagents/build
RUN cmake ..; make; ./exe.sh

WORKDIR /home/user/project/case6_all/build
RUN cmake ..; make; ./exe.sh

RUN echo "well done"
