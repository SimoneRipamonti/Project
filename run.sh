#!/bin/bash

#run doxygen
cd include; doxygen ../dconfig; cd ../
cd case6_all/src; doxygen ../../dconfig; cd ../../

#run docker file
docker build -t pacs_project .


## export the solutions
CONTAINER_ID=$(docker container create pacs_project)
PROJECT=/home/user/project
EXAMPLES=(case0_example case1_darcy case2_darcy case3_transport case4_linear_decay case5_2reagents case6_all)

for ex in "${EXAMPLES[@]}"; do
    #rm -rf ${ex}_solutions
    mkdir -p ${ex}_solutions
    docker cp ${CONTAINER_ID}:${PROJECT}/${ex}/build/solutions/ .
    mv solutions ${ex}_solutions
done


