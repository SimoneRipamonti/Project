#!/bin/bash

docker build -t pacs_project .
CONTAINER_ID=$(docker container create pacs_project)

PROJECT=/home/user/project
EXAMPLES=(case0_example case1_darcy case2_darcy case3_transport case4_linear_decay case5_2reagents case6_all)

# export the solutions
for ex in "${EXAMPLES[@]}"; do
    rm -rf ${ex}
    mkdir -p ${ex}
    docker cp ${CONTAINER_ID}:${PROJECT}/${ex}/build/solutions/ .
    mv solutions ${ex}
done
