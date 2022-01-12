#!/bin/bash

CONTAINER_ID=$1
PROJECT=/home/user/project
EXAMPLE=(case5_2reagents)
# export the solutions

rm -rf ./solutions
docker cp ${CONTAINER_ID}:${PROJECT}/${EXAMPLE}/build/solutions/ .


