#!/bin/sh
repo=${1}
docker build --rm=true -t ${repo}/ngseasy-ubuntu-base:v1.0 .
