#!/usr/bin/env bash
# Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved

cd 3rdparty/Pangolin
mkdir build
cd build
cmake ..
make -j

cd ../../nlohmann_json
mkdir build
cd build
cmake ..
make -j

cd ../../../
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j

cd ../
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
