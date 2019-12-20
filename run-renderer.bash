#!/usr/bin/env bash

ROOM="room_1"
INPATH="/media/eleanor/New\ Volume/replica-dataset-data/${ROOM}"
MESH="${INPATH}/mesh.ply"
TEXTURE="${INPATH}/textures/"
GLASS="${INPATH}/glass.sur"
CAMPATH="./data/camera_paths/room1-camera1.txt"
DYNMESH="./data/lpshead/head.obj"
CONFIG="./config.json"


eval ./build/ReplicaSDK/ReplicaRenderer ${MESH} ${TEXTURE} ${GLASS} ${CAMPATH} ${DYNMESH} ${CONFIG}
