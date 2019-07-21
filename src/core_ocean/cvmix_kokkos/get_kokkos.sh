#!/bin/bash

KOKKOS_GIT_HTTP_ADDRESS=https://github.com/kokkos/kokkos.git
KOKKOS_GIT_SSH_ADDRESS=git@github.com:kokkos/kokkos.git
KOKKOS_TAG=f4920c3

GIT=`which git`
PROTOCOL=""

if [ -d kokkos ]; then
  cd kokkos
  CURR_TAG=$(git rev-parse --short HEAD)
  if [ "${CURR_TAG}" == "${KOKKOS_TAG}" ]; then
    echo "kokkos version is current.  SKipping update"
  else
    rm -rf kokkos
  fi
else
  if [ "${GIT}" != "" ]; then
    echo " ** Using git to acquire kokkos. ** "
    PROTOCOL="git ssh"
    git clone ${KOKKOS_GIT_SSH_ADDRESS} &> /dev/null
    if [ -d kokkos ]; then
      cd kokkos
      git checkout ${KOKKOS_TAG} &> /dev/null
    fi
  fi
fi


