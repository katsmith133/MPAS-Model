#!/bin/bash
if [ -d cvmix_kokkos ]; then
	cd cvmix_kokkos
	./get_kokkos.sh
  echo " getting kokkos"
else
	echo " ****************************************************** "
	echo " ERROR: Build failed to acquire CVMix-Modernization source."
	echo ""
	echo " Please ensure your proxy information is setup properly for"
	echo " the protocol you use to acquire CVMix Modernization."
	echo ""
	echo " The automated script attempted to use: ${PROTOCOL}"
	echo ""
	if [ "${PROTOCOL}" == "git http" ]; then
		echo " This protocol requires setting up the http.proxy git config option."
	elif [ "${PROTOCOL}" == "git ssh" ]; then
		echo " This protocol requires having ssh-keys setup, and ssh access to git@github.com."
		echo " Please use 'ssh -vT git@github.com' to debug issues with ssh keys."
	elif [ "${PROTOCOL}" == "svn" ]; then
		echo " This protocol requires having svn proxys setup properly in ~/.subversion/servers."
	fi
	echo ""
	echo " ****************************************************** "
fi
