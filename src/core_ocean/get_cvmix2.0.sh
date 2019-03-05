#!/bin/bash

## CVMix Tag for build
CVMIX_TAG=v1.1
## Subdirectory in CVMix repo to use
CVMIX_SUBDIR=lib

## Available protocols for acquiring CVMix source code
CVMIX_GIT_HTTP_ADDRESS=https://gitlab.com/vanroekel/cvmix_modernization.git
CVMIX_GIT_SSH_ADDRESS=git@gitlab.com:CVMix/cvmix_modernization.git
CVMIX_SVN_ADDRESS=https://gitlab.com/vanroekel/cvmix_modernization/tags

GIT=`which git`
SVN=`which svn`
PROTOCOL=""

# CVMix exists. Check to see if it is the correct version.
# Otherwise, flush the directory to ensure it's updated.
if [ -d cvmix2 ]; then

	if [ -d .cvmix2_all/.git ]; then
		cd .cvmix2_all
		CURR_TAG=$(git describe --tags)
		cd ../
		if [ "${CURR_TAG}" == "${CVMIX_TAG}" ]; then
			echo "CVmix2 version is current. Skip update"
		else
			unlink cvmix2
			rm -rf .cvmix2_all
		fi
	else
		unlink cvmix2
		rm -rf .cvmix2_all
	fi
fi


# CVmix Doesn't exist, need to acquire souce code
# If might have been flushed from the above if, in the case where it was svn or wget that acquired the source.
if [ ! -d cvmix2 ]; then 
	if [ -d .cvmix2_all ]; then
		rm -rf .cvmix2_all
	fi

	if [ "${GIT}" != "" ]; then 
		echo " ** Using git to acquire cvmix2.0 source. ** "
		PROTOCOL="git ssh"
		git clone ${CVMIX_GIT_SSH_ADDRESS} .cvmix2_all &> /dev/null
		if [ -d .cvmix2_all ]; then 
			cd .cvmix2_all 
			git checkout ${CVMIX_TAG} &> /dev/null
			cd ../ 
			ln -sf .cvmix2_all/${CVMIX_SUBDIR} cvmix 
		else 
			git clone ${CVMIX_GIT_HTTP_ADDRESS} .cvmix2_all &> /dev/null
			PROTOCOL="git http"
			if [ -d .cvmix2_all ]; then 
				cd .cvmix2_all 
				git checkout ${CVMIX_TAG} &> /dev/null
				cd ../ 
				ln -sf .cvmix2_all/${CVMIX_SUBDIR} cvmix2 
			fi 
		fi 
	else 
		echo " ** Using svn to acquire cvmix source. ** "
		PROTOCOL="svn"
		svn co ${CVMIX_SVN_ADDRESS}/${CVMIX_TAG} .cvmix2_all &> /dev/null
		ln -sf .cvmix2_all/${CVMIX_SUBDIR} cvmix2
	fi 
fi

if [ -d cvmix2 ]; then
	cd cvmix2
	./get_kokkos.sh
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
