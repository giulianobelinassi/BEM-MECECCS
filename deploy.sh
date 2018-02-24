#!/bin/bash
# Deploy the BEM-MECECCS application.

# Color special caracters.
GREEN='\033[1;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

MAX_THREADS=4

function print_ok {
	echo -e "[${GREEN}OK${NC}]\n"
}

function print_failed {
	echo -e "[${RED}FALHOU${NC}]\n"
}

#Install OpenBLAS. version 0.20.0
function compile_openblas {
	rm -rf OpenBLAS-0.2.20
	wget http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz
	tar xf v0.2.20.tar.gz
	cd OpenBLAS-0.2.20

# USE_OPENMP=1 for parallel lapack, INTERFACE64=1 for 64-bits matrix indexing
	USE_OPENMP=1 INTERFACE64=1 make -j ${MAX_THREADS} || true

# Check if at least the static linking library was compiled.
	if [ ! -f libopenblas.a ]; then
		echo -e "[${RED}OpenBLAS Compilation Failed${NC}]\n"
		cd ..
		exit 1
	else
		echo -e "[${GREEN}OpenBLAS Succefully compiled${NC}]\n"
		cd ..
		return
	fi
}
#install libmagma, version 2.3.0
function compile_magma {
	rm -rf magma-2.3.0
	wget http://icl.utk.edu/projectsfiles/magma/downloads/magma-2.3.0.tar.gz
	tar xf magma-2.3.0.tar.gz

# Copy my custom configuration files. Also workarround bug #4
	cp -r magma-make.inc magma-2.3.0/make.inc
	cd magma-2.3.0
	make dense || true

#Check if at least the static linking library was compiled
	if [ ! -f lib/libmagma.a ]; then
		echo -e "[${RED}LibMAGMA Compilation Failed${NC}]\n"
		cd ..
		exit 1	
	else
		echo -e "[${GREEN}LibMAGMA Succefully compiled${NC}]\n"
		cd ..
		return
	fi
}

function compile_mececcs {
	CURR_DIR=$PWD
	cd src/

	echo -e "[${YELLOW}Running Tests...${NC}]\n"

	./tests.sh
	make clean

#extract meshes used in tests
	unzip -o ../malhas.zip
}

if [ ! -f OpenBLAS-0.2.20/libopenblas.a ]; then
	compile_openblas
else
	echo -e "[${GREEN}OpenBLAS already compiled${NC}]\n"
fi

if [ ! -f magma-2.3.0/lib/libmagma.a ]; then
	compile_magma
else
	echo -e "[${GREEN}libMAGMA already compiled${NC}]\n"
fi

compile_mececcs
