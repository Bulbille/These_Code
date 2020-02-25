#!/bin/bash

# This script removes arch specific files (symbolic links)
# Use ./select_arch.sh <arch> to re-generate these files

usage()
{
    echo "usage: ./clean_arch.sh"
}

# Check current directory (must be './')
if [ `dirname $0` != "." ]
then
    usage
    exit 1
fi

cd .. # root dir

# Clean arch/ directory
(
    rm -f arch/make_arch_inc
)

# Clean tp*/ directories
(
    rm -f tp*/parametres/make_tp_inc
    rm -f tp*/sans*/Makefile tp*/avec*/Makefile tp*/sol*/Makefile # do not delete JMFFT Makefile
    rm -f tp*/*/run.sh 
    rm -f tp*/*/batch.sh
    rm -f tp*/*/run_solutions.sh 
    rm -f tp*/*/batch_solutions.sh
)

# Remove JMFFT specific files (tp6)
(
    rm -f tp6/JMFFT-8.1/arch/Make.inc
)
