#!/bin/bash

# This script generates arch specific files (symbolic links)
# Use ./clean_arch.sh <arch> to delete these files

usage()
{
    echo "usage: ./select_arch.sh <arch>"
    echo -n " <arch>: "
    echo $(ls make_* | cut -d "_" -f2- | grep -v arch_inc)
}

# 
if [ $# == 0 ]
then
    usage
    exit 1
fi    

#
if [ -e make_$1 ]
then
    arch_cible=$1
else
    echo "ERREUR : argument invalide."
    echo
    usage
    exit 1
fi

echo "Vous avez choisi l'architecture :" $arch_cible
echo

./clean_arch.sh

# Traitement du repertoire arch
ln -s make_$arch_cible make_arch_inc

cd .. # root dir

# Traitement des repertoires tp1|...|tp9
for dir in `ls -d tp*`
do
    (
	echo "Traitement du repertoire $dir"
	cd $dir
	
	if [ ! -L parametres ]; then
	    (cd parametres; ln -s make_tp_$arch_cible make_tp_inc)
	fi

	for rep in avec_indications_openmp sans_indications_openmp solution
	do
	    if [ -d $rep ]; then
		(
		    cd $rep

		    if [[ $dir =~ tp.*_c ]]; then
			ln -s ../../arch/MakefileC Makefile
		    else
			ln -s ../../arch/Makefile Makefile
		    fi

		    if [ -e ../parametres/batch_$arch_cible.sh ]; then 
			ln -s ../parametres/batch_$arch_cible.sh batch.sh
		    fi

		    ln -s ../parametres/run_$arch_cible.sh run.sh

		    if [ $rep == "solution" ]
		    then
			if [ -e ../parametres/batch_solutions_$arch_cible.sh ]; then 
			    ln -s ../parametres/batch_solutions_$arch_cible.sh batch_solutions.sh
			fi
			
			if [ -e ../parametres/run_solutions_$arch_cible.sh ]; then 
			    ln -s ../parametres/run_solutions_$arch_cible.sh run_solutions.sh
			fi
		    fi

		)
	    fi
	done
    )
done

# Traitement specifique tp6 (JMFFT)
cd tp6/JMFFT-8.1
(cd arch; ln -s Make.inc.$arch_cible Make.inc)

echo
