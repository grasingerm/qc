#!/bin/bash
flag="$1"
destD="$2"
MY_PWD=$(pwd)

#
if [[ ! -d $destD ]]; then
	echo "input.sh : destination directory invalid."
	echo "arg passed are : "
	echo "1. $1"
	echo "2. $2"
	echo "exiting the script"
	touch error_script.txt
	exit
fi

checkDir="$(pwd)""/Input"
if [[ ! -d $checkDir ]]; then
	echo "input.sh : Input directory not found at $testDir"
	echo "exiting the script"
	touch error_script.txt
	exit
fi

cd Input
inputD=$(pwd)

# back to current directory
cd $MY_PWD

echo "flag = $flag"

if [[ "$flag" == "Build" ]]; then

	Source=$inputD"/"$flag
	Dest=$destD"/applications/nano-indentation"

	cd $Dest
	touch a.inp.gz a.dat a.ini a.out
	rm *inp.gz
	rm *dat
	rm *ini
	rm *out

	cp $Source/* .
	# cp $Source/quas* quasi.ini
elif [[ "$flag" == "Ar" ]]; then

	Source=$inputD"/"$flag
	Dest=$destD"/applications/nano-indentation"

	cd $Dest
	touch a.inp.gz a.dat a.ini a.out
	rm *inp.gz
	rm *dat
	rm *ini
	rm *out

	cp $Source/* .
	# cp $Source/quas* quasi.ini
else
	# find the "-" and index
	counter="0"

	size=${#flag}

	index="-1"
	while [ $counter -lt $size ]; do
		nn="${flag:$counter:1}"

		if [[ "$nn" == "-" ]]; then
			index=$counter
		fi

		counter=$[$counter+1]
	done

	# make sure index is not -1
	if [[ $index == -1 ]]; then
		echo "************error************"
		echo "wrong argument passed to input.sh"
		echo "argument = $flag"
		echo "************error************"
	fi

	nn="${flag:0:$index}"
	echo "copying input files for $nn"

	Source=$inputD"/"$nn
	Dest=$destD"/applications/nano-indentation"

	cd $Dest
	touch a.inp.gz a.dat a.ini a.out
	rm *inp.gz
	rm *dat
	rm *ini
	rm *out

	cp $Source/* .
	# cp $Source/quas* quasi.ini
fi

cp $inputD/materials.dat $destD/applications/nano-indentation/.
