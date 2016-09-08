#!/bin/bash +x
MY_PWD=$(pwd)
dateDay=`date +%Y-%m-%d`

################################################
#    read flags from quasi_script_parameters   #
################################################
#
#	runFlag : 0 - run quasi
#						 1 - run test
#						 other - do nothing
#
#	makeClean : 0 - no (default)
#							1 - yes
#
#	config : 0 - don't configure (default)
#					 1 - configure
#
#	libraryFlag : 0 - assuming that global library is installed in
#										 /usr/local
#								1 - assuming that global library is installed in
#										 /home/prashant/Softwares/local
#
#	testName : testName of directory where the QC code will be build.
#            <model>-<comment>
#						e.g.: NiAl-Electrostatics - NiAl is the model
#                 NiMn-NoElectrostatics - NiMn is the model
#
#						<model> - See Input folder for list of supported models
#           <comment> - it can be anything 
#
# *** testName also controls from where the input files will be copied.
#
# debugMode :  0 - no (default)
#							 1 - yes
#	*** debug_mode outputs the various folder names, exported paths,
#     value of flags used in this bash script, and proceeds only
#     if user types "yes" or "no".
#
# quasiBuildFlag : 0 (default)
#                    create quasiBuild at default location which is
#                    $(pwd)/quasiBuild
#                  1
#                    create quasiBuild at location specifed by
#                    quasiBuildHere
#
#	quasiRunOutFlag : 0 (default)
#											Outputs the make run data and quasi run data
#											to the shell
#										1
#											Outputs the make run data and quasi run data
#											to the corresponding .out files.
#
#

# read parameters
cd $MY_PWD
echo "******************************************"
echo "reading file : quasi_script_parameters.txt"

readParaCounter="0"
while IFS=$'\n' read line
do
	if [[ ! "$line" =~ \#.*  ]]; then
		# increment the counter
		readParaCounter=$[$readParaCounter+1]

		# reading uncommented files
		if [[ $readParaCounter == "1" ]]; then
			#let "runFlag = $line"
			runFlag="$line"
		elif [[ $readParaCounter == "2" ]]; then
			#let "makeClean = $line"
			makeClean="$line"
		elif [[ $readParaCounter == "3" ]]; then
			#let "config = $line"
			config="$line"
		elif [[ $readParaCounter == "4" ]]; then
			#let "libraryFlag = $line"
			libraryFlag="$line"
		elif [[ $readParaCounter == "5" ]]; then
			#let "libraryFlag = $line"
			libraryDir="$line"
		elif [[ $readParaCounter == "6" ]]; then
			#let "testName = $line"
			testName="$line"
		elif [[ $readParaCounter == "7" ]]; then
			#let "debugMode = $line"
			debugMode="$line"
		elif [[ $readParaCounter == "8" ]]; then
			#let "quasiBuildFlag = $line"
			quasiBuildFlag="$line"
		elif [[ $readParaCounter == "9" ]]; then
			#let "quasiBuildHere = $line"
			quasiBuildHere="$line"
		elif [[ $readParaCounter == "10" ]]; then
			#let "quasiRunOutFlag = $line"
			quasiRunOutFlag="$line"			
		fi
	fi
done < quasi_script_parameters.txt

echo "*************** done *********************"

# check if we have got all the flags
if [[ -z "$runFlag"  ]]; then
	echo "runFlag is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$makeClean"  ]]; then
	echo "makeClean is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$config"  ]]; then
	echo "config is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$libraryFlag"  ]]; then
	echo "libraryFlag is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$libraryDir"  ]]; then
	echo "libraryDir is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$testName"  ]]; then
	echo "testName is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$debugMode"  ]]; then
	echo "debugMode is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$quasiBuildFlag"  ]]; then
	echo "quasiBuildFlag is missing from quasi_script_parameters.txt"
	exit
fi

if [[ -z "$quasiBuildHere"  ]]; then
	echo "quasiBuildHere is missing from quasi_script_parameters.txt"
	exit
fi

# output flags if debugMode is 1
if [[ $debugMode == "1" ]]; then
	echo "      Debug Data : Flags"
	echo "        runFlag = $runFlag"
	echo "        makeClean = $makeClean"
	echo "        config = $config"
	echo "        libraryFlag = $libraryFlag"
	echo "        libraryDir = $libraryDir"
	echo "        testName = $testName"
	echo "        debugMode = $debugMode"
	echo "        quasiBuildFlag = $quasiBuildFlag"
	echo "        quasiBuildHere = $quasiBuildHere"
fi

##################################
#    QuasiBuild Folder     #
##################################
#
# if quasiBuildFlag = 0
#		 create QuasiBuild at pwd, if it does not exist
#    create testName folder at pwd, if it does not exist
# else
#    create QuasiBuild at specifed location, if it does not
#      exist there
#    create testName at specified folder, if it does not
#      exist there
#
cd $MY_PWD
if [[ $quasiBuildFlag == "0" ]]; then
	# check if it exist at default location, if not create it
	checkDir="$(pwd)/QuasiBuild"
	if [[ ! -d "$checkDir" ]]; then
		# directory does not exist. Creating one.
		mkdir QuasiBuild
	fi

	# check again
	if [[ ! -d "$checkDir" ]]; then
		# directory still does not exist. report error.
		echo "check runQuasi.sh script. QuasiBuild should be here $(pwd)"
		exit
	fi

	# create subdir inside QuasiBuild
	cd QuasiBuild
	quasiBuild=$(pwd)

	checkDir="$(pwd)/$testName"
	if [[ ! -d "$checkDir" ]]; then
		# directory does not exist. Creating one.
		mkdir $testName
	fi

	# check again
	if [[ ! -d "$checkDir" ]]; then
		# directory still does not exist. Report error.
		echo "check runQuasi.sh script. $testName should be here $(pwd)"
		exit
	fi

	cd $testName
	testDir=$(pwd)
elif [[ $quasiBuildFlag == "1" ]]; then
  # first check if $quasiBuildHere is valid directory
	if [[ ! -d "$quasiBuildHere" ]]; then
		# not valid directory
		echo "quasiBuildHere = $quasiBuildHere is not valid directory."
		echo "check quasi_script_parameters.txt"
		exit
	fi

	# go to quasiBuildHere
	cd $quasiBuildHere

	# check if it exist at specifed location, if not create it
	checkDir="$(pwd))/QuasiBuild"
	if [[ ! -d "$checkDir" ]]; then
		# directory does not exist. Creating one.
		mkdir QuasiBuild
	fi

	# check again
	if [[ ! -d "$checkDir" ]]; then
		# directory still does not exist. report error.
		echo "check runQuasi.sh script. QuasiBuild should be here $(pwd)"
		exit
	fi

	# create subdir inside QuasiBuild
	cd QuasiBuild
	quasiBuild=$(pwd)

	checkDir="$(pwd)/testName"
	if [[ ! -d "$checkDir" ]]; then
		# directory does not exist. Creating one.
		mkdir testName
	fi

	# check again
	if [[ ! -d "$checkDir" ]]; then
		# directory still does not exist. Report error.
		echo "check runQuasi.sh script. $testName should be here $(pwd)"
		exit
	fi

	cd $testName
	testDir=$(pwd)
else
	echo "quasiBuildFlag = $quasiBuildFlag is incorrect. check quasi_script_parameters.txt"
fi

#
# set QuasiSRC folder
#
cd $MY_PWD
checkDir="$(pwd)/QuasiSRC"

# check if folder exists
if [[ ! -d "$checkDir" ]]; then
	# directory still does not exist. Report error.
	echo "check runQuasi.sh script. QuasiSRC should be here $(pwd)"
	exit
fi

cd QuasiSRC
quasiSRC=$(pwd)

# get configure file and check if it exists in QuasiSRC
quasiConfigure=$quasiSRC"/configure"
if [[ ! -f "$quasiConfigure" ]] ; then
    echo "can not find configure=$quasiConfigure file in QuasiSRC=$quasiSRC, exiting bash"
    exit
fi

# output Data if debugMode is 1
if [[ $debugMode == "1" ]]; then
	echo "      Debug Data : Directories"
	echo "        quasiBuildFlag = $quasiBuildFlag"
	echo "        quasiBuild = $quasiBuild"
	echo "        testDir = $testDir"
	echo "        quasiSRC = $quasiSRC"
	echo "        quasiConfigure = $quasiConfigure"
fi

##################################
#    Build specific 					   #
##################################
cd $MY_PWD
if [[ $libraryFlag == "0" ]]; then
	# directory information
	localD="/usr/local"
	includeD="$localD""/include"
	libD="$localD""/lib"
	cgalIncludeD="$includeD""/CGAL"
	cgalLibD="$libD""/CGAL"
	directory_for_dependencies="-I""$includeD"
	directory_for_libraries="-L""$libD"

	# check if directories exist
	if [[ ! -d "$localD" ]]; then
		echo "local directory not in specifed location"
		echo "Check local directory = $localD in script file"
		exit
	fi

	if [[ ! -d "$includeD" ]]; then
		echo "local/inlcude directory not in specifed location"
		echo "Check local/inlcude directory = $includeD in script file"
		exit
	fi

	if [[ ! -d "$libD" ]]; then
		echo "local/lib directory not in specifed location"
		echo "Check local/lib directory = $libD in script file"
		exit
	fi

	if [[ ! -d "$cgalIncludeD" ]]; then
		echo "local/include/CGAL directory not in specifed location"
		echo "Check local/include/CGAL directory = $cgalIncludeD in script file"
		exit
	fi

	if [[ ! -d "$cgalLibD" ]]; then
		echo "local/lib/CGAL directory not in specifed location"
		echo "Check local/lib/CGAL directory = $cgalLibD in script file"
		exit
	fi

	# export paths for this flag
	export PATH=$PATH:"$localD"
	export PATH=$PATH:"$localD""/bin"
	export PATH=$PATH:"$localD""/include"
	export PATH=$PATH:"$cgalIncludeD"
	export PATH=$PATH:"$localD""/lib"
	export LD_LIBRARY_PATH="$localD""/include":"$localD""/lib"

	# output Data if debugMode is 1
	if [[ $debugMode == "1" ]]; then
		echo "      Debug Data : Build specific"
		echo "        library flag = $libraryFlag --> 0-/usr/local, 1-locally built"
		echo "        local dir = $localD"
		echo "        bin dir = $localD/bin"
		echo "        include dir = $includeD"
		echo "        lib dir = $libD"
		echo "        include CGAL dir = $cgalIncludeD"
		echo "        lib CGAL dir = $cgalLibD"
		echo "        *** path data ****"
		echo "        paths = $PATH"
		echo "        *** LD_LIBRARY_PATH data ****"
		echo "        LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
	fi
	# done
else
	# directory information
	localD="$libraryDir"
	includeD="$localD""/include"
	libD="$localD""/lib"
	cgalIncludeD="$includeD""/CGAL"
	cgalLibD="$libD""/CGAL"
	directory_for_dependencies="-I""$includeD"
	directory_for_libraries="-L""$libD"

	# check if directories exist
	if [[ ! -d "$localD" ]]; then
		echo "local directory not in specifed location"
		echo "Check local directory = $localD in script file"
		exit
	fi

	if [[ ! -d "$includeD" ]]; then
		echo "local/inlcude directory not in specifed location"
		echo "Check local/inlcude directory = $includeD in script file"
		exit
	fi

	if [[ ! -d "$libD" ]]; then
		echo "local/lib directory not in specifed location"
		echo "Check local/lib directory = $libD in script file"
		exit
	fi

	if [[ ! -d "$cgalIncludeD" ]]; then
		echo "local/include/CGAL directory not in specifed location"
		echo "Check local/include/CGAL directory = $cgalIncludeD in script file"
		exit
	fi

	if [[ ! -d "$cgalLibD" ]]; then
		echo "local/lib/CGAL directory not in specifed location"
		echo "Check local/lib/CGAL directory = $cgalLibD in script file"
		exit
	fi

	# export paths for this flag
	export PATH=$PATH:"$localD"
	export PATH=$PATH:"$localD""/bin"
	export PATH=$PATH:"$localD""/include"
	export PATH=$PATH:"$cgalIncludeD"
	export PATH=$PATH:"$localD""/lib"
	export LD_LIBRARY_PATH="$localD""/include":"$localD""/lib"

	# output Data if debugMode is 1
	if [[ $debugMode == "1" ]]; then
		echo "      Debug Data : Build specific"
		echo "        library flag = $libraryFlag --> 0-/usr/local, 1-locally built"
		echo "        local dir = $localD"
		echo "        bin dir = $localD/bin"
		echo "        include dir = $includeD"
		echo "        lib dir = $libD"
		echo "        include CGAL dir = $cgalIncludeD"
		echo "        lib CGAL dir = $cgalLibD"
		echo "        *** path data ****"
		echo "        paths = $PATH"
		echo "        *** LD_LIBRARY_PATH data ****"
		echo "        LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
	fi
	# done
fi

##################################
#    Compiler specific 					   #
##################################
cd $MY_PWD
Fortran_Compiler="gfortran"
C_Compiler="gcc"
CXX_Compiler="g++"
C_FLAGS="-O2 -g"
CXX_FLAGS="-O2 -g"
CPP_FLAGS="-std=c++11 -frounding-math ""$directory_for_dependencies"
LD_FLAGS="$directory_for_libraries"
numThreads="1"

#F77=gfortran CC=gcc CXX=g++ CFLAGS=-O2 CXXFLAGS=-O2 CPPFLAGS=-std=c++11 -frounding-math -I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure > $configOUT 2>&1


##################################
#    Quasi code Output 					   #
##################################
cd $MY_PWD

checkDir="$testDir""/OUTPUT"
if [[ ! -d $checkDir ]]; then
	cd $testDir
	mkdir OUTPUT
	cd $MY_PWD
fi

cd $testDir
cd OUTPUT
outDir=$(pwd)

# output file names
configOUT=$outDir"/configure.out"
makeOUT=$outDir"/make.out"
makeOUT2=$outDir"/main-make.out"
outputOUT=$outDir"/results.out"
testOUT=$outDir"/test.out"

# output Data if debugMode is 1
if [[ $debugMode == "1" ]]; then
	echo "      Debug Data : Output files"
	echo "        OUTPUT directory = $outDir"
	echo "        configure output = $configOUT"
	echo "        make output = $makeOUT"
	echo "        make main output = $makeOUT2"
	echo "        quasi output = $outputOUT"
	echo "        test output = $testOUT"
fi

##################################
#    Configure if needed  		   #
##################################
cd $MY_PWD

# #
# # 	generate config file if needed - autoreconf
# #
# cd $QuasiSRC
# autoreconf -fvi -I $QuasiSRC/system
# cd $MY_PWD

#
# configure if needed
#
if [[ $config == 1 ]]; then
	echo "********************** configure **********************"
	cd $testDir
	F77=$Fortran_Compiler CC=$C_Compiler CXX=$CXX_Compiler CFLAGS=$C_FLAGS CXXFLAGS=$CXX_FLAGS CPPFLAGS=$CPP_FLAGS LDFLAGS=$LD_FLAGS $quasiConfigure #> $configOUT 2>&1
	echo "************************* done ************************"
fi

#
# make in main directory
#
cd $MY_PWD
echo "************************ make *************************"
cd $testDir
if [[ $makeClean == 1 ]]; then
	make clean
fi
if [[ $quasiRunOutFlag == "0" ]]; then
	make
else
	make > $makeOUT 2>&1
fi
echo "************************* done ************************"

#
# make in nano-indentation directory
#
cd $MY_PWD
echo "************************ make *************************"
cd $testDir"/applications/nano-indentation/"
if [[ $makeClean == 1 ]]; then
	make clean
fi
if [[ $quasiRunOutFlag == "0" ]]; then
	make
else
	make > $makeOUT2 2>&1
fi
echo "************************* done ************************"

#
#	call another script file to copy the input file
# to corresponding build folder
#
cd $MY_PWD
if [[ $runFlag == 0 || $runFlag == 1 ]]; then
	echo "****************** copy input files *******************"
	cd $MY_PWD
	./input.sh $testName $testDir
	cd $MY_PWD
	echo "************************* done ************************"
fi

#
#	run the problem
#
cd $testDir"/applications/nano-indentation/"
if [[ $runFlag == 0 ]]; then
	#
	# run quasi on input data
	#
	echo "*******************************************************"
	echo "                     running quasi               "
	echo "*******************************************************"
	if [[ $quasiRunOutFlag == "0" ]]; then
		echo "quasi"
		./quasi
	else
		./quasi > $outputOUT 2>&1
	fi	
	echo "************************* done ************************"

	# ./quasi -r restart_00000.gz -n $numThreads > $outputOUT 2>&1
elif [[ $runFlag == 1 ]]; then
	#
	# run test on input data
	#
	echo "*******************************************************"
	echo "                     running test               "
	echo "*******************************************************"
	if [[ $quasiRunOutFlag == "0" ]]; then
		./test
	else
		./test > $testOUT 2>&1
	fi	
	echo "************************* done ************************"
fi

cd $MY_PWD
