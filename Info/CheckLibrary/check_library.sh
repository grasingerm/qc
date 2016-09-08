#!/bin/bash
MY_PWD=$(pwd)
dateDay=`date +%Y-%m-%d`

# if cgal and boost code can be compiled, cgal_success and boost_success
# files will be created.
#
# error or compiler output would be in error_cgal and error_boost

# if flag supplied is wrong, invalid_flag.txt will be created.

##################################
#	running in desktop or cluster	 #
##################################
# cd to home folder, irrespective of cluster or desktop
#
#	create desktop.txt at home folder in Desktop. If it is greenfield,
# create greenfield.txt. If it is any other clutser, add corresponding
# file and also add the else-if conditio here.
#
cd
if [[ -f desktop.txt ]]; then
	echo "        running on desktop machine"
	# set machine flag to 0
	machine="0"
	home="/home/prashant"
elif [[ -f greenfield.txt ]]; then
	echo "        running on greenfield machine"
	# set machine flag to 1
	machine="1"
	home="/home/pjha"
else
	echo "        can't find if it's desktop or greenfiel machine"
	echo "        if it is neither, modify the bash file to include it"
	echo "        exiting bash script"
	exit
fi

cd $MY_PWD
# remove file which indicates the success of cgal and boost comiplation
if [[ -f cgal_success.txt ]]; then
	rm cgal_success.txt
fi

if [[ -f boost_success.txt ]]; then
	rm boost_success.txt
fi

if [[ -f error_boost.txt ]]; then
	rm error_boost.txt
fi

if [[ -f error_cgal.txt ]]; then
	rm error_cgal.txt
fi

if [[ -f invalid_flag.txt ]]; then
	rm invalid_flag.txt
fi

touch a.out
rm a.out

# input parameters
# libraryFlag : 0 - libraries at /usr/local
#               1 - libraries locally built at ~/Softwares/local
libraryFlag="$1"

if [[ $machine == "1" ]]; then
  # greenfield, so put 1 to libraryFlag
  let "libraryFlag = 1"
fi

# check if there is no input argument
if [[ -z $libraryFlag ]]; then
  echo "        invalid or no argument supplied to libraryFlag"
  echo "        exiting script"
  touch invalid_flag.txt
  exit
fi

#
# set library directorios
#
if [[ $libraryFlag == "1" ]]; then
  localD="$home""/Softwares/local"
  includeD="$localD""/include"
  libD="$localD""/lib"
  cgalIncludeD="$includeD""/CGAL"
  cgalLibD="$libD""/CGAL"
else
  localD="/usr/local"
  includeD="$localD""/include"
  libD="$localD""/lib"
  cgalIncludeD="$includeD""/CGAL"
  cgalLibD="$libD""/CGAL"
fi

# export paths
export PATH=$PATH:"$localD"
export PATH=$PATH:"$localD""/bin"
export PATH=$PATH:"$localD""/include"
export PATH=$PATH:"$cgalIncludeD"
export PATH=$PATH:"$localD""/lib"
export LD_LIBRARY_PATH=$includeD:$libD

#
# check CGAL library
#
echo "        *** cgal check"

# remove object files
touch a.out
rm a.out

# run code
# echo "libraryFlag = $libraryFlag"
# echo "ld = $LD_LIBRARY_PATH"
g++ --std=c++11 cgal_check.cc -I"$includeD" -L"$libD" -lCGAL > error_cgal.txt

# check if it created a.out
if [[ ! -f a.out ]]; then
  echo "        *** CGAL code did not compile"
  echo "        *** check error_cgal.out in directory = $(pwd)"
else
  echo "        *** CGAL code compiled"
  touch cgal_success.txt
fi

#
# check BOOST library
#
echo "        *** boost check"

# remove object files
touch a.out
rm a.out

# run code
g++ --std=c++11 boost_check.cc -I"$includeD" -L"$libD" > error_boost.txt

# check if it created a.out
if [[ ! -f a.out ]]; then
  echo "        *** boost code did not compile"
  echo "        *** check error_boost.out in directory = $(pwd)"
else
  echo "        *** BOOST code compiled"
  touch boost_success.txt
fi

touch a.out
rm a.out
# done
