#!/bin/bash
MY_PWD=$(pwd)

#
# find what type of machine we are using
#
#	cd to home folder using machine independent command 
# and look for .txt file 
#
cd
if [[ -f desktop.txt ]]; then
	# local desktop machine
	home="/home/prashant"
elif [[ -f greenfield.txt ]]; then
	# greenfield
	home="/home/pjha"
else
	echo "can't verify if machine is local desktop machine or greenfield"
	echo "exiting script"
	exit 
fi

install_dir="$home""/Softwares/local"
packages="$home""/Softwares/Packages"

# check if folders exist
if [[ ! -d "$home" ]]; then
	echo "home folder = " $home " does not exist"
	exit
fi

if [[ ! -d "$install_dir" ]]; then
	echo "install directory = " $install_dir " does not exist"
	exit
fi

if [[ ! -d "$packages" ]]; then
	echo "packages directory = " $packages " does not exist"
	exit
fi

#
#	Cmake
#
cd $packages
echo "	"
echo "*************** CMAKE ********************"
if [ -f cmake*.tar.gz ]; then
	echo "cmake tar file exists"
	tar -zxf cmake*.tar.gz 
	cd cmake*
	./configure --prefix="$install_dir"
	make
	make install
else
	echo "cmake file does not exist"
	echo "skipping cmake installation"
fi

# for cmake also add path
export PATH=$PATH:"$install_dir"/bin

#
#	Boost
#
cd $packages
echo "	"
echo "*************** BOOST ********************"
if [ -f boost*.tar.gz ]; then
	echo "boost file exists"
else
	echo "boost file does not exist"
	echo "skipping boost library installation"
	echo "since other libraries depend on boost, exiting the bash"
	exit
fi
tar -zxf boost*.tar.gz
cd boost*
./bootstrap.sh --prefix="$install_dir"
./b2 install

# after installing boost export path of the library
export PATH=$PATH:"$install_dir"/include
export PATH=$PATH:"$install_dir"/lib

#
#	Gmp
#
cd $packages
echo "	"
echo "*************** GMP ********************"
if [ -f gmp*.tar.lz ]; then
	echo "gmp tar file exists"
	tar --lzip -xf gmp*.tar.lz
	cd gmp*
	./configure --prefix="$install_dir"
	make 
	make check
	make install
else
	echo "gmp tar file does not exist"
fi

#
#	Mpfr 
#
#	this requires gmp.h and hence need to provide the folder where
#	gmp is installed
#
cd $packages
echo "	"
echo "*************** MPFR ********************"
if [ -f mpfr*.tar.gz ]; then
	echo "mpfr tar file exists"
	# check if gmp.h is found
	if find "$install_dir"/include/gmp.h -maxdepth 0 -empty | read v; then 
		echo "gmp.h not found";
		echo "skipping mpfr installation"
	else
		tar -zxf mpfr*.tar.gz
		cd mpfr*
		./configure --prefix="$install_dir" --with-gmp="$install_dir"
		make
		make check
		make install
	fi
else
	echo "mpfr tar file does not exist"
	echo "skippping mpfr installation"
fi

#
#	Cgal
#
cd $packages
echo "	"
echo "*************** CGAL ********************"
if [ -f CGAL*.tar.gz ]; then
	echo "CGAL tar file exists"
	# check if cmake is found
	if [ -f "$install_dir"/bin/cmake ]; then
		echo "cmake found"
		echo "building cgal"
		tar -zxf CGAL*.tar.gz
		cd CGAL*
		cmake -DCMAKE_INSTALL_PREFIX="$install_dir" -DBUILD_SHARED_LIBS=TRUE
		make
		make install
	fi
else
	echo "CGAL tar file not found"
	echo "skipping CGAL installation"
fi
