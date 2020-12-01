#!/bin/bash
#
# Build and install the binaries of the spinal cord toolbox
#
# This script will compile (using CMake) the binaries of the spinal cord toolbox for the OS it is called.
# We assume here that this script is always in "install" folder of the toolbox

SCT_DIR_LOCAL=${PWD%/*}
echo ${SCT_DIR_LOCAL}

# program list
<<<<<<< HEAD
PROGRAM_LIST="sct_change_nifti_pixel_type sct_crop_image sct_detect_spinalcord sct_dice_coefficient sct_hausdorff_distance sct_modif_header sct_orientation sct_segmentation_propagation"
=======
PROGRAM_LIST="isct_vesselness isct_propseg sct_change_nifti_pixel_type sct_detect_spinalcord isct_dice_coefficient sct_hausdorff_distance sct_modif_header isct_bsplineapproximator" #
>>>>>>> upstream/master

PATH_BIN_SCT=osx
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
<<<<<<< HEAD
  PATH_BIN_SCT=debian
=======
  PATH_BIN_SCT=linux
>>>>>>> upstream/master
fi

# loop across programs
for program in $PROGRAM_LIST; do
  
  echo
  cmd="cd ${SCT_DIR_LOCAL}/dev/${program}/"
  echo ">> $cmd"; $cmd

  echo
<<<<<<< HEAD
  cmd="cmake ."
  echo ">> $cmd"; $cmd

  echo
  cmd="make"
  echo ">> $cmd"; $cmd

  echo
  cmd="cp ${program} ${SCT_DIR_LOCAL}/bin/${PATH_BIN_SCT}/"
  echo ">> $cmd"; $cmd

  echo
  cmd="make"
=======
  cmd="rm ._*"
  echo ">> $cmd"; $cmd

  echo
  cmd="mkdir build"
  echo ">> $cmd"; $cmd

  echo
  cmd="cd build"
  echo ">> $cmd"; $cmd

  echo
  cmd="cmake .."
  echo ">> $cmd"; $cmd

  echo
  cmd="make -j 4"
  echo ">> $cmd"; $cmd

  echo
  cmd="cp ${program} ${SCT_DIR_LOCAL}/bin/"
>>>>>>> upstream/master
  echo ">> $cmd"; $cmd

  echo
  cmd="make clean"
  echo ">> $cmd"; $cmd

<<<<<<< HEAD
=======
  echo
  cmd="cd .."
  echo ">> $cmd"; $cmd

  echo
  cmd="rm -rf build"
  echo ">> $cmd"; $cmd


>>>>>>> upstream/master
done
