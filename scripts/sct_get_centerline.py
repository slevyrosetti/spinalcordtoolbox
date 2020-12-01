#!/usr/bin/env python

<<<<<<< HEAD
## @package sct_get_centerline
#
# - run a gaussian-weighted slice-by-slice registration to roughtly estimate the spinal cord centerline.
# - fit the rought centerline to 3D spline function
#
#
# Description about how the function works:
#
# 1. slice-wise realignement
# ------------------------------------------------
# inputs:
# - spinal cord image
# - binary image with a point center of spinal cord.
#
# the algorithm iterates in the upwards and then downward direction (along Z). For each direction, it registers the slice i+1 on the slice i. It does that using a gaussian mask, centers on the spinal cord, and is applied on both the i and i+1 slices (inweight and refweight from flirt).
# NB: the code works on APRLIS
# output:
# - rough centerline
# - gaussian mask
# - transformation matrices Tx,Ty
#
# 2. smoothing of the centerline
# ------------------------------------------------
# input: centerline, transfo matrices
# it fits in 3d (using an independent decomposition of the planes XZ and YZ) a spline function of order 1.
# apply the smoothed transformation matrices to the image --> gives a "smoothly" straighted spinal cord (but stil with errors due to deformation of strucutre related to curvatyre)
# output:
# - centerline fitted
# - transformation matrices Tx,Ty smoothed
# - straighted spinal cord
#
#
# DEPENDENCIES
# ---------------------------------------------------------------------------------------
# EXTERNAL PYTHON PACKAGES
# - nibabel: <http://nipy.sourceforge.net/nibabel/>
# - numpy: <http://www.numpy.org>
#
# EXTERNAL SOFTWARE
# - FSL: <http://fsl.fmrib.ox.ac.uk/fsl/>
#
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 NeuroPoly, Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Julien Cohen-Adad, Geoffrey Leveque
#
# License: see the LICENSE.TXT
# ==========================================================================================

# TODO: apply previous transfo to mask (assuming spatial autocorrelation --> transformation is expected to be close to previous one)
# TODO: output data in short instead of float
# TODO: try to fit in 3D instead of 2D
# TODO: add option verbose for sct.run()
# TODO: include the schedule file in the python script to remove the dependence
# TODO: try applying N4 before
# TODO: 3d smooth of the input image (e.g., 3mm gaussian kernel) -> later
# TODO: subsample input data for faster processing (but might create problems of coordinate)-- alternatively, find appropriate schedule file without the 1mm at the end




## Default parameters
class param:
    ## The constructor
    def __init__(self):
        self.debug = 0
        self.schedule_file = 'schedule_TxTy.sch'
        self.gap = 1 # default gap between co-registered slices.
        self.gaussian_kernel = 4 # gaussian kernel for creating gaussian mask from center point.
        self.deg_poly = 10 # maximum degree of polynomial function for fitting centerline.
        self.remove_temp_files = 1 # remove temporary files

# check if needed Python libraries are already installed or not
import os
import getopt
import sys
import time
import sct_utils as sct
from sct_utils import fsloutput
try:
    import nibabel
except ImportError:
    print '--- nibabel not installed! Exit program. ---'
    sys.exit(2)
try:
    import numpy
except ImportError:
    print '--- numpy not installed! Exit program. ---'
    sys.exit(2)

# check if dependant software are installed
sct.check_if_installed('flirt -help','FSL')


#=======================================================================================================================
# main
#=======================================================================================================================
def main():

    # Initialization
    fname_anat = ''
    fname_point = ''
    slice_gap = param.gap
    remove_temp_files = param.remove_temp_files
    gaussian_kernel = param.gaussian_kernel
    start_time = time.time()

    # extract path of the script
    path_script = os.path.dirname(__file__)+'/'

    # Parameters for debug mode
    if param.debug == 1:
        print '\n*** WARNING: DEBUG MODE ON ***\n'
        fname_anat = path_script+'../dev/testing/data/errsm_22/t2/errsm_22_t2_cropped.nii.gz'
        fname_point = path_script+'../dev/testing/data/errsm_22/t2/center_cropped.nii.gz'
        slice_gap = 5
        import matplotlib.pyplot as plt

    # Check input param
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hi:p:g:r:k:')
    except getopt.GetoptError as err:
        print str(err)
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ('-i'):
            fname_anat = arg
        elif opt in ('-p'):
            fname_point = arg
        elif opt in ('-g'):
            slice_gap = int(arg)
        elif opt in ('-r'):
            remove_temp_files = int(arg)
        elif opt in ('-k'):
            gaussian_kernel = int(arg)

    # display usage if a mandatory argument is not provided
    if fname_anat == '' or fname_point == '':
        usage()

    # check existence of input files
    sct.check_file_exist(fname_anat)
    sct.check_file_exist(fname_point)

    # extract path/file/extension
    path_anat, file_anat, ext_anat = sct.extract_fname(fname_anat)
    path_point, file_point, ext_point = sct.extract_fname(fname_point)

    # extract path of schedule file
    # TODO: include schedule file in sct
    # TODO: check existence of schedule file
    file_schedule = path_script[0:-8]+'flirtsch/' + param.schedule_file

    # Get input image orientation
    status, output = sct.run('sct_orientation -i ' + fname_anat + ' -get')
    input_image_orientation = output[-3:]

    # Display arguments
    print '\nCheck input arguments...'
    print '  Anatomical image:     '+fname_anat
    print '  .... orientation:     '+input_image_orientation
    print '  Point in spinal cord: '+fname_point
    print '  Slice gap:            '+str(slice_gap)
    print '  Gaussian kernel:      '+str(gaussian_kernel)
    print '  Degree of polynomial: '+str(param.deg_poly)

    # convert to nii
    print '\nCopy input data...'
    sct.run('cp ' + fname_anat + ' tmp.anat'+ext_anat)
    sct.run('fslchfiletype NIFTI tmp.anat')
    sct.run('cp ' + fname_point + ' tmp.point'+ext_point)
    sct.run('fslchfiletype NIFTI tmp.point')

    # Reorient input anatomical volume into RL PA IS orientation
    print '\nReorient input volume to RL PA IS orientation...'
    sct.run(sct.fsloutput + 'fslswapdim tmp.anat RL PA IS tmp.anat_orient')

    # Reorient binary point into RL PA IS orientation
    print '\nReorient binary point into RL PA IS orientation...'
    sct.run(sct.fsloutput + 'fslswapdim tmp.point RL PA IS tmp.point_orient')

    # Get image dimensions
    print '\nGet image dimensions...'
    nx, ny, nz, nt, px, py, pz, pt = sct.get_dimension('tmp.anat_orient')
    print '.. matrix size: '+str(nx)+' x '+str(ny)+' x '+str(nz)
    print '.. voxel size:  '+str(px)+'mm x '+str(py)+'mm x '+str(pz)+'mm'

    # Split input volume
    print '\nSplit input volume...'
    sct.run(sct.fsloutput + 'fslsplit tmp.anat_orient tmp.anat_orient_z -z')
    file_anat_split = ['tmp.anat_orient_z'+str(z).zfill(4) for z in range(0,nz,1)]

    # Get the coordinates of the input point
    print '\nGet the coordinates of the input point...'
    file = nibabel.load('tmp.point_orient.nii')
    data = file.get_data()
    x_init, y_init, z_init = (data > 0).nonzero()
    x_init = x_init[0]
    y_init = y_init[0]
    z_init = z_init[0]
    print '('+str(x_init)+', '+str(y_init)+', '+str(z_init)+')'

    # Extract the slice corresponding to z=z_init
    print '\nExtract the slice corresponding to z='+str(z_init)+'...'
    file_point_split = ['tmp.point_orient_z'+str(z).zfill(4) for z in range(0,nz,1)]
    sct.run(sct.fsloutput+'fslroi tmp.point_orient '+file_point_split[z_init]+' 0 -1 0 -1 '+str(z_init)+' 1')

    # Create gaussian mask from point
    print '\nCreate gaussian mask from point...'
    file_mask_split = ['tmp.mask_orient_z'+str(z).zfill(4) for z in range(0,nz,1)]
    sct.run(sct.fsloutput+'fslmaths '+file_point_split[z_init]+' -s '+str(gaussian_kernel)+' '+file_mask_split[z_init])

    # Obtain max value from mask
    print '\nFind maximum value from mask...'
    file = nibabel.load(file_mask_split[z_init]+'.nii')
    data = file.get_data()
    max_value_mask = numpy.max(data)
    print '..'+str(max_value_mask)

    # Normalize mask beween 0 and 1
    print '\nNormalize mask beween 0 and 1...'
    sct.run(sct.fsloutput+'fslmaths '+file_mask_split[z_init]+' -div '+str(max_value_mask)+' '+file_mask_split[z_init])

    ## Take the square of the mask
    #print '\nCalculate the square of the mask...'
    #sct.run(sct.fsloutput+'fslmaths '+file_mask_split[z_init]+' -mul '+file_mask_split[z_init]+' '+file_mask_split[z_init])

    # initialize variables
    file_mat = ['tmp.mat_z'+str(z).zfill(4) for z in range(0,nz,1)]
    file_mat_inv = ['tmp.mat_inv_z'+str(z).zfill(4) for z in range(0,nz,1)]
    file_mat_inv_cumul = ['tmp.mat_inv_cumul_z'+str(z).zfill(4) for z in range(0,nz,1)]

    # create identity matrix for initial transformation matrix
    fid = open(file_mat_inv_cumul[z_init], 'w')
    fid.write('%i %i %i %i\n' %(1, 0, 0, 0) )
    fid.write('%i %i %i %i\n' %(0, 1, 0, 0) )
    fid.write('%i %i %i %i\n' %(0, 0, 1, 0) )
    fid.write('%i %i %i %i\n' %(0, 0, 0, 1) )
    fid.close()

    # initialize centerline: give value corresponding to initial point
    x_centerline = [x_init]
    y_centerline = [y_init]
    z_centerline = [z_init]
    warning_count = 0

    # go up (1), then down (2) in reference to the binary point
    for iUpDown in range(1, 3):

        if iUpDown == 1:
            # z increases
            slice_gap_signed = slice_gap
        elif iUpDown == 2:
            # z decreases
            slice_gap_signed = -slice_gap
            # reverse centerline (because values will be appended at the end)
            x_centerline.reverse()
            y_centerline.reverse()
            z_centerline.reverse()

        # initialization before looping
        z_dest = z_init # point given by user
        z_src = z_dest + slice_gap_signed

        # continue looping if 0 < z < nz
        while 0 <= z_src and z_src <= nz-1:

            # print current z:
            print 'z='+str(z_src)+':'

            # estimate transformation
            sct.run(fsloutput+'flirt -in '+file_anat_split[z_src]+' -ref '+file_anat_split[z_dest]+' -schedule '+file_schedule+ ' -verbose 0 -omat '+file_mat[z_src]+' -cost normcorr -forcescaling -inweight '+file_mask_split[z_dest]+' -refweight '+file_mask_split[z_dest])

            # display transfo
            status, output = sct.run('cat '+file_mat[z_src])
            print output

            # check if transformation is bigger than 1.5x slice_gap
            tx = float(output.split()[3])
            ty = float(output.split()[7])
            norm_txy = numpy.linalg.norm([tx, ty],ord=2)
            if norm_txy > 1.5*slice_gap:
                print 'WARNING: Transformation is too large --> using previous one.'
                warning_count = warning_count + 1
                # if previous transformation exists, replace current one with previous one
                if os.path.isfile(file_mat[z_dest]):
                    sct.run('cp '+file_mat[z_dest]+' '+file_mat[z_src])

            # estimate inverse transformation matrix
            sct.run('convert_xfm -omat '+file_mat_inv[z_src]+' -inverse '+file_mat[z_src])

            # compute cumulative transformation
            sct.run('convert_xfm -omat '+file_mat_inv_cumul[z_src]+' -concat '+file_mat_inv[z_src]+' '+file_mat_inv_cumul[z_dest])

            # apply inverse cumulative transformation to initial gaussian mask (to put it in src space)
            sct.run(fsloutput+'flirt -in '+file_mask_split[z_init]+' -ref '+file_mask_split[z_init]+' -applyxfm -init '+file_mat_inv_cumul[z_src]+' -out '+file_mask_split[z_src])

            # open inverse cumulative transformation file and generate centerline
            fid = open(file_mat_inv_cumul[z_src])
            mat = fid.read().split()
            x_centerline.append(x_init + float(mat[3]))
            y_centerline.append(y_init + float(mat[7]))
            z_centerline.append(z_src)
            #z_index = z_index+1

            # define new z_dest (target slice) and new z_src (moving slice)
            z_dest = z_dest + slice_gap_signed
            z_src = z_src + slice_gap_signed


    # Reconstruct centerline
    # ====================================================================================================

    # reverse back centerline (because it's been reversed once, so now all values are in the right order)
    x_centerline.reverse()
    y_centerline.reverse()
    z_centerline.reverse()

    # fit centerline in the Z-X plane using polynomial function
    print '\nFit centerline in the Z-X plane using polynomial function...'
    coeffsx = numpy.polyfit(z_centerline, x_centerline, deg=param.deg_poly)
    polyx = numpy.poly1d(coeffsx)
    x_centerline_fit = numpy.polyval(polyx, z_centerline)
    # calculate RMSE
    rmse = numpy.linalg.norm(x_centerline_fit-x_centerline)/numpy.sqrt( len(x_centerline) )
    # calculate max absolute error
    max_abs = numpy.max( numpy.abs(x_centerline_fit-x_centerline) )
    print '.. RMSE (in mm): '+str(rmse*px)
    print '.. Maximum absolute error (in mm): '+str(max_abs*px)

    # fit centerline in the Z-Y plane using polynomial function
    print '\nFit centerline in the Z-Y plane using polynomial function...'
    coeffsy = numpy.polyfit(z_centerline, y_centerline, deg=param.deg_poly)
    polyy = numpy.poly1d(coeffsy)
    y_centerline_fit = numpy.polyval(polyy, z_centerline)
    # calculate RMSE
    rmse = numpy.linalg.norm(y_centerline_fit-y_centerline)/numpy.sqrt( len(y_centerline) )
    # calculate max absolute error
    max_abs = numpy.max( numpy.abs(y_centerline_fit-y_centerline) )
    print '.. RMSE (in mm): '+str(rmse*py)
    print '.. Maximum absolute error (in mm): '+str(max_abs*py)

    # display
    if param.debug == 1:
        plt.figure()
        plt.plot(z_centerline,x_centerline,'.',z_centerline,x_centerline_fit,'r')
        plt.legend(['Data','Polynomial Fit'])
        plt.title('Z-X plane polynomial interpolation')
        plt.show()

        plt.figure()
        plt.plot(z_centerline,y_centerline,'.',z_centerline,y_centerline_fit,'r')
        plt.legend(['Data','Polynomial Fit'])
        plt.title('Z-Y plane polynomial interpolation')
        plt.show()

    # generate full range z-values for centerline
    z_centerline_full = [iz for iz in range(0, nz, 1)]

    # calculate X and Y values for the full centerline
    x_centerline_fit_full = numpy.polyval(polyx, z_centerline_full)
    y_centerline_fit_full = numpy.polyval(polyy, z_centerline_full)

    # Generate fitted transformation matrices and write centerline coordinates in text file
    print '\nGenerate fitted transformation matrices and write centerline coordinates in text file...'
    file_mat_inv_cumul_fit = ['tmp.mat_inv_cumul_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    file_mat_cumul_fit = ['tmp.mat_cumul_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    fid_centerline = open('tmp.centerline_coordinates.txt', 'w')
    for iz in range(0, nz, 1):
        # compute inverse cumulative fitted transformation matrix
        fid = open(file_mat_inv_cumul_fit[iz], 'w')
        fid.write('%i %i %i %f\n' %(1, 0, 0, x_centerline_fit_full[iz]-x_init) )
        fid.write('%i %i %i %f\n' %(0, 1, 0, y_centerline_fit_full[iz]-y_init) )
        fid.write('%i %i %i %i\n' %(0, 0, 1, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 0, 1) )
        fid.close()
        # compute forward cumulative fitted transformation matrix
        sct.run('convert_xfm -omat '+file_mat_cumul_fit[iz]+' -inverse '+file_mat_inv_cumul_fit[iz])
        # write centerline coordinates in x, y, z format
        fid_centerline.write('%f %f %f\n' %(x_centerline_fit_full[iz], y_centerline_fit_full[iz], z_centerline_full[iz]) )
    fid_centerline.close()


    # Prepare output data
    # ====================================================================================================

    # write centerline as text file
    for iz in range(0, nz, 1):
        # compute inverse cumulative fitted transformation matrix
        fid = open(file_mat_inv_cumul_fit[iz], 'w')
        fid.write('%i %i %i %f\n' %(1, 0, 0, x_centerline_fit_full[iz]-x_init) )
        fid.write('%i %i %i %f\n' %(0, 1, 0, y_centerline_fit_full[iz]-y_init) )
        fid.write('%i %i %i %i\n' %(0, 0, 1, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 0, 1) )
        fid.close()

    # write polynomial coefficients
    numpy.savetxt('tmp.centerline_polycoeffs_x.txt',coeffsx)
    numpy.savetxt('tmp.centerline_polycoeffs_y.txt',coeffsy)

    # apply transformations to data
    print '\nApply fitted transformation matrices...'
    file_anat_split_fit = ['tmp.anat_orient_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    file_mask_split_fit = ['tmp.mask_orient_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    file_point_split_fit = ['tmp.point_orient_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    for iz in range(0, nz, 1):
        # forward cumulative transformation to data
        sct.run(fsloutput+'flirt -in '+file_anat_split[iz]+' -ref '+file_anat_split[iz]+' -applyxfm -init '+file_mat_cumul_fit[iz]+' -out '+file_anat_split_fit[iz])
        # inverse cumulative transformation to mask
        sct.run(fsloutput+'flirt -in '+file_mask_split[z_init]+' -ref '+file_mask_split[z_init]+' -applyxfm -init '+file_mat_inv_cumul_fit[iz]+' -out '+file_mask_split_fit[iz])
        # inverse cumulative transformation to point
        sct.run(fsloutput+'flirt -in '+file_point_split[z_init]+' -ref '+file_point_split[z_init]+' -applyxfm -init '+file_mat_inv_cumul_fit[iz]+' -out '+file_point_split_fit[iz]+' -interp nearestneighbour')

    # Merge into 4D volume
    print '\nMerge into 4D volume...'
    sct.run(fsloutput+'fslmerge -z tmp.anat_orient_fit tmp.anat_orient_fit_z*')
    sct.run(fsloutput+'fslmerge -z tmp.mask_orient_fit tmp.mask_orient_fit_z*')
    sct.run(fsloutput+'fslmerge -z tmp.point_orient_fit tmp.point_orient_fit_z*')

    # Copy header geometry from input data
    print '\nCopy header geometry from input data...'
    sct.run(fsloutput+'fslcpgeom tmp.anat_orient.nii tmp.anat_orient_fit.nii ')
    sct.run(fsloutput+'fslcpgeom tmp.anat_orient.nii tmp.mask_orient_fit.nii ')
    sct.run(fsloutput+'fslcpgeom tmp.anat_orient.nii tmp.point_orient_fit.nii ')

    # Generate output file (in current folder)
    print '\nGenerate output file (in current folder)...'
    #sct.generate_output_file('tmp.centerline_polycoeffs_x.txt','./','centerline_polycoeffs_x','.txt')
    #sct.generate_output_file('tmp.centerline_polycoeffs_y.txt','./','centerline_polycoeffs_y','.txt')
    #sct.generate_output_file('tmp.centerline_coordinates.txt','./','centerline_coordinates','.txt')
    sct.generate_output_file('tmp.anat_orient.nii','./',file_anat+'_rpi',ext_anat)
    sct.generate_output_file('tmp.anat_orient_fit.nii','./',file_anat+'_rpi_align',ext_anat)
    sct.generate_output_file('tmp.mask_orient_fit.nii','./',file_anat+'_mask',ext_anat)
    fname_output_centerline = sct.generate_output_file('tmp.point_orient_fit.nii','./',file_anat+'_centerline',ext_anat)
    # Reorient the centerline into the initial orientation of the input image
    print '\nReorient the centerline into the initial orientation of the input image...'
    sct.run('sct_orientation -i '+fname_output_centerline+' -o ' +fname_output_centerline+' -orientation '+input_image_orientation)

    # Delete temporary files
    if remove_temp_files == 1:
        print '\nDelete temporary files...'
        sct.run('rm tmp.*')

    # print number of warnings
    print '\nNumber of warnings: '+str(warning_count)+' (if >10, you should probably reduce the gap and/or increase the kernel size'

    # display elapsed time
    elapsed_time = time.time() - start_time
    print '\nFinished! Elapsed time: '+str(int(round(elapsed_time)))+'s\n'


#=======================================================================================================================
# usage
#=======================================================================================================================
def usage():
    print '\n' \
        'sct_get_centerline\n' \
        '--------------------------------------------------------------------------------------------------------------\n' \
        'Part of the Spinal Cord Toolbox <https://sourceforge.net/projects/spinalcordtoolbox>\n' \
        '\n'\
        'DESCRIPTION\n' \
        '  This program estimates the spinal cord centerline from an anatomical image and from another image pointing\n' \
        '  anywhere in the spinal cord. It works by registering adjacent slices together (2DOF transformation). The\n' \
        '  transformation is constrained by a gaussian mask. After registering all slices, the centerline is approximated\n' \
        '  by a polynomial function to smooth it.\n' \
        '\n'\
        'USAGE\n' \
        '  sct_get_centerline.py -i <anat> -p <point>\n' \
        '\n'\
        'MANDATORY ARGUMENTS\n' \
        '  -i <anat>         anatomic nifti file. Image to extract centerline from.\n' \
        '  -p <point>        binary nifti file. Point in the spinal cord (e.g., created using fslview).\n' \
        '\n'\
        'OPTIONAL ARGUMENTS\n' \
        '  -g <gap>          gap between slices used for registration. Higher is faster but less robust. Default='+str(param.gap)+'.\n' \
        '  -k <kernel>       kernel size for gaussian mask. Higher is more robust but less accurate. Default='+str(param.gaussian_kernel)+'. \n' \
        '  -r <0,1>          remove temporary files. Default='+str(param.remove_temp_files)+'. \n' \
        '  -h                help. Show this message.\n'
    sys.exit(2)


#=======================================================================================================================
# Start program
#=======================================================================================================================
if __name__ == "__main__":
    # initialize parameters
    param = param()
    # call main function
    main()

=======
import os
import sys
import argparse

import numpy as np

from spinalcordtoolbox.image import Image
from spinalcordtoolbox.centerline.core import ParamCenterline, get_centerline, _call_viewer_centerline
from spinalcordtoolbox.reports.qc import generate_qc
from spinalcordtoolbox.utils.shell import Metavar, SmartFormatter, ActionCreateFolder, display_viewer_syntax
from spinalcordtoolbox.utils.sys import init_sct, printv
from spinalcordtoolbox.utils.fs import extract_fname


def get_parser():
    # Initialize the parser
    parser = argparse.ArgumentParser(
        description=(
            "This function extracts the spinal cord centerline. Three methods are available: 'optic' (automatic), "
            "'viewer' (manual), and 'fitseg' (applied on segmented image). These functions output (i) a NIFTI file "
            "with labels corresponding to the discrete centerline, and (ii) a csv file containing the float (more "
            "precise) coordinates of the centerline in the RPI orientation. \n"
            "\n"
            "Reference: C Gros, B De Leener, et al. Automatic spinal cord localization, robust to MRI contrast using "
            "global curve optimization (2017). doi.org/10.1016/j.media.2017.12.001"
        ),
        formatter_class=SmartFormatter,
        add_help=None,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        '-i',
        metavar=Metavar.file,
        required=True,
        help="Input image. Example: t1.nii.gz"
    )

    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )
    optional.add_argument(
        "-c",
        choices=['t1', 't2', 't2s', 'dwi'],
        help="Type of image contrast. Only with method=optic."
    )
    optional.add_argument(
        "-method",
        choices=['optic', 'viewer', 'fitseg'],
        default='optic',
        help="R|Method used for extracting the centerline.\n"
             "  - optic: automatic spinal cord detection method\n"
             "  - viewer: manual selection a few points followed by interpolation\n"
             "  - fitseg: fit a regularized centerline on an already-existing cord segmentation. It will "
             "interpolate if slices are missing and extrapolate beyond the segmentation boundaries (i.e., every "
             "axial slice will exhibit a centerline pixel)."
    )
    optional.add_argument(
        "-centerline-algo",
        choices=['polyfit', 'bspline', 'linear', 'nurbs'],
        default='bspline',
        help="Algorithm for centerline fitting. Only relevant with -method fitseg"
    )
    optional.add_argument(
        "-centerline-smooth",
        metavar=Metavar.int,
        type=int,
        default=30,
        help="Degree of smoothing for centerline fitting. Only for -centerline-algo {bspline, linear}."
    )
    optional.add_argument(
        "-o",
        metavar=Metavar.file,
        help="File name (without extension) for the centerline output files. By default, output file will be the "
             "input with suffix '_centerline'. Example: 'centerline_optic'"
    )
    optional.add_argument(
        "-gap",
        metavar=Metavar.float,
        type=float,
        default=20.0,
        help="Gap in mm between manually selected points. Only with method=viewer."
    )
    optional.add_argument(
        "-igt",
        metavar=Metavar.file,
        help="File name of ground-truth centerline or segmentation (binary nifti)."
    )
    optional.add_argument(
        "-v",
        choices=['0', '1'],
        default='1',
        help="Verbose. 1: display on, 0: display off (default)"
    )
    optional.add_argument(
        "-qc",
        metavar=Metavar.folder,
        action=ActionCreateFolder,
        help="The path where the quality control generated content will be saved."
    )
    optional.add_argument(
        "-qc-dataset",
        metavar=Metavar.str,
        help="If provided, this string will be mentioned in the QC report as the dataset the process was run on."
    )
    optional.add_argument(
        "-qc-subject",
        metavar=Metavar.str,
        help="If provided, this string will be mentioned in the QC report as the subject the process was run on."
    )
    return parser


def run_main():
    init_sct()
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    # Input filename
    fname_input_data = arguments.i
    fname_data = os.path.abspath(fname_input_data)

    # Method used
    method = arguments.method

    # Contrast type
    contrast_type = arguments.c
    if method == 'optic' and not contrast_type:
        # Contrast must be
        error = "ERROR: -c is a mandatory argument when using 'optic' method."
        printv(error, type='error')
        return

    # Gap between slices
    interslice_gap = arguments.gap

    param_centerline = ParamCenterline(
        algo_fitting=arguments.centerline_algo,
        smooth=arguments.centerline_smooth,
        minmax=True)

    # Output folder
    if arguments.o is not None:
        file_output = arguments.o
    else:
        path_data, file_data, ext_data = extract_fname(fname_data)
        file_output = os.path.join(path_data, file_data + '_centerline')

    verbose = int(arguments.v)
    init_sct(log_level=verbose, update=True)  # Update log level

    if method == 'viewer':
        # Manual labeling of cord centerline
        im_labels = _call_viewer_centerline(Image(fname_data), interslice_gap=interslice_gap)
    elif method == 'fitseg':
        im_labels = Image(fname_data)
    elif method == 'optic':
        # Automatic detection of cord centerline
        im_labels = Image(fname_data)
        param_centerline.algo_fitting = 'optic'
        param_centerline.contrast = contrast_type
    else:
        printv("ERROR: The selected method is not available: {}. Please look at the help.".format(method), type='error')
        return

    # Extrapolate and regularize (or detect if optic) cord centerline
    im_centerline, arr_centerline, _, _ = get_centerline(im_labels,
                                                         param=param_centerline,
                                                         verbose=verbose)

    # save centerline as nifti (discrete) and csv (continuous) files
    im_centerline.save(file_output + '.nii.gz')
    np.savetxt(file_output + '.csv', arr_centerline.transpose(), delimiter=",")

    display_viewer_syntax([fname_input_data, file_output + '.nii.gz'], colormaps=['gray', 'red'], opacities=['', '1'])

    path_qc = arguments.qc
    qc_dataset = arguments.qc_dataset
    qc_subject = arguments.qc_subject

    # Generate QC report
    if path_qc is not None:
        generate_qc(fname_input_data, fname_seg=file_output + '.nii.gz', args=sys.argv[1:], path_qc=os.path.abspath(path_qc),
                    dataset=qc_dataset, subject=qc_subject, process='sct_get_centerline')
    display_viewer_syntax([fname_input_data, file_output + '.nii.gz'], colormaps=['gray', 'red'], opacities=['', '0.7'])


if __name__ == '__main__':
    run_main()
>>>>>>> upstream/master
