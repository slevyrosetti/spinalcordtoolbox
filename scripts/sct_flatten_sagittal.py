#!/usr/bin/env python
#########################################################################################
#
# Flatten spinal cord in sagittal plane.
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Benjamin De Leener, Julien Cohen-Adad
# Modified: 2014-06-02
#
# About the license: see the file LICENSE.TXT
#########################################################################################

<<<<<<< HEAD
# TODO: check input param with -s flag

## Default parameters
class param:
    ## The constructor
    def __init__(self):
        self.debug = 0
        self.interp = 'sinc' # final interpolation
        self.deg_poly = 10 # maximum degree of polynomial function for fitting centerline.
        self.remove_temp_files = 1 # remove temporary files

# check if needed Python libraries are already installed or not
import os
import getopt
import commands
import sys
import sct_utils as sct
from sct_nurbs import NURBS
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
sct.check_if_installed('WarpImageMultiTransform -h','ANTS')

#=======================================================================================================================
# main
#=======================================================================================================================
def main():
    
    # Initialization
    fname_anat = ''
    fname_centerline = ''
    centerline_fitting = 'splines'
    remove_temp_files = param.remove_temp_files
    interp = param.interp
    degree_poly = param.deg_poly
    
    # extract path of the script
    path_script = os.path.dirname(__file__)+'/'
    
    # Parameters for debug mode
    if param.debug == 1:
        print '\n*** WARNING: DEBUG MODE ON ***\n'
        fname_anat = path_script+'../testing/sct_straighten_spinalcord/data/errsm_22_t2_cropped_rpi.nii.gz'
        fname_centerline = path_script+'../testing/sct_straighten_spinalcord/data/errsm_22_t2_cropped_centerline.nii.gz'
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
    
    # Check input param
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hi:c:r:d:f:s:')
    except getopt.GetoptError as err:
        print str(err)
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ('-i'):
            fname_anat = arg
        elif opt in ('-c'):
            fname_centerline = arg
        elif opt in ('-r'):
            remove_temp_files = int(arg)
        elif opt in ('-d'):
            degree_poly = int(arg)
        elif opt in ('-f'):
            centerline_fitting = str(arg)
        elif opt in ('-s'):
            interp = str(arg)
    
    # display usage if a mandatory argument is not provided
    if fname_anat == '' or fname_centerline == '':
        usage()
    
    # check existence of input files
    sct.check_file_exist(fname_anat)
    sct.check_file_exist(fname_centerline)
    
    # extract path/file/extension
    path_anat, file_anat, ext_anat = sct.extract_fname(fname_anat)
    
    # Display arguments
    print '\nCheck input arguments...'
    print '  Input volume ...................... '+fname_anat
    print '  Centerline ........................ '+fname_centerline
    print ''
    
    # Get input image orientation
    status, output = sct.run('sct_orientation -i ' + fname_anat + ' -get')
    input_image_orientation = output[-3:]

    # Reorient input data into RL PA IS orientation
    sct.run('sct_orientation -i '+fname_anat+' -o tmp.anat_orient.nii -orientation RPI')
    sct.run('sct_orientation -i '+fname_centerline+' -o tmp.centerline_orient.nii -orientation RPI')

    # Open centerline
    #==========================================================================================
    print '\nGet dimensions of input centerline...'
    nx, ny, nz, nt, px, py, pz, pt = sct.get_dimension('tmp.centerline_orient.nii')
    print '.. matrix size: '+str(nx)+' x '+str(ny)+' x '+str(nz)
    print '.. voxel size:  '+str(px)+'mm x '+str(py)+'mm x '+str(pz)+'mm'
    
    print '\nOpen centerline volume...'
    file = nibabel.load('tmp.centerline_orient.nii')
    data = file.get_data()

    X, Y, Z = (data>0).nonzero()
    min_z_index, max_z_index = min(Z), max(Z)
    
    
    # loop across z and associate x,y coordinate with the point having maximum intensity
    x_centerline = [0 for iz in range(min_z_index, max_z_index+1, 1)]
    y_centerline = [0 for iz in range(min_z_index, max_z_index+1, 1)]
    z_centerline = [iz for iz in range(min_z_index, max_z_index+1, 1)]

    # Two possible scenario:
    # 1. the centerline is probabilistic: each slices contains voxels with the probability of containing the centerline [0:...:1]
    # We only take the maximum value of the image to aproximate the centerline.
    # 2. The centerline/segmentation image contains many pixels per slice with values {0,1}.
    # We take all the points and approximate the centerline on all these points.

    X, Y, Z = ((data<1)*(data>0)).nonzero() # X is empty if binary image
    if (len(X) > 0): # Scenario 1
        for iz in range(min_z_index, max_z_index+1, 1):
            x_centerline[iz-min_z_index], y_centerline[iz-min_z_index] = numpy.unravel_index(data[:,:,iz].argmax(), data[:,:,iz].shape)
    else: # Scenario 2
        for iz in range(min_z_index, max_z_index+1, 1):
            x_seg, y_seg = (data[:,:,iz]>0).nonzero()
            if len(x_seg) > 0:
                x_centerline[iz-min_z_index] = numpy.mean(x_seg)
                y_centerline[iz-min_z_index] = numpy.mean(y_seg)

    # TODO: find a way to do the previous loop with this, which is more neat:
    # [numpy.unravel_index(data[:,:,iz].argmax(), data[:,:,iz].shape) for iz in range(0,nz,1)]
    
    # clear variable
    del data
    
    # Fit the centerline points with the kind of curve given as argument of the script and return the new smoothed coordinates
    if centerline_fitting == 'splines':
        x_centerline_fit, y_centerline_fit = b_spline_centerline(x_centerline,y_centerline,z_centerline)
    elif centerline_fitting == 'polynome':
        x_centerline_fit, y_centerline_fit = polynome_centerline(x_centerline,y_centerline,z_centerline)

    #==========================================================================================
    # Split input volume
    print '\nSplit input volume...'
    sct.run(sct.fsloutput + 'fslsplit tmp.anat_orient.nii tmp.anat_z -z')
    file_anat_split = ['tmp.anat_z'+str(z).zfill(4) for z in range(0,nz,1)]

    # initialize variables
    file_mat_inv_cumul = ['tmp.mat_inv_cumul_z'+str(z).zfill(4) for z in range(0,nz,1)]
    z_init = min_z_index
    displacement_max_z_index = x_centerline_fit[z_init-min_z_index]-x_centerline_fit[max_z_index-min_z_index]

    # write centerline as text file
    print '\nGenerate fitted transformation matrices...'
    file_mat_inv_cumul_fit = ['tmp.mat_inv_cumul_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    for iz in range(min_z_index, max_z_index+1, 1):
        # compute inverse cumulative fitted transformation matrix
        fid = open(file_mat_inv_cumul_fit[iz], 'w')
        if (x_centerline[iz-min_z_index] == 0 and y_centerline[iz-min_z_index] == 0):
            displacement = 0
        else:
            displacement = x_centerline_fit[z_init-min_z_index]-x_centerline_fit[iz-min_z_index]
        fid.write('%i %i %i %f\n' %(1, 0, 0, displacement) )
        fid.write('%i %i %i %f\n' %(0, 1, 0, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 1, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 0, 1) )
        fid.close()

    # we complete the displacement matrix in z direction
    for iz in range(0, min_z_index, 1):
        fid = open(file_mat_inv_cumul_fit[iz], 'w')
        fid.write('%i %i %i %f\n' %(1, 0, 0, 0) )
        fid.write('%i %i %i %f\n' %(0, 1, 0, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 1, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 0, 1) )
        fid.close()
    for iz in range(max_z_index+1, nz, 1):
        fid = open(file_mat_inv_cumul_fit[iz], 'w')
        fid.write('%i %i %i %f\n' %(1, 0, 0, displacement_max_z_index) )
        fid.write('%i %i %i %f\n' %(0, 1, 0, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 1, 0) )
        fid.write('%i %i %i %i\n' %(0, 0, 0, 1) )
        fid.close()

    # apply transformations to data
    print '\nApply fitted transformation matrices...'
    file_anat_split_fit = ['tmp.anat_orient_fit_z'+str(z).zfill(4) for z in range(0,nz,1)]
    for iz in range(0, nz, 1):
        # forward cumulative transformation to data
        sct.run(fsloutput+'flirt -in '+file_anat_split[iz]+' -ref '+file_anat_split[iz]+' -applyxfm -init '+file_mat_inv_cumul_fit[iz]+' -out '+file_anat_split_fit[iz]+' -interp '+interp)

    # Merge into 4D volume
    print '\nMerge into 4D volume...'
    sct.run(fsloutput+'fslmerge -z tmp.anat_orient_fit tmp.anat_orient_fit_z*')

    # Reorient data as it was before
    print '\nReorient data back into native orientation...'
    sct.run('sct_orientation -i tmp.anat_orient_fit.nii -o tmp.anat_orient_fit_reorient.nii -orientation '+input_image_orientation)

    # Generate output file (in current folder)
    print '\nGenerate output file (in current folder)...'
    sct.generate_output_file('tmp.anat_orient_fit_reorient.nii','./',file_anat+'_flatten',ext_anat)

    # Delete temporary files
    if remove_temp_files == 1:
        print '\nDelete temporary files...'
        sct.run('rm tmp.*')

    # to view results
    print '\nDone! To view results, type:'
    print 'fslview '+file_anat+ext_anat+' '+file_anat+'_flatten'+ext_anat+' &\n'


#=======================================================================================================================
# usage
#=======================================================================================================================
def usage():
    print 'USAGE: \n' \
        ''+os.path.basename(__file__)+'\n' \
        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' \
        'Part of the Spinal Cord Toolbox <https://sourceforge.net/projects/spinalcordtoolbox>\n' \
        '\n'\
        'DESCRIPTION\n' \
        '  Flatten the spinal cord in the sagittal plane (to make nice pictures).\n' \
        '\n' \
        'USAGE\n' \
        '  '+os.path.basename(__file__)+' -i <source> -c <centerline>\n' \
        '\n' \
        'MANDATORY ARGUMENTS\n' \
        '  -i <source>       input volume.\n' \
        '  -c                centerline.\n' \
        '\n'\
        'OPTIONAL ARGUMENTS\n' \
        '  -s {nearestneighbour, trilinear, sinc}       final interpolation. Default='+str(param.interp)+'\n' \
        '  -d <deg>          degree of fitting polynome. Default='+str(param.deg_poly)+'\n' \
        '  -r {0, 1}         remove temporary files. Default='+str(param.remove_temp_files)+'\n' \
        '  -h                help. Show this message.\n' \
        '\n'\
        'EXAMPLE:\n' \
        '  sct_flatten_sagittal.py -i t2.nii.gz -c centerline.nii.gz\n'
    sys.exit(2)

def b_spline_centerline(x_centerline,y_centerline,z_centerline):
    """Give a better fitting of the centerline than the method 'spline_centerline' using b-splines"""
    
    
    points = [[x_centerline[n],y_centerline[n],z_centerline[n]] for n in range(len(x_centerline))]
    
    nurbs = NURBS(3,len(z_centerline)*3,points) # for the third argument (number of points), give at least len(z_centerline)
    # (len(z_centerline)+500 or 1000 is ok)
    P = nurbs.getCourbe3D()
    x_centerline_fit=P[0]
    y_centerline_fit=P[1]
    
    return x_centerline_fit, y_centerline_fit



def polynome_centerline(x_centerline,y_centerline,z_centerline):
    """Fit polynomial function through centerline"""
    
    # Fit centerline in the Z-X plane using polynomial function
    print '\nFit centerline in the Z-X plane using polynomial function...'
    coeffsx = np.polyfit(z_centerline, x_centerline, deg=5)
    polyx = np.poly1d(coeffsx)
    x_centerline_fit = np.polyval(polyx, z_centerline)
    
    #Fit centerline in the Z-Y plane using polynomial function
    print '\nFit centerline in the Z-Y plane using polynomial function...'
    coeffsy = np.polyfit(z_centerline, y_centerline, deg=5)
    polyy = np.poly1d(coeffsy)
    y_centerline_fit = np.polyval(polyy, z_centerline)
    
    
    return x_centerline_fit,y_centerline_fit

#=======================================================================================================================
# Start program
#=======================================================================================================================
if __name__ == "__main__":
    # initialize parameters
    param = param()
    # call main function
    main()
=======
import sys
import os
import argparse

import numpy as np
from skimage import transform, img_as_float

from spinalcordtoolbox.image import Image, add_suffix, change_type
from spinalcordtoolbox.centerline.core import ParamCenterline, get_centerline
from spinalcordtoolbox.utils import Metavar, SmartFormatter, init_sct, display_viewer_syntax


# Default parameters
class Param:
    # The constructor
    def __init__(self):
        self.debug = 0
        self.interp = 'sinc'  # final interpolation
        self.remove_temp_files = 1  # remove temporary files
        self.verbose = 1


def flatten_sagittal(im_anat, im_centerline, verbose):
    """
    Flatten a 3D volume using the segmentation, such that the spinal cord is centered in the R-L medial plane.
    :param im_anat:
    :param im_centerline:
    :param verbose:
    :return:
    """
    # re-oriente to RPI
    orientation_native = im_anat.orientation
    im_anat.change_orientation("RPI")
    im_centerline.change_orientation("RPI")
    nx, ny, nz, nt, px, py, pz, pt = im_anat.dim

    # smooth centerline and return fitted coordinates in voxel space
    _, arr_ctl, _, _ = get_centerline(im_centerline, param=ParamCenterline(), verbose=verbose)
    x_centerline_fit, y_centerline_fit, z_centerline = arr_ctl

    # Extend the centerline by copying values below zmin and above zmax to avoid discontinuities
    zmin, zmax = z_centerline.min().astype(int), z_centerline.max().astype(int)
    x_centerline_extended = np.concatenate([np.ones(zmin) * x_centerline_fit[0],
                                            x_centerline_fit,
                                            np.ones(nz - zmax) * x_centerline_fit[-1]])

    # change type to float32 and scale between -1 and 1 as requested by img_as_float(). See #1790, #2069
    im_anat_flattened = change_type(im_anat, np.float32)
    min_data, max_data = np.min(im_anat_flattened.data), np.max(im_anat_flattened.data)
    im_anat_flattened.data = 2 * im_anat_flattened.data / (max_data - min_data) - 1

    # loop and translate each axial slice, such that the flattened centerline is centered in the medial plane (R-L)
    for iz in range(nz):
        # compute translation along x (R-L)
        translation_x = x_centerline_extended[iz] - np.round(nx / 2.0)
        # apply transformation to 2D image with linear interpolation
        # tform = tf.SimilarityTransform(scale=1, rotation=0, translation=(translation_x, 0))
        tform = transform.SimilarityTransform(translation=(0, translation_x))
        # important to force input in float to skikit image, because it will output float values
        img = img_as_float(im_anat_flattened.data[:, :, iz])
        img_reg = transform.warp(img, tform)
        im_anat_flattened.data[:, :, iz] = img_reg

    # change back to native orientation
    im_anat_flattened.change_orientation(orientation_native)

    return im_anat_flattened


def main(fname_anat, fname_centerline, verbose):
    """
    Main function
    :param fname_anat:
    :param fname_centerline:
    :param verbose:
    :return:
    """
    # load input images
    im_anat = Image(fname_anat)
    im_centerline = Image(fname_centerline)

    # flatten sagittal
    im_anat_flattened = flatten_sagittal(im_anat, im_centerline, verbose)

    # save output
    fname_out = add_suffix(fname_anat, '_flatten')
    im_anat_flattened.save(fname_out)

    display_viewer_syntax([fname_anat, fname_out])


def get_parser():
    parser = argparse.ArgumentParser(
        description="Flatten the spinal cord such within the medial sagittal plane. Useful to make nice pictures. "
                    "Output data has suffix _flatten. Output type is float32 (regardless of input type) to minimize "
                    "loss of precision during conversion.",
        formatter_class=SmartFormatter,
        add_help=None,
        prog=os.path.basename(__file__).strip(".py")
    )
    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        '-i',
        metavar=Metavar.file,
        required=True,
        help="Input volume. Example: t2.nii.gz"
    )
    mandatory.add_argument(
        '-s',
        metavar=Metavar.file,
        required=True,
        help="Spinal cord segmentation or centerline. Example: t2_seg.nii.gz"
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )
    optional.add_argument(
        '-v',
        choices=['0', '1', '2'],
        default=str(param_default.verbose),
        help="Verbosity. 0: no verbose (default), 1: min verbose, 2: verbose + figures"
    )

    return parser


if __name__ == "__main__":
    init_sct()
    # initialize parameters
    param = Param()
    param_default = Param()

    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    fname_anat = arguments.i
    fname_centerline = arguments.s
    verbose = int(arguments.v)
    init_sct(log_level=verbose, update=True)  # Update log level

    # call main function
    main(fname_anat, fname_centerline, verbose)
>>>>>>> upstream/master
