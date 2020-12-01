#!/usr/bin/env python
#########################################################################################
# Register a volume (e.g., EPI from fMRI or DTI scan) to an anatomical image.
#
# See Usage() below for more information.
#
<<<<<<< HEAD
#
# DEPENDENCIES
# ---------------------------------------------------------------------------------------
# EXTERNAL PYTHON PACKAGES
# none
#
# EXTERNAL SOFTWARE
# - itksnap/c3d <http://www.itksnap.org/pmwiki/pmwiki.php?n=Main.HomePage>
# - ants <http://stnava.github.io/ANTs/>
#
#
=======
>>>>>>> upstream/master
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Julien Cohen-Adad
# Modified: 2014-06-03
#
# About the license: see the file LICENSE.TXT
#########################################################################################

<<<<<<< HEAD
# TODO: testing script for all cases
# TODO: try to combine seg and image based for 2nd stage
# TODO: output name file for warp using "src" and "dest" file name, i.e. warp_filesrc2filedest.nii.gz
# TODO: set gradient-step-length in mm instead of vox size.

# Note for the developer: DO NOT use --collapse-output-transforms 1, otherise inverse warping field is not output


# DEFAULT PARAMETERS
class param:
    ## The constructor
    def __init__(self):
        self.debug               = 0
        self.remove_temp_files   = 1 # remove temporary files
        self.outSuffix           = "_reg"
        self.padding             = 3 # add 'padding' slices at the top and bottom of the volumes if deformation at the edge is not good. Default=5. Put 0 for no padding.
#        self.convertDeformation  = 0 # Convert deformation field to 4D volume (readable by fslview)
        self.numberIterations    = "15x3" # number of iterations
        self.numberIterationsStep2 = "10" # number of iterations at step 2
        self.verbose             = 0 # verbose
        self.compute_dest2sr     = 0 # compute dest2src warping field
        self.gradientStep = ['0.2', '0.5']  # gradientStep in SyN transformation. First value is for image-based, second is for segmentation-based (if exist)

import sys
import getopt
import os
import commands
import time
import sct_utils as sct

# MAIN
# ==========================================================================================
def main():

    # Initialization
    fname_src = ''
    fname_dest = ''
    fname_src_seg = ''
    fname_dest_seg = ''
    fname_output = ''
    padding = param.padding
    numberIterations = param.numberIterations
    numberIterationsStep2 = param.numberIterationsStep2
    remove_temp_files = param.remove_temp_files
    verbose = param.verbose
    use_segmentation = 0 # use spinal cord segmentation to improve robustness
    fname_init_transfo = ''
    fname_init_transfo_inv = ''
    use_init_transfo = ''
    gradientStep = param.gradientStep
    gradientStep_input = ''
    compute_dest2src = param.compute_dest2sr
    start_time = time.time()
    print ''

    # get path of the toolbox
    status, path_sct = commands.getstatusoutput('echo $SCT_DIR')

    # Parameters for debug mode
    if param.debug:
        print '\n*** WARNING: DEBUG MODE ON ***\n'
        #fname_src = path_sct+'/data/template/MNI-Poly-AMU_T2.nii.gz'
        #fname_src = path_sct+'/testing/data/errsm_23/mt/mtc0.nii.gz'
        #fname_dest = path_sct+'/testing/data/errsm_23/mt/mtc1.nii.gz'
        fname_src = '/Users/julien/data/errsm/errsm_30/t2/template2anat.nii.gz'
        fname_dest = '/Users/julien/data/errsm/errsm_30/t1/t1.nii.gz'
        #fname_src_seg = path_sct+'/data/template/MNI-Poly-AMU_cord.nii.gz'
        #fname_dest_seg = path_sct+'/testing/data/errsm_23/mt/segmentation_binary.nii.gz'
        #fname_init_transfo = path_sct+'/testing/data/errsm_23/template/warp_template2anat.nii.gz'
        #fname_init_transfo_inv = path_sct+'/testing/data/errsm_23/template/warp_anat2template.nii.gz'
        numberIterations = '3x0'
        numberIterationsStep2 = "1"
        gradientStep_input = '0.2'
        compute_dest2src = 1
        verbose = 1

    # Check input parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:g:i:m:n:o:p:q:r:s:t:v:x:y:z:')
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ("-d"):
            fname_dest = arg
        elif opt in ('-g'):
            gradientStep_input = arg
        elif opt in ("-i"):
            fname_src = arg
        elif opt in ("-m"):
            fname_mask = arg
        elif opt in ("-n"):
            numberIterations = arg
        elif opt in ("-o"):
            fname_output = arg
        elif opt in ('-p'):
            padding = arg
        elif opt in ('-q'):
            fname_init_transfo = arg
        elif opt in ('-r'):
            remove_temp_files = int(arg)
        elif opt in ("-s"):
            fname_src_seg = arg
        elif opt in ("-t"):
            fname_dest_seg = arg
        elif opt in ('-v'):
            verbose = int(arg)
        elif opt in ('-x'):
            compute_dest2src = int(arg)
        elif opt in ('-y'):
            numberIterationsStep2 = arg
        elif opt in ('-z'):
            fname_init_transfo_inv = arg

    # display usage if a mandatory argument is not provided
    if fname_src == '' or fname_dest == '':
        print "ERROR: Input file missing. Exit program."
        usage()

    # check segmentation data
    if (fname_src_seg != '' and fname_dest_seg == '') or (fname_src_seg == '' and fname_dest_seg != ''):
        print "ERROR: You need to select a segmentation file for the source AND the destination image. Exit program."
        usage()
    elif fname_src_seg != '' and fname_dest_seg != '':
        use_segmentation = 1

    # Parse gradient step
    print '\nParse gradient step...'
    gradientStep_input = gradientStep_input.replace(' ', '')  # remove spaces
    gradientStep_input = gradientStep_input.split(",")  # parse with comma
    for i in range(len(gradientStep_input)):
        try:
            float(gradientStep_input[i])
            gradientStep[i] = gradientStep_input[i]
        except:
            print '  WARNING: Not a float. Use default value for gradientStep['+str(i)+']'

    # print arguments
    print '\nInput parameters:'
    print '  Source .............. '+fname_src
    print '  Destinationf ........ '+fname_dest
    print '  Segmentation source . '+fname_src_seg
    print '  Segmentation dest ... '+fname_dest_seg
    print '  Init transfo ........ '+fname_init_transfo
    print '  Output name ......... '+fname_output
    print '  Iterations at step1 (seg) .... '+str(numberIterations)
    print '  Iterations at step2 (image) .. '+str(numberIterationsStep2)
    print '  Gradient step ....... '+str(gradientStep)
    print '  Remove temp files ... '+str(remove_temp_files)
    print '  Verbose ............. '+str(verbose)
    #print '.. gradient step:    '+str(gradientStepLength)
    #print '.. metric type:      '+metricType

    # check existence of input files
    print '\nCheck if files exist...'
    sct.check_file_exist(fname_src)
    sct.check_file_exist(fname_dest)
    if use_segmentation:
        sct.check_file_exist(fname_src_seg)
        sct.check_file_exist(fname_dest_seg)

    # get full path
    fname_src = os.path.abspath(fname_src)
    fname_dest = os.path.abspath(fname_dest)
    fname_src_seg = os.path.abspath(fname_src_seg)
    fname_dest_seg = os.path.abspath(fname_dest_seg)
    if not fname_init_transfo == '':
        fname_init_transfo = os.path.abspath(fname_init_transfo)  # test if not empty, otherwise it will transform the empty string into a string with path, which is a problem because the emptiness of the string is tested later.
    if not fname_init_transfo_inv == '':
        fname_init_transfo_inv = os.path.abspath(fname_init_transfo_inv)
    #fname_output = os.path.abspath(fname_output)

    # Extract path, file and extension
    path_src, file_src, ext_src = sct.extract_fname(fname_src)
    path_dest, file_dest, ext_dest = sct.extract_fname(fname_dest)
    if use_segmentation:
        path_src_seg, file_src_seg, ext_src_seg = sct.extract_fname(fname_src_seg)
        path_dest_seg, file_dest_seg, ext_dest_seg = sct.extract_fname(fname_dest_seg)

    # define output folder and file name
    if fname_output == '':
        path_out = ''  # output in user's current directory
        file_out = file_src+"_reg"
        ext_out = ext_src
    else:
        path_out, file_out, ext_out = sct.extract_fname(fname_output)

    # create temporary folder
    print('\nCreate temporary folder...')
    path_tmp = 'tmp.'+time.strftime("%y%m%d%H%M%S")
    status, output = sct.run('mkdir '+path_tmp)

    # copy files to temporary folder
    print('\nCopy files...')
    file_src_tmp = 'src'
    file_dest_tmp = 'dest'
    status, output = sct.run('c3d '+fname_src+' -o '+path_tmp+'/'+file_src_tmp+'.nii')
    status, output = sct.run('c3d '+fname_dest+' -o '+path_tmp+'/'+file_dest_tmp+'.nii')
    if use_segmentation:
        file_src_seg_tmp = 'src_seg'
        file_dest_seg_tmp = 'dest_seg'
        status, output = sct.run('c3d '+fname_src_seg+' -o '+path_tmp+'/'+file_src_seg_tmp+'.nii')
        status, output = sct.run('c3d '+fname_dest_seg+' -o '+path_tmp+'/'+file_dest_seg_tmp+'.nii')

    # go to tmp folder
    os.chdir(path_tmp)

    # Find orientation of source data
    print('\nFind orientation of source data...')
    orientation = sct.get_orientation('src.nii')
    sct.printv('.. '+orientation, verbose)

    # Find the dimension corresponding to z
    sct.printv('\nFind the dimension corresponding to the superior-inferior direction...', verbose)
    dimension_si = 0
    if orientation.find('I') != '-1':
        dimension_si = orientation.find('I')
    elif orientation.find('S') != '-1':
        dimension_si = orientation.find('S')
    else:
        print "ERROR: Cannot find proper dimension. Exit program.\n"
        sys.exit(2)
    sct.printv('.. '+str(dimension_si), verbose)

    # Adjust ANTs variable so that the deformation is restricted in the slice plane
    restrict_deformation = '1x1x1'
    if dimension_si == 0:
        restrict_deformation = '0x1x1'
    elif dimension_si == 1:
        restrict_deformation = '1x0x1'
    elif dimension_si == 2:
        restrict_deformation = '1x1x0'
    else:
        print "ERROR: Cannot adjust variable: restrict_deformation. Exit program.\n"
        sys.exit(2)

    # if use initial transformation (!! needs to be inserted before the --transform field in antsRegistration)
    if fname_init_transfo != '':
        file_src_reg_tmp = file_src_tmp+'_reg'
        # apply initial transformation to moving image, and then estimate transformation between this output and
        # destination image. This approach was chosen instead of inputting the transfo into ANTs, because if the transfo
        # does not bring the image to the same space as the destination image, then warping fields cannot be concatenated at the end.
        print('\nApply initial transformation to moving image...')
        sct.run('WarpImageMultiTransform 3 '+file_src_tmp+'.nii '+file_src_reg_tmp+'.nii -R '+file_dest_tmp+'.nii '+fname_init_transfo+' --use-BSpline')
        file_src_tmp = file_src_reg_tmp
        if use_segmentation:
            file_src_seg_reg_tmp = file_src_seg_tmp+'_reg'
            sct.run('WarpImageMultiTransform 3 '+file_src_seg_tmp+'.nii '+file_src_seg_reg_tmp+'.nii -R '+file_dest_seg_tmp+'.nii '+fname_init_transfo+' --use-BSpline')
            file_src_seg_tmp = file_src_seg_reg_tmp

    # Pad the target and source image (because ants doesn't deform the extremities)
    if padding:
        # Pad source image
        print('\nPad source...')
        pad_image(file_src_tmp, file_src_tmp+'_pad.nii', padding)
        file_src_tmp += '_pad'  # update file name
        # Pad destination image
        print('\nPad destination...')
        pad_image(file_dest_tmp, file_dest_tmp+'_pad.nii', padding)
        file_dest_tmp += '_pad'  # update file name
        if use_segmentation:
            # Pad source image
            print('\nPad source segmentation...')
            file_src_seg_tmp = file_src_seg_tmp
            pad_image(file_src_seg_tmp, file_src_seg_tmp+'_pad.nii', padding)
            file_src_seg_tmp += '_pad'  # update file name
            # Pad destination image
            print('\nPad destination segmentation...')
            pad_image(file_dest_seg_tmp, file_dest_seg_tmp+'_pad.nii', padding)
            file_dest_seg_tmp += '_pad'  # update file name

    # don't use spinal cord segmentation
    if use_segmentation == 0:

        # Estimate transformation using ANTS
        print('\nEstimate transformation using ANTS (might take a couple of minutes)...')

        cmd = 'antsRegistration \
--dimensionality 3 \
'+use_init_transfo+' \
--transform SyN['+str(gradientStep[0])+',3,0] \
--metric MI['+file_dest_tmp+'.nii,'+file_src_tmp+'.nii,1,32] \
--convergence '+numberIterations+' \
--shrink-factors 2x1 \
--smoothing-sigmas 0x0mm \
--Restrict-Deformation '+restrict_deformation+' \
--output [reg,'+file_src_tmp+'_reg.nii] \
--collapse-output-transforms 1 \
--interpolation BSpline[3] \
--winsorize-image-intensities [0.005,0.995]'

        status, output = sct.run(cmd)
        if verbose:
            print output

    # use spinal cord segmentation
    elif use_segmentation == 1:

        # Estimate transformation using ANTS
        print('\nStep #1: Estimate transformation using spinal cord segmentations...')

        cmd = 'antsRegistration \
--dimensionality 3 \
--transform SyN['+str(gradientStep[1])+',3,0] \
--metric MI['+file_dest_seg_tmp+'.nii,'+file_src_seg_tmp+'.nii,1,32] \
--convergence '+numberIterations+' \
--shrink-factors 4x1 \
--smoothing-sigmas 1x1mm \
--Restrict-Deformation '+restrict_deformation+' \
--output [regSeg,regSeg.nii]'

        status, output = sct.run(cmd)
        if verbose:
            print output

        print('\nStep #2: Improve local deformation using images (start from previous transformation)...')

        cmd = 'antsRegistration \
--dimensionality 3 \
--initial-moving-transform regSeg0Warp.nii.gz \
--transform SyN['+str(gradientStep[0])+',3,0] \
--metric MI['+file_dest_tmp+'.nii,'+file_src_tmp+'.nii,1,32] \
--convergence '+numberIterationsStep2+' \
--shrink-factors 1 \
--smoothing-sigmas 0mm \
--Restrict-Deformation '+restrict_deformation+' \
--output [reg,'+file_src_tmp+'_reg.nii] \
--collapse-output-transforms 0 \
--interpolation BSpline[3]'

        status, output = sct.run(cmd)
        if verbose:
            print output

    # Concatenate transformations
    print('\nConcatenate transformations...')

    # if user has initial transfo:
    if fname_init_transfo != '':
        if use_segmentation == 0:
            # src --> dest
            cmd1 = 'ComposeMultiTransform 3 warp_src2dest.nii.gz -R dest.nii reg0Warp.nii.gz '+fname_init_transfo
            # dest --> src
            if compute_dest2src:
                cmd2 = 'ComposeMultiTransform 3 warp_dest2src.nii.gz -R src.nii '+fname_init_transfo_inv+' reg0InverseWarp.nii.gz'

        elif use_segmentation == 1:
            # src --> dest
            cmd1 = 'ComposeMultiTransform 3 warp_src2dest.nii.gz -R dest.nii reg1Warp.nii.gz regSeg0Warp.nii.gz '+fname_init_transfo
            # dest --> src
            if compute_dest2src:
                cmd2 = 'ComposeMultiTransform 3 warp_dest2src.nii.gz -R src.nii '+fname_init_transfo_inv+' regSeg0InverseWarp.nii.gz reg1InverseWarp.nii.gz'

    # if user does not have initial transfo:
    else:
        if use_segmentation == 0:
            # src --> dest
            cmd1 = 'ComposeMultiTransform 3 warp_src2dest.nii.gz -R dest.nii reg0Warp.nii.gz'
            # dest --> src
            if compute_dest2src:
                cmd2 = 'ComposeMultiTransform 3 warp_dest2src.nii.gz -R src.nii reg0InverseWarp.nii.gz'

        elif use_segmentation == 1:
            # src --> dest
            cmd1 = 'ComposeMultiTransform 3 warp_src2dest.nii.gz -R dest.nii reg1Warp.nii.gz regSeg0Warp.nii.gz'
            # dest --> src
            if compute_dest2src:
                cmd2 = 'ComposeMultiTransform 3 warp_dest2src.nii.gz -R src.nii regSeg0InverseWarp.nii.gz reg1InverseWarp.nii.gz'

    print('>> ' + cmd1)
    commands.getstatusoutput(cmd1)  # here cannot use sct.run() because of wrong output status in ComposeMultiTransform
    if compute_dest2src:
        print('>> ' + cmd2)
        commands.getstatusoutput(cmd2)  # here cannot use sct.run() because of wrong output status in ComposeMultiTransform

    # Apply warping field to src data
    print('\nApply transfo source --> dest...')
    status, output = sct.run('WarpImageMultiTransform 3 src.nii src_reg.nii -R dest.nii warp_src2dest.nii.gz --use-BSpline')
    if compute_dest2src:
        print('\nApply transfo dest --> source...')
        status, output = sct.run('WarpImageMultiTransform 3 dest.nii dest_reg.nii -R src.nii warp_dest2src.nii.gz --use-BSpline')

    # come back to parent folder
    os.chdir('..')

    # Generate output files
    print('\nGenerate output files...')
    fname_src2dest = sct.generate_output_file(path_tmp+'/src_reg.nii', path_out, file_out, ext_out)
    sct.generate_output_file(path_tmp+'/warp_src2dest.nii.gz', path_out, 'warp_src2dest', '.nii.gz')
    if compute_dest2src:
        fname_dest2src = sct.generate_output_file(path_tmp+'/dest_reg.nii', path_out, file_dest+'_reg', ext_dest)
        sct.generate_output_file(path_tmp+'/warp_dest2src.nii.gz', path_out, 'warp_dest2src', '.nii.gz')

    # Delete temporary files
    if remove_temp_files == 1:
        print '\nRemove temporary files...'
        sct.run('rm -rf '+path_tmp)

    # display elapsed time
    elapsed_time = time.time() - start_time
    print '\nFinished! Elapsed time: '+str(int(round(elapsed_time)))+'s'

    # to view results
    print '\nTo view results, type:'
    print 'fslview '+fname_dest+' '+fname_src2dest+' &'
    if compute_dest2src:
        print 'fslview '+fname_src+' '+fname_dest2src+' &'
    print ''


# Print usage
# ==========================================================================================
def usage():
    print """
"""+os.path.basename(__file__)+"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Part of the Spinal Cord Toolbox <https://sourceforge.net/projects/spinalcordtoolbox>

DESCRIPTION
  This program co-registers two spinal cord volumes. The deformation is non-rigid and is constrained
  in the Z direction (i.e., axial plane). Hence, this function assumes that orientation of the DEST
  image is axial. If you need to register two volumes with large deformations and/or different
  contrasts, it is recommended to input spinal cord segmentations (binary mask) in order to achieve
  maximum robustness. To do so, you can use sct_segmentation_propagation.
  The program outputs a warping field that can be used to register other images to the destination
  image. To apply the warping field to another image, type this:
    WarpImageMultiTransform 3 another_image another_image_reg -R dest_image warp_src2dest.

USAGE
  """+os.path.basename(__file__)+""" -i <source> -d <dest>

MANDATORY ARGUMENTS
  -i <source>                  source image
  -d <dest>                    destination image

OPTIONAL ARGUMENTS
  -s <source_seg>              segmentation for source image (mandatory if -t is used)
  -t <dest_seg>                segmentation for destination image (mandatory if -s is used)
  -q <init_transfo>            transformation file (ITK-based) to apply to source image before
                               registration. Default=none
  -x {0,1}                     compute inverse transformation (dest --> source)
  -z <init_transfo_inv>        inverse transformation file to obtain warp_dest2src deformation field
                               N.B. Only use this flag with -q and -x
  -o <output>                  name of output file. Default=source_reg
  -n <N1xN2>                   number of iterations for first and second stage. Default="""+param.numberIterations+"""
  -y <N>                       number of iterations at step 2 (if using segmentation). Default="""+param.numberIterationsStep2+"""
  -g <gradientStep>            gradientStep for SyN transformation. The larger the more deformation.
                               If you use a segmentation, you can specify gradientStep for each
                               step as follow: val1,val2 (val1: image, val2: seg).
                               Default="""+param.gradientStep[0]+""","""+param.gradientStep[1]+"""
  -p <padding>                 size of padding (top & bottom), to enable deformation at edges.
                               Default="""+str(param.padding)+"""
  -r {0,1}                     remove temporary files. Default='+str(param.remove_temp_files)+'
  -v {0,1}                     verbose. Default="""+str(param.verbose)+"""

EXAMPLE
  Register mean DWI data to the T1 volume using segmentations:
  """+os.path.basename(__file__)+"""
        -i dwi_mean.nii.gz -d t1.nii.gz -s dwi_mean_seg.nii.gz -t t1_seg.nii.gz

  Register another volume to the template using previously-estimated transformations:
  """+os.path.basename(__file__)+"""
        -i $SCT_DIR/data/template/MNI-Poly-AMU_T2.nii.gz
        -d t1.nii.gz
        -s $SCT_DIR/data/template/MNI-Poly-AMU_cord.nii.gz
        -t segmentation_binary.nii.gz
        -q ../t2/warp_template2anat.nii.gz
        -x 1
        -z ../t2/warp_anat2template.nii.gz \n"""

    # exit program
    sys.exit(2)



# pad an image
# ==========================================================================================
def pad_image(fname_in,file_out,padding):
    cmd = 'c3d '+fname_in+' -pad 0x0x'+str(padding)+'vox 0x0x'+str(padding)+'vox 0 -o '+file_out
    print(">> "+cmd)
    os.system(cmd)
    return



# remove padding
# ==========================================================================================
def remove_padding(file_ref,file_in,file_out):
    # remove padding by reslicing padded data into unpadded space
    cmd = 'c3d '+file_ref+' '+file_in+' -reslice-identity -o '+file_out
    print(">> "+cmd)    
    os.system(cmd)
    return

=======
# TODO: add flag -owarpinv
# TODO: if user specified -param, then ignore the default paramreg
# TODO: check syn with shrink=4
# TODO: output name file for warp using "src" and "dest" file name, i.e. warp_filesrc2filedest.nii.gz
# TODO: testing script for all cases
# TODO: add following feature:
# -r of isct_antsRegistration at the initial step (step 0).
# -r [' dest ',' src ',0] --> align the geometric center of the two images
# -r [' dest ',' src ',1] --> align the maximum intensities of the two images I use that quite often...
# TODO: output reg for ants2d and centermass (2016-02-25)

# Note for the developer: DO NOT use --collapse-output-transforms 1, otherwise inverse warping field is not output

# TODO: make three possibilities:
# - one-step registration, using only image registration (by sliceReg or antsRegistration)
# - two-step registration, using first segmentation-based registration (based on sliceReg or antsRegistration) and
# second the image registration (and allow the choice of algo, metric, etc.)
# - two-step registration, using only segmentation-based registration

import sys
import os
import time
import argparse

import numpy as np

from spinalcordtoolbox.reports.qc import generate_qc
from spinalcordtoolbox.registration.register import Paramreg, ParamregMultiStep
from spinalcordtoolbox.utils.shell import Metavar, SmartFormatter, ActionCreateFolder, list_type, display_viewer_syntax
from spinalcordtoolbox.utils.sys import init_sct, printv
from spinalcordtoolbox.utils.fs import extract_fname
from spinalcordtoolbox.image import check_dim

from sct_register_to_template import register_wrapper


def get_parser(paramregmulti=None):
    # Initialize the parser

    if paramregmulti is None:
        step0 = Paramreg(step='0', type='im', algo='syn', metric='MI', iter='0', shrink='1', smooth='0', gradStep='0.5',
                         slicewise='0', dof='Tx_Ty_Tz_Rx_Ry_Rz')  # only used to put src into dest space
        step1 = Paramreg(step='1', type='im')
        paramregmulti = ParamregMultiStep([step0, step1])

    parser = argparse.ArgumentParser(
        description="This program co-registers two 3D volumes. The deformation is non-rigid and is constrained along "
                    "Z direction (i.e., axial plane). Hence, this function assumes that orientation of the destination "
                    "image is axial (RPI). If you need to register two volumes with large deformations and/or "
                    "different contrasts, it is recommended to input spinal cord segmentations (binary mask) in order "
                    "to achieve maximum robustness. The program outputs a warping field that can be used to register "
                    "other images to the destination image. To apply the warping field to another image, use "
                    "'sct_apply_transfo'\n"
                    "\n"
                    "Tips:\n"
                    " - For a registration step using segmentations, use the MeanSquares metric. Also, simple "
                    "algorithm will be very efficient, for example centermass as a 'preregistration'.\n"
                    " - For a registration step using images of different contrast, use the Mutual Information (MI) "
                    "metric.\n"
                    " - Combine the steps by increasing the complexity of the transformation performed in each step, "
                    "for example: -param step=1,type=seg,algo=slicereg,metric=MeanSquares:"
                    "step=2,type=seg,algo=affine,metric=MeanSquares,gradStep=0.2:"
                    "step=3,type=im,algo=syn,metric=MI,iter=5,shrink=2\n"
                    " - When image contrast is low, a good option is to perform registration only based on the image "
                    "segmentation, i.e. using type=seg\n"
                    " - Columnwise algorithm needs to be applied after a translation and rotation such as centermassrot "
                    "algorithm. For example: -param step=1,type=seg,algo=centermassrot,metric=MeanSquares:"
                    "step=2,type=seg,algo=columnwise,metric=MeanSquares",
        formatter_class=SmartFormatter,
        add_help=None,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        '-i',
        metavar=Metavar.file,
        required=True,
        help="Image source. Example: src.nii.gz"
    )
    mandatory.add_argument(
        '-d',
        metavar=Metavar.file,
        required=True,
        help="Image destination. Example: dest.nii.gz"
    )

    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )
    optional.add_argument(
        '-iseg',
        metavar=Metavar.file,
        help="Segmentation source. Example: src_seg.nii.gz"
    )
    optional.add_argument(
        '-dseg',
        metavar=Metavar.file,
        help="Segmentation destination. Example: dest_seg.nii.gz"
    )
    optional.add_argument(
        '-ilabel',
        metavar=Metavar.file,
        help="Labels source."
    )
    optional.add_argument(
        '-dlabel',
        metavar=Metavar.file,
        help="Labels destination."
    )
    optional.add_argument(
        '-initwarp',
        metavar=Metavar.file,
        help="Initial warping field to apply to the source image."
    )
    optional.add_argument(
        '-initwarpinv',
        metavar=Metavar.file,
        help="Initial inverse warping field to apply to the destination image (only use if you wish to generate the "
             "dest->src warping field)"
    )
    optional.add_argument(
        '-m',
        metavar=Metavar.file,
        help="Mask that can be created with sct_create_mask to improve accuracy over region of interest. This mask "
             "will be used on the destination image. Example: mask.nii.gz"
    )
    optional.add_argument(
        '-o',
        metavar=Metavar.file,
        help="Name of output file. Example: src_reg.nii.gz"
    )
    optional.add_argument(
        '-owarp',
        metavar=Metavar.file,
        help="Name of output forward warping field."
    )
    optional.add_argument(
        '-param',
        metavar=Metavar.list,
        type=list_type(':', str),
        help=(f"R|Parameters for registration. Separate arguments with \",\". Separate steps with \":\".\n"
              f"Example: step=1,type=seg,algo=slicereg,metric=MeanSquares:step=2,type=im,algo=syn,metric=MI,iter=5,"
              f"shrink=2\n"
              f"  - step: <int> Step number (starts at 1, except for type=label).\n"
              f"  - type: {{im, seg, imseg, label}} type of data used for registration. Use type=label only at "
              f"step=0.\n"
              f"  - algo: The algorithm used to compute the transformation. Default={paramregmulti.steps['1'].algo}\n"
              f"    * translation: translation in X-Y plane (2dof)\n"
              f"    * rigid: translation + rotation in X-Y plane (4dof)\n"
              f"    * affine: translation + rotation + scaling in X-Y plane (6dof)\n"
              f"    * syn: non-linear symmetric normalization\n"
              f"    * bsplinesyn: syn regularized with b-splines\n"
              f"    * slicereg: regularized translations (see: goo.gl/Sj3ZeU)\n"
              f"    * centermass: slicewise center of mass alignment (seg only).\n"
              f"    * centermassrot: slicewise center of mass and rotation alignment using method specified in "
              f"'rot_method'\n"
              f"    * columnwise: R-L scaling followed by A-P columnwise alignment (seg only).\n"
              f"  - slicewise: <int> Slice-by-slice 2d transformation. "
              f"Default={paramregmulti.steps['1'].slicewise}.\n"
              f"  - metric: {{CC, MI, MeanSquares}}. Default={paramregmulti.steps['1'].metric}.\n"
              f"    * CC: The cross correlation metric compares the images based on their intensities but with a small "
              f"normalization. It can be used with images with the same contrast (for ex. T2-w with T2-w). In this "
              f"case it is very efficient but the computation time can be very long.\n"
              f"    * MI: the mutual information metric compares the images based on their entropy, therefore the "
              f"images need to be big enough to have enough information. It works well for images with different "
              f"contrasts (for example T2-w with T1-w) but not on segmentations.\n"
              f"    * MeanSquares: The mean squares metric compares the images based on their intensities. It can be "
              f"used only with images that have exactly the same contrast (with the same intensity range) or with "
              f"segmentations.\n"
              f"  - iter: <int> Number of iterations. Default={paramregmulti.steps['1'].iter}.\n"
              f"  - shrink: <int> Shrink factor. A shrink factor of 2 will down sample the images by a factor of 2 to "
              f"do the registration, and thus allow bigger deformations (and be faster to compute). It is usually "
              f"combined with a smoothing. (only for syn/bsplinesyn). Default={paramregmulti.steps['1'].shrink}.\n"
              f"  - smooth: <int> Smooth factor (in mm). Note: if algo={{centermassrot,columnwise}} the smoothing "
              f"kernel is: SxSx0. Otherwise it is SxSxS. Default={paramregmulti.steps['1'].smooth}.\n"
              f"  - laplacian: <int> Laplace filter using Gaussian second derivatives, applied before registration. "
              f"The input number correspond to the standard deviation of the Gaussian filter. "
              f"Default={paramregmulti.steps['1'].laplacian}.\n"
              f"  - gradStep: <float> The gradient step used by the function opitmizer. A small gradient step can lead "
              f"to a more accurate registration but will take longer to compute, with the risk to not reach "
              f"convergence. A bigger gradient step will make the registration faster but the result can be far from "
              f"an optimum. Default={paramregmulti.steps['1'].gradStep}.\n"
              f"  - deformation: ?x?x?: Restrict deformation (for ANTs algo). Replace ? by 0 (no deformation) or 1 "
              f"(deformation). Default={paramregmulti.steps['1'].deformation}.\n"
              f"  - init: Initial translation alignment based on:\n"
              f"    * geometric: Geometric center of images\n"
              f"    * centermass: Center of mass of images\n"
              f"    * origin: Physical origin of images\n"
              f"  - poly: <int> Polynomial degree of regularization (only for algo=slicereg). "
              f"Default={paramregmulti.steps['1'].poly}.\n"
              f"  - filter_size: <float> Filter size for regularization (only for algo=centermassrot). "
              f"Default={paramregmulti.steps['1'].filter_size}.\n"
              f"  - smoothWarpXY: <int> Smooth XY warping field (only for algo=columnwize). "
              f"Default={paramregmulti.steps['1'].smoothWarpXY}.\n"
              f"  - pca_eigenratio_th: <int> Min ratio between the two eigenvalues for PCA-based angular adjustment "
              f"(only for algo=centermassrot and rot_method=pca). "
              f"Default={paramregmulti.steps['1'].pca_eigenratio_th}.\n"
              f"  - dof: <str> Degree of freedom for type=label. Separate with '_'. T stands for translation and R "
              f"stands for rotation, x, y, and z indicating the direction. For example, Tx_Ty_Tz_Rx_Ry_Rz would allow "
              f"translation on x, y and z axes and rotation on x, y and z axes. "
              f"Default={paramregmulti.steps['0'].dof}.\n"
              f"  - rot_method {{pca, hog, pcahog}}: rotation method to be used with algo=centermassrot. If using hog "
              f"or pcahog, type should be set to imseg. Default={paramregmulti.steps['1'].rot_method}\n"
              f"    * pca: approximate cord segmentation by an ellipse and finds it orientation using PCA's "
              f"eigenvectors\n"
              f"    * hog: finds the orientation using the symmetry of the image\n"
              f"    * pcahog: tries method pca and if it fails, uses method hog.\n")
    )
    optional.add_argument(
        '-identity',
        metavar=Metavar.int,
        type=int,
        choices=[0, 1],
        default=0,
        help="Just put source into destination (no optimization)."
    )
    optional.add_argument(
        '-z',
        metavar=Metavar.int,
        type=int,
        default=Param().padding,
        help="Size of z-padding to enable deformation at edges when using SyN."
    )
    optional.add_argument(
        '-x',
        choices=['nn', 'linear', 'spline'],
        default='linear',
        help="Final interpolation."
    )
    optional.add_argument(
        '-ofolder',
        metavar=Metavar.folder,
        action=ActionCreateFolder,
        help="Output folder. Example: reg_results/"
    )
    optional.add_argument(
        '-qc',
        metavar=Metavar.folder,
        action=ActionCreateFolder,
        help="The path where the quality control generated content will be saved."
    )
    optional.add_argument(
        '-qc-dataset',
        metavar=Metavar.str,
        help="If provided, this string will be mentioned in the QC report as the dataset the process was run on."
    )
    optional.add_argument(
        '-qc-subject',
        metavar=Metavar.str,
        help="If provided, this string will be mentioned in the QC report as the subject the process was run on."
    )
    optional.add_argument(
        '-r',
        metavar=Metavar.int,
        type=int,
        choices=[0, 1],
        default=1,
        help="Whether to remove temporary files. 0 = no, 1 = yes"
    )
    optional.add_argument(
        '-v',
        choices=['0', '1', '2'],
        default='1',
        help="Verbose. 0: nothing, 1: basic, 2: extended."
    )
    return parser


# DEFAULT PARAMETERS

class Param:
    # The constructor
    def __init__(self):
        self.debug = 0
        self.outSuffix = "_reg"
        self.padding = 5
        self.remove_temp_files = 1


# MAIN
# ==========================================================================================
def main(args=None):
    if args is None:
        args = sys.argv[1:]

    # initialize parameters
    param = Param()

    # Initialization
    fname_output = ''
    path_out = ''
    fname_src_seg = ''
    fname_dest_seg = ''
    fname_src_label = ''
    fname_dest_label = ''
    generate_warpinv = 1

    start_time = time.time()

    # get default registration parameters
    # step1 = Paramreg(step='1', type='im', algo='syn', metric='MI', iter='5', shrink='1', smooth='0', gradStep='0.5')
    step0 = Paramreg(step='0', type='im', algo='syn', metric='MI', iter='0', shrink='1', smooth='0', gradStep='0.5',
                     slicewise='0', dof='Tx_Ty_Tz_Rx_Ry_Rz')  # only used to put src into dest space
    step1 = Paramreg(step='1', type='im')
    paramregmulti = ParamregMultiStep([step0, step1])

    parser = get_parser(paramregmulti=paramregmulti)

    arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    # get arguments
    fname_src = arguments.i
    fname_dest = arguments.d
    if arguments.iseg is not None:
        fname_src_seg = arguments.iseg
    if arguments.dseg is not None:
        fname_dest_seg = arguments.dseg
    if arguments.ilabel is not None:
        fname_src_label = arguments.ilabel
    if arguments.dlabel is not None:
        fname_dest_label = arguments.dlabel
    if arguments.o is not None:
        fname_output = arguments.o
    if arguments.ofolder is not None:
        path_out = arguments.ofolder
    if arguments.owarp is not None:
        fname_output_warp = arguments.owarp
    else:
        fname_output_warp = ''
    if arguments.initwarp is not None:
        fname_initwarp = os.path.abspath(arguments.initwarp)
    else:
        fname_initwarp = ''
    if arguments.initwarpinv is not None:
        fname_initwarpinv = os.path.abspath(arguments.initwarpinv)
    else:
        fname_initwarpinv = ''
    if arguments.m is not None:
        fname_mask = arguments.m
    else:
        fname_mask = ''
    padding = arguments.z
    if arguments.param is not None:
        paramregmulti_user = arguments.param
        # update registration parameters
        for paramStep in paramregmulti_user:
            paramregmulti.addStep(paramStep)
    path_qc = arguments.qc
    qc_dataset = arguments.qc_dataset
    qc_subject = arguments.qc_subject

    identity = arguments.identity
    interp = arguments.x
    remove_temp_files = arguments.r
    verbose = int(arguments.v)
    init_sct(log_level=verbose, update=True)  # Update log level

    # printv(arguments)
    printv('\nInput parameters:')
    printv('  Source .............. ' + fname_src)
    printv('  Destination ......... ' + fname_dest)
    printv('  Init transfo ........ ' + fname_initwarp)
    printv('  Mask ................ ' + fname_mask)
    printv('  Output name ......... ' + fname_output)
    # printv('  Algorithm ........... '+paramregmulti.algo)
    # printv('  Number of iterations  '+paramregmulti.iter)
    # printv('  Metric .............. '+paramregmulti.metric)
    printv('  Remove temp files ... ' + str(remove_temp_files))
    printv('  Verbose ............. ' + str(verbose))

    # update param
    param.verbose = verbose
    param.padding = padding
    param.fname_mask = fname_mask
    param.remove_temp_files = remove_temp_files

    # Get if input is 3D
    printv('\nCheck if input data are 3D...', verbose)
    check_dim(fname_src, dim_lst=[3])
    check_dim(fname_dest, dim_lst=[3])

    # Check if user selected type=seg, but did not input segmentation data
    if 'paramregmulti_user' in locals():
        if True in ['type=seg' in paramregmulti_user[i] for i in range(len(paramregmulti_user))]:
            if fname_src_seg == '' or fname_dest_seg == '':
                printv('\nERROR: if you select type=seg you must specify -iseg and -dseg flags.\n', 1, 'error')

    # Put source into destination space using header (no estimation -- purely based on header)
    # TODO: Check if necessary to do that
    # TODO: use that as step=0
    # printv('\nPut source into destination space using header...', verbose)
    # run_proc('isct_antsRegistration -d 3 -t Translation[0] -m MI[dest_pad.nii,src.nii,1,16] -c 0 -f 1 -s 0 -o
    # [regAffine,src_regAffine.nii] -n BSpline[3]', verbose)
    # if segmentation, also do it for seg

    fname_src2dest, fname_dest2src, _, _ = \
        register_wrapper(fname_src, fname_dest, param, paramregmulti, fname_src_seg=fname_src_seg,
                         fname_dest_seg=fname_dest_seg, fname_src_label=fname_src_label,
                         fname_dest_label=fname_dest_label, fname_mask=fname_mask, fname_initwarp=fname_initwarp,
                         fname_initwarpinv=fname_initwarpinv, identity=identity, interp=interp,
                         fname_output=fname_output,
                         fname_output_warp=fname_output_warp,
                         path_out=path_out)

    # display elapsed time
    elapsed_time = time.time() - start_time
    printv('\nFinished! Elapsed time: ' + str(int(np.round(elapsed_time))) + 's', verbose)

    if path_qc is not None:
        if fname_dest_seg:
            generate_qc(fname_src2dest, fname_in2=fname_dest, fname_seg=fname_dest_seg, args=args,
                        path_qc=os.path.abspath(path_qc), dataset=qc_dataset, subject=qc_subject,
                        process='sct_register_multimodal')
        else:
            printv('WARNING: Cannot generate QC because it requires destination segmentation.', 1, 'warning')

    if generate_warpinv:
        display_viewer_syntax([fname_src, fname_dest2src], verbose=verbose)
    display_viewer_syntax([fname_dest, fname_src2dest], verbose=verbose)
>>>>>>> upstream/master


# START PROGRAM
# ==========================================================================================
if __name__ == "__main__":
<<<<<<< HEAD
    # initialize parameters
    param = param()
    # call main function
    main()



    # Convert deformation field to 4D volume (readable by fslview)
    # DONE: clean code below-- right now it does not work
    #===========
    #if convertDeformation:
    #    print('\nConvert deformation field...')
    #    cmd = 'c3d -mcs tmp.regWarp.nii -oo tmp.regWarp_x.nii tmp.regWarp_y.nii tmp.regWarp_z.nii'
    #    print(">> "+cmd)
    #    os.system(cmd)
    #    cmd = 'fslmerge -t '+path_out+'warp_comp.nii tmp.regWarp_x.nii tmp.regWarp_y.nii tmp.regWarp_z.nii'
    #    print(">> "+cmd)
    #    os.system(cmd)
    #===========
=======
    init_sct()
    # call main function
    main()
>>>>>>> upstream/master
