#!/usr/bin/env python
#########################################################################################
#
# Separate b=0 and DW images from diffusion dataset.
#
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Julien Cohen-Adad
<<<<<<< HEAD
# Modified: 2014-05-22
=======
# Modified: 2014-08-14
>>>>>>> upstream/master
#
# About the license: see the file LICENSE.TXT
#########################################################################################

<<<<<<< HEAD


# DEFAULT PARAMETERS
class param:
    ## The constructor
    def __init__(self):
        self.debug              = 0
        self.verbose            = 0 # verbose

import re
import sys
import getopt
import os
import math
import time
import sct_utils as sct

=======
import sys
import math
import time
import argparse
import os

import numpy as np

from spinalcordtoolbox.image import Image, generate_output_file
from spinalcordtoolbox.utils.shell import Metavar, SmartFormatter, ActionCreateFolder
from spinalcordtoolbox.utils.sys import init_sct, run_proc, printv
from spinalcordtoolbox.utils.fs import tmp_create, copy, extract_fname, rmtree

from sct_image import split_data
from sct_convert import convert


class Param:
    def __init__(self):
        self.debug = 0
        self.average = 1
        self.remove_temp_files = 1
        self.verbose = 1
        self.bval_min = 100  # in case user does not have min bvalues at 0, set threshold.


def get_parser():
    # Initialize parser
    param_default = Param()
    parser = argparse.ArgumentParser(
        description="Separate b=0 and DW images from diffusion dataset. The output files will have a suffix "
                    "(_b0 and _dwi) appended to the input file name.",
        formatter_class=SmartFormatter,
        add_help=None,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        '-i',
        metavar=Metavar.file,
        required=True,
        help="Diffusion data. Example: dmri.nii.gz"
    )
    mandatory.add_argument(
        '-bvec',
        metavar=Metavar.file,
        required=True,
        help="Bvecs file. Example: bvecs.txt"
    )

    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )
    optional.add_argument(
        '-a',
        type=int,
        choices=(0, 1),
        default=param_default.average,
        help="Average b=0 and DWI data. 0 = no, 1 = yes"
    )
    optional.add_argument(
        '-bval',
        metavar=Metavar.file,
        default="",
        help='bvals file. Used to identify low b-values (in case different from 0). Example: bvals.nii.gz',
    )
    optional.add_argument(
        '-bvalmin',
        type=float,
        metavar=Metavar.float,
        help='B-value threshold (in s/mm2) below which data is considered as b=0. Example: 50.0',
    )
    optional.add_argument(
        '-ofolder',
        metavar=Metavar.folder,
        action=ActionCreateFolder,
        default='./',
        help='Output folder. Example: dmri_separate_results/',
    )
    optional.add_argument(
        "-r",
        choices=('0', '1'),
        default=param_default.remove_temp_files,
        help="Remove temporary files. 0 = no, 1 = yes"
    )
    optional.add_argument(
        '-v',
        choices=('0', '1'),
        default=param_default.verbose,
        help='Verbose. 0 = nothing, 1 = expanded',
    )
    return parser
>>>>>>> upstream/master


# MAIN
# ==========================================================================================
<<<<<<< HEAD
def main():

    # Initialization
    path_script = os.path.dirname(__file__)
    fsloutput = 'export FSLOUTPUTTYPE=NIFTI; ' # for faster processing, all outputs are in NIFTI
    # THIS DOES NOT WORK IN MY LAPTOP: path_sct = os.environ['SCT_DIR'] # path to spinal cord toolbox
    path_sct = path_script[:-8] # TODO: make it cleaner!
    fname_data = ''
    fname_bvecs = ''
    verbose = param.verbose
    start_time = time.time()

    # Parameters for debug mode
    if param.debug:
        fname_data = os.path.expanduser("~")+'/code/spinalcordtoolbox_dev/testing/data/errsm_22/dmri/dmri.nii.gz'
        fname_bvecs = os.path.expanduser("~")+'/code/spinalcordtoolbox_dev/testing/data/errsm_22/dmri/bvecs.txt'
        verbose = 1

    # Check input parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hb:i:v:')
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ("-b"):
            fname_bvecs = arg
        elif opt in ("-i"):
            fname_data = arg
        elif opt in ('-v'):
            verbose = int(arg)

    # display usage if a mandatory argument is not provided
    if fname_data == '' or fname_bvecs == '':
        usage()

    # check existence of input files
    sct.check_file_exist(fname_data)
    sct.check_file_exist(fname_bvecs)

    # print arguments
    print '\nCheck parameters:'
    print '.. DWI data:             '+fname_data
    print '.. bvecs file:           '+fname_bvecs

    # Extract path, file and extension
    path_data, file_data, ext_data = sct.extract_fname(fname_data)

    # create temporary folder
    path_tmp = 'tmp.'+time.strftime("%y%m%d%H%M%S")
    sct.run('mkdir '+path_tmp)

    # copy files into tmp folder
    sct.run('cp '+fname_data+' '+path_tmp)
    sct.run('cp '+fname_bvecs+' '+path_tmp)

    # go to tmp folder
    os.chdir(path_tmp)

    # Get size of data
    print '\nGet dimensions data...'
    nx, ny, nz, nt, px, py, pz, pt = sct.get_dimension(fname_data)
    print '.. '+str(nx)+' x '+str(ny)+' x '+str(nz)+' x '+str(nt)

    # Open bvecs file
    bvecs = []
    with open(fname_bvecs) as f:
        for line in f:
            bvecs_new = map(float, line.split())
            bvecs.append(bvecs_new)

    # Check if bvecs file is nx3
    if not len(bvecs[0][:]) == 3:
        print 'WARNING: bvecs file is 3xn instead of nx3. Consider using sct_dmri_transpose_bvecs'
        # transpose bvecs
        bvecs = zip(*bvecs)

    # Identify b=0 and DW images
    print '\nIdentify b=0 and DW images...'
    index_b0 = []
    index_dwi = []
    for it in xrange(0,nt):
        if math.sqrt(math.fsum([i**2 for i in bvecs[it]])) < 0.01:
            index_b0.append(it)
        else:
            index_dwi.append(it)
    nb_b0 = len(index_b0)
    nb_dwi = len(index_dwi)
    print '.. Number of b=0: '+str(nb_b0)+' '+str(index_b0)
    print '.. Number of DWI: '+str(nb_dwi)+' '+str(index_dwi)

    #TODO: check if number of bvecs and nt match

    # Split into T dimension
    print '\nSplit along T dimension...'
    sct.run(fsloutput+' fslsplit '+fname_data+' data_splitT')

    # retrieve output names
    status, output = sct.run('ls data_splitT*.*')
    file_data_split = output.split()
    # Remove .nii extension
    file_data_split = [file_data_split[i].replace('.nii','') for i in xrange (0,len(file_data_split))]

    # Merge b=0 images
    print '\nMerge b=0...'
    cmd = fsloutput+'fslmerge -t b0'
    for it in xrange(0,nb_b0):
        cmd += ' '+file_data_split[index_b0[it]]
    sct.run(cmd)

    # Merge DWI images
    print '\nMerge DWI...'
    cmd = fsloutput+'fslmerge -t dwi'
    for it in xrange(0,nb_dwi):
        cmd += ' '+file_data_split[index_dwi[it]]
    sct.run(cmd)

    # come back to parent folder
    os.chdir('..')

    # Generate output files
    print('\nGenerate output files...')
    sct.generate_output_file(path_tmp+'/b0.nii',path_data,'b0',ext_data)
    sct.generate_output_file(path_tmp+'/dwi.nii',path_data,'dwi',ext_data)

    # Remove temporary files
    print('\nRemove temporary files...')
    sct.run('rm -rf '+path_tmp)

    # display elapsed time
    elapsed_time = time.time() - start_time
    print '\nFinished! Elapsed time: '+str(int(round(elapsed_time)))+'s'

    # to view results
    print '\nTo view results, type:'
    print 'fslview b0 dwi &\n'

    # End of Main


# Print usage
# ==========================================================================================
def usage():
    print '\n' \
        ''+os.path.basename(__file__)+'\n' \
        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' \
        'Part of the Spinal Cord Toolbox <https://sourceforge.net/projects/spinalcordtoolbox>\n' \
        '\n'\
        'DESCRIPTION\n' \
        '  Separate b=0 and DW images from diffusion dataset.\n' \
        '\n' \
        'USAGE\n' \
        '  '+os.path.basename(__file__)+' -i <dmri> -w <bvecs>\n' \
        '\n' \
        'MANDATORY ARGUMENTS\n' \
        '  -i <dmri>                  diffusion data\n' \
        '  -b <bvecs>                 bvecs file\n' \
        '\n' \
        'OPTIONAL ARGUMENTS\n' \
        '  -v <0,1>                   verbose. Default='+str(param.verbose)+'.\n'

    # exit program
    sys.exit(2)
=======
def main(args=None):
    # initialize parameters
    param = Param()
    # call main function
    parser = get_parser()
    if args:
        arguments = parser.parse_args(args)
    else:
        arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    fname_data = arguments.i
    fname_bvecs = arguments.bvec
    average = arguments.a
    verbose = int(arguments.v)
    init_sct(log_level=verbose, update=True)  # Update log level
    remove_temp_files = arguments.r
    path_out = arguments.ofolder

    fname_bvals = arguments.bval
    if arguments.bvalmin:
        param.bval_min = arguments.bvalmin

    # Initialization
    start_time = time.time()

    # printv(arguments)
    printv('\nInput parameters:', verbose)
    printv('  input file ............' + fname_data, verbose)
    printv('  bvecs file ............' + fname_bvecs, verbose)
    printv('  bvals file ............' + fname_bvals, verbose)
    printv('  average ...............' + str(average), verbose)

    # Get full path
    fname_data = os.path.abspath(fname_data)
    fname_bvecs = os.path.abspath(fname_bvecs)
    if fname_bvals:
        fname_bvals = os.path.abspath(fname_bvals)

    # Extract path, file and extension
    path_data, file_data, ext_data = extract_fname(fname_data)

    # create temporary folder
    path_tmp = tmp_create(basename="dmri_separate")

    # copy files into tmp folder and convert to nifti
    printv('\nCopy files into temporary folder...', verbose)
    ext = '.nii'
    dmri_name = 'dmri'
    b0_name = file_data + '_b0'
    b0_mean_name = b0_name + '_mean'
    dwi_name = file_data + '_dwi'
    dwi_mean_name = dwi_name + '_mean'

    if not convert(fname_data, os.path.join(path_tmp, dmri_name + ext)):
        printv('ERROR in convert.', 1, 'error')
    copy(fname_bvecs, os.path.join(path_tmp, "bvecs"), verbose=verbose)

    # go to tmp folder
    curdir = os.getcwd()
    os.chdir(path_tmp)

    # Get size of data
    im_dmri = Image(dmri_name + ext)
    printv('\nGet dimensions data...', verbose)
    nx, ny, nz, nt, px, py, pz, pt = im_dmri.dim
    printv('.. ' + str(nx) + ' x ' + str(ny) + ' x ' + str(nz) + ' x ' + str(nt), verbose)

    # Identify b=0 and DWI images
    printv(fname_bvals)
    index_b0, index_dwi, nb_b0, nb_dwi = identify_b0(fname_bvecs, fname_bvals, param.bval_min, verbose)

    # Split into T dimension
    printv('\nSplit along T dimension...', verbose)
    im_dmri_split_list = split_data(im_dmri, 3)
    for im_d in im_dmri_split_list:
        im_d.save()

    # Merge b=0 images
    printv('\nMerge b=0...', verbose)
    from sct_image import concat_data
    l = []
    for it in range(nb_b0):
        l.append(dmri_name + '_T' + str(index_b0[it]).zfill(4) + ext)
    im_out = concat_data(l, 3).save(b0_name + ext)

    # Average b=0 images
    if average:
        printv('\nAverage b=0...', verbose)
        run_proc(['sct_maths', '-i', b0_name + ext, '-o', b0_mean_name + ext, '-mean', 't'], verbose)

    # Merge DWI
    l = []
    for it in range(nb_dwi):
        l.append(dmri_name + '_T' + str(index_dwi[it]).zfill(4) + ext)
    im_out = concat_data(l, 3).save(dwi_name + ext)

    # Average DWI images
    if average:
        printv('\nAverage DWI...', verbose)
        run_proc(['sct_maths', '-i', dwi_name + ext, '-o', dwi_mean_name + ext, '-mean', 't'], verbose)

    # come back
    os.chdir(curdir)

    # Generate output files
    fname_b0 = os.path.abspath(os.path.join(path_out, b0_name + ext_data))
    fname_dwi = os.path.abspath(os.path.join(path_out, dwi_name + ext_data))
    fname_b0_mean = os.path.abspath(os.path.join(path_out, b0_mean_name + ext_data))
    fname_dwi_mean = os.path.abspath(os.path.join(path_out, dwi_mean_name + ext_data))
    printv('\nGenerate output files...', verbose)
    generate_output_file(os.path.join(path_tmp, b0_name + ext), fname_b0, verbose=verbose)
    generate_output_file(os.path.join(path_tmp, dwi_name + ext), fname_dwi, verbose=verbose)
    if average:
        generate_output_file(os.path.join(path_tmp, b0_mean_name + ext), fname_b0_mean, verbose=verbose)
        generate_output_file(os.path.join(path_tmp, dwi_mean_name + ext), fname_dwi_mean, verbose=verbose)

    # Remove temporary files
    if remove_temp_files == 1:
        printv('\nRemove temporary files...', verbose)
        rmtree(path_tmp, verbose=verbose)

    # display elapsed time
    elapsed_time = time.time() - start_time
    printv('\nFinished! Elapsed time: ' + str(int(np.round(elapsed_time))) + 's', verbose)

    return fname_b0, fname_b0_mean, fname_dwi, fname_dwi_mean


# ==========================================================================================
# identify b=0 and DW images
# ==========================================================================================
def identify_b0(fname_bvecs, fname_bvals, bval_min, verbose):

    # Identify b=0 and DWI images
    printv('\nIdentify b=0 and DWI images...', verbose)
    index_b0 = []
    index_dwi = []

    # if bval is not provided
    if not fname_bvals:
        # Open bvecs file
        bvecs = []
        with open(fname_bvecs) as f:
            for line in f:
                bvecs_new = [x for x in map(float, line.split())]
                bvecs.append(bvecs_new)

        # Check if bvecs file is nx3
        if not len(bvecs[0][:]) == 3:
            printv('  WARNING: bvecs file is 3xn instead of nx3. Consider using sct_dmri_transpose_bvecs.', verbose, 'warning')
            printv('  Transpose bvecs...', verbose)
            # transpose bvecs
            bvecs = list(zip(*bvecs))

        # get number of lines
        nt = len(bvecs)

        # identify b=0 and dwi
        for it in range(0, nt):
            if math.sqrt(math.fsum([i**2 for i in bvecs[it]])) < 0.01:
                index_b0.append(it)
            else:
                index_dwi.append(it)

    # if bval is provided
    else:

        # Open bvals file
        from dipy.io import read_bvals_bvecs
        bvals, bvecs = read_bvals_bvecs(fname_bvals, fname_bvecs)

        # get number of lines
        nt = len(bvals)

        # Identify b=0 and DWI images
        printv('\nIdentify b=0 and DWI images...', verbose)
        for it in range(0, nt):
            if bvals[it] < bval_min:
                index_b0.append(it)
            else:
                index_dwi.append(it)

    # check if no b=0 images were detected
    if index_b0 == []:
        printv('ERROR: no b=0 images detected. Maybe you are using non-null low bvals? in that case use flag -bvalmin. Exit program.', 1, 'error')
        sys.exit(2)

    # display stuff
    nb_b0 = len(index_b0)
    nb_dwi = len(index_dwi)
    printv('  Number of b=0: ' + str(nb_b0) + ' ' + str(index_b0), verbose)
    printv('  Number of DWI: ' + str(nb_dwi) + ' ' + str(index_dwi), verbose)

    # return
    return index_b0, index_dwi, nb_b0, nb_dwi
>>>>>>> upstream/master


# START PROGRAM
# ==========================================================================================
if __name__ == "__main__":
<<<<<<< HEAD
    # initialize parameters
    param = param()
    # call main function
    main()
=======
    init_sct()
    main()
>>>>>>> upstream/master
