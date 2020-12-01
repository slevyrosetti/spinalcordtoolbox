#!/usr/bin/env python
<<<<<<< HEAD
#########################################################################################
# Convert bvecs file to column, in case they are in line.
#
#
# USAGE
# ---------------------------------------------------------------------------------------
#   sct_dmri_transpose_bvecs.py <bvecs>
#
#
# INPUT
# ---------------------------------------------------------------------------------------
# bvecs             bvecs ASCII file (FSL format).
#
#
# OUTPUT
# ---------------------------------------------------------------------------------------
# bvecs_col         bvecs in column.
#
#
# DEPENDENCIES
# ---------------------------------------------------------------------------------------
# EXTERNAL PYTHON PACKAGES
# none
#
# EXTERNAL SOFTWARE
# none
#
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 NeuroPoly, Polytechnique Montreal <www.neuropoly.info>
# Author: Julien Cohen-Adad
# Modified: 2013-10-19
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#########################################################################################

import sys
import os


# init
#fname_in = '/Users/julien/MRI/david_cadotte/2013-10-19_multiparametric/bvecs.txt'


# Extracts path, file and extension
#########################################################################################
def extract_fname(fname):
    # extract path
    path_fname = os.path.dirname(fname)+'/'
    # check if only single file was entered (without path)
    if path_fname == '/':
        path_fname = ''
    # extract file and extension
    file_fname = fname
    file_fname = file_fname.replace(path_fname,'')
    file_fname, ext_fname = os.path.splitext(file_fname)
    # check if .nii.gz file
    if ext_fname == '.gz':
        file_fname = file_fname[0:len(file_fname)-4]
        ext_fname = ".nii.gz"
    return path_fname, file_fname, ext_fname
#########################################################################################


# MAIN
#########################################################################################

# Check inputs
path_func, file_func, ext_func = extract_fname(sys.argv[0])
if len(sys.argv) < 2:
    print 'Usage: '+file_func+ext_func+' <bvecs>'
    sys.exit(1)
fname_in = sys.argv[1]

# Extracts path, file and extension
path_in, file_in, ext_in = extract_fname(fname_in)

# read ASCII file
print('Read file...')
text_file = open(fname_in, 'r')
bvecs = text_file.readlines()
text_file.close()

# Parse each line
# TODO: find a better way to do it, maybe with string or numpy...
lin0 = bvecs[0].split()
lin1 = bvecs[1].split()
lin2 = bvecs[2].split()

# Write new file
print('Transpose bvecs...')
fname_out = path_in+file_in+'_t'+ext_in
fid = open(fname_out,'w')
for iCol in xrange(0, len(lin0)):
    fid.write(lin0[iCol]+' '+lin1[iCol]+' '+lin2[iCol]+'\n')
fid.close()

# Display
print('File created: '+fname_out)

=======
# =======================================================================================================================
#
# Transpose bvecs file (if necessary) to get nx3 structure
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Julien Cohen-Adad
#
# About the license: see the file LICENSE.TXT
#########################################################################################

#!/usr/bin/env python
#########################################################################################
#
# Compute DTI.
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Julien Cohen-Adad
#
# About the license: see the file LICENSE.TXT
#########################################################################################

import os
import sys
import argparse

from spinalcordtoolbox.utils import Metavar, SmartFormatter, init_sct, extract_fname, printv


def get_parser():
    # Initialize the parser
    parser = argparse.ArgumentParser(
        description='Transpose bvecs file (if necessary) to get nx3 structure.',
        formatter_class=SmartFormatter,
        add_help=None,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        '-bvec',
        metavar=Metavar.file,
        required=True,
        help="Input bvecs file. Example: bvecs.txt"
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )
    optional.add_argument(
        '-o',
        metavar=Metavar.file,
        default='',
        help="Output bvecs file. By default, input file is overwritten. Example: bvecs_t.txt"
    )
    optional.add_argument(
        '-v',
        choices=['0', '1', '2'],
        default='1',
        help="Verbose: 0 = nothing, 1 = basic, 2 = extended."
    )

    return parser


# MAIN
# ==========================================================================================
def main(args=None):

    parser = get_parser()
    if args:
        arguments = parser.parse_args(args)
    else:
        arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    fname_in = arguments.bvec
    fname_out = arguments.o
    verbose = int(arguments.v)
    init_sct(log_level=verbose, update=True)  # Update log level

    # get bvecs in proper orientation
    from dipy.io import read_bvals_bvecs
    bvals, bvecs = read_bvals_bvecs(None, fname_in)

    # # Transpose bvecs
    # printv('Transpose bvecs...', verbose)
    # # from numpy import transpose
    # bvecs = bvecs.transpose()

    # Write new file
    if fname_out == '':
        path_in, file_in, ext_in = extract_fname(fname_in)
        fname_out = path_in + file_in + ext_in
    fid = open(fname_out, 'w')
    for iLine in range(bvecs.shape[0]):
        fid.write(' '.join(str(i) for i in bvecs[iLine, :]) + '\n')
    fid.close()

    # display message
    printv('Created file:\n--> ' + fname_out + '\n', verbose, 'info')


# Start program
# =======================================================================================================================
if __name__ == "__main__":
    init_sct()
    # call main function
    main()
>>>>>>> upstream/master
