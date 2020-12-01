#!/usr/bin/env python
#########################################################################################
#
<<<<<<< HEAD
# Register anatomical image to the template using the spinal cord centerline/segmentation.
# we assume here that we have a RPI orientation, where Z axis is inferior-superior direction
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2013 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Benjamin De Leener, Julien Cohen-Adad
# Modified: 2014-06-03
=======
# All sort of utilities for labels.
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Benjamin De Leener, Julien Cohen-Adad
# Modified: 2015-02-11
>>>>>>> upstream/master
#
# About the license: see the file LICENSE.TXT
#########################################################################################

<<<<<<< HEAD
# TODO: currently it seems like cross_radius is given in pixel instead of mm

import os, sys
import getopt
import commands
import sys
import sct_utils as sct
import nibabel
import numpy as np

# DEFAULT PARAMETERS
class param:
    ## The constructor
    def __init__(self):
        self.debug               = 0


#=======================================================================================================================
# main
#=======================================================================================================================
def main():

    # Initialization
    fname_label = ''
    fname_label_output = ''
    cross_radius = 5
    dilate = False;
    fname_ref = ''
    type_process = ''
    output_level = 0 # 0 for image with point ; 1 for txt file

    # get path of the toolbox
    status, path_sct = commands.getstatusoutput('echo $SCT_DIR')

    # Parameters for debug mode
    if param.debug:
        print '\n*** WARNING: DEBUG MODE ON ***\n'
        fname_label = path_sct+'/testing/data/errsm_23/t2/landmarks_rpi.nii.gz'
        fname_label_output = 'landmarks_rpi_output.nii.gz'
        type_process = 'cross'
        cross_radius = 5
        dilate = True

    # extract path of the script
    path_script = os.path.dirname(__file__)+'/'
    
    # Check input param
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hi:o:c:r:t:l:d')
    except getopt.GetoptError as err:
        print str(err)
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ('-i'):
            fname_label = arg
        elif opt in ('-o'):
            fname_label_output = arg
        elif opt in ('-c'):
            cross_radius = int(arg)
        elif opt in ('-d'):
            dilate = True
        elif opt in ('-r'):
            fname_ref = arg
        elif opt in ('-t'):
            type_process = arg
        elif opt in ('-l'):
            output_level = int(arg)

    # display usage if a mandatory argument is not provided
    if fname_label == '':
        usage()
        
    # check existence of input files
    sct.check_file_exist(fname_label)
    if fname_ref != '':
        sct.check_file_exist(fname_ref)
    
    # extract path/file/extension
    path_label, file_label, ext_label = sct.extract_fname(fname_label)
    path_label_output, file_label_output, ext_label_output = sct.extract_fname(fname_label_output)

    # read nifti input file
    img = nibabel.load(fname_label)
    # 3d array for each x y z voxel values for the input nifti image
    data = img.get_data()
    hdr = img.get_header()

    #print '\nGet dimensions of input centerline...'
    nx, ny, nz, nt, px, py, pz, pt = sct.get_dimension(fname_label)
    #print '.. matrix size: '+str(nx)+' x '+str(ny)+' x '+str(nz)
    #print '.. voxel size:  '+str(px)+'mm x '+str(py)+'mm x '+str(pz)+'mm'


    if type_process == 'cross':
        data = cross(data, cross_radius, fname_ref, dilate, px, py)
    elif type_process == 'remove':
        data = remove_label(data, fname_ref)
    elif type_process == 'disk':
        extract_disk_position(data, fname_ref, output_level, fname_label_output)
    elif type_process == 'centerline':
        extract_centerline(data, fname_label_output)
        output_level = 1
    elif type_process == 'segmentation':
        extract_segmentation(data, fname_label_output)
        output_level = 1
    elif type_process == 'fraction-volume':
        fraction_volume(data,fname_ref,fname_label_output)
        output_level = 1
    elif type_process == 'write-vert-levels':
        write_vertebral_levels(data,fname_ref)
    elif type_process == 'display-voxel':
        display_voxel(data)
        output_level = 1

    if (output_level == 0):
        hdr.set_data_dtype('int32') # set imagetype to uint8, previous: int32. 
        print '\nWrite NIFTI volumes...'
        data.astype('int')
        img = nibabel.Nifti1Image(data, None, hdr)
        nibabel.save(img, 'tmp.'+file_label_output+'.nii.gz')
        sct.generate_output_file('tmp.'+file_label_output+'.nii.gz','./',file_label_output,ext_label_output)


#=======================================================================================================================
def cross(data, cross_radius, fname_ref, dilate, px, py):
    X, Y, Z = (data > 0).nonzero()
    a = len(X)
    d = cross_radius # cross radius in pixel
    dx = d/px # cross radius in mm
    dy = d/py

    # for all points with non-zeros neighbors, force the neighbors to 0
    for i in range(0,a):
        value = data[X[i]][Y[i]][Z[i]]
        data[X[i]][Y[i]][Z[i]] = 0 # remove point on the center of the spinal cord
        if fname_ref == '':
            data[X[i]][Y[i]+dy][Z[i]] = value*10+1 # add point at distance from center of spinal cord
            data[X[i]+dx][Y[i]][Z[i]] = value*10+2
            data[X[i]][Y[i]-dy][Z[i]] = value*10+3
            data[X[i]-dx][Y[i]][Z[i]] = value*10+4

            # dilate cross to 3x3
            if dilate:
                data[X[i]-1][Y[i]+dy-1][Z[i]] = data[X[i]][Y[i]+dy-1][Z[i]] = data[X[i]+1][Y[i]+dy-1][Z[i]] = data[X[i]+1][Y[i]+dy][Z[i]] = data[X[i]+1][Y[i]+dy+1][Z[i]] = data[X[i]][Y[i]+dy+1][Z[i]] = data[X[i]-1][Y[i]+dy+1][Z[i]] = data[X[i]-1][Y[i]+dy][Z[i]] = data[X[i]][Y[i]+dy][Z[i]]
                data[X[i]+dx-1][Y[i]-1][Z[i]] = data[X[i]+dx][Y[i]-1][Z[i]] = data[X[i]+dx+1][Y[i]-1][Z[i]] = data[X[i]+dx+1][Y[i]][Z[i]] = data[X[i]+dx+1][Y[i]+1][Z[i]] = data[X[i]+dx][Y[i]+1][Z[i]] = data[X[i]+dx-1][Y[i]+1][Z[i]] = data[X[i]+dx-1][Y[i]][Z[i]] = data[X[i]+dx][Y[i]][Z[i]]
                data[X[i]-1][Y[i]-dy-1][Z[i]] = data[X[i]][Y[i]-dy-1][Z[i]] = data[X[i]+1][Y[i]-dy-1][Z[i]] = data[X[i]+1][Y[i]-dy][Z[i]] = data[X[i]+1][Y[i]-dy+1][Z[i]] = data[X[i]][Y[i]-dy+1][Z[i]] = data[X[i]-1][Y[i]-dy+1][Z[i]] = data[X[i]-1][Y[i]-dy][Z[i]] = data[X[i]][Y[i]-dy][Z[i]]
                data[X[i]-dx-1][Y[i]-1][Z[i]] = data[X[i]-dx][Y[i]-1][Z[i]] = data[X[i]-dx+1][Y[i]-1][Z[i]] = data[X[i]-dx+1][Y[i]][Z[i]] = data[X[i]-dx+1][Y[i]+1][Z[i]] = data[X[i]-dx][Y[i]+1][Z[i]] = data[X[i]-dx-1][Y[i]+1][Z[i]] = data[X[i]-dx-1][Y[i]][Z[i]] = data[X[i]-dx][Y[i]][Z[i]]
        else:
            # read nifti input file
            img_ref = nibabel.load(fname_ref)
            # 3d array for each x y z voxel values for the input nifti image
            data_ref = img_ref.get_data()
            profile = []; p_median = []; p_gradient = []
            for j in range(0,d+1):
                profile.append(data_ref[X[i]][Y[i]+j][Z[i]])
            for k in range(1,d):
                a = np.array([profile[k-1],profile[k],profile[k+1]])
                p_median.append(np.median(a))
            for l in range(0,d-2):
                p_gradient.append(p_median[l+1]-p_median[l])
            d1 = p_gradient.index(max(p_gradient))

            profile = []; p_median = []; p_gradient = []
            for j in range(0,d+1):
                profile.append(data_ref[X[i]+j][Y[i]][Z[i]])
            for k in range(1,d):
                a = np.array([profile[k-1],profile[k],profile[k+1]])
                p_median.append(np.median(a))
            for l in range(0,d-2):
                p_gradient.append(p_median[l+1]-p_median[l])
            d2 = p_gradient.index(max(p_gradient))

            profile = []; p_median = []; p_gradient = []
            for j in range(0,d+1):
                profile.append(data_ref[X[i]][Y[i]-j][Z[i]])
            for k in range(1,d):
                a = np.array([profile[k-1],profile[k],profile[k+1]])
                p_median.append(np.median(a))
            for l in range(0,d-2):
                p_gradient.append(p_median[l+1]-p_median[l])
            d3 = p_gradient.index(max(p_gradient))

            profile = []; p_median = []; p_gradient = []
            for j in range(0,d+1):
                profile.append(data_ref[X[i]-j][Y[i]][Z[i]])
            for k in range(1,d):
                a = np.array([profile[k-1],profile[k],profile[k+1]])
                p_median.append(np.median(a))
            for l in range(0,d-2):
                p_gradient.append(p_median[l+1]-p_median[l])
            d4 = p_gradient.index(max(p_gradient))

            data[X[i]][Y[i]+d1][Z[i]] = value*10+1 # add point at distance from center of spinal cord
            data[X[i]+d2][Y[i]][Z[i]] = value*10+2
            data[X[i]][Y[i]-d3][Z[i]] = value*10+3
            data[X[i]-d4][Y[i]][Z[i]] = value*10+4

            # dilate cross to 3x3
            if dilate:
                data[X[i]-1][Y[i]+d1-1][Z[i]] = data[X[i]][Y[i]+d1-1][Z[i]] = data[X[i]+1][Y[i]+d1-1][Z[i]] = data[X[i]+1][Y[i]+d1][Z[i]] = data[X[i]+1][Y[i]+d1+1][Z[i]] = data[X[i]][Y[i]+d1+1][Z[i]] = data[X[i]-1][Y[i]+d1+1][Z[i]] = data[X[i]-1][Y[i]+d1][Z[i]] = data[X[i]][Y[i]+d1][Z[i]]
                data[X[i]+d2-1][Y[i]-1][Z[i]] = data[X[i]+d2][Y[i]-1][Z[i]] = data[X[i]+d2+1][Y[i]-1][Z[i]] = data[X[i]+d2+1][Y[i]][Z[i]] = data[X[i]+d2+1][Y[i]+1][Z[i]] = data[X[i]+d2][Y[i]+1][Z[i]] = data[X[i]+d2-1][Y[i]+1][Z[i]] = data[X[i]+d2-1][Y[i]][Z[i]] = data[X[i]+d2][Y[i]][Z[i]]
                data[X[i]-1][Y[i]-d3-1][Z[i]] = data[X[i]][Y[i]-d3-1][Z[i]] = data[X[i]+1][Y[i]-d3-1][Z[i]] = data[X[i]+1][Y[i]-d3][Z[i]] = data[X[i]+1][Y[i]-d3+1][Z[i]] = data[X[i]][Y[i]-d3+1][Z[i]] = data[X[i]-1][Y[i]-d3+1][Z[i]] = data[X[i]-1][Y[i]-d3][Z[i]] = data[X[i]][Y[i]-d3][Z[i]]
                data[X[i]-d4-1][Y[i]-1][Z[i]] = data[X[i]-d4][Y[i]-1][Z[i]] = data[X[i]-d4+1][Y[i]-1][Z[i]] = data[X[i]-d4+1][Y[i]][Z[i]] = data[X[i]-d4+1][Y[i]+1][Z[i]] = data[X[i]-d4][Y[i]+1][Z[i]] = data[X[i]-d4-1][Y[i]+1][Z[i]] = data[X[i]-d4-1][Y[i]][Z[i]] = data[X[i]-d4][Y[i]][Z[i]]

    return data


#=======================================================================================================================
def remove_label(data, fname_ref):
    X, Y, Z = (data > 0).nonzero()

    img_ref = nibabel.load(fname_ref)
    # 3d array for each x y z voxel values for the input nifti image
    data_ref = img_ref.get_data()
    X_ref, Y_ref, Z_ref = (data_ref > 0).nonzero()

    nbLabel = len(X)
    nbLabel_ref = len(X_ref)
    for i in range(0,nbLabel):
        value = data[X[i]][Y[i]][Z[i]]
        isInRef = False
        for j in range(0,nbLabel_ref):
            value_ref = data_ref[X_ref[j]][Y_ref[j]][Z_ref[j]]
            if value_ref == value:
                isInRef = True
        if isInRef == False:
            data[X[i]][Y[i]][Z[i]] = 0

    return data

# need binary centerline and segmentation with vertebral level. output_level=1 -> write .txt file. output_level=1 -> write centerline with vertebral levels
#=======================================================================================================================
def extract_disk_position(data_level, fname_centerline, output_level, fname_label_output):
    X, Y, Z = (data_level > 0).nonzero()
    
    img_centerline = nibabel.load(fname_centerline)
    # 3d array for each x y z voxel values for the input nifti image
    data_centerline = img_centerline.get_data()
    Xc, Yc, Zc = (data_centerline > 0).nonzero()
    nbLabel = len(X)
    nbLabel_centerline = len(Xc)
    # sort Xc, Yc, and Zc depending on Yc
    cent = [Xc, Yc, Zc]
    indices = range(nbLabel_centerline)
    indices.sort(key = cent[1].__getitem__)
    for i, sublist in enumerate(cent):
        cent[i] = [sublist[j] for j in indices]
    Xc = []
    Yc = []
    Zc = []
    # remove double values
    for i in range(0,len(cent[1])):
        if Yc.count(cent[1][i])==0:
            Xc.append(cent[0][i])
            Yc.append(cent[1][i])
            Zc.append(cent[2][i])
    nbLabel_centerline = len(Xc)
    
    centerline_level = [0 for a in range(nbLabel_centerline)]
    for i in range(0,nbLabel_centerline):
        centerline_level[i] = data_level[Xc[i]][Yc[i]][Zc[i]]
        data_centerline[Xc[i]][Yc[i]][Zc[i]] = 0
    for i in range(0,nbLabel_centerline-1):
        centerline_level[i] = abs(centerline_level[i+1]-centerline_level[i])
    centerline_level[-1] = 0

    C = [i for i, e in enumerate(centerline_level) if e != 0]
    nb_disks = len(C)
    
    if output_level==0:
        for i in range(0,nb_disks):
            data_centerline[Xc[C[i]]][Yc[C[i]]][Zc[C[i]]] = data_level[Xc[C[i]]][Yc[C[i]]][Zc[C[i]]]
    elif output_level==1:
        fo = open(fname_label_output, "wb")
        for i in range(0,nb_disks):
            line = (data_level[Xc[C[i]]][Yc[C[i]]][Zc[C[i]]],Xc[C[i]],Yc[C[i]],Zc[C[i]])
            fo.write("%i %i %i %i\n" %line)
        fo.close()

    return data_centerline

#=======================================================================================================================
def extract_centerline(data,fname_label_output):
    # the Z image is assume to be in second dimension
    X, Y, Z = (data > 0).nonzero()
    cent = [X, Y, Z]
    indices = range(0,len(X))
    indices.sort(key = cent[1].__getitem__)
    for i, sublist in enumerate(cent):
        cent[i] = [sublist[j] for j in indices]
    X = []; Y = []; Z = []
    # remove double values
    for i in range(0,len(cent[1])):
        if Y.count(cent[1][i])==0:
            X.append(cent[0][i])
            Y.append(cent[1][i])
            Z.append(cent[2][i])
    
    fo = open(fname_label_output, "wb")
    for i in range(0,len(X)):
        line = (X[i],Y[i],Z[i])
        fo.write("%i %i %i\n" %line)
    fo.close()

#=======================================================================================================================
def extract_segmentation(data,fname_label_output):
    # the Z image is assume to be in second dimension
    X, Y, Z = (data > 0).nonzero()
    cent = [X, Y, Z]
    indices = range(0,len(X))
    indices.sort(key = cent[1].__getitem__)
    for i, sublist in enumerate(cent):
        cent[i] = [sublist[j] for j in indices]
    X = []; Y = []; Z = []
    # remove double values
    for i in range(0,len(cent[1])):
        X.append(cent[0][i])
        Y.append(cent[1][i])
        Z.append(cent[2][i])
    
    fo = open(fname_label_output, "wb")
    for i in range(0,len(X)):
        line = (X[i],Y[i],Z[i])
        fo.write("%i %i %i\n" %line)
    fo.close()

#=======================================================================================================================
def fraction_volume(data,fname_ref,fname_label_output):
    nx, ny, nz, nt, px, py, pz, pt = sct.get_dimension(fname_ref)
    img_ref = nibabel.load(fname_ref)
    # 3d array for each x y z voxel values for the input nifti image
    data_ref = img_ref.get_data()
    Xr, Yr, Zr = (data_ref > 0).nonzero()
    ref_matrix = [Xr, Yr, Zr]
    indices = range(0,len(Xr))
    indices.sort(key = ref_matrix[1].__getitem__)
    for i, sublist in enumerate(ref_matrix):
        ref_matrix[i] = [sublist[j] for j in indices]
    Xr = []; Yr = []; Zr = []
    for i in range(0,len(ref_matrix[1])):
        Xr.append(ref_matrix[0][i])
        Yr.append(ref_matrix[1][i])
        Zr.append(ref_matrix[2][i])
    
    X, Y, Z = (data > 0).nonzero()
    data_matrix = [X, Y, Z]
    indices = range(0,len(X))
    indices.sort(key = data_matrix[1].__getitem__)
    for i, sublist in enumerate(data_matrix):
        data_matrix[i] = [sublist[j] for j in indices]
    X = []; Y = []; Z = []
    for i in range(0,len(data_matrix[1])):
        X.append(data_matrix[0][i])
        Y.append(data_matrix[1][i])
        Z.append(data_matrix[2][i])

    volume_fraction = []
    for i in range(ny):
        r = []
        for j,p in enumerate(Yr):
            if p == i:
                r.append(j)
        d = []
        for j,p in enumerate(Y):
            if p == i:
                d.append(j)
        volume_ref = 0.0
        for k in range(len(r)):
            value = data_ref[Xr[r[k]]][Yr[r[k]]][Zr[r[k]]]
            if value > 0.5:
                volume_ref = volume_ref + value #suppose 1mm isotropic resolution
        volume_data = 0.0
        for k in range(len(d)):
            value = data[X[d[k]]][Y[d[k]]][Z[d[k]]]
            if value > 0.5:
                volume_data = volume_data + value #suppose 1mm isotropic resolution
        if volume_ref!=0:
            volume_fraction.append(volume_data/volume_ref)
        else:
            volume_fraction.append(0)

    fo = open(fname_label_output, "wb")
    for i in range(ny):
        fo.write("%i %f\n" %(i,volume_fraction[i]))
    fo.close()

#=======================================================================================================================
def write_vertebral_levels(data,fname_vert_level_input):
    fo = open(fname_vert_level_input)
    vertebral_levels = fo.readlines()
    vert = [int(n[2]) for n in [line.strip().split() for line in vertebral_levels]]
    vert.reverse()
    fo.close()
    
    X, Y, Z = (data > 0).nonzero()
    length_points = len(X)
    
    for i in range(0,length_points):
        if Y[i] > vert[0]:
            data[X[i]][Y[i]][Z[i]] = 0
        elif Y[i] < vert[-1]:
            data[X[i]][Y[i]][Z[i]] = 0
        else:
            for k in range(0,len(vert)-1):
                if vert[k+1] < Y[i] <= vert[k]:
                    data[X[i]][Y[i]][Z[i]] = k+1

#=======================================================================================================================
def display_voxel(data):
    # the Z image is assume to be in second dimension
    X, Y, Z = (data > 0).nonzero()
    for k in range(0,len(X)):
        print 'Position=('+str(X[k])+','+str(Y[k])+','+str(Z[k])+') -- Value= '+str(data[X[k],Y[k],Z[k]])

#=======================================================================================================================
# usage
#=======================================================================================================================
def usage():
    print 'USAGE: \n' \
        '  sct_label_utils.py -i <inputdata> -o <outputdata> -c <crossradius>\n' \
        '\n'\
        'MANDATORY ARGUMENTS\n' \
        '  -i           input volume.\n' \
        '  -o           output volume.\n' \
        '  -t           process: cross, remove, display-voxel\n' \
        '  -c           cross radius in mm (default=5mm).\n' \
        '  -r           reference image for label removing' \
        '\n'\
        'OPTIONAL ARGUMENTS\n' \
        '  -h           help. Show this message.\n' \
        '\n'\
        'EXAMPLE:\n' \
        '  sct_label_utils.py -i t2.nii.gz -c 5\n'
    sys.exit(2)
    
    
#=======================================================================================================================
# Start program
#=======================================================================================================================
if __name__ == "__main__":
    # initialize parameters
    param = param()
=======
# TODO: for vert-disc: make it faster! currently the module display-voxel is very long (esp. when ran on PAM50). We can find an alternative approach by sweeping through centerline voxels.
# TODO: label_disc: for top vertebrae, make label at the center of the cord (currently it's at the tip)
# TODO: check if use specified several processes.
# TODO: currently it seems like cross_radius is given in pixel instead of mm

import os
import sys
import argparse
from typing import Sequence

import numpy as np

import spinalcordtoolbox.labels as sct_labels
from spinalcordtoolbox.image import Image, zeros_like
from spinalcordtoolbox.types import Coordinate
from spinalcordtoolbox.reports.qc import generate_qc
from spinalcordtoolbox.utils import Metavar, SmartFormatter, ActionCreateFolder, list_type, init_sct, printv
from spinalcordtoolbox.utils.shell import display_viewer_syntax


def get_parser():
    parser = argparse.ArgumentParser(
        description="Utility functions for label images.",
        formatter_class=SmartFormatter,
        add_help=None,
        prog=os.path.basename(__file__).strip(".py")
    )

    req_group = parser.add_argument_group("\nREQUIRED I/O")
    req_group.add_argument(
        '-i',
        metavar=Metavar.file,
        required=True,
        help="Input image (Required) Example: t2_labels.nii.gz"
    )

    io_group = parser.add_argument_group("\nOPTIONAL I/O")

    io_group.add_argument(
        '-o',
        metavar=Metavar.file,
        default='labels.nii.gz',
        help=("Output image. Note: Only some label utilities create an output image. Example: t2_labels.nii.gz")
    )

    io_group.add_argument(
        '-ilabel',
        metavar=Metavar.file,
        help="File that contain labels that you want to correct. It is possible to add new points with this option. "
             "Use with -create-viewer. Example: t2_labels_auto.nii.gz"
    )

    functions = parser.add_argument_group("\nLABEL FUNCTIONS")
    func_group = functions.add_mutually_exclusive_group(required=True)

    func_group.add_argument(
        '-add',
        metavar=Metavar.int,
        type=int,
        help="Add value to all labels. Value can be negative."
    )

    func_group.add_argument(
        '-create',
        metavar=Metavar.list,
        type=list_type(':', Coordinate),
        help="Create labels in a new image. List labels as: x1,y1,z1,value1:x2,y2,z2,value2. "
             "Example: 12,34,32,1:12,35,33,2"
    )

    func_group.add_argument(
        '-create-add',
        metavar=Metavar.list,
        type=list_type(':', Coordinate),
        help="Same as '-create', but add labels to the input image instead of creating a new image. "
             "Example: 12,34,32,1:12,35,33,2"
    )

    func_group.add_argument(
        '-create-seg',
        metavar=Metavar.list,
        type=list_type(':', list_type(',', int)),
        help="R|Create labels on a cord segmentation (or centerline) image defined by '-i'. Each label should be "
             "specified using the form 'v1,v2' where 'v1' is value of the slice index along the inferior-superior "
             "axis, and 'v2' is the value of the label. Separate each label with ':'. \n"
             "Example: '-create-seg 5,1:14,2:23,3' adds three labels at the axial slices 5, 14, and 23 (starting from the most inferior slice).\n"
             "You can also choose a slice value of '-1' to automatically select the mid-point in the "
             "inferior-superior direction. For example, if you know that the C2-C3 disc is centered in the I-S "
             "direction, then you can enter '-1,3' for that label instead."
    )
    func_group.add_argument(
        '-create-viewer',
        metavar=Metavar.list,
        type=list_type(',', int),
        help="Manually label from a GUI a list of labels IDs, separated with ','. Example: 2,3,4,5"
    )

    func_group.add_argument(
        '-cubic-to-point',
        action="store_true",
        help="Compute the center-of-mass for each label value."
    )

    func_group.add_argument(
        '-disc',
        metavar=Metavar.file,
        help="Create an image with regions labelized depending on values from reference"
    )

    func_group.add_argument(
        '-display',
        action="store_true",
        help="Display all labels (i.e. non-zero values)."
    )
    func_group.add_argument(
        '-increment',
        action="store_true",
        help="Takes all non-zero values, sort them along the inverse z direction, and attributes the values "
             "1, 2, 3, etc."
    )
    func_group.add_argument(
        '-vert-body',
        metavar=Metavar.list,
        type=list_type(',', int),
        help="R|From vertebral labeling, create points that are centered at the mid-vertebral levels. Separate "
             "desired levels with ','. Example: 3,8\n"
             "To get all levels, enter 0."
    )

    func_group.add_argument(
        '-vert-continuous',
        action="store_true",
        help="Convert discrete vertebral labeling to continuous vertebral labeling.",
    )
    func_group.add_argument(
        '-MSE',
        metavar=Metavar.file,
        help="Compute Mean Square Error between labels from input and reference image. Specify reference image here."
    )
    func_group.add_argument(
        '-remove-reference',
        metavar=Metavar.file,
        help="Remove labels from input image (-i) that are not in reference image (specified here)."
    )
    func_group.add_argument(
        '-remove-sym',
        metavar=Metavar.file,
        help="Remove labels from input image (-i) and reference image (specified here) that don't match. You must "
             "provide two output names separated by ','."
    )
    func_group.add_argument(
        '-remove',
        metavar=Metavar.list,
        type=list_type(',', int),
        help="Remove labels of specific value (specified here) from reference image."
    )
    func_group.add_argument(
        '-keep',
        metavar=Metavar.list,
        type=list_type(',', int),
        help="Keep labels of specific value (specified here) from reference image."
    )

    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit."
    )

    optional.add_argument(
        '-msg',
        metavar=Metavar.str,
        help="Display a message to explain the labeling task. Use with -create-viewer"
    )

    optional.add_argument(
        '-v',
        choices=[0, 1, 2],
        default=1,
        metavar=Metavar.int,
        type=int,
        help="Verbose. 0: nothing. 1: basic. 2: extended."
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

    return parser


# MAIN
# ==========================================================================================
def main(args=None):
    parser = get_parser()
    if args:
        arguments = parser.parse_args(args)
    else:
        arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    verbosity = arguments.v
    init_sct(log_level=verbosity, update=True)  # Update log level

    input_filename = arguments.i
    output_fname = arguments.o

    img = Image(input_filename)
    dtype = None

    if arguments.add is not None:
        value = arguments.add
        out = sct_labels.add(img, value)
    elif arguments.create is not None:
        labels = arguments.create
        out = sct_labels.create_labels_empty(img, labels)
    elif arguments.create_add is not None:
        labels = arguments.create_add
        out = sct_labels.create_labels(img, labels)
    elif arguments.create_seg is not None:
        labels = arguments.create_seg
        out = sct_labels.create_labels_along_segmentation(img, labels)
    elif arguments.cubic_to_point:
        out = sct_labels.cubic_to_point(img)
    elif arguments.display:
        display_voxel(img, verbosity)
        return
    elif arguments.increment:
        out = sct_labels.increment_z_inverse(img)
    elif arguments.disc is not None:
        ref = Image(arguments.disc)
        out = sct_labels.labelize_from_discs(img, ref)
    elif arguments.vert_body is not None:
        levels = arguments.vert_body
        if len(levels) == 1 and levels[0] == 0:
            levels = None  # all levels
        out = sct_labels.label_vertebrae(img, levels)
    elif arguments.vert_continuous:
        out = sct_labels.continuous_vertebral_levels(img)
        dtype = 'float32'
    elif arguments.MSE is not None:
        ref = Image(arguments.MSE)
        mse = sct_labels.compute_mean_squared_error(img, ref)
        printv(f"Computed MSE: {mse}")
        return
    elif arguments.remove_reference is not None:
        ref = Image(arguments.remove_reference)
        out = sct_labels.remove_missing_labels(img, ref)
    elif arguments.remove_sym is not None:
        # first pass use img as source
        ref = Image(arguments.remove_reference)
        out = sct_labels.remove_missing_labels(img, ref)

        # second pass use previous pass result as reference
        ref_out = sct_labels.remove_missing_labels(ref, out)
        ref_out.save(path=ref.absolutepath)
    elif arguments.remove is not None:
        labels = arguments.remove
        out = sct_labels.remove_labels_from_image(img, labels)
    elif arguments.keep is not None:
        labels = arguments.keep
        out = sct_labels.remove_other_labels_from_image(img, labels)
    elif arguments.create_viewer is not None:
        msg = "" if arguments.msg is None else f"{arguments.msg}\n"
        if arguments.ilabel is not None:
            input_labels_img = Image(arguments.ilabel)
            out = launch_manual_label_gui(img, input_labels_img, arguments.create_viewer, msg)
        else:
            out = launch_sagittal_viewer(img, arguments.create_viewer, msg)

    printv("Generating output files...")
    out.save(path=output_fname, dtype=dtype)
    display_viewer_syntax([input_filename, output_fname])

    if arguments.qc is not None:
        generate_qc(fname_in1=input_filename, fname_seg=output_fname, args=args,
                    path_qc=os.path.abspath(arguments.qc), dataset=arguments.qc_dataset,
                    subject=arguments.qc_subject, process='sct_label_utils')


def display_voxel(img: Image, verbose: int = 1) -> Sequence[Coordinate]:
    """
    Display all the labels that are contained in the input image.
    :param img: source image
    :param verbose: verbosity level
    """

    coordinates_input = img.getNonZeroCoordinates(sorting='value')
    useful_notation = ''

    for coord in coordinates_input:
        printv('Position=(' + str(coord.x) + ',' + str(coord.y) + ',' + str(coord.z) + ') -- Value= ' + str(coord.value), verbose=verbose)
        if useful_notation:
            useful_notation = useful_notation + ':'
        useful_notation += str(coord)

    printv('All labels (useful syntax):', verbose=verbose)
    printv(useful_notation, verbose=verbose)


def launch_sagittal_viewer(img: Image, labels: Sequence[int], msg: str, previous_points: Sequence[Coordinate] = None, output_img: Image = None) -> Image:
    from spinalcordtoolbox.gui import base
    from spinalcordtoolbox.gui.sagittal import launch_sagittal_dialog
    params = base.AnatomicalParams()
    params.vertebraes = labels
    params.input_file_name = img.absolutepath

    if output_img is not None:
        params.output_file_name = output_img.absolutepath
    else:
        params.output_file_name = img.absolutepath

    params.subtitle = msg

    if previous_points is not None:
        params.message_warn = 'Please select the label you want to add \nor correct in the list below before clicking \non the image'

    out = zeros_like(img)
    out.absolutepath = params.output_file_name
    launch_sagittal_dialog(img, out, params, previous_points)

    return out


def launch_manual_label_gui(img: Image, input_labels_img: Image, labels: Sequence[int], msg):
    # the input image is reoriented to 'SAL' when open by the GUI
    input_labels_img.change_orientation('SAL')
    mid = int(np.round(input_labels_img.data.shape[2] / 2))
    previous_points = input_labels_img.getNonZeroCoordinates()

    # boolean used to mark first element to initiate the list.
    first = True

    previous_label = None

    for i in range(len(previous_points)):
        if int(previous_points[i].value) in labels:
            pass
        else:
            labels.append(int(previous_points[i].value))
        if first:
            points = np.array([previous_points[i]. x, previous_points[i].y, previous_points[i].z, previous_points[i].value])
            points = np.reshape(points, (1, 4))
            previous_label = points
            first = False
        else:
            points = np.array([previous_points[i].x, previous_points[i].y, previous_points[i].z, previous_points[i].value])
            points = np.reshape(points, (1, 4))
            previous_label = np.append(previous_label, points, axis=0)
        labels.sort()

    # check if variable was created which means the file was not empty and contains some points asked in labels
    if previous_label is not None:
        # project onto mid sagittal plane
        for i in range(len(previous_label)):
            previous_label[i][2] = mid

    out = launch_sagittal_viewer(img, labels, msg, previous_points=previous_label)

    return out


if __name__ == "__main__":
    init_sct()
>>>>>>> upstream/master
    # call main function
    main()
