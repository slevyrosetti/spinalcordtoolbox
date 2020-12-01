#!/usr/bin/env python
<<<<<<< HEAD
#
# Create package with appropriate version number.
#
# Author: Julien Cohen-Adad, Benjamin De Leener
=======
# -*- coding: utf-8 -*-
#
# Create package with appropriate version number.
#
# Author: Julien Cohen-Adad, Benjamin De Leener, P-O Quirion
>>>>>>> upstream/master
#

# TODO: remove quick fix with folder_sct_temp

<<<<<<< HEAD
import os
import getopt
import sys
from numpy import loadtxt
=======
DEBUG = False

import getopt
import os
import platform
import shutil
import sys
import tarfile
import tempfile
import urllib2


>>>>>>> upstream/master
sys.path.append('../scripts')
import sct_utils as sct

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
        '  Create a package of the Spinal Cord Toolbox.\n' \
<<<<<<< HEAD
        '\n' \
        'USAGE\n' \
        '  '+os.path.basename(__file__)+'\n' \
        '\n' \
        'MANDATORY ARGUMENTS\n' \
        '  -s <OS name>      name of the OS {osx,debian}.\n' \
        '\n'\
        'EXAMPLE:\n' \
        '  create_package.py -s debian\n'
    sys.exit(2)

listOS = ['osx', 'debian']
OSname = ''
# Check input param
try:
    opts, args = getopt.getopt(sys.argv[1:],'hs:')
=======
        '\n'\
    # sys.exit(2)


if platform.mac_ver()[0]:
    local_os = 'osx'
elif platform.linux_distribution()[0]:
    local_os = 'linux'

listOS = ["linux", "osx"]
deb_fsl = []
tar_fsl = []



# OSname = ''
# Check input param
try:
    opts, args = getopt.getopt(sys.argv[1:], 'h')
>>>>>>> upstream/master
except getopt.GetoptError as err:
    print str(err)
    usage()
for opt, arg in opts:
    if opt == '-h':
        usage()
<<<<<<< HEAD
    elif opt in ('-s'):
        OSname = str(arg)

if OSname not in listOS:
    print 'ERROR: OS name should be one of the following: '+'[%s]' % ', '.join(map(str,listOS))
    usage()
=======
    # elif opt in ('-s'):
    #     OSname = str(arg)

# if OSname not in listOS:
#     print 'ERROR: OS name should be one of the following: '+'[%s]' % ', '.join(map(str,listOS))
#     usage()
>>>>>>> upstream/master

# get version
with open ("../version.txt", "r") as myfile:
    version = myfile.read().replace('\n', '')

# create output folder
<<<<<<< HEAD
folder_sct = '../spinalcordtoolbox_v'+version+'_'+OSname+'/'
if os.path.exists(folder_sct):
    sct.run('rm -rf '+folder_sct)
sct.run('mkdir '+folder_sct)

# copy folders
sct.run('mkdir '+folder_sct+'spinalcordtoolbox')
sct.run('cp installer.sh '+folder_sct)
sct.run('cp ../README.md '+folder_sct+'spinalcordtoolbox/')
sct.run('cp ../LICENSE '+folder_sct+'spinalcordtoolbox/')
sct.run('cp ../version.txt '+folder_sct+'spinalcordtoolbox/')
sct.run('cp -r ../flirtsch '+folder_sct+'spinalcordtoolbox/')
sct.run('cp -r ../scripts '+folder_sct+'spinalcordtoolbox/')
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/bin')
if OSname == 'osx':
    sct.run('cp -r ../bin/osx/* '+folder_sct+'spinalcordtoolbox/bin/')
    sct.run('cp -r osx/ants/* '+folder_sct+'spinalcordtoolbox/bin/')
    sct.run('cp -r osx/c3d/* '+folder_sct+'spinalcordtoolbox/bin/')
elif OSname == 'debian':
    sct.run('cp -r ../bin/debian/* '+folder_sct+'spinalcordtoolbox/bin/')
    sct.run('cp -r debian/ants/* '+folder_sct+'spinalcordtoolbox/bin/')
    sct.run('cp -r debian/c3d/* '+folder_sct+'spinalcordtoolbox/bin/')

# copy colormap
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/data/')
sct.run('cp -rf ../data/colormap '+folder_sct+'spinalcordtoolbox/data/')

# copy atlas
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/data/atlas')
sct.run('ls ../data/atlas/vol*.* | while read F; do cp $F '+folder_sct+'spinalcordtoolbox/data/atlas/; done')
sct.run('cp ../data/atlas/list.txt '+folder_sct+'spinalcordtoolbox/data/atlas/')

# copy spinal_level
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/data/spinal_level')
sct.run('cp -rf ../data/spinal_level '+folder_sct+'spinalcordtoolbox/data/')

# copy template
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/data/template')
sct.run('cp ../data/template/MNI-Poly-AMU_T2.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/landmarks_center.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/mask_gaussian_templatespace_sigma20.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/MNI-Poly-AMU_cord.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/MNI-Poly-AMU_CSF.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/MNI-Poly-AMU_GM.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/MNI-Poly-AMU_level.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/MNI-Poly-AMU_seg.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/MNI-Poly-AMU_WM.nii.gz '+folder_sct+'spinalcordtoolbox/data/template/')
sct.run('cp ../data/template/version.txt '+folder_sct+'spinalcordtoolbox/data/template/')

# testing
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing')

# testing - sct_segmentation_propagation
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing/sct_segmentation_propagation')
sct.run('cp ../testing/sct_segmentation_propagation/test_sct_segmentation_propagation.sh '+folder_sct+'spinalcordtoolbox/testing/sct_segmentation_propagation/')
sct.run('cp -r ../testing/sct_segmentation_propagation/snapshots '+folder_sct+'spinalcordtoolbox/testing/sct_segmentation_propagation/')

# testing - sct_register_to_template
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing/sct_register_to_template')
sct.run('cp ../testing/sct_register_to_template/test_sct_register_to_template.sh '+folder_sct+'spinalcordtoolbox/testing/sct_register_to_template/')

# testing - sct_register_multimodal
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing/sct_register_multimodal')
sct.run('cp ../testing/sct_register_multimodal/test_sct_register_multimodal.sh '+folder_sct+'spinalcordtoolbox/testing/sct_register_multimodal/')
sct.run('cp -r ../testing/sct_register_multimodal/snapshots '+folder_sct+'spinalcordtoolbox/testing/sct_register_multimodal/')
sct.run('cp -r ../testing/sct_register_multimodal/check_integrity '+folder_sct+'spinalcordtoolbox/testing/sct_register_multimodal/')

# testing - sct_warp_atlas2metric
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing/sct_warp_atlas2metric')
sct.run('cp ../testing/sct_warp_atlas2metric/test_sct_warp_atlas2metric.sh '+folder_sct+'spinalcordtoolbox/testing/sct_warp_atlas2metric/')
sct.run('cp -r ../testing/sct_warp_atlas2metric/snapshots '+folder_sct+'spinalcordtoolbox/testing/sct_warp_atlas2metric/')

# testing - sct_estimate_MAP_tracts
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing/sct_estimate_MAP_tracts')
sct.run('cp ../testing/sct_estimate_MAP_tracts/test_sct_estimate_MAP_tracts.sh '+folder_sct+'spinalcordtoolbox/testing/sct_estimate_MAP_tracts/')

# testing - data
sct.run('mkdir '+folder_sct+'spinalcordtoolbox/testing/data')
sct.run('cp -r ../testing/data/errsm_23 '+folder_sct+'spinalcordtoolbox/testing/data/')

# remove .DS_Store files
sct.run('find '+folder_sct+' -type f -name .DS_Store -delete')

# remove AppleDouble files - doesn't work on Linux
if OSname == 'osx':
    sct.run('find '+folder_sct+' -type d | xargs dot_clean -m')

# remove pycharm files
sct.run('rm '+folder_sct+'spinalcordtoolbox/scripts/*.pyc')
if os.path.exists(folder_sct+'spinalcordtoolbox/scripts/.idea'):
    sct.run('rm -rf '+folder_sct+'spinalcordtoolbox/scripts/.idea')

# go to parent directory to be able to tar without de-tarring error in mac OSX
os.chdir('../')

# compress folder
folder_sct_temp = 'spinalcordtoolbox_v'+version+'_'+OSname+'/'

sct.run('tar -cvzf spinalcordtoolbox_v'+version+'_'+OSname+'.tar.gz '+folder_sct_temp)

print "done!\n"


=======


folder_sct = '{0}/sct_v{1}'.format(tempfile.gettempdir(), version)
print("output dir creation {0}".format(folder_sct))
if os.path.exists(folder_sct) and not DEBUG:
    shutil.rmtree(folder_sct)


offline = "{0}/offline".format(folder_sct)

try:
    os.mkdir(folder_sct)
except OSError:
    pass
try:
    os.mkdir(offline)
except OSError:
    pass


def download(url, dest_dir):

    if os.path.isfile(dest_dir) and DEBUG:
        return

    file_name = url.split('/')[-1]
    u = urllib2.urlopen(url)
    f = open(dest_dir, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print "Downloading: %s Bytes: %s" % (file_name, file_size)

    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8)*(len(status)+1)
        print status,

    f.close()

for OS in listOS:

    download("https://dl.dropboxusercontent.com/u/20592661/sct/bin/{OS}/sct_binaries.tar.gz".format(OS=OS),
             "{sct}/sct_binaries_{OS}.tar.gz".format(sct=offline, OS=OS))
    download("https://dl.dropboxusercontent.com/u/20592661/sct/miniconda/{OS}/miniconda.sh".format(OS=OS),
             "{sct}/miniconda_{OS}.sh".format(sct=offline, OS=OS))
    download("https://dl.dropboxusercontent.com/u/20592661/sct/wheel/{OS}/wheels.tar.gz".format(OS=OS),
             "{sct}/wheels_{OS}.tar.gz".format(sct=offline, OS=OS))

download("https://dl.dropboxusercontent.com/u/20592661/sct/MNI-Poly-AMU.tar.gz",
         "{sct}/data.tar.gz".format(sct=offline))
download("https://pypi.python.org/packages/b7/12/71ff08d828319e80ad575762a5e0d7d07db"
         "72a51e3d648a388ef1aa2e0ae/nibabel-2.0.1-py2-none-any.whl#md5=17c20bde23d9c7e5cbcf48df47c27d63",
         "{sct}/nibabel-2.0.1-py2-none-any.whl".format(sct=offline))

print("sct web dependent package downloaded")


# all needed package
src_pkgs_dir = '../python/pkgs'
src_pkgs = os.listdir(src_pkgs_dir)
dest_pkgs = '{0}/pkgs_{1}'.format(offline, local_os)
try:
    os.mkdir(dest_pkgs)
except OSError:
    pass
for p in src_pkgs:
    if p.endswith('.bz2'):
        src_bz = '{0}/{1}'.format(src_pkgs_dir, p)
        shutil.copy2(src_bz, dest_pkgs)

print('Python package copied')

# copy following folders and file

cp_list = [os.path.abspath(e.strip()) for e in ['../install_sct', '../README.md', '../LICENSE',
                                                '../version.txt', '../commit.txt ', '../batch_processing.sh',
                                                '../batch_processing.sh', '../scripts', '../install', '../testing']]
for elem in cp_list:
    dest = "{0}/{1}".format(folder_sct, os.path.basename(elem))
    if os.path.isfile(dest):
        os.remove(dest)
    elif os.path.isdir(dest):
        shutil.rmtree(dest)
        # continue
    if os.path.isfile(elem):
        shutil.copyfile(elem, dest)
    elif os.path.isdir(elem):
        shutil.copytree(elem, dest)

# remove .DS_Store files
sct.run('find '+folder_sct+' -type f -name .DS_Store -delete')
# remove Pycharm-related files
sct.run('find '+folder_sct+' -type f -name *.pyc -delete')
sct.run('find '+folder_sct+' -type f -name *.pyo -delete')
sct.run('find '+folder_sct+' -type f -name *.idea -delete')

# remove AppleDouble files - doesn't work on Linux
# if OSname == 'osx':
#     sct.run('find '+folder_sct+' -type d | xargs dot_clean -m')


def make_tarfile(output_filename, source_dir):
    print(" building {0} file".format(output_filename))
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


make_tarfile("sct_v{0}_offline_{1}.tar.gz".format(version, local_os), folder_sct)

if not DEBUG:
    shutil.rmtree(folder_sct)

print "done!\n"
>>>>>>> upstream/master
