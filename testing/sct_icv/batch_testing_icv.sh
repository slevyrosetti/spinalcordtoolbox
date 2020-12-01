#!/bin/bash
#
<<<<<<< HEAD
# testing sct_icv.py

sct_icv.py -i data/errsm_20_t1.nii.gz -c t1 -o tmp.output_sienax -d sienax
=======
# testing sct_icv

sct_icv.py data/errsm_20_t1.nii.gz -c t1 -o tmp.output_sienax -d sienax
>>>>>>> upstream/master

sct_icv.py -i data/errsm_20_t1.nii.gz -c t1 -o tmp.output_rbm -d rbm

