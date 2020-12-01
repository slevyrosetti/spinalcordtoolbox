#!/bin/sh

<<<<<<< HEAD
matlab_exec=/Applications/MATLAB_R2012a.app/bin/matlab
=======
matlab_exec=/Applications/MATLAB_R2014b.app/bin/matlab
>>>>>>> upstream/master
X="${1}(${2}); exit"
echo ${X} > matlab_command.m
cat matlab_command.m
${matlab_exec} -nodesktop -nosplash < matlab_command.m
