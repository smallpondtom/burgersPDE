#! /bin/bash
#
matlab -nodisplay -nosplash -nodesktop -batch \
  "run('burgers_solution_test.m');exit;" > burgers_solution_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
