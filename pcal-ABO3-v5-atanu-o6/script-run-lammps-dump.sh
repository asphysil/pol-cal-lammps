#!/bin/bash

echo "*****Please change the path of the executable*****"

for i in $(seq 1 1); do
rm dump.xyz
if test -f ../dump.xyz-${i}
then
cp -vr ../dump.xyz-${i} dump.xyz
nl=`grep "ITEM: BOX BOUNDS pp pp pp" dump.xyz | wc -l`
echo ${nl}
# Run Fortran program with input redirection
/global/homes/a/atanu/fe-tool-box/pol-cal-lammps/pcal-ABO3-v5-atanu-o6/exeO6O12.x <<EOF
${nl}
EOF
fi
done
