#!/bin/bash
#Usage: Data_transform.sh <directory> <Group_ID>

cd $1
sed 's:cat/::g' $1/$2 > traj_filelist
ls traj_file* | sed 's/traj_file_//g' > traj_filelist
echo "File list made"
cat traj_filelist | while read i; do 
echo $i
tail -n40000 $i | awk -v FS=" " '$3>=2246400{print $0}' | awk -v FS=" " '$8==0{print $0}' > $i.txt
done
