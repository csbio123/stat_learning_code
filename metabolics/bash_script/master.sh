#!/bin/sh


files=`ls /users/spjtcoi/git/stat_learning_code/metabolics/bash_script/combat_scripts/*.sh`
for i in $files;
do
 qsub $i
done


