#!/bin/bash


#-------------------------------------------------------
# This scrips pulls the newest version of RAMSES in the
# "clean" repo and copies the new stuff into the
# masterarbeit repo.
#-------------------------------------------------------


srcdir=$OR
thisdir=$PWD

echo "========================================================"
echo "updating ramses files to newest version."
echo "skipping files in patch/mergertree and patch/unbinding"
echo "========================================================"


cd $srcdir
git checkout master
git pull upstream master 


echo
echo "========================================================"
echo "copying new files to " $ma/project/ramses
echo "skipping files in patch/mergertree and patch/unbinding"
echo "========================================================"

for dir in amr aton bin doc hydro io mhd namelist pario pm poisson rhd rt utils patch/challenge patch/hydro patch/init patch/insta   patch/mhd patch/mom2 patch/phantom patch/rt patch/test_suite;
    do 
        cp -rvu $srcdir/$dir/* $MA/project/ramses/$dir;
    done

for file in README README.md;
    do
        cp -rvu $srcdir/$file $MA/project/ramses/$file;
    done




echo
echo "========================================================"
echo "pushing new commits to my ramses repo"
echo "========================================================"
# Origin: my repo
# upstream: Romain's repo
git push origin master


echo "Done updating files"







cd $thisdir


exit 0
