#!/bin/bash


if [ $# -ne 1 ]; then
    echo 'Supply commit message, please!'
    exit
fi

find . -name *~ -exec rm {} \;

git init
git add -A
git commit -m "$1"
#git remote add origin https://github.com/tagoreas/Lensing-code.git
git push -u origin master

