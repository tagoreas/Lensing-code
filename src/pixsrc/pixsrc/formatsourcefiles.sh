#!/bin/bash

for file in $(find -maxdepth 1 -name '*.?pp' -o -name '*.c' -o -name '*.cu'); do
    emacs -batch $file -l emacs-format-file -f emacs-format-function
done
