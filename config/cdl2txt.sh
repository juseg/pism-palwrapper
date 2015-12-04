#!/bin/bash

# Convert cdl file to format needed by palwrapper

sed -r \
    -e '1,3d' \
    -e '$d' \
    -e 's|^[[:space:]]*//|\n//|' \
    -e 's|^[[:space:]]*pism_overrides:(.*) = (.*) ;$|\1: \2|' \
    -e 's|"||g' \
    $1
