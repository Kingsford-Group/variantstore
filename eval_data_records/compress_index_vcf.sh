#!/bin/bash

file=$1
echo $file
#gunzip $file

gz=".gz"
gzfile=$1$gz

if [ -s "$file" ]
then
        bgzip -f $file
        tabix -f $gzfile
      fi
