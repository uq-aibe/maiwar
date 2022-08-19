#!/bin/bash

#script to parse templates to replace text ampl model

#sed= stream editor
#$1 is the first parameter on the command line
cat template.mod \
  | sed -e "s/__INVSEC__/$1/"
