#!/bin/bash

for((i=1;i<42;i++));
do
python rewriteWMsmall.py;
qsub runWMsmall.sh;
done