#!/bin/bash
#let's run a test
steve='13'
steven=`printf '%3.3d' $steve`
echo $steven
balrog=$(($steven + $1))
balrogn=`printf '%3.3d' $balrog`
echo $balrog
echo $balrogn
