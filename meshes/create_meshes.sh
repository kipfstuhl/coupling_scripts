#!/bin/bash

PYTHON=/usr/bin/python2.7
GMSH=/Applications/Gmsh.app/Contents/MacOS/gmsh
THIS=$(pwd)

rm *.msh
rm -r -f refinement*

$GMSH -3 bifurcation_smooth_distorted_fine.geo

ref=1
mult=1
while [ $ref -le 12 ]
do
	DIR=refinement$ref
	mkdir $DIR
	cd $DIR

    # inflow
    FILE=inflow_distorted
    cp ../$FILE.geo .
    meshsize=$(grep -n "lc = " $FILE.geo)
    meshsize=${meshsize:7:5}
    echo $meshsize
    mymult=$(echo "$meshsize / $mult" | bc -l)
    # linux
    # sed -i "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    # mac
    sed -i '.original' "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    $GMSH -3 $FILE.geo

    # outflow 1
    FILE=outflow1_distorted
    cp ../$FILE.geo .
    meshsize=$(grep -n "lc = " $FILE.geo)
    meshsize=${meshsize:7:4}
    echo $meshsize
    mymult=$(echo "$meshsize / $mult" | bc -l)
    # linux
    # sed -i "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    # mac
    sed -i '.original' "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    $GMSH -3 $FILE.geo

    # outflow 2
    FILE=outflow2_distorted
    cp ../$FILE.geo .
    meshsize=$(grep -n "lc = " $FILE.geo)
    meshsize=${meshsize:7:5}
    echo $meshsize
    mymult=$(echo "$meshsize / $mult" | bc -l)
    # linux
    # sed -i "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    # mac
    sed -i '.original' "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    $GMSH -3 $FILE.geo

	# smooth bifurcation
    FILE=bifurcation_smooth_distorted
    cp ../$FILE.geo .
    meshsize=$(grep -n "lc = " $FILE.geo)
    meshsize=${meshsize:7:3}
    echo $meshsize
    mymult=$(echo "$meshsize / $mult" | bc -l)
    # linux
    # sed -i "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    # mac
    sed -i '.original' "s/lc = $meshsize;/lc = $mymult;/g" $FILE.geo
    $GMSH -3 $FILE.geo

    rm *.original

	cd $THIS
	mult=$(echo "$mult * 1.189207115002721" | bc)
	((ref++))
done
