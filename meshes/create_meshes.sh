#!/bin/bash

GMSH="/Applications/Gmsh.app/Contents/MacOS/gmsh"

$GMSH -3 inflow_symmetric.geo
$GMSH -3 outflow1_symmetric.geo
$GMSH -3 outflow2_symmetric.geo
