# neupy
Neupy is a collection of scripts for use in sorting nuclear data and generating new data. 
It's main goal is to generate neutrino flux and spectrum for use in toy reactor models and for 
experimental data. There is however other things neupy can do.

Neupy can be seperated into 2 parts, the parsing of nuclear data and generation of sorted data
and then the neutrino emission part. These two parts should be run seperately to save computing
time. 

This script also has a few dependencies including pandas, numpy, matplotlib, math, os, scipy.

# nubaseSorter

Starting with it's parsing of nuclear data, you can access the raw nubase and fission yeild data
using

fission yeild data = nubaseSorter.fiss_data
nuclear data = nubaseSorter.nubase

The original databases should be kept in a folder ./databases/ with file headings fyu235thermal.txt 
and Parsed_nubase2016.xlsx for the fission data and nubase data respectively. 

~~~
import nubaseSorter as ns
fissionData = ns.fiss_data           
nubaseData = ns.nubase               
fiss_column_Z = ns.Z                 
nubase_column_Z = ns.nuZ              
~~~