#!/bin/bash

usage() { echo "Visualize ACV restraints in VMD
Usage: vmd_restraints.sh -c <structure file (.gro/.pdb)> -x <xtc trajectory> (optional) -v <vmd visualization file> (optional) -s <start,stop,stride (default 1,-1,1)" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

while getopts ":c:x:v:s:h" opt; do
    case $opt in
        c) 
            gro_file=$OPTARG
            ;;
        x)
            xtc_file=$OPTARG
            ;;
        v)
            vis_file=$OPTARG
            ;;
        s)
            sss=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            invalidOpt
            ;;
        :)
            missingArg
            ;;
        *)
            usage
            ;;
    esac
done


# no cmd line arguments given
if [ -z "$gro_file" ]; then
    usage
fi

echo "mol new $gro_file first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" > tmp_vis.vmd


# start,stop,stride specified
if [ ! -z "$sss" ]; then
    start=`echo "$sss" | cut -d',' -f1`
    stop=`echo "$sss" | cut -d',' -f2`
    stride=`echo "$sss" | cut -d',' -f3`
else
    start=1
    stop=-1
    stride=1
fi


# xtc argument specified
if [ ! -z "$xtc_file" ]; then
    echo "mol addfile $xtc_file type xtc first $start last $stop step $stride filebonds 1 autobonds 1 waitfor all" >> tmp_vis.vmd
fi

# vmd file specified
if [ ! -z "$vis_file" ]; then
    cat "$vis_file" >> tmp_vis.vmd
fi

vmd.exe -e tmp_vis.vmd || vmd -e tmp_vis.vmd || echo "-> Error: VMD can not be found. Please add it to your path."
rm tmp_vis.vmd
