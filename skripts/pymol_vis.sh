#!/bin/bash

usage() { echo "Visualize ACV restraints in PyMOL
Usage: pymol_restraints.sh -c <structure file (.gro/.pdb)> -x <xtc trajectory> (optional) -v <pymol visualization file> (optional), -s <start,stop,stride (default 1,-1,1)>" 1>&2; exit 1; }
invalidOpt() { echo "Invalid option: -$OPTARG" 1>&2; exit 1; }
missingArg() { echo "Option -$OPTARG requires an argument" 1>&2; exit 1; }
cleanup() { if ls -f $1/\#* 1> /dev/null 2>&1 ; then rm $1/\#* ; fi ; }

while getopts ":c:x:v:s:h" opt; do
    case $opt in
        c) 
            gro_file=$OPTARG
            gro_base=${gro_file##*/}
            gro_path=${gro_file%$gro_base}
            gro_name=${gro_base%.*}
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

echo "$gro_path"
echo "cmd.load('$gro_file')" > tmp_vis.py

pwd

# start,stop,stride specified
if [ ! -z "$sss" ]; then
    start=`echo $sss | cut -d',' -f1`
    stop=`echo $sss | cut -d',' -f2`
    stride=`echo $sss | cut -d',' -f3`
else
    start=1
    stop=-1
    stride=1
fi

# xtc argument specified
if [ ! -z "$xtc_file" ]; then
    echo "cmd.load_traj('$xtc_file', '$gro_name', 0, start=$start, stop=$stop, interval=$stride)" >> tmp_vis.py
fi

# PyMOL file specified
if [ ! -z "$vis_file" ]; then
    cat "$vis_file" >> tmp_vis.py
fi

pymol tmp_vis.py || (PyMOLWin.exe tmp_vis.py & sleep 10 ) || echo "-> Error: PyMOL can not be found. Please add it to your path."
rm tmp_vis.py
