#!bin/bash

# Command line input
usage="Calculate the dependence of transmission and reflection rate on momentum expectations

$(basename "$0") [-h] exe input

positional arguments:
  exe            the path to OneDimDVR-scatter.exe
  input          the path to input directory (containing all input files except OneDimDVR-scatter.in)

optional arguments:
  -h             show this help message and exit"

while getopts ':hq:' option; do
  case "$option" in
     h) echo "$usage"
        exit;;
     :) printf "missing argument for -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit;;
  esac
done
shift $((OPTIND - 1))

exe=`realpath $1`
input=`realpath $2`

# Go through each momentum
for i in 4 4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.5 8 8.5 9 9.5 10 11 12 13 14 15 16 17 18 19 20; do
    # Create work directory
    if [ -d $i ]; then
        rm -r $i
    fi
    mkdir $i
    cd $i
    # Create input
    cp -r $input/* .
    python ~/Software/Mine/OneDimDVR/example/super-exchange/generate_input.py $i
    # Run
    $exe > log
    # Clean up
    cd ..
done

# Collect
echo K$'\t'T1$'\t'R1$'\t'T2$'\t'R2$'\t'T3$'\t'R3 > log
for i in 4 4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.5 8 8.5 9 9.5 10 11 12 13 14 15 16 17 18 19 20; do
    echo -n $i $'\t' >> log
    readarray -t lines < $i/transmission.txt
    for ((state=1; state<4; state++)); do
        IFS=' ' read -r -a strs <<< "${lines[$state]}"
        for ((j=1; j<${#strs[@]}; j++)); do
            echo -n ${strs[$j]}$'\t' >> log
        done
    done
    echo >> log
done

# Clean up
for i in 4 4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.5 8 8.5 9 9.5 10 11 12 13 14 15 16 17 18 19 20; do
    rm -r $i
done
