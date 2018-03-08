#!/bin/bash

file=$1
folder=$2

if [ -d "$folder" ]; then
    ./BetaScatt "$file" "$folder" > "${file}.stdout"
else
    ./BetaScatt "$file" "$PWD" > "${file}.stdout"
fi


