#!/bin/bash

FILE="missing_run16_1.list"

runnumbers=()

while IFS= read -r line; do
   numbers+=("$line")
done < "$FILE"

for number in "${numbers[@]}"; do
    mkdir -p /star/data01/pwg/aspovarov/for_Maxim/Polarization/200GeV/log/$number
    mkdir -p /star/data01/pwg/aspovarov/for_Maxim/Polarization/200GeV/data/$number

    echo "${number} was created"
done
