#!/bin/bash

FILE="missing_run16_1.list"

runnumbers=()

while IFS= read -r line; do
   numbers+=("$line")
done < "$FILE"

for number in "${numbers[@]:1}"; do

	star-submit-beta submiter_${number}.xml
done

