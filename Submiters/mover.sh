#!/bin/bash

INPUT_FILE="run16_1_runnum.list"
OUTPUT_FILE="run16_2_runnum.list"

START_LINE=1381

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "$INPUT_FILE dont found"
    exit 1
fi

tail -n +"$START_LINE" "$INPUT_FILE" > "$OUTPUT_FILE"

echo "String from $INPUT_FILE moved to $OUTPUT_FILE"
