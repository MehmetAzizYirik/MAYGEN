#!/bin/bash
MAYGEN="java -jar /Users/steinbeck/develop/MAYGEN/target/MAYGEN-jar-with-dependencies.jar -c -f "
while IFS= read -r line; do
    $MAYGEN $line
done < "$1"
