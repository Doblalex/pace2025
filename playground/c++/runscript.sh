#!/usr/bin/env bash
cd $3
file=$(basename $2)
$1 < $2 > "./results/${file}.out"