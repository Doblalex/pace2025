#!/usr/bin/env bash
cd $3
file=$(basename $2)
(time $1 < $2) &> "./$4/${file}.out"