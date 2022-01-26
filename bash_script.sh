#!bin/bash

for f in ./$1/*
do
    python3 clusters_on_time.py $f
done
