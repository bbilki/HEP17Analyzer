#!/bin/sh

if [ ! -d ../Plots/$1 ]
then
	mkdir ../Plots/$1
fi

if [ ! -d ../Histos/$1 ]
then
	mkdir ../Histos/$1
fi

./plotter $2 $1
wait


