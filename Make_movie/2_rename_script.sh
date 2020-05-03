#!/bin/bash
f=0
while [ $f -le 9 ]
do
	mv yo_$f.jpg yo_0$f.jpg
	f=$(( f + 1 ))
	echo $f
done
