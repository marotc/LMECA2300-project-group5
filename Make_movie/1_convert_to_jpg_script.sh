#!/bin/bash
i=0
f=0
while [ $f -le 44500 ]
do
	convert -quality 100 myEllipse_$f  yo_$i.jpg;
	i=$(( i + 1 ))
	f=$(( f + 500 ))
	echo $f
done
# convert -delay 20 -quality 100 yo_*.jpg movie.mpg
