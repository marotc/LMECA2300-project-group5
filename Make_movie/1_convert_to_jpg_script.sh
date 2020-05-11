#!/bin/bash
for f in ../screenshots/*; do 
	echo `basename $f`;
	convert -quality 100 $f `basename $f`.jpg;
done
#do convert quality 100 $f $f.jpg; done
