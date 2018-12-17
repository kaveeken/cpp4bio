#!/bin/bash
#19ff00 groen 
#ff0000 rood 

#list of images
imageR=$(ls -1 *jpg | head -n 21)
imageS=$(ls -1 *jpg | tail -n 21)

# convert
for file in $imageR; do
	convert file +level-colors blue,white ${file}
done

for file in $imageS; do
	convert file +level-colors red,white ${file}
done
