
# These are rules which don't correspond to an actual file
.PHONY: calibrate flatten cleanup

# We maybe don't want to keep all these huge bitmaps lying around
.INTERMEDIATE: out/%.ppm orig/%.ppm

# There is a tool ... which make knows how to compile ...
ppmunwarp: ppmunwarp.cc

# Calibration requires the creation of an even and an odd deformation
# map.  This indirectly requires the files `even.jpg' and `odd.jpg'.
calibrate: even.deformation odd.deformation ppmunwarp

# To create a deformation map, look for an image of the same name.
# Also creates a [thing]-check.ppm file, in case you're interested in
# what the machine sees.
%.deformation: %.ppm
	./ppmunwarp -cc $< -d $@ -m $*-check.ppm

# We can create a ppm from a jpeg.
%.ppm: %.jpg
	convert $< $@

# We can create png images too, it's not that complicated
%.png: %.ppm
	convert $< $@

# Determining whether a number is even or odd, in the year 2017,
# requires spawning two (three?) new processes and using a regular
# expression.
NUMMAP = $(shell echo $(1) | sed -e 's/.*(.)$/\1/' -e 's/[02468]/even/' -e 's/[13579]/odd/' )

# Here's a rule to unwarp a single image.
out/%.ppm: orig/%.ppm
	./ppmunwarp -i $< -o $@ -d $(call NUMMAP $*).deformation

RENAME = out/$(basename $(1)).png

flatten: calibrate $(call rename orig/*)

all: calibrate flatten cleanup
