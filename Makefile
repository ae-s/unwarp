
# These are rules which don't correspond to an actual file
.PHONY: calibrate flatten

export MAGICK_TMPDIR=/tmp/3

# We maybe don't want to keep all these huge bitmaps lying around
#.INTERMEDIATE: out/%.ppm orig/%.ppm

# There is a tool ... which make knows how to compile ...
ppmunwarp: ppmunwarp.cc

# Calibration requires the creation of an even and an odd deformation
# map.  This indirectly requires the files `even.jpg' and `odd.jpg'.
calibrate: cal-even.deformation cal-odd.deformation ppmunwarp

# To create a deformation map, look for an image of the same name.
# Also creates a [thing]-check.ppm file, in case you're interested in
# what the machine sees.
%.deformation: %.ppm
	/idk/astrid/unwarp/ppmunwarp -cc $< -d $@ -m $*-check.ppm

# We can create a ppm from an otherwise identically named jpeg.
%.ppm: %.jpg
	convert $< $@

# We can create png images too, it's not that complicated
#%.flat.png: %.flat.ppm
#	convert $< $@


# Here's the ruleset to unwarp a single image.
#
# Determining whether a number is even or odd, in the year 2017,
# requires spawning two (three?) new processes and using a regular
# expression.

%.in.ppm: in/%.jpg
	convert $< $@

%.dew.ppm: %.in.ppm
	/idk/astrid/unwarp/ppmunwarp -i $< -o $@ \
		-d cal-`echo $* | sed -e 's/.*\(.\)/\1/' -e 's/[02468]/even/' -e 's/[13579]/odd/' `.deformation

%.out.ppm: %.dew.ppm
	pnmflip `echo $* | sed -e 's/.*\(.\)/\1/' -e 's/[02468]/-ccw/' -e 's/[13579]/-cw/' ` $< > $@

out/%.png: %.out.ppm
	convert $< $@

flatten: calibrate $(addprefix out/,$(notdir $(addsuffix .png,$(basename $(wildcard in/[0-9]*.jpg)))))

all: calibrate flatten
