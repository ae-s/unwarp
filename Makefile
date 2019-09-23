
# These are rules which don't correspond to an actual file
.PHONY: calibrate cal-whitepoint flatten all

export MAGICK_TMPDIR=/tmp/3
TOOLDIR=/idk/astrid/unwarp

# We maybe don't want to keep all these huge bitmaps lying around
#.INTERMEDIATE: out/%.ppm orig/%.ppm

# There is a tool ... which make knows how to compile ...
$(TOOLDIR)/ppmunwarp: $(TOOLDIR)/ppmunwarp.cc

$(TOOLDIR)/ppmwhitebalance: $(TOOLDIR)/ppmwhitebalance.cc

# Calibration requires the creation of an even and an odd deformation
# map.  This indirectly requires the files `even.jpg' and `odd.jpg'.
calibrate: cal-even.deformation cal-odd.deformation $(TOOLDIR)/ppmunwarp

cal-whitepoint: cal-even.white cal-odd.white $(TOOLDIR)/ppmwhitebalance

# To create a deformation map, look for an image of the same name.
# Also creates a [thing]-check.ppm file, in case you're interested in
# what the machine sees.
%.deformation: %.ppm $(TOOLDIR)/ppmunwarp
	$(TOOLDIR)/ppmunwarp -cc $< -d $@ -m $*-check.ppm 2>&1 | tee $*-deflog
	cat $*-deflog | grep 'PPI: ' | cut -d\  -f2 > $*.dpi

%.white: %-grey.ppm %.deformation $(TOOLDIR)/ppmwhitebalance
	$(TOOLDIR)/ppmunwarp -i $< -o $*.dew.ppm -d $*.deformation
	$(TOOLDIR)/ppmwhitebalance -gc 90 $*.dew.ppm > $@
#	$(TOOLDIR)/ppmwhitebalance -c -i $*.dew.ppm -o $@ -m $*-check.ppm
	rm $*.dew.ppm

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

%.dew.ppm: %.in.ppm cal-even.deformation cal-odd.deformation
	$(TOOLDIR)/ppmunwarp -i $< -o $@ \
		-d cal-`echo $* | sed -e 's/.*\(.\)/\1/' -e 's/[02468]/even/' -e 's/[13579]/odd/' `.deformation
	cp cal-`echo $* | sed -e 's/.*\(.\)/\1/' -e 's/[02468]/even/' -e 's/[13579]/odd/' `.dpi $*.dpi

%.wht.ppm: %.dew.ppm cal-even.white cal-odd.white cal-even.dpi cal-odd.dpi
	$(TOOLDIR)/ppmwhitebalance -i $< -o $@ \
		-d cal-`echo $* | sed -e 's/.*\(.\)/\1/' -e 's/[02468]/even/' -e 's/[13579]/odd/' `.white

%.cropinfo: %.wht.ppm
	convert $< -virtual-pixel edge -blur 0x30 -fuzz 20% -trim -format '%wx%h%O' info: > $*.cropinfo

%.out.jpg: %.wht.ppm %.cropinfo
	convert $< -crop `cat $*.cropinfo` +repage \
		-rotate `echo $* | sed -e 's/.*\(.\)/\1/' -e '/^[02468]$$/{s/.*/270/;q}' -e '/^[13579]$$/{s/.*/90/;q}' ` \
		-density `cat $*.dpi` \
		-quality 95 $@

out/%.pdf: %.out.jpg
	convert $< $@

flatten: $(addprefix out/,$(notdir $(addsuffix .pdf,$(basename $(wildcard in/[0-9]*.jpg)))))

out.pdf: flatten calibrate cal-whitepoint
	pdfunite out/*.pdf $@

all: calibrate cal-whitepoint out.pdf
