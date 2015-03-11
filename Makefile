MAKE = make
SRCDIR = src
EXEDIR = bin
OBJDIR = src

.PHONY: build clean dist-clean

all: build

build:
	$(MAKE) -C $(SRCDIR)

solver:
	$(MAKE) -C $(SRCDIR) solve_limid

ve:
	$(MAKE) -C $(SRCDIR) run_ve

clean:
	$(MAKE) -C $(SRCDIR) clean

dist-clean:
	$(MAKE) -C $(SRCDIR) dist-clean
