# 2026 10 30, GPT-6 Astra (OpenAI): Initial version (very slight modification by DK)

ifeq ($(origin FC),default)
FC = gfortran
endif
FC ?= gfortran

FFLAGS ?= -O2 -g
PICFLAGS ?= -fPIC
SHLIBFLAGS ?= -shared

# Local Windows default; BinaryBuilder supplies the target extension.
dlext ?= dll
prefix ?= /usr/local
libdir ?= $(prefix)/lib

LIBRARY = libdqawse.$(dlext)

.PHONY: all install clean
.DELETE_ON_ERROR:

all: $(LIBRARY)

$(LIBRARY): dqawse.f Makefile
	$(FC) $(FFLAGS) $(PICFLAGS) $(SHLIBFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

install: $(LIBRARY)
	install -d "$(DESTDIR)$(libdir)"
	install -m 755 "$(LIBRARY)" "$(DESTDIR)$(libdir)/"

clean:
	$(RM) libdqawse.dll libdqawse.so libdqawse.dylib