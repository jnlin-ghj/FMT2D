FC = gfortran
FFLAGS = -O2 -fbacktrace -cpp
LDFLAGS =

SRCDIR = src
OBJDIR = obj
BINDIR = bin

FWD_OBJS = $(OBJDIR)/mod_fmt2d.o \
           $(OBJDIR)/propagate2d.o \
           $(OBJDIR)/rays2d.o \
           $(OBJDIR)/read_tomo2d.o \
           $(OBJDIR)/fmt2d_main.o

INV_OBJS = $(OBJDIR)/mod_fmt2d.o \
           $(OBJDIR)/propagate2d.o \
           $(OBJDIR)/rays2d.o \
           $(OBJDIR)/read_tomo2d.o \
           $(OBJDIR)/frechet2d.o \
           $(OBJDIR)/inversion2d.o \
           $(OBJDIR)/invert2d_main.o

TARGET_FWD = $(BINDIR)/fmt2d
TARGET_INV = $(BINDIR)/invert2d

all: $(OBJDIR) $(BINDIR) $(TARGET_FWD) $(TARGET_INV)

fmt2d: $(OBJDIR) $(BINDIR) $(TARGET_FWD)

invert2d: $(OBJDIR) $(BINDIR) $(TARGET_INV)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(TARGET_FWD): $(FWD_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(TARGET_INV): $(INV_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/mod_fmt2d.o: $(SRCDIR)/mod_fmt2d.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/propagate2d.o: $(SRCDIR)/propagate2d.f90 $(OBJDIR)/mod_fmt2d.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/rays2d.o: $(SRCDIR)/rays2d.f90 $(OBJDIR)/mod_fmt2d.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/read_tomo2d.o: $(SRCDIR)/read_tomo2d.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/frechet2d.o: $(SRCDIR)/frechet2d.f90 $(OBJDIR)/mod_fmt2d.o $(OBJDIR)/rays2d.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/inversion2d.o: $(SRCDIR)/inversion2d.f90 $(OBJDIR)/mod_fmt2d.o $(OBJDIR)/propagate2d.o $(OBJDIR)/rays2d.o $(OBJDIR)/frechet2d.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/fmt2d_main.o: $(SRCDIR)/fmt2d_main.f90 $(OBJDIR)/mod_fmt2d.o $(OBJDIR)/rays2d.o $(OBJDIR)/read_tomo2d.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/invert2d_main.o: $(SRCDIR)/invert2d_main.f90 $(OBJDIR)/mod_fmt2d.o $(OBJDIR)/rays2d.o $(OBJDIR)/read_tomo2d.o $(OBJDIR)/frechet2d.o $(OBJDIR)/inversion2d.o
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET_FWD) $(TARGET_INV)

.PHONY: all fmt2d invert2d clean