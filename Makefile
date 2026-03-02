# zblas — z/OS Native BLAS Library
# Makefile for z390 (off-platform) and HLASM (on z/OS)
#
# "Make, like COBOL, will outlive us all." — anon

# ---- Config ----

Z390     ?= z390
ASMA     ?= as
PYTHON   ?= python

SRCDIR   = src
TESTDIR  = test
BLAS1    = $(SRCDIR)/blas1
BLAS2    = $(SRCDIR)/blas2
BLAS3    = $(SRCDIR)/blas3
MACDIR   = $(SRCDIR)/macros

# BLAS Level 1 kernels
B1_SRC   = daxpy ddot dscal dnrm2 dcopy dswap idamax dasum
B1_MLC   = $(addprefix $(BLAS1)/,$(addsuffix .mlc,$(B1_SRC)))

# ---- z390 targets (off-platform testing) ----

.PHONY: all test clean verify

all:
	@echo "zblas: use 'make test' for z390, or assemble on z/OS"

test: test_daxpy test_ddot test_dgemm test_dgemv test_dtrsv test_dtrsm test_dsyrk test_dger
	@echo "All z390 tests passed. The mainframe is pleased."

test_daxpy:
	$(Z390) $(TESTDIR)/test_daxpy.mlc

test_ddot:
	$(Z390) $(TESTDIR)/test_ddot.mlc

test_dgemm:
	$(Z390) $(TESTDIR)/z390_dgemm.mlc

test_dgemv:
	$(Z390) $(TESTDIR)/z390_dgemv.mlc

test_dtrsv:
	$(Z390) $(TESTDIR)/z390_dtrsv.mlc

test_dtrsm:
	$(Z390) $(TESTDIR)/z390_dtrsm.mlc

test_dsyrk:
	$(Z390) $(TESTDIR)/z390_dsyrk.mlc

test_dger:
	$(Z390) $(TESTDIR)/z390_dger.mlc

verify:
	$(PYTHON) $(TESTDIR)/verify.py

clean:
	@echo "Cleaning is beneath us. We are mainframe programmers."
	rm -f $(SRCDIR)/**/*.obj $(TESTDIR)/*.obj
	rm -f $(SRCDIR)/**/*.lst $(TESTDIR)/*.lst
	rm -f $(SRCDIR)/**/*.prn $(TESTDIR)/*.prn

# ---- z/OS targets (on-platform) ----
# These assume HLASM and the z/OS binder are available.
# Uncomment and adjust for your shop's JCL preferences.

# zos: $(B1_MLC:.mlc=.o)
#	$(ASMA) -o libzblas.a $^

# %.o: %.mlc
#	$(ASMA) -I $(MACDIR) -o $@ $<
