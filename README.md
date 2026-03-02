# zblas

A native BLAS library for z/OS, written in HLASM, because the mainframe deserves better than calling out to Fortran and pretending that's fine.

## What It Is

zblas is a double-precision BLAS (Basic Linear Algebra Subprograms) library targeting z/Architecture. The hot kernels exploit the z13+ Vector Facility for 128-bit SIMD with two doubles per register and fused multiply-add. If VF isn't available it falls back to scalar BFP loops, so it won't crash on older iron, you'll just be waiting longer.

14 routines across all three BLAS levels, about 3,500 lines of hand-written assembler. All IEEE 754 binary floating point throughout, never hexadecimal. Column-major because that's what BLAS expects and arguing with a 45-year-old API convention is not a productive use of anyone's time.

The verification suite cross-checks every routine against reference implementations I wrote in Python and a Fortran driver built on [slatec](https://github.com/Zaneham/slatec), the classic numerical library I modernised. Both layers have to agree for a test to pass. 130 tests, zero tolerance.

## Why

Most z/OS shops that need linear algebra either call IBM's ESSL (proprietary, expensive, opaque) or link against a Fortran BLAS and hope for the best. There's no open-source z/OS-native option that actually uses the hardware. zblas is that option.

Also I wanted to write a BLAS in assembler and nobody could stop me.

## What's Implemented

All eight Level 1 routines, three Level 2, and three Level 3. This covers the full LU factorisation path since DGETRF needs DGER + DSWAP + IDAMAX + DTRSM + DGEMM, all of which are present.

**Level 1** (vector operations):
DAXPY, DDOT, DSCAL, DNRM2, DCOPY, DSWAP, IDAMAX, DASUM

**Level 2** (matrix-vector operations):
DGEMV (general multiply), DTRSV (triangular solve), DGER (rank-1 update)

**Level 3** (matrix-matrix operations):
DGEMM (general multiply), DTRSM (triangular solve), DSYRK (symmetric rank-k)

## Calling Convention

Every routine supports dual entry points. Standard OS linkage for COBOL, PL/I, and legacy assembler where parameters go by reference via the R1 address list, and XPLINK for XL C/C++ under Language Environment where they go in registers and on the stack. Both use Fortran-style pass-by-reference so C callers write `daxpy_(&n, &alpha, x, &incx, y, &incy)` with the trailing underscore, same as every other BLAS on earth.

See [docs/CALLING_CONVENTION.md](docs/CALLING_CONVENTION.md) for the full details including COBOL and PL/I examples. Yes, you can call DGEMM from COBOL. No, I'm not sorry.

## Building

You need HLASM on z/OS for production, or z390 for off-platform testing:

```
make test      # z390 test harnesses
make verify    # Python + Fortran cross-check (needs python3, optionally gfortran)
```

## Project Layout

```
src/blas1/     Level 1 kernels (8 routines)
src/blas2/     Level 2 kernels (3 routines)
src/blas3/     Level 3 kernels (3 routines)
src/macros/    Shared HLASM macros (prolog/epilog, VF helpers, linkage)
src/fortran/   Slatec routines (linpack, slap, special functions)
include/       C header (zblas.h)
test/          z390 harnesses, verify.py, slatec_driver.f90
docs/          Calling convention docs
```

## License

Apache 2.0. See [LICENSE](LICENSE).
