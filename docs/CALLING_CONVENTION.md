# Calling Convention — zblas

How to call zblas routines from C, COBOL, and PL/I.

Every routine has dual entry points:
- **Standard OS linkage** — for COBOL, PL/I, legacy assembler
- **XPLINK** — for XL C/C++, Language Environment

## OS Linkage (Standard)

Traditional z/OS calling convention. Parameters passed by reference
via an address list pointed to by R1.

```
R1  → +0  addr of parm1
       +4  addr of parm2
       ...
       +n  addr of parmN  (high bit set)
R13 = save area (72 bytes, caller-provided)
R14 = return address
R15 = entry point
```

Return values:
- Integer results in R15
- Float results in FPR0

### From COBOL

```cobol
       IDENTIFICATION DIVISION.
       PROGRAM-ID. TESTAXPY.
       DATA DIVISION.
       WORKING-STORAGE SECTION.
       01 WS-N      PIC S9(9) BINARY VALUE 100.
       01 WS-ALPHA  COMP-2 VALUE 2.5.
       01 WS-INCX   PIC S9(9) BINARY VALUE 1.
       01 WS-INCY   PIC S9(9) BINARY VALUE 1.
       01 WS-X.
          05 WS-X-ELEM COMP-2 OCCURS 100 TIMES.
       01 WS-Y.
          05 WS-Y-ELEM COMP-2 OCCURS 100 TIMES.
       PROCEDURE DIVISION.
           CALL 'DAXPY' USING WS-N WS-ALPHA WS-X WS-INCX
                               WS-Y WS-INCY.
           STOP RUN.
```

COBOL passes everything by reference and sets the high bit on the
last parameter address automatically. This maps directly to OS linkage.

### From PL/I

```pli
DCL DAXPY ENTRY OPTIONS(ASM);
DCL N     FIXED BIN(31) INIT(100);
DCL ALPHA FLOAT BIN(53) INIT(2.5);
DCL X(100) FLOAT BIN(53);
DCL Y(100) FLOAT BIN(53);
DCL INCX  FIXED BIN(31) INIT(1);
DCL INCY  FIXED BIN(31) INIT(1);

CALL DAXPY(N, ALPHA, X, INCX, Y, INCY);
```

PL/I with `OPTIONS(ASM)` uses standard OS linkage.

### From Assembler

```hlasm
         LA    R1,PLIST            Point to parameter list
         L     R15,=V(DAXPY)       Load entry point
         BALR  R14,R15             Call DAXPY
         ...
PLIST    DC    A(N)                Addr of N
         DC    A(ALPHA)            Addr of alpha
         DC    A(DX)               Addr of X array
         DC    A(INCX)             Addr of INCX
         DC    A(DY)               Addr of Y array
         DC    X'80',AL3(INCY)     Addr of INCY (high bit set)
```

## XPLINK

For XL C/C++ compiled with `-Wc,XPLINK`. Parameters in registers
and on the downward-growing stack.

```
R1-R3  = first three pointer parameters
FPR0/2/4/6 = first four float parameters
R4     = DSA (stack pointer)
R7     = return address
```

Return values:
- Integer in R3
- Float in FPR0

### From C (XPLINK)

```c
#include "zblas.h"

int main(void) {
    int n = 100, incx = 1, incy = 1;
    double alpha = 2.5;
    double x[100], y[100];

    /* Standard BLAS call — all by reference */
    daxpy_(&n, &alpha, x, &incx, y, &incy);

    /* Dot product — returns double */
    double result = ddot_(&n, x, &incx, y, &incy);

    return 0;
}
```

Compile and link:
```
xlc -Wc,XPLINK -o testaxpy testaxpy.c -lzblas
```

## Notes

- All parameters are passed **by reference** (pointers), following
  the Fortran BLAS convention. This is why C callers use `&n` etc.

- The trailing underscore on C names (`daxpy_`, `ddot_`) matches
  the Fortran name-mangling convention used by most BLAS libraries.
  This ensures drop-in compatibility with existing BLAS-consuming code.

- Integer parameters are 32-bit signed (`FIXED BIN(31)` / `int`).

- Float parameters are 64-bit IEEE 754 (`COMP-2` / `FLOAT BIN(53)` /
  `double`). zblas uses BFP (Binary Floating Point) instructions
  throughout, not the old hexadecimal floating point.

- Array parameters point to the first element. The routines compute
  offsets from the base using the increment parameters.

- Negative increments are supported (Phase 2). The starting position
  is adjusted per BLAS convention: `offset = (1-N)*INC*8`.
