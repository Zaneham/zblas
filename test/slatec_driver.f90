! slatec_driver.f90 — Test driver that calls slatec-modern BLAS
! and prints results for cross-verification against zblas.
!
! Compile:
!   gfortran -o slatec_driver slatec_driver.f90 \
!     -I<slatec_mod_dir> <service.f90> <blas.f90>
!
! Output: one line per test, "ROUTINE TEST_NAME VALUE"
! verify.py parses this and compares against expected values.

program slatec_driver
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  call test_daxpy()
  call test_ddot()
  call test_dscal()
  call test_dcopy()
  call test_dswap()
  call test_dasum()
  call test_dnrm2()
  call test_idamax()
  call test_dgemm()
  call test_dgemv()
  call test_dtrsv()
  call test_dtrsm()
  call test_dsyrk()
  call test_dger()

contains

  ! ---- Standalone DAXPY (extracted from slatec-modern) ----
  ! We inline it here to avoid module dependency headaches.
  pure subroutine my_daxpy(n, da, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: da, dx(*)
    real(dp), intent(inout) :: dy(*)
    integer :: i, ix, iy
    if (n <= 0 .or. da == 0.0_dp) return
    if (incx == 1 .and. incy == 1) then
      do i = 1, n
        dy(i) = dy(i) + da * dx(i)
      end do
    else
      ix = 1; iy = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      if (incy < 0) iy = (-n + 1) * incy + 1
      do i = 1, n
        dy(iy) = dy(iy) + da * dx(ix)
        ix = ix + incx; iy = iy + incy
      end do
    end if
  end subroutine my_daxpy

  ! ---- Standalone DDOT (from reference BLAS) ----
  pure function my_ddot(n, dx, incx, dy, incy) result(res)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: dx(*), dy(*)
    real(dp) :: res
    integer :: i, ix, iy
    res = 0.0_dp
    if (n <= 0) return
    if (incx == 1 .and. incy == 1) then
      do i = 1, n
        res = res + dx(i) * dy(i)
      end do
    else
      ix = 1; iy = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      if (incy < 0) iy = (-n + 1) * incy + 1
      do i = 1, n
        res = res + dx(ix) * dy(iy)
        ix = ix + incx; iy = iy + incy
      end do
    end if
  end function my_ddot

  ! ---- Standalone DSCAL ----
  pure subroutine my_dscal(n, da, dx, incx)
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: da
    real(dp), intent(inout) :: dx(*)
    integer :: i, ix
    if (n <= 0) return
    if (incx == 1) then
      do i = 1, n
        dx(i) = da * dx(i)
      end do
    else
      ix = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      do i = 1, n
        dx(ix) = da * dx(ix)
        ix = ix + incx
      end do
    end if
  end subroutine my_dscal

  ! ---- Standalone DCOPY ----
  pure subroutine my_dcopy(n, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: dx(*)
    real(dp), intent(inout) :: dy(*)
    integer :: i, ix, iy
    if (n <= 0) return
    if (incx == 1 .and. incy == 1) then
      do i = 1, n
        dy(i) = dx(i)
      end do
    else
      ix = 1; iy = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      if (incy < 0) iy = (-n + 1) * incy + 1
      do i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx; iy = iy + incy
      end do
    end if
  end subroutine my_dcopy

  ! ---- Standalone DSWAP ----
  pure subroutine my_dswap(n, dx, incx, dy, incy)
    integer, intent(in) :: n, incx, incy
    real(dp), intent(inout) :: dx(*), dy(*)
    real(dp) :: tmp
    integer :: i, ix, iy
    if (n <= 0) return
    if (incx == 1 .and. incy == 1) then
      do i = 1, n
        tmp = dx(i); dx(i) = dy(i); dy(i) = tmp
      end do
    else
      ix = 1; iy = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      if (incy < 0) iy = (-n + 1) * incy + 1
      do i = 1, n
        tmp = dx(ix); dx(ix) = dy(iy); dy(iy) = tmp
        ix = ix + incx; iy = iy + incy
      end do
    end if
  end subroutine my_dswap

  ! ---- Standalone DASUM ----
  pure function my_dasum(n, dx, incx) result(res)
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: dx(*)
    real(dp) :: res
    integer :: i, ix
    res = 0.0_dp
    if (n <= 0) return
    if (incx == 1) then
      do i = 1, n
        res = res + abs(dx(i))
      end do
    else
      ix = 1
      if (incx < 0) ix = (-n + 1) * incx + 1
      do i = 1, n
        res = res + abs(dx(ix))
        ix = ix + incx
      end do
    end if
  end function my_dasum

  ! ---- Standalone DNRM2 (naive) ----
  pure function my_dnrm2(n, dx, incx) result(res)
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: dx(*)
    real(dp) :: res, ssq
    integer :: i, ix
    res = 0.0_dp
    if (n <= 0) return
    if (n == 1) then
      ix = 1
      if (incx < 0) ix = 1
      res = abs(dx(ix))
      return
    end if
    ssq = 0.0_dp
    ix = 1
    if (incx < 0) ix = (-n + 1) * incx + 1
    do i = 1, n
      ssq = ssq + dx(ix) * dx(ix)
      ix = ix + incx
    end do
    res = sqrt(ssq)
  end function my_dnrm2

  ! ---- Standalone DGEMM (reference BLAS Level 3) ----
  pure subroutine my_dgemm(transa, transb, m, n, k, alpha, a, lda, &
                           b, ldb, beta, c, ldc)
    character, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real(dp), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
    real(dp), intent(inout) :: c(ldc,*)
    real(dp) :: temp
    integer :: i, j, l
    logical :: nota, notb

    nota = (transa == 'N' .or. transa == 'n')
    notb = (transb == 'N' .or. transb == 'n')

    ! Quick return
    if (m <= 0 .or. n <= 0) return
    ! Alpha = 0: just scale C
    if (alpha == 0.0_dp) then
      do j = 1, n
        if (beta == 0.0_dp) then
          do i = 1, m
            c(i,j) = 0.0_dp
          end do
        else
          do i = 1, m
            c(i,j) = beta * c(i,j)
          end do
        end if
      end do
      return
    end if

    if (nota) then
      ! ---- A not transposed (NN or NT) — DAXPY form ----
      do j = 1, n
        if (beta == 0.0_dp) then
          do i = 1, m
            c(i,j) = 0.0_dp
          end do
        else if (beta /= 1.0_dp) then
          do i = 1, m
            c(i,j) = beta * c(i,j)
          end do
        end if
        do l = 1, k
          if (notb) then
            temp = alpha * b(l,j)
          else
            temp = alpha * b(j,l)
          end if
          if (temp /= 0.0_dp) then
            do i = 1, m
              c(i,j) = c(i,j) + temp * a(i,l)
            end do
          end if
        end do
      end do
    else
      ! ---- A transposed (TN or TT) — dot product form ----
      do j = 1, n
        do i = 1, m
          temp = 0.0_dp
          do l = 1, k
            if (notb) then
              temp = temp + a(l,i) * b(l,j)
            else
              temp = temp + a(l,i) * b(j,l)
            end if
          end do
          if (beta == 0.0_dp) then
            c(i,j) = alpha * temp
          else
            c(i,j) = alpha * temp + beta * c(i,j)
          end if
        end do
      end do
    end if
  end subroutine my_dgemm

  ! ---- Standalone IDAMAX ----
  pure function my_idamax(n, dx, incx) result(res)
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: dx(*)
    integer :: res
    real(dp) :: best
    integer :: i, ix
    res = 0
    if (n <= 0) return
    res = 1
    if (n == 1) return
    ix = 1
    if (incx < 0) ix = (-n + 1) * incx + 1
    best = abs(dx(ix))
    ix = ix + incx
    do i = 2, n
      if (abs(dx(ix)) > best) then
        best = abs(dx(ix))
        res = i
      end if
      ix = ix + incx
    end do
  end function my_idamax

  subroutine test_daxpy()
    real(dp) :: x(10), y(10), alpha
    integer :: n

    ! Test 1: basic, alpha=2, N=4
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    alpha = 2.0_dp; n = 4
    call my_daxpy(n, alpha, x, 1, y, 1)
    write(*,'(A,4(F20.15))') 'DAXPY basic ', y(1), y(2), y(3), y(4)

    ! Test 2: alpha=0 (should not modify Y)
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_daxpy(4, 0.0_dp, x, 1, y, 1)
    write(*,'(A,4(F20.15))') 'DAXPY zero_a ', y(1), y(2), y(3), y(4)

    ! Test 3: N=0
    y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_daxpy(0, 2.0_dp, x, 1, y, 1)
    write(*,'(A,4(F20.15))') 'DAXPY n_zero ', y(1), y(2), y(3), y(4)

    ! Test 4: N=1
    y(1) = 7.0_dp; x(1) = 5.0_dp
    call my_daxpy(1, 3.0_dp, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DAXPY n_one  ', y(1)

    ! Test 5: INCX=2
    x = [1.0_dp, 99.0_dp, 2.0_dp, 99.0_dp, 3.0_dp, &
         99.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_daxpy(3, 1.0_dp, x, 2, y, 1)
    write(*,'(A,3(F20.15))') 'DAXPY incx2  ', y(1), y(2), y(3)

    ! Test 6: negative INCX=-1
    x = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_daxpy(3, 1.0_dp, x, -1, y, 1)
    write(*,'(A,3(F20.15))') 'DAXPY neg_inc', y(1), y(2), y(3)

  end subroutine test_daxpy

  subroutine test_ddot()
    real(dp) :: x(10), y(10), result

    ! Test 1: basic [1,2,3,4].[5,6,7,8] = 70
    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_ddot(4, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DDOT basic   ', result

    ! Test 2: N=0
    result = my_ddot(0, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DDOT n_zero  ', result

    ! Test 3: orthogonal
    x = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_ddot(2, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DDOT ortho   ', result

    ! Test 4: N=1
    x(1) = 3.0_dp; y(1) = 7.0_dp
    result = my_ddot(1, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DDOT n_one   ', result

    ! Test 5: INCX=2
    x = [1.0_dp, 99.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_ddot(2, x, 2, y, 1)
    write(*,'(A,F20.15)') 'DDOT incx2   ', result

    ! Test 6: negative INCX=-1
    x = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [10.0_dp, 20.0_dp, 30.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_ddot(3, x, -1, y, 1)
    write(*,'(A,F20.15)') 'DDOT neg_inc ', result

  end subroutine test_ddot

  subroutine test_dscal()
    real(dp) :: x(10)

    ! basic: alpha=2, X=[1,2,3]
    x = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dscal(3, 2.0_dp, x, 1)
    write(*,'(A,3(F20.15))') 'DSCAL basic  ', x(1), x(2), x(3)

    ! alpha=0
    x = [5.0_dp, 10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dscal(2, 0.0_dp, x, 1)
    write(*,'(A,2(F20.15))') 'DSCAL zero_a ', x(1), x(2)

    ! N=0
    x(1) = 1.0_dp
    call my_dscal(0, 2.0_dp, x, 1)
    write(*,'(A,F20.15)') 'DSCAL n_zero ', x(1)

    ! N=1
    x(1) = 5.0_dp
    call my_dscal(1, 3.0_dp, x, 1)
    write(*,'(A,F20.15)') 'DSCAL n_one  ', x(1)
  end subroutine test_dscal

  subroutine test_dcopy()
    real(dp) :: x(10), y(10)

    ! basic
    x = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = 0.0_dp
    call my_dcopy(3, x, 1, y, 1)
    write(*,'(A,3(F20.15))') 'DCOPY basic  ', y(1), y(2), y(3)

    ! N=0
    y(1) = 99.0_dp
    call my_dcopy(0, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DCOPY n_zero ', y(1)

    ! INCX=2
    x = [1.0_dp, 99.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = 0.0_dp
    call my_dcopy(2, x, 2, y, 1)
    write(*,'(A,2(F20.15))') 'DCOPY incx2  ', y(1), y(2)
  end subroutine test_dcopy

  subroutine test_dswap()
    real(dp) :: x(10), y(10)

    ! basic
    x = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    y = [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dswap(2, x, 1, y, 1)
    write(*,'(A,2(F20.15))') 'DSWAP basic_x', x(1), x(2)
    write(*,'(A,2(F20.15))') 'DSWAP basic_y', y(1), y(2)

    ! N=0
    x(1) = 1.0_dp; y(1) = 2.0_dp
    call my_dswap(0, x, 1, y, 1)
    write(*,'(A,F20.15)') 'DSWAP n_zero_x', x(1)
    write(*,'(A,F20.15)') 'DSWAP n_zero_y', y(1)
  end subroutine test_dswap

  subroutine test_dasum()
    real(dp) :: x(10), result

    ! basic: |1|+|-2|+|3|+|-4| = 10
    x = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_dasum(4, x, 1)
    write(*,'(A,F20.15)') 'DASUM basic  ', result

    ! N=0
    result = my_dasum(0, x, 1)
    write(*,'(A,F20.15)') 'DASUM n_zero ', result

    ! all positive
    x = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_dasum(3, x, 1)
    write(*,'(A,F20.15)') 'DASUM allpos ', result

    ! N=1
    x(1) = -5.0_dp
    result = my_dasum(1, x, 1)
    write(*,'(A,F20.15)') 'DASUM n_one  ', result
  end subroutine test_dasum

  subroutine test_dnrm2()
    real(dp) :: x(10), result

    ! basic: sqrt(9+16) = 5
    x = [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_dnrm2(2, x, 1)
    write(*,'(A,F20.15)') 'DNRM2 basic  ', result

    ! N=0
    result = my_dnrm2(0, x, 1)
    write(*,'(A,F20.15)') 'DNRM2 n_zero ', result

    ! N=1
    x(1) = 3.0_dp
    result = my_dnrm2(1, x, 1)
    write(*,'(A,F20.15)') 'DNRM2 n_one  ', result
  end subroutine test_dnrm2

  subroutine test_idamax()
    real(dp) :: x(10)
    integer :: result

    ! basic: [1,-5,3] -> index 2
    x = [1.0_dp, -5.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_idamax(3, x, 1)
    write(*,'(A,F20.15)') 'IDAMAX basic ', real(result, dp)

    ! N=0
    result = my_idamax(0, x, 1)
    write(*,'(A,F20.15)') 'IDAMAX n_zero', real(result, dp)

    ! N=1
    x(1) = 3.0_dp
    result = my_idamax(1, x, 1)
    write(*,'(A,F20.15)') 'IDAMAX n_one ', real(result, dp)

    ! first is max
    x = [10.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_idamax(4, x, 1)
    write(*,'(A,F20.15)') 'IDAMAX first ', real(result, dp)

    ! last is max
    x = [1.0_dp, 2.0_dp, 3.0_dp, 10.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    result = my_idamax(4, x, 1)
    write(*,'(A,F20.15)') 'IDAMAX last  ', real(result, dp)
  end subroutine test_idamax

  subroutine test_dgemm()
    real(dp) :: a(4), b(4), c(4)

    ! T1: NN, alpha=1, beta=0
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    b = [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
    c = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemm('N','N', 2, 2, 2, 1.0_dp, a, 2, b, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DGEMM nn_basic', c(1), c(2), c(3), c(4)

    ! T2: NN, alpha=2, beta=1, C=[1,1,1,1]
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    b = [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
    c = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    call my_dgemm('N','N', 2, 2, 2, 2.0_dp, a, 2, b, 2, 1.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DGEMM nn_ab   ', c(1), c(2), c(3), c(4)

    ! T3: alpha=0, beta=2
    c = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call my_dgemm('N','N', 2, 2, 2, 0.0_dp, a, 2, b, 2, 2.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DGEMM alpha0  ', c(1), c(2), c(3), c(4)

    ! T5: TN
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    b = [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
    c = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemm('T','N', 2, 2, 2, 1.0_dp, a, 2, b, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DGEMM tn_basic', c(1), c(2), c(3), c(4)

    ! T6: NT
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    b = [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
    c = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemm('N','T', 2, 2, 2, 1.0_dp, a, 2, b, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DGEMM nt_basic', c(1), c(2), c(3), c(4)

    ! T7: TT
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    b = [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
    c = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemm('T','T', 2, 2, 2, 1.0_dp, a, 2, b, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DGEMM tt_basic', c(1), c(2), c(3), c(4)

  end subroutine test_dgemm

  ! ---- Standalone DGEMV ----
  pure subroutine my_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    character, intent(in) :: trans
    integer, intent(in) :: m, n, lda, incx, incy
    real(dp), intent(in) :: alpha, beta, a(lda,*), x(*)
    real(dp), intent(inout) :: y(*)
    real(dp) :: temp
    integer :: i, j, ix, iy, jx, jy, kx, ky, lenx, leny
    logical :: nota

    nota = (trans == 'N' .or. trans == 'n')
    if (nota) then; leny = m; lenx = n
    else;           leny = n; lenx = m; end if

    if (m <= 0 .or. n <= 0) return

    kx = 1; ky = 1
    if (incx < 0) kx = (-lenx + 1) * incx + 1
    if (incy < 0) ky = (-leny + 1) * incy + 1

    ! Scale y by beta
    iy = ky
    if (beta == 0.0_dp) then
      do i = 1, leny; y(iy) = 0.0_dp; iy = iy + incy; end do
    else if (beta /= 1.0_dp) then
      do i = 1, leny; y(iy) = beta * y(iy); iy = iy + incy; end do
    end if
    if (alpha == 0.0_dp) return

    if (nota) then
      jx = kx
      do j = 1, n
        temp = alpha * x(jx)
        if (temp /= 0.0_dp) then
          iy = ky
          do i = 1, m
            y(iy) = y(iy) + temp * a(i,j)
            iy = iy + incy
          end do
        end if
        jx = jx + incx
      end do
    else
      jy = ky
      do j = 1, n
        temp = 0.0_dp
        ix = kx
        do i = 1, m
          temp = temp + a(i,j) * x(ix)
          ix = ix + incx
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy
      end do
    end if
  end subroutine my_dgemv

  ! ---- Standalone DTRSV ----
  pure subroutine my_dtrsv(uplo, trans, diag, n, a, lda, x, incx)
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: n, lda, incx
    real(dp), intent(in) :: a(lda,*)
    real(dp), intent(inout) :: x(*)
    real(dp) :: temp
    integer :: i, j, ix, jx, kx
    logical :: nounit

    if (n <= 0) return
    nounit = (diag == 'N' .or. diag == 'n')
    kx = 1
    if (incx < 0) kx = (-n + 1) * incx + 1

    if (trans == 'N' .or. trans == 'n') then
      if (uplo == 'U' .or. uplo == 'u') then
        jx = kx + (n - 1) * incx
        do j = n, 1, -1
          if (nounit) x(jx) = x(jx) / a(j,j)
          temp = x(jx)
          ix = jx
          do i = j - 1, 1, -1
            ix = ix - incx
            x(ix) = x(ix) - temp * a(i,j)
          end do
          jx = jx - incx
        end do
      else
        jx = kx
        do j = 1, n
          if (nounit) x(jx) = x(jx) / a(j,j)
          temp = x(jx)
          ix = jx
          do i = j + 1, n
            ix = ix + incx
            x(ix) = x(ix) - temp * a(i,j)
          end do
          jx = jx + incx
        end do
      end if
    else
      if (uplo == 'U' .or. uplo == 'u') then
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = kx
          do i = 1, j - 1
            temp = temp - a(i,j) * x(ix)
            ix = ix + incx
          end do
          if (nounit) temp = temp / a(j,j)
          x(jx) = temp
          jx = jx + incx
        end do
      else
        jx = kx + (n - 1) * incx
        do j = n, 1, -1
          temp = x(jx)
          ix = kx + (n - 1) * incx
          do i = n, j + 1, -1
            temp = temp - a(i,j) * x(ix)
            ix = ix - incx
          end do
          if (nounit) temp = temp / a(j,j)
          x(jx) = temp
          jx = jx - incx
        end do
      end if
    end if
  end subroutine my_dtrsv

  ! ---- Standalone DTRSM ----
  pure subroutine my_dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    real(dp), intent(in) :: alpha, a(lda,*)
    real(dp), intent(inout) :: b(ldb,*)
    real(dp) :: temp
    integer :: i, j, k
    logical :: lside, upper, nota, nounit

    if (m <= 0 .or. n <= 0) return
    lside = (side == 'L' .or. side == 'l')
    upper = (uplo == 'U' .or. uplo == 'u')
    nota = (transa == 'N' .or. transa == 'n')
    nounit = (diag == 'N' .or. diag == 'n')

    if (alpha == 0.0_dp) then
      do j = 1, n
        do i = 1, m
          b(i,j) = 0.0_dp
        end do
      end do
      return
    end if

    if (lside) then
      if (nota) then
        if (upper) then
          do j = 1, n
            if (alpha /= 1.0_dp) then
              do i = 1, m; b(i,j) = alpha * b(i,j); end do
            end if
            do k = m, 1, -1
              if (b(k,j) /= 0.0_dp) then
                if (nounit) b(k,j) = b(k,j) / a(k,k)
                do i = 1, k - 1
                  b(i,j) = b(i,j) - b(k,j) * a(i,k)
                end do
              end if
            end do
          end do
        else
          do j = 1, n
            if (alpha /= 1.0_dp) then
              do i = 1, m; b(i,j) = alpha * b(i,j); end do
            end if
            do k = 1, m
              if (b(k,j) /= 0.0_dp) then
                if (nounit) b(k,j) = b(k,j) / a(k,k)
                do i = k + 1, m
                  b(i,j) = b(i,j) - b(k,j) * a(i,k)
                end do
              end if
            end do
          end do
        end if
      else
        if (upper) then
          do j = 1, n
            do i = 1, m
              temp = alpha * b(i,j)
              do k = 1, i - 1
                temp = temp - a(k,i) * b(k,j)
              end do
              if (nounit) temp = temp / a(i,i)
              b(i,j) = temp
            end do
          end do
        else
          do j = 1, n
            do i = m, 1, -1
              temp = alpha * b(i,j)
              do k = i + 1, m
                temp = temp - a(k,i) * b(k,j)
              end do
              if (nounit) temp = temp / a(i,i)
              b(i,j) = temp
            end do
          end do
        end if
      end if
    else
      if (nota) then
        if (upper) then
          do j = 1, n
            if (alpha /= 1.0_dp) then
              do i = 1, m; b(i,j) = alpha * b(i,j); end do
            end if
            do k = 1, j - 1
              if (a(k,j) /= 0.0_dp) then
                do i = 1, m
                  b(i,j) = b(i,j) - a(k,j) * b(i,k)
                end do
              end if
            end do
            if (nounit) then
              temp = 1.0_dp / a(j,j)
              do i = 1, m; b(i,j) = temp * b(i,j); end do
            end if
          end do
        else
          do j = n, 1, -1
            if (alpha /= 1.0_dp) then
              do i = 1, m; b(i,j) = alpha * b(i,j); end do
            end if
            do k = j + 1, n
              if (a(k,j) /= 0.0_dp) then
                do i = 1, m
                  b(i,j) = b(i,j) - a(k,j) * b(i,k)
                end do
              end if
            end do
            if (nounit) then
              temp = 1.0_dp / a(j,j)
              do i = 1, m; b(i,j) = temp * b(i,j); end do
            end if
          end do
        end if
      else
        if (upper) then
          do k = n, 1, -1
            if (nounit) then
              temp = 1.0_dp / a(k,k)
              do i = 1, m; b(i,k) = temp * b(i,k); end do
            end if
            do j = 1, k - 1
              if (a(j,k) /= 0.0_dp) then
                temp = a(j,k)
                do i = 1, m
                  b(i,j) = b(i,j) - temp * b(i,k)
                end do
              end if
            end do
            if (alpha /= 1.0_dp) then
              do i = 1, m; b(i,k) = alpha * b(i,k); end do
            end if
          end do
        else
          do k = 1, n
            if (nounit) then
              temp = 1.0_dp / a(k,k)
              do i = 1, m; b(i,k) = temp * b(i,k); end do
            end if
            do j = k + 1, n
              if (a(j,k) /= 0.0_dp) then
                temp = a(j,k)
                do i = 1, m
                  b(i,j) = b(i,j) - temp * b(i,k)
                end do
              end if
            end do
            if (alpha /= 1.0_dp) then
              do i = 1, m; b(i,k) = alpha * b(i,k); end do
            end if
          end do
        end if
      end if
    end if
  end subroutine my_dtrsm

  ! ---- Standalone DSYRK ----
  pure subroutine my_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    character, intent(in) :: uplo, trans
    integer, intent(in) :: n, k, lda, ldc
    real(dp), intent(in) :: alpha, beta, a(lda,*)
    real(dp), intent(inout) :: c(ldc,*)
    real(dp) :: temp
    integer :: i, j, l
    logical :: upper, nota

    if (n <= 0) return
    upper = (uplo == 'U' .or. uplo == 'u')
    nota = (trans == 'N' .or. trans == 'n')
    if (alpha == 0.0_dp .and. beta == 1.0_dp) return

    if (nota) then
      if (upper) then
        do j = 1, n
          if (beta == 0.0_dp) then
            do i = 1, j; c(i,j) = 0.0_dp; end do
          else if (beta /= 1.0_dp) then
            do i = 1, j; c(i,j) = beta * c(i,j); end do
          end if
          do l = 1, k
            if (a(j,l) /= 0.0_dp) then
              temp = alpha * a(j,l)
              do i = 1, j; c(i,j) = c(i,j) + temp * a(i,l); end do
            end if
          end do
        end do
      else
        do j = 1, n
          if (beta == 0.0_dp) then
            do i = j, n; c(i,j) = 0.0_dp; end do
          else if (beta /= 1.0_dp) then
            do i = j, n; c(i,j) = beta * c(i,j); end do
          end if
          do l = 1, k
            if (a(j,l) /= 0.0_dp) then
              temp = alpha * a(j,l)
              do i = j, n; c(i,j) = c(i,j) + temp * a(i,l); end do
            end if
          end do
        end do
      end if
    else
      if (upper) then
        do j = 1, n
          do i = 1, j
            temp = 0.0_dp
            do l = 1, k; temp = temp + a(l,i) * a(l,j); end do
            if (beta == 0.0_dp) then
              c(i,j) = alpha * temp
            else
              c(i,j) = alpha * temp + beta * c(i,j)
            end if
          end do
        end do
      else
        do j = 1, n
          do i = j, n
            temp = 0.0_dp
            do l = 1, k; temp = temp + a(l,i) * a(l,j); end do
            if (beta == 0.0_dp) then
              c(i,j) = alpha * temp
            else
              c(i,j) = alpha * temp + beta * c(i,j)
            end if
          end do
        end do
      end if
    end if
  end subroutine my_dsyrk

  ! ---- Standalone DGER ----
  pure subroutine my_dger(m, n, alpha, x, incx, y, incy, a, lda)
    integer, intent(in) :: m, n, lda, incx, incy
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda,*)
    real(dp) :: temp
    integer :: i, j, ix, jy, kx

    if (m <= 0 .or. n <= 0 .or. alpha == 0.0_dp) return
    jy = 1
    if (incy < 0) jy = (-n + 1) * incy + 1
    kx = 1
    if (incx < 0) kx = (-m + 1) * incx + 1

    do j = 1, n
      if (y(jy) /= 0.0_dp) then
        temp = alpha * y(jy)
        ix = kx
        do i = 1, m
          a(i,j) = a(i,j) + temp * x(ix)
          ix = ix + incx
        end do
      end if
      jy = jy + incy
    end do
  end subroutine my_dger

  subroutine test_dgemv()
    real(dp) :: a(4), x(4), y(4)

    ! T1: N, a=1,b=0
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    x = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemv('N', 2, 2, 1.0_dp, a, 2, x, 1, 0.0_dp, y, 1)
    write(*,'(A,2(F20.15))') 'DGEMV n_basic ', y(1), y(2)

    ! T2: T, a=1,b=0
    x = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemv('T', 2, 2, 1.0_dp, a, 2, x, 1, 0.0_dp, y, 1)
    write(*,'(A,2(F20.15))') 'DGEMV t_basic ', y(1), y(2)

    ! T3: a=0,b=2
    y = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemv('N', 2, 2, 0.0_dp, a, 2, x, 1, 2.0_dp, y, 1)
    write(*,'(A,2(F20.15))') 'DGEMV alpha0  ', y(1), y(2)

    ! T5: N, a=2,b=1
    x = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    y = [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    call my_dgemv('N', 2, 2, 2.0_dp, a, 2, x, 1, 1.0_dp, y, 1)
    write(*,'(A,2(F20.15))') 'DGEMV n_ab    ', y(1), y(2)

  end subroutine test_dgemv

  subroutine test_dtrsv()
    real(dp) :: a(4), x(4)

    ! T1: U,N,nonunit
    a = [1.0_dp, 0.0_dp, 3.0_dp, 2.0_dp]
    x = [7.0_dp, 4.0_dp, 0.0_dp, 0.0_dp]
    call my_dtrsv('U', 'N', 'N', 2, a, 2, x, 1)
    write(*,'(A,2(F20.15))') 'DTRSV un_nonu ', x(1), x(2)

    ! T2: L,N,nonunit
    a = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp]
    x = [3.0_dp, 12.0_dp, 0.0_dp, 0.0_dp]
    call my_dtrsv('L', 'N', 'N', 2, a, 2, x, 1)
    write(*,'(A,2(F20.15))') 'DTRSV ln_nonu ', x(1), x(2)

    ! T3: U,N,unit
    a = [1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp]
    x = [7.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    call my_dtrsv('U', 'N', 'U', 2, a, 2, x, 1)
    write(*,'(A,2(F20.15))') 'DTRSV un_unit ', x(1), x(2)

  end subroutine test_dtrsv

  subroutine test_dtrsm()
    real(dp) :: a(4), b(4)

    ! T1: L,U,N,nonunit,a=1
    a = [1.0_dp, 0.0_dp, 3.0_dp, 2.0_dp]
    b = [7.0_dp, 4.0_dp, 13.0_dp, 6.0_dp]
    call my_dtrsm('L', 'U', 'N', 'N', 2, 2, 1.0_dp, a, 2, b, 2)
    write(*,'(A,4(F20.15))') 'DTRSM lun_basic', b(1), b(2), b(3), b(4)

    ! T2: L,L,N,nonunit,a=1
    a = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp]
    b = [3.0_dp, 12.0_dp, 5.0_dp, 16.0_dp]
    call my_dtrsm('L', 'L', 'N', 'N', 2, 2, 1.0_dp, a, 2, b, 2)
    write(*,'(A,4(F20.15))') 'DTRSM lln_basic', b(1), b(2), b(3), b(4)

    ! T3: alpha=0
    b = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call my_dtrsm('L', 'U', 'N', 'N', 2, 2, 0.0_dp, a, 2, b, 2)
    write(*,'(A,4(F20.15))') 'DTRSM alpha0  ', b(1), b(2), b(3), b(4)

  end subroutine test_dtrsm

  subroutine test_dsyrk()
    real(dp) :: a(4), c(4)

    ! T1: U,N,a=1,b=0
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    c = [0.0_dp, 99.0_dp, 0.0_dp, 0.0_dp]
    call my_dsyrk('U', 'N', 2, 2, 1.0_dp, a, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DSYRK un_basic', c(1), c(2), c(3), c(4)

    ! T2: L,N,a=1,b=0
    c = [0.0_dp, 0.0_dp, 99.0_dp, 0.0_dp]
    call my_dsyrk('L', 'N', 2, 2, 1.0_dp, a, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DSYRK ln_basic', c(1), c(2), c(3), c(4)

    ! T3: a=0,b=2 (scale upper triangle)
    c = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call my_dsyrk('U', 'N', 2, 2, 0.0_dp, a, 2, 2.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DSYRK alpha0  ', c(1), c(2), c(3), c(4)

    ! T5: U,T
    c = [0.0_dp, 99.0_dp, 0.0_dp, 0.0_dp]
    call my_dsyrk('U', 'T', 2, 2, 1.0_dp, a, 2, 0.0_dp, c, 2)
    write(*,'(A,4(F20.15))') 'DSYRK ut_basic', c(1), c(2), c(3), c(4)

  end subroutine test_dsyrk

  subroutine test_dger()
    real(dp) :: x(4), y(4), a(4)

    ! T1: basic, alpha=1
    x = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    y = [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp]
    a = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call my_dger(2, 2, 1.0_dp, x, 1, y, 1, a, 2)
    write(*,'(A,4(F20.15))') 'DGER basic   ', a(1), a(2), a(3), a(4)

    ! T2: accumulate, alpha=2
    x = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
    y = [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call my_dger(2, 2, 2.0_dp, x, 1, y, 1, a, 2)
    write(*,'(A,4(F20.15))') 'DGER accum   ', a(1), a(2), a(3), a(4)

    ! T3: alpha=0
    a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call my_dger(2, 2, 0.0_dp, x, 1, y, 1, a, 2)
    write(*,'(A,4(F20.15))') 'DGER alpha0  ', a(1), a(2), a(3), a(4)

  end subroutine test_dger

end program slatec_driver
