module newsplinemod

! this is newspline, a fast, means-preserving spline for interval data
! Leo O Lai and Jed O. Kaplan
! 2020

use parametersmod, only : i4,sp
use utilitiesmod,  only : matsol,findloc

implicit none

public  :: newspline_all        ! Wrapper subroutine with all optional arguments
private :: newspline            ! Subroutine of newspline base function with no optional arguments
private :: newspline_bound      ! Subroutine with absolute upper and lower limit options
private :: newspline_pbound     ! Subroutine with percentage limit option
private :: newspline_bound_all  ! Subroutine with both absolute and percentage limit options
private :: monocheck            ! Monotonicity check subroutine for each control point
private :: mono_adjust          ! Adjustment subroutine for un-monotonic control points
private :: days_even            ! calculates Hermite cubic values with even number of partitions
private :: days_odd             ! calculates Hermite cubic values with odd number of partitions
private :: days_osc_even        ! calculates Hermite cubic values with even number of partitions in intervals adjusted for oscillation
private :: days_osc_odd         ! calculates Hermite cubic values with odd number of partitions in intervals adjusted for oscillation

contains

! ------------------------------------------------------------------------------------------------------------------

subroutine newspline_all(monthdata,nk,daydata,llim,ulim,plim,prec)

! COMPOSITE WRAPPER SUBROUTINE : newspline_all
!---
! newspline_all(monthdata, nk, daydata, llim, ulim, plim, prec)
!---
! monthdata (input) : mean series of N length (e.g. monthly temperature)
! nk        (input) : series of smaller parititions of N length (e.g. number of days per month)
! daydata   (output): series of interpolated data of sum(nk) length (e.g. daily temperature)
!---
! llim (optional): absolute lower limit (real number)
! ulim (optional): absolute upper limit (real number)
! plim (optional): percentage limit (real number > 0.)
! prec (optional): decimal places of output (integer(i4) => 0)
!---
! All options can be applied simultaneously.

implicit none

real(sp),    dimension(:), intent(in)  :: monthdata
integer(i4), dimension(:), intent(in)  :: nk
real(sp),    dimension(:), intent(out) :: daydata

real(sp),    optional, intent(in)  :: llim
real(sp),    optional, intent(in)  :: ulim
real(sp),    optional, intent(in)  :: plim
integer(i4), optional, intent(in)  :: prec

integer(i4) :: precision
integer(i4), dimension(:), allocatable :: int_daydata

!---

if (present(plim)) then

  if (present(llim)) then

    if (present(ulim)) then

      call newspline_bound_all(monthdata, nk, daydata, llim, ulim, plim)

    else if (.not. present(ulim)) then

      call newspline_bound_all(monthdata, nk, daydata, llim, 999999., plim)

    end if

  else if (.not. present(llim)) then

    if (present(ulim)) then

      call newspline_bound_all(monthdata, nk, daydata, -999999., ulim, plim)

    else if (.not. present(ulim)) then

      call newspline_pbound(monthdata, nk, daydata, plim)

    end if

  end if

else if (.not. present(plim)) then

  if (present(llim)) then

    if (present(ulim)) then

      call newspline_bound(monthdata, nk, daydata, llim, ulim)

    else if (.not. present(ulim)) then

      call newspline_bound(monthdata, nk, daydata, llim, 999999.)

    end if

  else if (.not. present(llim)) then

    if (present(ulim)) then

      call newspline_bound(monthdata, nk, daydata, -999999., ulim)

    else if (.not. present(ulim)) then

      call newspline(monthdata, nk, daydata)

    end if

  end if

end if

!------

if(present(prec)) then

  allocate(int_daydata(size(daydata)))

  precision = 10 ** prec

  int_daydata = nint(daydata * real(precision))

  daydata = real(int_daydata) / precision

end if

end subroutine newspline_all

! ------------------------------------------------------------------------------------------------------------------

subroutine newspline(monthdata,nk,daydata)

implicit none

real(sp),    dimension(:), intent(in)  :: monthdata
integer(i4), dimension(:), intent(in)  :: nk
real(sp),    dimension(:), intent(out) :: daydata

! Local variables for first mid-control points
real(sp), dimension(:), allocatable :: fmc
real(sp), dimension(:), allocatable :: d
real(sp), dimension(:), allocatable :: m

! Local variables for wall control points
real(sp), dimension(:), allocatable :: swc

! Local variables for linear system of mid control adjustments
real(sp), dimension(:,:), allocatable :: mat
real(sp), dimension(:),   allocatable :: solution
real(sp), dimension(:),   allocatable :: all_solution

! Final vector of all control points
real(sp), dimension(:), allocatable :: all_cont

! Local variables for generating daily values
real(sp), dimension(:), allocatable :: d_cont
real(sp), dimension(:), allocatable :: m_cont

! Hermite cubic and quartic spline basis functions
real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11
real(sp) :: G_00
real(sp) :: G_01
real(sp) :: G_10
real(sp) :: G_11
real(sp) :: u

integer(i4) :: len, len_cont
integer(i4) :: i, j, k, n

!----------------------------------------------------------------
! Start of the newspline routine

len = size(monthdata)

allocate(fmc(len+2))
allocate(d(len+1))
allocate(m(len+2))

!------
! Define first mid-control point as equal to original monthdata
fmc(1) = monthdata(1)
fmc(2:(len+1)) = monthdata(1:len)
fmc(len+2) = monthdata(len)

!---

do i = 1, (len+1)

  d(i) = fmc(i+1) - fmc(i)

end do

!---

do i = 2, (len+1)

  m(i) = (d(i-1) + d(i)) / 2

  !---
  ! Monotonic adjustment to slope
  if(d(i-1) > 0 .and. d(i) < 0) then

    m(i) = 0

  else if(d(i-1) < 0 .and. d(i) > 0) then

    m(i) = 0

  end if

end do

m(1)     = (d(1) + 0.) / 2.
m(len+2) = (d(len+1) + 0.) / 2.

do i = 1, (len+1)

  if(d(i) == 0) then

    m(i) = 0.
    m(i+1) = 0.

  end if

end do

!------
! Calculate wall control based on interception of Hermite functions
allocate(swc(len+1))

u = 0.5

H_00 = 1. + (u**2) * (2.*u - 3.)
H_10 = u * (u - 1.) * (u - 1.)
H_01 = (u**2) * (3. - 2.*u)
H_11 = (u**2) * (u - 1.)

do i = 1, (len+1)

  swc(i) = (fmc(i) * H_00) + (m(i) * H_10) + (fmc(i+1) * H_01) + (m(i+1) * H_11)

end do

!------
! Generate matrix for final adjustments to mid-control points
allocate(mat(len,len))
allocate(solution(len))

!---

u = 1.

G_00 = u - (u**3) - (u**4) / 2.
G_10 = (u**2) * (3.*(u**2) - 8.*u + 6.) / 12.
G_01 = (u**3) * (1. - u / 2.)
G_11 = (u**3) * (3.*u - 4.) / 12.

!---

mat = 0.

! Consider two "buffer midpoints" outside of first and last interval
mat(1,1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.
mat(1,2) = 0.5 * G_11

mat(len,len-1) = -0.5 * G_10
mat(len,len)   = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

n = 1

do i = 2, (len-1)

  mat(i,n)   = -0.5 * G_10

  mat(i,n+1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

  mat(i,n+2) = 0.5 * G_11

  n = n + 1

end do

!---

solution(1)   = 2. * monthdata(1) - swc(1) * G_00 - 0.5 * (swc(2) - swc(1)) * G_11 -    &
               0.5 * (swc(2) - swc(1)) * G_10 - swc(2) * G_01 - swc(1)

solution(1)   = solution(1) + 0.5 * G_10 * monthdata(1)

solution(len) =  2. * monthdata(len) - swc(len) * G_00 - 0.5 * (swc(len+1) - swc(len)) * G_11 -   &
                0.5 * (swc(len+1) - swc(len)) * G_10 - swc(len+1) * G_01 - swc(len)

solution(len) = solution(len) - 0.5 * G_11 * monthdata(len)

do i = 2, (len-1)

  solution(i) = 2. * monthdata(i) - swc(i) * G_00 - 0.5 * (swc(i+1) - swc(i)) * G_11 - 0.5 * (swc(i+1) - swc(i)) *   &
                     G_10 - swc(i+1) * G_01 - swc(i)

end do

! Solve linear system to get final mid control points

call matsol(mat,solution)

!------
! Compile wall control with newly adjusted mid control points (all_cont)

allocate(all_solution(len+2))

all_solution(1)       = monthdata(1)
all_solution(2:len+1) = solution
all_solution(len+2)   = monthdata(len)

!---
allocate(all_cont(size(swc)+size(all_solution)))

all_cont(1) = all_solution(1)

n = 2
do i = 1, size(swc)

  all_cont(n)   = swc(i)
  all_cont(n+1) = all_solution(i+1)

  n = n + 2

end do

!------
! Construct the spline for daily values based on all_cont

len_cont = size(all_cont)

allocate(d_cont(len_cont-1))
allocate(m_cont(len_cont-2))

do i = 1, (len_cont-1)

  d_cont(i) = all_cont(i+1) - all_cont(i)

end do

!---

do i = 1, (len_cont-2)

  m_cont(i) = (d_cont(i) + d_cont(i+1)) / 2

end do

do i = 2, (len_cont-3)

  if(d_cont(i) == 0) then

    m_cont(i)   = 0.
    m_cont(i+1) = 0.

  end if

end do

!------
!Assign the daily value based on sequence

n = 1
k = 1

do i = 1, len !outer loop start, for all monthly intervals N

  if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

    u = 1. / nk(i)

    do j = 1, (nk(i) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k+1) * H_00 + m_cont(k) * H_10 + all_cont(k+2) * H_01 + m_cont(k+1) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

    !---

    u = 1. / nk(i)

    do j = 1, (nk(i) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k+2) * H_00 + m_cont(k+1) * H_10 + all_cont(k+3) * H_01 + m_cont(k+2) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

  else ! if odd months (or odd number of smaller time step)

    u = 1. / nk(i)

    do j = 1, ((nk(i)+1) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k+1) * H_00 + m_cont(k) * H_10 + all_cont(k+2) * H_01 + m_cont(k+1) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

    !---

    u = 2. / nk(i)

    do j = 1, ((nk(i)-1) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k+2) * H_00 + m_cont(k+1) * H_10 + all_cont(k+3) * H_01 + m_cont(k+2) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

  end if

  k = k + 2

end do !end of outer loop

end subroutine newspline

! ------------------------------------------------------------------------------------------------------------------

subroutine newspline_bound(monthdata,nk,daydata,llim,ulim)

implicit none

real(sp),     dimension(:), intent(in)  :: monthdata
integer(i4),  dimension(:), intent(in)  :: nk
real(sp),     dimension(:), intent(out) :: daydata
real(sp),                   intent(in)  :: llim
real(sp),                   intent(in)  :: ulim

! Local variables for first mid control points
real(sp), dimension(:), allocatable :: fmc
real(sp), dimension(:), allocatable :: d
real(sp), dimension(:), allocatable :: m

! Local variables for wall control points
real(sp), dimension(:), allocatable :: swc

! Local variables for linear system of mid control adjustments
real(sp), dimension(:,:), allocatable :: mat
real(sp), dimension(:),   allocatable :: solution
real(sp), dimension(:),   allocatable :: all_solution

integer(i4) :: len
integer(i4) :: i, n

! Final vector of all control points
real(sp), dimension(:), allocatable :: all_cont

! Hermite cubic quartic spline basis functions
real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11

real(sp)    :: u
integer(i4) :: l

! Local variables for generating daily values after minmax bound adjustment of all_cont
real(sp), dimension(:), allocatable :: d_new
real(sp), dimension(:), allocatable :: m_new

integer(i4) :: len_new
integer(i4) :: slpe_l, slpe_r
integer(i4) :: k

! Local variables for max and min bound adjustments
integer(i4), dimension(:), allocatable :: d_orig
logical,     dimension(:), allocatable :: osc_check
real(sp),    dimension(:), allocatable :: c2
real(sp),    dimension(:), allocatable :: root
integer(i4), dimension(:), allocatable :: root_days

! Local variables for calculating root of quadratic approximation
real(sp)    :: diff_yi1
real(sp)    :: diff_yi
real(sp)    :: top, bot
real(sp)    :: root_adj
integer(i4) :: count

real(sp) :: del
real(sp) :: G_00, G_10, G_01, G_11
real(sp) :: yi, yi1, y2i, y2i1
real(sp) :: mi, m2i1
real(sp) :: area_total
real(sp) :: area_int

! Local variables for x_new and y_new
real(sp), dimension(:), allocatable :: x_new
real(sp), dimension(:), allocatable :: y_new
real(sp)    :: kk
integer(i4) :: nn, mm

!------------------
! PART 1: Determine all wall and mid control points
!------------------

! Start of the spline routine
len = size(monthdata)

allocate(fmc(len+2))
allocate(d(len+1))
allocate(m(len+2))

!------
! Define first mid-control point as equal to original monthdata
fmc(1) = monthdata(1)
fmc(2:(len+1)) = monthdata(1:len)
fmc(len+2) = monthdata(len)

!---

do i = 1, (len+1)

  d(i) = fmc(i+1) - fmc(i)

end do

!---

do i = 2, (len+1)

  m(i) = (d(i-1) + d(i)) / 2

  !---
  ! Monotonic adjustment to slope

  if(d(i-1) > 0 .and. d(i) < 0) then

    m(i) = 0

  else if(d(i-1) < 0 .and. d(i) > 0) then

    m(i) = 0

  end if

end do

m(1)     = d(1)
m(len+2) = d(len+1)

do i = 1, (len+1)

  if(d(i) == 0) then

    m(i) = 0.
    m(i+1) = 0.

  end if

end do

!------
! Calculate "second" wall control based on interception of Hermite functions
allocate(swc(len+1))

u = 0.5

H_00 = 1. + (u**2) * (2*u - 3)
H_10 = u * (u - 1) * (u - 1)
H_01 = (u**2) * (3 - 2*u)
H_11 = (u**2) * (u - 1)

do i = 1, (len+1)

  swc(i) = (fmc(i) * H_00) + (m(i) * H_10) + (fmc(i+1) * H_01) + (m(i+1) * H_11)

end do

!------
! Compile the NxN linear matrix to adjust mid-control points
allocate(mat(len,len))
allocate(solution(len))

!---

u = 1.

G_00 = u - (u**3) - (u**4) / 2.
G_10 = (u**2) * (3.*(u**2) - 8.*u + 6.) / 12.
G_01 = (u**3) * (1. - u / 2.)
G_11 = (u**3) * (3*u - 4.) / 12.

!---

mat = 0.

! Consider two "midpoints" outside of first and last interval
mat(1,1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.
mat(1,2) = 0.5 * G_11

mat(len,len-1) = -0.5 * G_10
mat(len,len)   = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

n = 1
do i = 2, (len-1)

  mat(i,n)   = -0.5 * G_10

  mat(i,n+1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

  mat(i,n+2) = 0.5 * G_11

  n = n + 1

end do

!---

solution(1) = 2. * monthdata(1) - swc(1) * G_00 - 0.5 * (swc(2) - swc(1)) * G_11 - 0.5 * (swc(2) - swc(1)) * G_10 - swc(2) *    &
                  G_01 - swc(1) + 0.5 * G_10 * monthdata(1)

solution(len) = 2. * monthdata(len) - swc(len) * G_00 - 0.5 * (swc(len+1) - swc(len)) * G_11 - 0.5 * (swc(len+1) - swc(len)) *   &
                   G_10 - swc(len+1) * G_01 - swc(len) - 0.5 * G_11 * monthdata(len)

do i = 2, (len-1)

  solution(i) = 2. * monthdata(i) - swc(i) * G_00 - 0.5 * (swc(i+1) - swc(i)) * G_11 - 0.5 * (swc(i+1) - swc(i)) *    &
                  G_10 - swc(i+1) * G_01 - swc(i)

end do

! Solve linear system to get final mid control points
call matsol(mat, solution)

!------
! Compile "second" wall control with newly adjusted mid control points (all_cont)
allocate(all_solution(len+2))

all_solution(1)       = monthdata(1)
all_solution(2:len+1) = solution
all_solution(len+2)   = monthdata(len)

!---
allocate(all_cont(size(swc)+size(all_solution)))

n = 1
do i = 1, size(swc)

  all_cont(n)   = swc(i)
  all_cont(n+1) = all_solution(i+1)

  n = n + 2

end do

!---

all_cont(size(all_cont)) = swc(size(swc))

!------------------
! PART 2: Adjust monotonicty using the tridiagonal equation
!------------------

call mono_adjust(monthdata, all_cont)

!------------------
! PART 3: Adjustment to maximum and minimum bound
!------------------

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))
allocate(root(len))
allocate(root_days(len))

d_orig = -9999.
osc_check = .FALSE.
c2 = -9999.
root = -9999.
root_days = -9999

!------
! Assign -1 for negative slope, 1 for postive and 0 for turning point

do i = 2, (len-1)

  if ((monthdata(i+1) - monthdata(i)) < 0) then

    d_orig(i) = -1

  else if ((monthdata(i+1) - monthdata(i)) > 0) then

    d_orig(i) = 1

  end if

  !---

  if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

    d_orig(i) = 0

  else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

    d_orig(i) = 0

  end if

end do

!------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold

do i = 2, (len-1)

  !j = 2 * i

  if (d_orig(i) == 0 .and. monthdata(i) > 0) then

    if (all_cont(2*i) > ulim) then

      osc_check(i) = .TRUE.

    end if

  else if (d_orig(i) == 0 .and. monthdata(i) < 0) then

    if (all_cont(2*i) < llim) then

      osc_check(i) = .TRUE.

    end if

  end if

end do

!------
! Calculate the amount of adjustment required and insert into the c2 variable

do i = 2, (len-1)

  if (osc_check(i) .and. monthdata(i) > 0) then

    c2(i) = ulim - monthdata(i)

  else if (osc_check(i) .and. monthdata(i) < 0) then

    c2(i) = llim - monthdata(i)

  end if

end do

!------
! Calculate the root deviation based on c2 and triangular approximation

count = 0 ! Count the number of extra Hermite spline segments required

do i = 2, (len-1)

  if(osc_check(i)) then

    diff_yi = monthdata(i) - all_cont((2*i)-1)
    diff_yi1 = monthdata(i) - all_cont((2*i)+1)

    !---

    top = 3. * (diff_yi + diff_yi1)
    bot = 8. * c2(i) + 3. * (diff_yi + diff_yi1)

    !---

    root(i) = top / bot

    root(i) = root(i) + 1. / real(nk(i))

    !---

    if (mod(nk(i),2) == 0) then ! EVEN month partitions

      u = 1. / nk(i)

      root_adj = 9999.

        l = 0
        do while (root_adj > 0.)

          root_adj = root(i) - u

          u = u + (2. / nk(i))

          l = l + 1

        end do

      !---

      root_days(i) = l - 1

      !---

    else if (mod(nk(i),2) /= 0) then ! ODD month partitions

      u = 0. !1. / (real(nk(i)))

      root_adj = 9999.

        l = 0
        do while (root_adj > 0.)

          root_adj = root(i) - u

          u = u + (2. / real(nk(i)))

          l = l + 1

        end do

      !---

      root_days(i) = l - 1

      !---

    end if

    !---

    count = count + 2

  end if

end do

!------
! Re-estimate c2 based on integral area under Hermite spline

do i = 1, (len-1)

  if(osc_check(i) ) then

    !--- Construct fourth degree Hermite at u = 1 (integral [0,1])

    del = 1 - root(i)

    u = 1.

    G_00 = u - (u**3) - (u**4) / 2.
    G_01 = (u**3) * (1. - u/2.)
    G_10 = (U**2) * (3. * (u**2) - 8. *u + 6.) / 12.
    G_11 = (u**3) * (3.*u - 4.) / 12.

    !--- Assign local control points

    yi   = all_cont((2*i)-1)
    yi1  = monthdata(i)
    y2i  = monthdata(i)
    y2i1 = all_cont((2*i)+1)

    !--- Assign local slope

    mi   = ((yi - all_cont((2*i)-2)) / 1. + (all_cont(2*i) - yi) / 1.) / 2.
    m2i1 = ((all_cont((2*i)+2) - y2i1) / 1. + (y2i1 - all_cont(2*i)) / 1.) / 2.

    !--- Calculate new area approximation based on Hermite intergral

    area_total = 2. * del * monthdata(i)

    top = (G_00 * yi) + (G_10 * del * mi) + (G_01 * yi1) + (G_00 * y2i) + (G_01 * y2i1) + (G_11 * del * m2i1)
    bot = 1. + ((3. * del**2) / (2. * root(i)**2)) * G_11 - ((3. * del**2) / (2. * root(i)**2)) * G_10

    area_int = (area_total - del * (top + yi + y2i)) / bot

    !--- Re-assign c2 as the integral-estimated value

    c2(i) = (3. * area_int) / (4. * root(i))

  end if

end do

!------
! Generate x_new and y_new series that contains the extra quadratic adjusted segments

allocate(x_new((2*len)+2+count))
allocate(y_new((2*len)+2+count))

x_new(1:3) = [1.,2.,3.]
y_new(1:3) = all_cont(1:3)

nn = 4
mm = 4
kk = 4.

do i = 2, len

  if (.not.osc_check(i)) then

    x_new(nn)   = kk
    x_new(nn+1) = kk + 1.

    y_new(nn)   = all_cont(mm)
    y_new(nn+1) = all_cont(mm+1)

    nn = nn + 2

    mm = mm + 2

    kk = kk + 2.

  else if (osc_check(i) ) then

    x_new(nn)   = kk - root(i)
    x_new(nn+1) = kk
    x_new(nn+2) = kk + root(i)
    x_new(nn+3) = kk + 1.

    y_new(nn)   = monthdata(i)
    y_new(nn+1) = monthdata(i) + c2(i)
    y_new(nn+2) = monthdata(i)
    y_new(nn+3) = all_cont(mm+1)

    nn = nn + 4

    mm = mm + 2

    kk = kk + 2.

  end if

end do

x_new(nn) = kk
y_new(nn) = all_cont(mm)

!------
! Construct the spline for daily values based on all_cont
len_new = size(x_new)

allocate(d_new(len_new-1))
allocate(m_new(len_new))

do i = 1, (len_new-1)

  d_new(i) = (y_new(i+1) - y_new(i)) / (x_new(i+1) - x_new(i))

end do

!---

do i = 2, (len_new-1)

  m_new(i) = (d_new(i-1) + d_new(i)) / 2.

  !---
  ! Monotonic adjustment to slope
  if(d_new(i-1) > 0 .and. d_new(i) < 0) then

    m_new(i) = 0

  else if(d_new(i-1) < 0 .and. d_new(i) > 0) then

    m_new(i) = 0

  end if

end do

m_new(1)       = d_new(1) / 2.
m_new(len_new) = d_new(len_new-1) / 2.

!---

do i = 1, (len_new-1)

  if(d_new(i) == 0) then

    m_new(i)   = 0.
    m_new(i+1) = 0.

  end if

end do

!------
! Reassign quadratic approximation slopes and original slopes to adjusted intervals
do i = 2, (len-1)

  if (osc_check(i) ) then

    call findloc(x_new, (real(2*i)-root(i)), slpe_l)

    slpe_r = slpe_l + 2

    !---

    m_new(slpe_l) = 2. * c2(i) / root(i)

    m_new(slpe_r) = -2. * c2(i) / root(i)

    !---

    m_new(slpe_l-1) = ((all_cont(2*i-1) - all_cont((2*i)-2)) / 1. + (all_cont(2*i) - all_cont(2*i-1)) / 1.) / 2.

    m_new(slpe_r+1) = ((all_cont((2*i)+2) - all_cont(2*i+1)) / 1. + (all_cont(2*i+1) - all_cont(2*i)) / 1.) / 2.

  end if

end do

!---
! Ensure smooth monotonic interpolation between control points

call monocheck(d_new, m_new)

!------
!Assign the daily value based on sequence

n = 1
k = 1

do i = 1, len !outer loop start, for all monthly intervals N

  if (.not.osc_check(i)) then ! For intervals that did not required adjustment

    if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

      l = n + nk(i) - 1

      call days_even(nk(i), y_new(k:k+2), m_new(k:k+2), daydata(n:l))

      n = n + nk(i)

    else ! if odd months (or odd number of smaller time step)

      l = n + nk(i) - 1

      call days_odd(nk(i), y_new(k:k+2), m_new(k:k+2), daydata(n:l))

      n = n + nk(i)

    end if

    k = k + 2

  else if (osc_check(i) ) then ! for intervals adjusted by minmax bound

    if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

      l = n + nk(i) - 1

      call days_osc_even(nk(i), root_days(i), root(i), y_new(k:k+4), m_new(k:k+4), daydata(n:l))

      n = n + nk(i)

    else ! if odd months (or odd number of smaller time step)

      l = n + nk(i) - 1

      call days_osc_odd(nk(i), root_days(i), root(i), y_new(k:k+4), m_new(k:k+4), daydata(n:l))

      n = n + nk(i)

    end if

    k = k + 4

  end if

end do !end of outer loop

end subroutine newspline_bound

! ------------------------------------------------------------------------------------------------------------------

subroutine newspline_pbound(monthdata,nk,daydata,plim)

implicit none

real(sp),    dimension(:), intent(in)  :: monthdata
integer(i4), dimension(:), intent(in)  :: nk
real(sp),    dimension(:), intent(out) :: daydata
real(sp),                  intent(in)  :: plim      !taken in as 0-100%

! Local variables for first mid control points
real(sp), dimension(:), allocatable :: fmc
real(sp), dimension(:), allocatable :: d
real(sp), dimension(:), allocatable :: m

! Local variables for second wall control points
real(sp), dimension(:), allocatable :: swc

! Local variables for linear system of mid control adjustments
real(sp), dimension(:,:), allocatable :: mat
real(sp), dimension(:),   allocatable :: solution
real(sp), dimension(:),   allocatable :: all_solution

integer(i4) :: len
integer(i4) :: i, n

! Hermite cubic spline basis functions
real(sp)     :: H_00
real(sp)     :: H_01
real(sp)     :: H_10
real(sp)     :: H_11

real(sp)     :: u
integer(i4) :: l

! Final vector of all control points
real(sp), dimension(:), allocatable :: all_cont

! Local variables for generating daily values after minmax bound adjustment of all_cont
real(sp), dimension(:), allocatable :: d_new
real(sp), dimension(:), allocatable :: m_new

integer(i4) :: len_new
integer(i4) :: slpe_l, slpe_r
integer(i4) :: k

! Local variables for max and min bound adjustments
integer(i4), dimension(:), allocatable :: d_orig
logical,     dimension(:), allocatable :: osc_check
real(sp),    dimension(:), allocatable :: c2
real(sp),    dimension(:), allocatable :: root
integer(i4), dimension(:), allocatable :: root_days
real(sp)                               :: perc

! Local variables for calculating root of quadratic approximation
real(sp)     :: diff_yi1
real(sp)     :: diff_yi
real(sp)     :: top, bot
real(sp)     :: root_adj
integer(i4)  :: count

real(sp) :: del

real(sp) :: G_00
real(sp) :: G_01
real(sp) :: G_10
real(sp) :: G_11

real(sp) :: yi, yi1
real(sp) :: y2i, y2i1

real(sp) :: mi, m2i1
real(sp) :: area_total
real(sp) :: area_int

! Local variables for x_new and y_new
real(sp), dimension(:), allocatable :: x_new
real(sp), dimension(:), allocatable :: y_new

real(sp)     :: kk
integer(i4)  :: nn, mm

!------------------
! PART 1: Determine all wall and mid control points
!------------------

! Start of the spline routine
len = size(monthdata)

allocate(fmc(len+2))
allocate(d(len+1))
allocate(m(len+2))

!------
! Define first mid-control point as equal to original monthdata
fmc(1) = monthdata(1)
fmc(2:(len+1)) = monthdata(1:len)
fmc(len+2) = monthdata(len)

!---

do i = 1, (len+1)

  d(i) = fmc(i+1) - fmc(i)

end do

!---

do i = 2, (len+1)

  m(i) = (d(i-1) + d(i)) / 2

  !---
  ! Monotonic adjustment to slope

  if(d(i-1) > 0 .and. d(i) < 0) then

    m(i) = 0

  else if(d(i-1) < 0 .and. d(i) > 0) then

    m(i) = 0

  end if

end do

m(1)     = d(1)
m(len+2) = d(len+1)

do i = 1, (len+1)

  if(d(i) == 0) then

    m(i) = 0.
    m(i+1) = 0.

  end if

end do

!------
! Calculate "second" wall control based on interception of Hermite functions

allocate(swc(len+1))

u = 0.5

H_00 = 1. + (u**2) * (2*u - 3)
H_10 = u * (u - 1) * (u - 1)
H_01 = (u**2) * (3 - 2*u)
H_11 = (u**2) * (u - 1)

do i = 1, (len+1)

  swc(i) = (fmc(i) * H_00) + (m(i) * H_10) + (fmc(i+1) * H_01) + (m(i+1) * H_11)

end do

!------
! Compile the NxN linear matrix to adjust mid-control points
allocate(mat(len,len))
allocate(solution(len))

!---

u = 1.

G_00 = u - (u**3) - (u**4) / 2.
G_10 = (u**2) * (3.*(u**2) - 8.*u + 6.) / 12.
G_01 = (u**3) * (1. - u / 2.)
G_11 = (u**3) * (3*u - 4.) / 12.

!---

mat = 0.

! Consider two "midpoints" outside of first and last interval

mat(1,1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.
mat(1,2) = 0.5 * G_11

mat(len,len-1) = -0.5 * G_10
mat(len,len)   = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

n = 1
do i = 2, (len-1)

  mat(i,n)   = -0.5 * G_10

  mat(i,n+1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

  mat(i,n+2) = 0.5 * G_11

  n = n + 1

end do

!---

solution(1) = 2. * monthdata(1) - swc(1) * G_00 - 0.5 * (swc(2) - swc(1)) * G_11 - 0.5 * (swc(2) - swc(1)) * G_10 - swc(2) *    &
                G_01 - swc(1) + 0.5 * G_10 * monthdata(1)

solution(len) =  2. * monthdata(len) - swc(len) * G_00 - 0.5 * (swc(len+1) - swc(len)) * G_11 - 0.5 * (swc(len+1) - swc(len)) *   &
                   G_10 - swc(len+1) * G_01 - swc(len) - 0.5 * G_11 * monthdata(len)

do i = 2, (len-1)

  solution(i) = 2. * monthdata(i) - swc(i) * G_00 - 0.5 * (swc(i+1) - swc(i)) * G_11 - 0.5 * (swc(i+1) - swc(i)) *   &
                 G_10 - swc(i+1) * G_01 - swc(i)

end do

! Solve linear system to get final mid control points
call matsol(mat,solution)

!------
! Compile "second" wall control with newly adjusted mid control points (all_cont)

allocate(all_solution(len+2))

all_solution(1)       = monthdata(1)
all_solution(2:len+1) = solution
all_solution(len+2)   = monthdata(len)

!---
allocate(all_cont(size(swc)+size(all_solution)))

n = 1
do i = 1, size(swc)

  all_cont(n)   = swc(i)
  all_cont(n+1) = all_solution(i+1)

  n = n + 2

end do

!---

all_cont(size(all_cont)) = swc(size(swc))

!------------------
! PART 2: Adjust monotonicty using the tridiagonal equation
!------------------

call mono_adjust(monthdata, all_cont)

!------------------
! PART 3: Adjustment to maximum and minimum bound
!------------------

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))
allocate(root(len))
allocate(root_days(len))

d_orig = -9999.
osc_check = .FALSE.
c2 = -9999.
root = -9999.
root_days = -9999

!------
! Assign -1 for negative slope, 1 for postive and 0 for turning point
do i = 2, (len-1)

  if ((monthdata(i+1) - monthdata(i)) < 0) then

    d_orig(i) = -1

  else if ((monthdata(i+1) - monthdata(i)) > 0) then

    d_orig(i) = 1

  end if

  !---

  if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

    d_orig(i) = 0

  else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

    d_orig(i) = 0

  end if

end do

!------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold
perc = plim / 100.

do i = 2, (len-1)

  if (d_orig(i) == 0 .and. monthdata(i) > 0) then

    if (all_cont(2*i) > (1+perc) * monthdata(i) .OR. all_cont(2*i) < (1-perc) * monthdata(i)) then

      osc_check(i) = .TRUE.

    end if

  else if (d_orig(i) == 0 .and. monthdata(i) < 0) then

    if (all_cont(2*i) < (1+perc) * monthdata(i) .OR. all_cont(2*i) > (1-perc) * monthdata(i)) then

      osc_check(i) = .TRUE.

    end if

  end if

end do

!------
! Calculate the amount of adjustment required and insert into the c2 variable
do i = 2, (len-1)

  if (osc_check(i)  .and. monthdata(i) > 0) then

    !---

    if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

      c2(i) = perc * monthdata(i)

    else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

      c2(i) = -perc * monthdata(i)

    end if

    !---

  else if (osc_check(i)  .and. monthdata(i) < 0) then

    !---

    if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

      c2(i) = -perc * monthdata(i)

    else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

      c2(i) = perc * monthdata(i)

    end if

  end if

end do

!------
! Calculate the root deviation based on c2 and triangular approximation

count = 0 ! Count the number of extra Hermite spline segments required

do i = 2, (len-1)

  if(osc_check(i) ) then

    diff_yi = monthdata(i) - all_cont((2*i)-1)
    diff_yi1 = monthdata(i) - all_cont((2*i)+1)

    !---

    top = 3. * (diff_yi + diff_yi1)
    bot = 8. * c2(i) + 3. * (diff_yi + diff_yi1)

    !---

    root(i) = top / bot

    root(i) = root(i) + 1. / real(nk(i))

    !---

    if (mod(nk(i),2) == 0) then ! EVEN month partitions

      u = 1. / nk(i)

      root_adj = 9999.

      !---

        l = 0
        do while (root_adj > 0.)

          root_adj = root(i) - u

          u = u + (2. / nk(i))

          l = l + 1

        end do

      root_days(i) = l - 1

    else if (mod(nk(i),2) /= 0) then ! ODD month partitions

      u = 0.

      root_adj = 9999.

      !---

        l = 0
        do while (root_adj > 0.)

          root_adj = root(i) - u

          u = u + (2. / real(nk(i)))

          l = l + 1

        end do

      root_days(i) = l - 1

    end if

    !---

    count = count + 2

  end if

end do

!------
! Re-estimate c2 based on integral area under Hermite spline

do i = 1, (len-1)

  if(osc_check(i) ) then

    !--- Construct fourth degree Hermite at u = 1 (integral [0,1])

    del = 1. - root(i)

    u = 1.

    G_00 = u - (u**3) - (u**4) / 2.
    G_01 = (u**3) * (1. - u/2.)
    G_10 = (U**2) * (3. * (u**2) - 8. *u + 6.) / 12.
    G_11 = (u**3) * (3.*u - 4.) / 12.

    !--- Assign local control points

    yi   = all_cont((2*i)-1)
    yi1  = monthdata(i)
    y2i  = monthdata(i)
    y2i1 = all_cont((2*i)+1)

    !--- Assign local slope

    mi   = ((yi - all_cont((2*i)-2)) / 1. + (all_cont(2*i) - yi) / 1.) / 2.
    m2i1 = ((all_cont((2*i)+2) - y2i1) / 1. + (y2i1 - all_cont(2*i)) / 1.) / 2.

    !--- Calculate new area approximation based on Hermite intergral

    area_total = 2. * del * monthdata(i)

    top = (G_00 * yi) + (G_10 * del * mi) + (G_01 * yi1) + (G_00 * y2i) + (G_01 * y2i1) + (G_11 * del * m2i1)
    bot = 1. + ((3. * del**2) / (2. * root(i)**2)) * G_11 - ((3. * del**2) / (2. * root(i)**2)) * G_10

    area_int = (area_total - del * (top + yi + y2i)) / bot

    !--- Re-assign c2 as the integral-estimated value

    c2(i) = (3. * area_int) / (4 * root(i))

  end if

end do

!------
! Generate x_new and y_new series that contains the extra quadratic adjusted segments

allocate(x_new((2*len)+2+count))
allocate(y_new((2*len)+2+count))

x_new(1:3) = [1.,2.,3.]
y_new(1:3) = all_cont(1:3)

nn = 4
mm = 4
kk = 4.

do i = 2, len

  if (.not.osc_check(i)) then

    x_new(nn)   = kk
    x_new(nn+1) = kk + 1.

    y_new(nn)   = all_cont(mm)
    y_new(nn+1) = all_cont(mm+1)

    nn = nn + 2

    mm = mm + 2

    kk = kk + 2.

  else if (osc_check(i) ) then

    x_new(nn)   = kk - root(i)
    x_new(nn+1) = kk
    x_new(nn+2) = kk + root(i)
    x_new(nn+3) = kk + 1.

    y_new(nn)   = monthdata(i)
    y_new(nn+1) = monthdata(i) + c2(i)
    y_new(nn+2) = monthdata(i)
    y_new(nn+3) = all_cont(mm+1)

    nn = nn + 4

    mm = mm + 2

    kk = kk + 2.

  end if

end do

x_new(nn) = kk
y_new(nn) = all_cont(mm)

!------
! Construct the spline for daily values based on all_cont

len_new = size(x_new)

allocate(d_new(len_new-1))
allocate(m_new(len_new))

do i = 1, (len_new-1)

  d_new(i) = (y_new(i+1) - y_new(i)) / (x_new(i+1) - x_new(i))

end do

!---

do i = 2, (len_new-1)

  m_new(i) = (d_new(i-1) + d_new(i)) / 2

  !---
  ! Monotonic adjustment to slope

  if(d_new(i-1) > 0 .and. d_new(i) < 0) then

    m_new(i) = 0

  else if(d_new(i-1) < 0 .and. d_new(i) > 0) then

    m_new(i) = 0

  end if

end do

m_new(1)       = d_new(1) / 2.
m_new(len_new) = d_new(len_new-1) / 2.

!---

do i = 1, (len_new-1)

  if(d_new(i) == 0) then

    m_new(i)   = 0.
    m_new(i+1) = 0.

  end if

end do

!------
! Reassign quadratic approximation slopes and original slopes to adjusted intervals

do i = 2, (len-1)

  if (osc_check(i) ) then

    ! slpe_l = findloc(x_new, (real(2*i)-root(i)), dim = 1) ! Find index in x_new = 2*i ('which' function in R)

    call findloc(x_new, (real(2*i)-root(i)), slpe_l)

    slpe_r = slpe_l + 2

    !---

    m_new(slpe_l) = 2. * c2(i) / root(i)

    m_new(slpe_r) = -2. * c2(i) / root(i)

    !---

    m_new(slpe_l-1) = ((all_cont(2*i-1) - all_cont((2*i)-2)) / 1. + (all_cont(2*i) - all_cont(2*i-1)) / 1.) / 2.

    m_new(slpe_r+1) = ((all_cont((2*i)+2) - all_cont(2*i+1)) / 1. + (all_cont(2*i+1) - all_cont(2*i)) / 1.) / 2.

    !---

  end if

end do

!---
! Ensure smooth monotonic interpolation between control points

call monocheck(d_new, m_new)

!------
!Assign the daily value based on sequence

n = 1
k = 1

do i = 1, len !outer loop start, for all monthly intervals N

  if (.not.osc_check(i)) then ! For intervals that did not required adjustment

    !---

    if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

      l = n + nk(i) - 1

      call days_even(nk(i), y_new(k:k+2), m_new(k:k+2), daydata(n:l))

      n = n + nk(i)

    else ! if odd months (or odd number of smaller time step)

      l = n + nk(i) - 1

      call days_odd(nk(i), y_new(k:k+2), m_new(k:k+2), daydata(n:l))

      n = n + nk(i)

    end if

    k = k + 2

  else if (osc_check(i) ) then ! for intervals adjusted by minmax bound

    if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

      l = n + nk(i) - 1

      call days_osc_even(nk(i), root_days(i), root(i), y_new(k:k+4), m_new(k:k+4), daydata(n:l))

      n = n + nk(i)

    else ! if odd months (or odd number of smaller time step)

      l = n + nk(i) - 1

      call days_osc_odd(nk(i), root_days(i), root(i), y_new(k:k+4), m_new(k:k+4), daydata(n:l))

      n = n + nk(i)

    end if

    k = k + 4

  end if

end do !end of outer loop

end subroutine newspline_pbound

! ------------------------------------------------------------------------------------------------------------------

subroutine newspline_bound_all(monthdata,nk,daydata,llim,ulim,plim)

implicit none

real(sp),    dimension(:), intent(in)  :: monthdata
integer(i4), dimension(:), intent(in)  :: nk
real(sp),    dimension(:), intent(out) :: daydata
real(sp),                  intent(in)  :: llim
real(sp),                  intent(in)  :: ulim
real(sp),                  intent(in)  :: plim ! taken in a 0-100

! Local variables for first mid control points
real(sp), dimension(:), allocatable :: fmc
real(sp), dimension(:), allocatable :: d
real(sp), dimension(:), allocatable :: m

! Local variables for wall control points
real(sp), dimension(:), allocatable :: swc

! Local variables for linear system of mid control adjustments
real(sp),    dimension(:,:), allocatable :: mat
real(sp),    dimension(:),   allocatable :: solution
real(sp),    dimension(:),   allocatable :: all_solution

integer(i4) :: len
integer(i4) :: i,n

! Final vector of all control points
real(sp), dimension(:), allocatable :: all_cont

! Hermite cubic quartic spline basis functions
real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11
real(sp) :: u

integer(i4) :: l

! Local variables for generating daily values after minmax bound adjustment of all_cont
real(sp), dimension(:), allocatable :: d_new
real(sp), dimension(:), allocatable :: m_new

integer(i4) :: len_new
integer(i4) :: slpe_l, slpe_r
integer(i4) :: k

! Local variables for max and min bound adjustments
integer(i4), dimension(:), allocatable :: d_orig
logical,     dimension(:), allocatable :: osc_check
real(sp),    dimension(:), allocatable :: c2
real(sp),    dimension(:), allocatable :: root
integer(i4), dimension(:), allocatable :: root_days
real(sp)                               :: perc

! Local variables for calculating root of quadratic approximation
real(sp)    :: diff_yi1
real(sp)    :: diff_yi
real(sp)    :: top, bot
real(sp)    :: root_adj
integer(i4) :: count

real(sp) :: del

real(sp) :: G_00
real(sp) :: G_01
real(sp) :: G_10
real(sp) :: G_11

real(sp) :: yi, yi1
real(sp) :: y2i, y2i1
real(sp) :: mi, m2i1
real(sp) :: area_total
real(sp) :: area_int

! Local variables for x_new and y_new
real(sp), dimension(:), allocatable :: x_new
real(sp), dimension(:), allocatable :: y_new

real(sp)    :: kk
integer(i4) :: nn, mm

!------------------
! PART 1: Determine all wall and mid control points
!------------------

! Start of the spline routine
len = size(monthdata)

allocate(fmc(len+2))
allocate(d(len+1))
allocate(m(len+2))

!------
! Define first mid-control point as equal to original monthdata
fmc(1) = monthdata(1)
fmc(2:(len+1)) = monthdata(1:len)
fmc(len+2) = monthdata(len)

!---

do i = 1, (len+1)

  d(i) = fmc(i+1) - fmc(i)

end do

!---

do i = 2, (len+1)

  m(i) = (d(i-1) + d(i)) / 2

  !---
  ! Monotonic adjustment to slope
  if(d(i-1) > 0 .and. d(i) < 0) then

    m(i) = 0

  else if(d(i-1) < 0 .and. d(i) > 0) then

    m(i) = 0

  end if

end do

m(1)     = d(1)
m(len+2) = d(len+1)

do i = 1, (len+1)

  if(d(i) == 0) then

    m(i) = 0.
    m(i+1) = 0.

  end if

end do

!------
! Calculate "second" wall control based on interception of Hermite functions
allocate(swc(len+1))

u = 0.5

H_00 = 1. + (u**2) * (2*u - 3)
H_10 = u * (u - 1) * (u - 1)
H_01 = (u**2) * (3 - 2*u)
H_11 = (u**2) * (u - 1)

do i = 1, (len+1)

  swc(i) = (fmc(i) * H_00) + (m(i) * H_10) + (fmc(i+1) * H_01) + (m(i+1) * H_11)

end do

!------
! Compile the NxN linear matrix to adjust mid-control points
allocate(mat(len,len))
allocate(solution(len))

!---

u = 1.

G_00 = u - (u**3) - (u**4) / 2.
G_10 = (u**2) * (3.*(u**2) - 8.*u + 6.) / 12.
G_01 = (u**3) * (1. - u / 2.)
G_11 = (u**3) * (3*u - 4.) / 12.

!---

mat = 0.

! Consider two "midpoints" outside of first and last interval
mat(1,1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.
mat(1,2) = 0.5 * G_11

mat(len,len-1) = -0.5 * G_10
mat(len,len)   = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

n = 1

do i = 2, (len-1)

  mat(i,n)   = -0.5 * G_10

  mat(i,n+1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11 + 1.

  mat(i,n+2) = 0.5 * G_11

  n = n + 1

end do

!---

solution(1) = 2. * monthdata(1) - swc(1) * G_00 - 0.5 * (swc(2) - swc(1)) * G_11 - 0.5 * (swc(2) - swc(1)) * G_10 - swc(2) *   &
                G_01 - swc(1) + 0.5 * G_10 * monthdata(1)

solution(len) = 2. * monthdata(len) - swc(len) * G_00 - 0.5 * (swc(len+1) - swc(len)) * G_11 - 0.5 * (swc(len+1) - swc(len)) *  &
                  G_10 - swc(len+1) * G_01 - swc(len) - 0.5 * G_11 * monthdata(len)

do i = 2, (len-1)

  solution(i) = 2. * monthdata(i) - swc(i) * G_00 - 0.5 * (swc(i+1) - swc(i)) * G_11 - 0.5 * (swc(i+1) - swc(i)) *    &
                  G_10 - swc(i+1) * G_01 - swc(i)

end do

! Solve linear system to get final mid control points

call matsol(mat, solution)

!------
! Compile "second" wall control with newly adjusted mid control points (all_cont)

allocate(all_solution(len+2))

all_solution(1)       = monthdata(1)
all_solution(2:len+1) = solution
all_solution(len+2)   = monthdata(len)

!---

allocate(all_cont(size(swc)+size(all_solution)))

n = 1
do i = 1, size(swc)

  all_cont(n)   = swc(i)
  all_cont(n+1) = all_solution(i+1)

  n = n + 2

end do

!---

all_cont(size(all_cont)) = swc(size(swc))

!------------------
! PART 2: Adjust monotonicty using the tridiagonal equation
!------------------

call mono_adjust(monthdata, all_cont)

!------------------
! PART 3: Adjustment to maximum and minimum bound
!------------------

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))
allocate(root(len))
allocate(root_days(len))

d_orig = -9999.
osc_check = .FALSE.
c2 = -9999.
root = -9999.
root_days = -9999

!------
! Assign -1 for negative slope, 1 for postive and 0 for turning point

do i = 2, (len-1)

  if ((monthdata(i+1) - monthdata(i)) < 0) then

    d_orig(i) = -1

  else if ((monthdata(i+1) - monthdata(i)) > 0) then

    d_orig(i) = 1

  end if

  !---

  if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

    d_orig(i) = 0

  else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

    d_orig(i) = 0

  end if

end do

!------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold

perc = plim / 100.

do i = 2, (len-1)

  !j = 2 * i

  if (d_orig(i) == 0 .and. monthdata(i) > 0) then

    if (all_cont(2*i) > (1+perc) * monthdata(i) .OR. all_cont(2*i) < (1-perc) * monthdata(i)) then

      osc_check(i) = .TRUE.

    end if

  else if (d_orig(i) == 0 .and. monthdata(i) < 0) then

    if (all_cont(2*i) < (1+perc) * monthdata(i) .OR. all_cont(2*i) > (1-perc) * monthdata(i)) then

      osc_check(i) = .TRUE.

    end if

  end if

end do

!------
! Calculate the amount of adjustment required and insert into the c2 variable

do i = 2, (len-1)

  if (osc_check(i)  .and. monthdata(i) > 0) then

    !---

    if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

      c2(i) = perc * monthdata(i)

    else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

      c2(i) = -perc * monthdata(i)

    end if

    !---

    if (all_cont(2*i) > ulim) then

      c2(i) = ulim - monthdata(i)

    end if

    !---

  else if (osc_check(i)  .and. monthdata(i) < 0) then

    !---

    if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

      c2(i) = -perc * monthdata(i)

    else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

      c2(i) = perc * monthdata(i)

    end if

    !---

    if (all_cont(2*i) < llim) then

      c2(i) = llim - monthdata(i)

    end if

    !---

  end if

end do

!------
! Calculate the root deviation based on c2 and triangular approximation

count = 0 ! Count the number of extra Hermite spline segments required

do i = 2, (len-1)

  if(osc_check(i)) then

    diff_yi = monthdata(i) - all_cont((2*i)-1)
    diff_yi1 = monthdata(i) - all_cont((2*i)+1)

    !---

    top = 3. * (diff_yi + diff_yi1)
    bot = 8. * c2(i) + 3. * (diff_yi + diff_yi1)

    !---

    root(i) = top / bot

    root(i) = root(i) + 1. / real(nk(i))

    !---

    if (mod(nk(i),2) == 0) then ! EVEN month partitions

      u = 1. / nk(i)

      root_adj = 9999.

      !---

        l = 0
        do while (root_adj > 0.)

          root_adj = root(i) - u

          u = u + (2. / nk(i))

          l = l + 1

        end do

      root_days(i) = l - 1

      !---

    else if (mod(nk(i),2) /= 0) then ! ODD month partitions

      u = 0.

      root_adj = 9999.

        l = 0
        do while (root_adj > 0.)

          root_adj = root(i) - u

          u = u + (2. / real(nk(i)))

          l = l + 1

        end do

      root_days(i) = l - 1

    end if

    count = count + 2

  end if

end do

!------
! Re-estimate c2 based on integral area under Hermite spline

do i = 1, (len-1)

  if(osc_check(i) ) then

    !--- Construct fourth degree Hermite at u = 1 (integral [0,1])

    del = 1 - root(i)

    u = 1.

    G_00 = u - (u**3) - (u**4) / 2.
    G_01 = (u**3) * (1. - u/2.)
    G_10 = (U**2) * (3. * (u**2) - 8. *u + 6.) / 12.
    G_11 = (u**3) * (3.*u - 4.) / 12.

    !--- Assign local control points

    yi   = all_cont((2*i)-1)
    yi1  = monthdata(i)
    y2i  = monthdata(i)
    y2i1 = all_cont((2*i)+1)

    !--- Assign local slope

    mi   = ((yi - all_cont((2*i)-2)) / 1. + (all_cont(2*i) - yi) / 1.) / 2.
    m2i1 = ((all_cont((2*i)+2) - y2i1) / 1. + (y2i1 - all_cont(2*i)) / 1.) / 2.

    !--- Calculate new area approximation based on Hermite intergral

    area_total = 2. * del * monthdata(i)

    top = (G_00 * yi) + (G_10 * del * mi) + (G_01 * yi1) + (G_00 * y2i) + (G_01 * y2i1) + (G_11 * del * m2i1)
    bot = 1. + ((3. * del**2) / (2. * root(i)**2)) * G_11 - ((3. * del**2) / (2. * root(i)**2)) * G_10

    area_int = (area_total - del * (top + yi + y2i)) / bot

    !--- Re-assign c2 as the integral-estimated value

    c2(i) = (3. * area_int) / (4. * root(i))

  end if

end do

!------
! Generate x_new and y_new series that contains the extra quadratic adjusted segments

allocate(x_new((2*len)+2+count))
allocate(y_new((2*len)+2+count))

x_new(1:3) = [1.,2.,3.]
y_new(1:3) = all_cont(1:3)

nn = 4
mm = 4
kk = 4.

do i = 2, len

  if (.not.osc_check(i)) then

    x_new(nn)   = kk
    x_new(nn+1) = kk + 1.

    y_new(nn)   = all_cont(mm)
    y_new(nn+1) = all_cont(mm+1)

    nn = nn + 2

    mm = mm + 2

    kk = kk + 2.

  else if (osc_check(i) ) then

    x_new(nn)   = kk - root(i)
    x_new(nn+1) = kk
    x_new(nn+2) = kk + root(i)
    x_new(nn+3) = kk + 1.

    y_new(nn)   = monthdata(i)
    y_new(nn+1) = monthdata(i) + c2(i)
    y_new(nn+2) = monthdata(i)
    y_new(nn+3) = all_cont(mm+1)

    nn = nn + 4

    mm = mm + 2

    kk = kk + 2.

  end if

end do

x_new(nn) = kk
y_new(nn) = all_cont(mm)

!------
! Construct the spline for daily values based on all_cont

len_new = size(x_new)

allocate(d_new(len_new-1))
allocate(m_new(len_new))

do i = 1, (len_new-1)

  d_new(i) = (y_new(i+1) - y_new(i)) / (x_new(i+1) - x_new(i))

end do

!---

do i = 2, (len_new-1)

  m_new(i) = (d_new(i-1) + d_new(i)) / 2.

  !---
  ! Monotonic adjustment to slope

  if(d_new(i-1) > 0 .and. d_new(i) < 0) then

    m_new(i) = 0

  else if(d_new(i-1) < 0 .and. d_new(i) > 0) then

    m_new(i) = 0

  end if

end do

m_new(1)       = d_new(1) / 2.
m_new(len_new) = d_new(len_new-1) / 2.

!---

do i = 1, (len_new-1)

  if(d_new(i) == 0) then

    m_new(i)   = 0.
    m_new(i+1) = 0.

  end if

end do

!------
! Reassign quadratic approximation slopes and original slopes to adjusted intervals

do i = 2, (len-1)

  if (osc_check(i) ) then

    call findloc(x_new, (real(2*i)-root(i)), slpe_l)

    slpe_r = slpe_l + 2

    !---

    m_new(slpe_l) = 2. * c2(i) / root(i)

    m_new(slpe_r) = -2. * c2(i) / root(i)

    !---

    m_new(slpe_l-1) = ((all_cont(2*i-1) - all_cont((2*i)-2)) / 1. + (all_cont(2*i) - all_cont(2*i-1)) / 1.) / 2.

    m_new(slpe_r+1) = ((all_cont((2*i)+2) - all_cont(2*i+1)) / 1. + (all_cont(2*i+1) - all_cont(2*i)) / 1.) / 2.

  end if

end do

!---
! Ensure smooth monotonic interpolation between control points

call monocheck(d_new, m_new)

!------
!Assign the daily value based on sequence

n = 1
k = 1

do i = 1, len !outer loop start, for all monthly intervals N

  if (.not.osc_check(i)) then ! For intervals that did not required adjustment

    if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

      l = n + nk(i) - 1

      call days_even(nk(i), y_new(k:k+2), m_new(k:k+2), daydata(n:l))

      n = n + nk(i)

    else ! if odd months (or odd number of smaller time step)

      l = n + nk(i) - 1

      call days_odd(nk(i), y_new(k:k+2), m_new(k:k+2), daydata(n:l))

      n = n + nk(i)

    end if

    k = k + 2

  else if (osc_check(i) ) then ! for intervals adjusted by minmax bound

    if(mod(nk(i),2) == 0) then !seperate into even or odd months, starting with EVEN

      l = n + nk(i) - 1

      call days_osc_even(nk(i), root_days(i), root(i), y_new(k:k+4), m_new(k:k+4), daydata(n:l))

      n = n + nk(i)

    else ! if odd months (or odd number of smaller time step)

      l = n + nk(i) - 1

      call days_osc_odd(nk(i), root_days(i), root(i), y_new(k:k+4), m_new(k:k+4), daydata(n:l))

      n = n + nk(i)

    end if

    k = k + 4

  end if

end do !end of outer loop

end subroutine newspline_bound_all

! ------------------------------------------------------------------------------------------------------------------

subroutine monocheck(d,m)

implicit none

real(sp), dimension(:), intent(in)    :: d
real(sp), dimension(:), intent(inout) :: m

!Local variables for adjusting monotonicity

real(sp), dimension(:), allocatable :: a,b,r
real(sp), dimension(:), allocatable :: check

integer(i4) :: len, i

!---

len = size(m)

allocate(a(len-1))
allocate(b(len-1))
allocate(r(len-1))
allocate(check(len-1))

a = 0.
b = 0.
r = 0.
check = 0.

!------

do i = 1, (len-1)

  if(d(i) /= 0) then

    a(i) = m(i) / d(i)

    b(i) = m(i+1) / d(i)

  end if

end do

!------

do i = 1, (len-1)

  if((a(i)**2 + b(i)**2) /= 0) then

    r(i) = 3. / ( (a(i)**2 + b(i)**2) ** 0.5 )

  end if

  check(i) = a(i)**2 + b(i)**2

end do

!------

do i = 1, (len-1)

  if(check(i) > 9.) then

    m(i)   = r(i) * a(i) * d(i)
    m(i+1) = r(i) * b(i) * d(i)

  end if

end do

end subroutine monocheck

! ------------------------------------------------------------------------------------------------------------------

subroutine mono_adjust(monthdata, all_cont)

implicit none

real(sp), dimension(:), intent(in)    :: monthdata
real(sp), dimension(:), intent(inout) :: all_cont

! Local variables for checking monotonicity
integer(i4), dimension(:), allocatable :: d_orig
logical    , dimension(:), allocatable :: slope_check

! Local 2x2 matrix system to determine new wall / mid-control points
real(sp)   , dimension(2,2) :: s_mat
real(sp)   , dimension(2)   :: s_sol

real(sp) :: G_00, G_10, G_01, G_11
real(sp) :: u

integer(i4) :: i, len

!------

len =  size(monthdata)

allocate(d_orig(len))
allocate(slope_check(len))

d_orig = -9999
slope_check = .TRUE.

s_mat = 0.
s_sol = 0.

!------

do i = 2, (len-1)

  if (monthdata(i+1) - monthdata(i) < 0) then

    d_orig(i) = -1

  else if (monthdata(i+1) - monthdata(i) > 0) then

    d_orig(i) = 1

  end if

end do

!---

do i = 2, (len-1)

  if (monthdata(i) > monthdata(i-1) .and. monthdata(i) > monthdata(i+1)) then

    d_orig(i) = 0

  else if (monthdata(i) < monthdata(i-1) .and. monthdata(i) < monthdata(i+1)) then

    d_orig(i) = 0

  end if

end do

!------
! Procedures to check monotonicity and adjust mid-control + left/right wall-control accordingly
! all_cont WILL BE ALTERED after this procedure
do i = 2, (len-1)

  if (d_orig(i) == 1) then

    if (all_cont(2*i) - all_cont(2*i-1) < 0. .OR. all_cont(2*i-1) - all_cont(2*i) < 0.) then

      slope_check(i) = .FALSE.

    end if

  else if (d_orig(i) == -1) then

    if (all_cont(2*i) - all_cont(2*i-1) > 0. .OR. all_cont(2*i-1) - all_cont(2*i) > 0.) then

      slope_check(i) = .FALSE.

    end if

  end if

end do

!---

u = 1.

G_00 = u - (u**3) - (u**4) / 2.
G_10 = (u**2) * (3.*(u**2) - 8.*u + 6.) / 12.
G_01 = (u**3) * (1. - u / 2.)
G_11 = (u**3) * (3*u - 4.) / 12.

!---

do i = 2, (len-1)

  if (.not.slope_check(i)) then

    if ((d_orig(i) == 1 .and. (all_cont(2*i) - all_cont(2*i-1)) < 0.) .or.   &
       (d_orig(i) == -1 .and. (all_cont(2*i) - all_cont(2*i-1)) > 0))  then

      all_cont(2*i) = monthdata(i)

      !---

      s_mat(1,1) = G_01 + 0.5 * G_10 + 0.5 * G_11
      s_mat(1,2) = 0.5 * G_11

      s_mat(2,1) = G_00 - 0.5 * G_10 - 0.5 * G_11 + 1.
      s_mat(2,2) = G_00 + G_01 + 0.5 * G_10 - 0.5 * G_11 + 1.

      !---

      s_sol(1) = 2. * monthdata(i) - (G_00 - 0.5 * G_10 - 0.5 * G_11 + 1) * all_cont(2*i-1) + 0.5 * G_10 * all_cont(2*i-2) -   &
                  (G_00 + G_01 + 0.5 * G_10 - 0.5 * G_11 + 1.) * all_cont(2*i)

      s_sol(2) = 2. * monthdata(i+1) - (G_01 + 0.5 * G_10 + 0.5 * G_11) * all_cont(2*i+3) + 0.5 * G_10 * all_cont(2*i) -   &
                   0.5 * G_11 * all_cont(2*i+4)

      !---

      call matsol(s_mat, s_sol)

      !---

      all_cont(2*i+1) = s_sol(1)

      all_cont(2*i+2) = s_sol(2)

      !---

    else if ((d_orig(i) == 1 .and. (all_cont(2*i+1) - all_cont(2*i)) < 0.) .OR. (d_orig(i) == -1 .and.   &
           (all_cont(2*i+1) - all_cont(2*i)) > 0))  then

      all_cont(2*i) = monthdata(i)

      !---

      s_mat(1,1) = G_00 + G_01 + 0.5 * G_10 - 0.5 * G_11 + 1.
      s_mat(1,2) = G_01 + 0.5 * G_10 + 0.5 * G_11

      s_mat(2,1) = -0.5 * G_10
      s_mat(2,2) = G_00 - 0.5 * G_10 - 0.5 * G_11 + 1.

      !---

      s_sol(1) = 2 * monthdata(i-1) - (G_00 - 0.5 * G_10 - 0.5 * G_11 + 1.) * all_cont(2*i-3) + 0.5 * G_10 * all_cont(2*i-4) -   &
                                       0.5 * G_11 * all_cont(2*i)

      s_sol(2) = 2 * monthdata(i) - (G_01 + 0.5 * G_10 + 0.5 * G_11) * all_cont(2*i+1) -   &
                                    (G_00 + G_01 + 0.5 * G_10 - 0.5 * G_11 + 1.) * all_cont(2*i) - 0.5 * G_11 * all_cont(2*i+2)

      !---

      call matsol(s_mat, s_sol)

      !---

      all_cont(2*i-2) = s_sol(1)

      all_cont(2*i-1) = s_sol(2)

      !---

    end if

  end if

end do

end subroutine mono_adjust

! ------------------------------------------------------------------------------------------------------------------

subroutine days_even(nk,y_val,m_val,daydata)

implicit none

integer(i4),               intent(in)    :: nk
real(sp),    dimension(:), intent(in)    :: y_val, m_val ! Take in three control points
real(sp),    dimension(:), intent(inout) :: daydata

! Local variable for generating day fraction from [0,1] for Hermite curve
integer(i4) :: i,n

real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11

real(sp) :: u

!------
u = 1. / nk

n = 1

do i = 1, (nk / 2)

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(1) * H_00 + m_val(1) * H_10 + y_val(2) * H_01 + m_val(2) * H_11

  u = u + (2. / nk)

  n = n + 1

end do

!---

u = 1. / nk

do i = 1, (nk / 2)

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(2) * H_00 + m_val(2) * H_10 + y_val(3) * H_01 + m_val(3) * H_11

  u = u + (2. / nk)

  n = n + 1

end do

end subroutine days_even

! ------------------------------------------------------------------------------------------------------------------

subroutine days_odd(nk,y_val,m_val,daydata)

implicit none

integer(i4),               intent(in)    :: nk
real(sp),    dimension(:), intent(in)    :: y_val, m_val ! Take in three control points
real(sp),    dimension(:), intent(inout) :: daydata

! Local variable for generating day fraction from [0,1] for Hermite curve

integer(i4) :: i,n

real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11

real(sp) :: u

!------

u = 1. / nk

n = 1

do i = 1, ((nk+1) / 2)

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(1) * H_00 + m_val(1) * H_10 + y_val(2) * H_01 + m_val(2) * H_11

  u = u + (2. / nk)

  n = n + 1

end do

!---

u = 2. / nk

do i = 1, ((nk-1) / 2)

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(2) * H_00 + m_val(2) * H_10 + y_val(3) * H_01 + m_val(3) * H_11

  u = u + (2. / nk)

  n = n + 1

end do

end subroutine days_odd

! ------------------------------------------------------------------------------------------------------------------

subroutine days_osc_even(nk,root_days,root,y_val,m_val,daydata)

implicit none

integer(i4),               intent(in)    :: nk
integer(i4),               intent(in)    :: root_days
real(sp),                  intent(in)    :: root
real(sp),    dimension(:), intent(in)    :: y_val, m_val ! Take in five control points
real(sp),    dimension(:), intent(inout) :: daydata

! Local variable for generating day fraction from [0,1] for Hermite curve
integer(i4) :: day_insd
integer(i4) :: day_outsd
integer(i4) :: i, n, num

real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11
real(sp) :: del
real(sp) :: a,b
real(sp) :: u

!------

day_insd = root_days

day_outsd = (nk / 2) - day_insd

del = 1 - root

!------
! LEFT INTERVAL (to the midpoint)

a = del - ((real(day_outsd) * 2. - 1.) / real(nk))

b = root - ((real(day_insd) * 2. - 1.) / real(nk))

!------
! Start with days outside quadratic partition on LEFT interval

n = 1

num = (2 * day_outsd) - 1

u = (1. / nk) / del

do i = 1, day_outsd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(1) * H_00 + del * m_val(1) * H_10 + y_val(2) * H_01 + del * m_val(2) * H_11

  u = u + (2. / nk) / del

  n = n + 1

end do

!---
! Days inside the quadratic partition on LEFT interval

num = (2 * day_insd) - 1

u = b / root

do i = 1, day_insd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(2) * H_00 + root * m_val(2) * H_10 + y_val(3) * H_01 + root * m_val(3) * H_11

  u = u + (2. / nk) / root

  n = n + 1

end do

!------
! RIGHT INTERVAL (to the midpoint)

a = root - ((real(day_insd) * 2. - 1.) / real(nk))

b = del - ((real(day_outsd) * 2. - 1.) / real(nk))

! Days inside quadratic partition on RIGHT interval

num = (2 * day_insd) - 1

u = (1. / nk) / root

do i = 1, day_insd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(3) * H_00 + root * m_val(3) * H_10 + y_val(4) * H_01 + root * m_val(4) * H_11

  u = u + (2. / nk) / root

  n = n + 1

end do

!---
! Days outside the quadratic partition on RIGHT interval

num = (2 * day_outsd) - 1

u = b / del

do i = 1, day_outsd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(4) * H_00 + del * m_val(4) * H_10 + y_val(5) * H_01 + del * m_val(5) * H_11

  u = u + (2. / nk) / del

  n = n + 1

end do

end subroutine days_osc_even

! ------------------------------------------------------------------------------------------------------------------

subroutine days_osc_odd(nk,root_days,root,y_val,m_val,daydata)

implicit none

integer(i4)           , intent(in)    :: nk
integer(i4)           , intent(in)    :: root_days
real(sp)              , intent(in)    :: root
real(sp), dimension(:), intent(in)    :: y_val, m_val ! Take in five control points
real(sp), dimension(:), intent(inout) :: daydata

! Local variable for generating day fraction from [0,1] for Hermite curve

integer(i4) :: day_insd, day_outsd
integer(i4) :: i, n, num

real(sp) :: H_00
real(sp) :: H_01
real(sp) :: H_10
real(sp) :: H_11
real(sp) :: del
real(sp) :: a,b
real(sp) :: u

!------
! LEFT INTERVAL (to the midpoint)

day_insd = root_days !+ 1

day_outsd = ((nk+1) / 2) - day_insd

del = 1. - root

!---

a = del - ((real(day_outsd) * 2. - 1.) / real(nk))

b = root - ((real(day_insd) * 2. - 2.) / real(nk))

!------
! Start with days outside quadratic partition on LEFT interval

n = 1

num = (2 * day_outsd) - 1

u = (1. / real(nk)) / del

do i = 1, day_outsd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(1) * H_00 + del * m_val(1) * H_10 + y_val(2) * H_01 + del * m_val(2) * H_11

  u = u + (2. / real(nk)) / del

  n = n + 1

end do

!---
! Days inside the quadratic partition on LEFT interval

num = (2 * day_insd) - 2

u = b / root

do i = 1, day_insd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(2) * H_00 + root * m_val(2) * H_10 + y_val(3) * H_01 + root * m_val(3) * H_11

  u = u + (2. / real(nk)) / root

  n = n + 1

end do

!------
! RIGHT INTERVAL (to the midpoint)

day_insd = root_days - 1

day_outsd = ((nk-1) / 2) - day_insd

!------

a = root - ((real(day_insd) * 2.) / real(nk))

b = del - ((real(day_outsd) * 2. - 1.) / real(nk))

! Days inside quadratic partition on RIGHT interval

num = (2 * day_insd)

u = (2. / real(nk)) / root

do i = 1, day_insd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(3) * H_00 + root * m_val(3) * H_10 + y_val(4) * H_01 + root * m_val(4) * H_11

  u = u + (2. / real(nk)) / root

  n = n + 1

end do

!---
! Days outside the quadratic partition on RIGHT interval

num = (2 * day_outsd) - 1

u = b / del

do i = 1, day_outsd

  H_00 = 1. + (u**2) * (2*u - 3)
  H_01 = (u**2) * (3 - 2*u)
  H_10 = u * (u - 1) * (u - 1)
  H_11 = (u**2) * (u - 1)

  daydata(n) = y_val(4) * H_00 + del * m_val(4) * H_10 + y_val(5) * H_01 + del * m_val(5) * H_11

  u = u + (2. / real(nk)) / del

  n = n + 1

end do

end subroutine days_osc_odd

! ------------------------------------------------------------------------------------------------------------------

end module newsplinemod
