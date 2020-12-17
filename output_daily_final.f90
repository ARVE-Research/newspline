program output_daily_final

use parametersmod
use newsplinemod,   only : newspline_all

implicit none

real(sp)    , dimension(53) :: temp_F, temp_C
integer(i4),  dimension(53) :: date
integer(i4),  dimension(38) :: nk

real(sp), dimension(38)             :: temp
real(sp), dimension(:), allocatable :: daydata

integer(i4)  :: summ, i, j, n
character(6) :: header

!------
! Read temperature data file
open(10, file = "captaincook_monthly_temp_NOAA.txt", status = 'old')

read(10,*) header

do i = 1, 53

  read(10,*) date(i), temp_F(i), temp_C(i)

end do

!---
! Define number of days in each month
nk = [31,31,28,31,30,31,30,31,31,30,31,30,31,  &
      31,29,31,30,31,30,31,31,30,31,30,31,  &
      31,28,31,30,31,30,31,31,30,31,30,31,31]

!---
! Initiate temperature array. Copy first and last month to start and end.
summ = sum(nk)

allocate(daydata(summ))

daydata = 0.

temp(1) = temp_C(1)
temp(2:37) = temp_C(1:36)
temp(38) = temp_C(36)

!---
! Interpolation without any limit options
call newspline_all(temp, nk, daydata)

!---
! Write 36 months (1096 days) of daily values
do i = 32, summ-31

  write(*,*) daydata(i)

end do

end program output_daily_final
