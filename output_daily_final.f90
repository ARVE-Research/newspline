program output_daily_final

use nrtype
use nrutil
use newsplinemod

implicit none

real(sp)    , dimension(53) :: temp_F, temp_C
integer(I4B), dimension(53) :: date
integer(I4B), dimension(36) :: nk

real(sp), dimension(:), allocatable :: daydata

integer(I4B) :: summ, i, j, n

character(6) :: header

!------

open(10, file = "captaincook_monthly_temp_NOAA.txt", status = 'old')

read(10,*) header

do i = 1, 53

  read(10,*) date(i), temp_F(i), temp_C(i)

end do

!---

nk = [31,28,31,30,31,30,31,31,30,31,30,31,  &
      31,29,31,30,31,30,31,31,30,31,30,31,  &
      31,28,31,30,31,30,31,31,30,31,30,31]

!---

summ = sum(nk)

allocate(daydata(summ))

daydata = 0.

!---

call newspline_all(temp_C(1:36), nk, daydata)

!---

 do i = 1, summ

    print *, daydata(i)

 end do


end program output_daily_final
