program output_daily_final

use nrtype
use nrutil
use mcsplinemod
use rmsmoothmod
use newsplinemodfinal

implicit none

real(sp)    , dimension(53) :: temp_F, temp_C
integer(I4B), dimension(53) :: date
real(sp)    , dimension(2)  :: bcond
integer(I4B), dimension(36) :: nk

real(sp), dimension(:), allocatable :: daydata_ns, daydata_mc, daydata_rm
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

nk = nk * 1

!---

bcond(1) = temp_C(1)
bcond(2) = temp_C(36)

summ = sum(nk)
allocate(daydata_mc(summ))
allocate(daydata_ns(summ))
allocate(daydata_rm(summ))

daydata_mc = 0.
daydata_ns = 0.
daydata_rm = 0.

!---

call mcspline(temp_C(1:36), nk, daydata_mc)

call newspline(temp_C(1:36), nk, daydata_ns)

call rmsmooth(temp_C(1:36), nk, bcond, daydata_rm)

!---

 do i = 1, summ
   print *, daydata_mc(i)
 end do

!---

 do i = 1, summ
    print *, daydata_ns(i)
 end do

!---

 do i = 1, summ
   print *, daydata_rm(i)
 end do


end program output_daily_final
