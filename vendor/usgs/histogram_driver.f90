! Independent harness for the SELMVS histogram statistics used in calibration.
program histogram
  implicit none
  integer*2 :: values(3720),dummy(3720),n,i,idum,ib,il
  real*4 :: se,av,mode,median
  common /easy/ idum,se,ib,il,av,mode,median
  read(*,*) n
  do i=1,n
    read(*,*) values(i)
  end do
  call mstat(values,n,1_2,10_2,0_2,0_2,0_2,dummy)
  write(*,'(3ES24.16)') av,se,mode
end program
