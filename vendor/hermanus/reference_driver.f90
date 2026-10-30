! Standalone reference harness: does not link the Python adapter or C ABI wrapper.
program reference
  implicit none
  real :: signal(1440), raw(1440), smooth(1440), mu, qfactor
  integer :: qs, qe, seed, i
  integer(kind=8) :: state
  common /reference_rng/ state
  character :: qi
  read(*,*) mu, qi, qfactor, qs, qe, seed
  state = seed
  do i=1,1440
    read(*,*) signal(i)
  end do
  call kfilter(signal, raw, smooth, mu, qi, qfactor, qs, qe)
  do i=1,1440
    write(*,'(2ES25.16)') raw(i), smooth(i)
  end do
end program

function random() result(value)
  implicit none
  real(kind=8) :: value
  integer(kind=8) :: state
  common /reference_rng/ state
  state = mod(16807_8 * state, 2147483647_8)
  value = 2.0_8 * dble(state) / 2147483647.0_8 - 1.0_8
end function
