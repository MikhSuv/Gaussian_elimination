program main
  use precision_mod
  use Gaussian_elimination
  implicit none
  real(dp), allocatable :: A(:,:), B(:)
  integer :: n

  call read_linear_system("data.dat", A, B, n)
  print *, A
  print *, B
  print *, n

end program main
