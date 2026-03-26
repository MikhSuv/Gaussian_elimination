program main
  use precision_mod
  use Gaussian_elimination
  implicit none
  real(dp), allocatable :: A(:,:), B(:), T(:)
  integer :: n

  call read_linear_system("data.dat", A, B, n)
  call write_result("result.dat", B)
  T = gaussian_solve(A, b) 
  print *, T

end program main
