program check
  use precision_mod
  use Gaussian_elimination
implicit none

integer :: n = 3000
integer :: i, j
real(dp), allocatable :: A(:,:), B(:), X(:)

allocate(A(n, n), B(n))
call random_number(A)
call random_number(B)

X = solve_linear_system(A, B, "gauss")
print *, residual(A, B, X)
Deallocate(X)
X = solve_linear_system(A, B, "jordan")
print *, residual(A, B, X)
Deallocate( X)
X = solve_linear_system(A, B, "choice")
print *, residual(A, B, X)
Deallocate(A, B, X)

allocate(A(n, n), B(n))
forall (i = 1:n, j = 1:n) A(i,j) = 1.0_dp/ (real(i + j - 1, kind = dp))
call RANDOM_NUMBER(B)
X = solve_linear_system(A, B, "gauss")
print *, residual(A, B, X)
Deallocate(X)
X = solve_linear_system(A, B, "jordan")
print *, residual(A, B, X)
Deallocate(X)
X = solve_linear_system(A, B, "choice")
print *, residual(A, B, X)
Deallocate(A, B, X)
end program check
