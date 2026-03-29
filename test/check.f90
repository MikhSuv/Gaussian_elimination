program check
  use precision_mod
  use Gaussian_elimination
implicit none

integer :: n = 1000
integer :: i, j
real(dp), allocatable :: A(:,:), B(:), X(:)
real(dp) :: r(6)

allocate(A(n, n), B(n))
call random_number(A)
call random_number(B)

X = solve_linear_system(A, B, "gauss")
r(1)= residual(A, B, X)
Deallocate(X)
X = solve_linear_system(A, B, "jordan")
r(2)= residual(A, B, X)
Deallocate( X)
X = solve_linear_system(A, B, "choice")
r(3)= residual(A, B, X)
Deallocate(A, B, X)

allocate(A(n, n), B(n))
forall (i = 1:n, j = 1:n) A(i,j) = 1.0_dp/ (real(i + j - 1, kind = dp))
call RANDOM_NUMBER(B)
X = solve_linear_system(A, B, "gauss")
r(4)= residual(A, B, X)
Deallocate(X)
X = solve_linear_system(A, B, "jordan")
r(5)= residual(A, B, X)
Deallocate(X)
X = solve_linear_system(A, B, "choice")
r(6)= residual(A, B, X)
Deallocate(A, B, X)

print *, "Модули векторов невзяки"
print *, "Случайная матрица"
print *, "Метод Гаусса: ", r(1)
print *, "Метод Йордана: ", r(2)
print *, "Метод Гаусса с выбором элемента: ", r(3)
print *, "Матрица Гильберта"
print *, "Метод Гаусса: ", r(4)
print *, "Метод Йордана: ", r(5)
print *, "Метод Гаусса с выбором элемента: ", r(6)
end program check
