module Gaussian_elimination
   use precision_mod
   implicit none
   ! private

   ! public :: read_linear_system, write_result

contains

   subroutine read_linear_system(filename, A, B, n)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: A(:, :), B(:) ! СЛУ (A|B)
      integer, intent(out), optional :: n

      integer :: iunit, iostatus, i
      character(len=256) :: line ! для чтения первой строки

      open (newunit=iunit, file=filename, status='old', &
            action='read', iostat=iostatus)
      if (iostatus /= 0) then
         error stop 'Error occured while opening file'
      end if

      read (iunit, '(a)', iostat=iostatus) line

      if (iostatus /= 0) then
         error stop 'Error occured while reading line'
      end if
      read (line(2:), *) n

      ! Выделение памяти под матрицу A
      allocate (A(n, n))
      ! Выделение памяти под столбец B
      allocate (B(n))
      do i = 1, n
         read (iunit, *) A(i, :)
      end do

      do i = 1, n
         read (iunit, *) B(i)
      end do

      close (iunit)

   end subroutine read_linear_system

   subroutine write_result(filename, X)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) ::  X(:) ! Столбец результата

      integer :: n
      integer :: ounit, i, iostatus

      n = size(X)
      open (newunit=ounit, file=filename, action='write', iostat=iostatus)
      if (iostatus /= 0) then
         error stop 'Error occured while opening file'
      end if

      write (ounit, '("# ", i0)') n

      do i = 1, n
         write (ounit, '(*(e23.15, 1x))') X(i)
      end do

      close (ounit)

   end subroutine write_result

   ! function solve_linear_system(A, B, methjod) result(X)
   !   real(dp), intent(in) :: A(:,:), B(:)
   !   real(dp), allocatable :: X(:)
   !   character(len=*), intent(in) :: method
   !
   !   integer :: n
   !   n = size(B)
   !
   ! end function gaussian_solve

   function gaussian_solve(AA, B) result(X)
      real(dp), intent(in) :: AA(:, :), B(:)
      real(dp), allocatable :: X(:)

      real(dp), allocatable :: A(:, :) ! расширенная матрица системы
      real(dp) :: leading_entry, factor

      integer :: n
      integer :: i, k
      n = size(B)
      allocate (A(n, n + 1))
      A(:, 1:n) = AA
      A(:, n + 1) = B
      ! Прамой ход
      do k = 1, n
         leading_entry = a(k, k)
         if (abs(leading_entry) <= eps) then
            print *, "WARNING: the leading_entry in step", k, " is too small: ", leading_entry
         end if

         a(k, k:n + 1) = a(k, k:n + 1)/leading_entry

         do concurrent(i=k + 1:n)
            a(i, k:n + 1) = a(i, k:n + 1) - a(i, k)*a(k, k:n + 1)
         end do
      end do

      ! Обратный ход
      allocate(X(n))
      do i = n, 1, -1
      x(i) = a(i, n+1) - dot_product(a(i, i+1:n), x(i+1:n))
      end do
      
   end function gaussian_solve
end module Gaussian_elimination
