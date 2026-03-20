module Gaussian_elimination
  use precision_mod
  implicit none
  private

  public :: read_linear_system, write_result

contains

  subroutine read_linear_system(filename, A, B, n)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: A(:, :), B(:) ! СЛУ (A|B)
    integer, intent(out), optional :: n

    integer :: iunit, iostatus, i
    character(len=256) :: line ! для чтения первой строки

    open(newunit=iunit, file=filename, status='old', &
    action = 'read', iostat=iostatus)
    if (iostatus /= 0) then
      error stop 'Error occured while opening file'
    end if

    read(iunit, '(a)', iostat=iostatus) line

    if (iostatus /= 0) then
      error stop 'Error occured while reading line'
    end if
    read(line(2:), *) n 

    ! Выделение памяти под матрицу A
    allocate(A(n,n))
    ! Выделение памяти под столбец B
    allocate(B(n))
    do i = 1, n
      read(iunit, *) A(i, :)
    end do

    do i = 1, n
      read(iunit, *) B(i)
    end do

    close(iunit)

  end subroutine read_linear_system

  subroutine write_result(filename, X)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) ::  X(:) ! Столбец результата
    
    integer :: n
    integer :: ounit, i, iostatus

    n = size(X)
    open(newunit=ounit, file=filename, action='write', iostat=iostatus)
    if (iostatus /= 0) then
      error stop 'Error occured while opening file'
    end if
    
    write(ounit, '("# ", i0)') n

    do i = 1, n
      write(ounit, '(*(e23.15, 1x))') X(i)
    end do

    close(ounit)
    
  end subroutine write_result
end module Gaussian_elimination
