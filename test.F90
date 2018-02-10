module mod1
  implicit none
!  private ::gen_rand, madd_, msub_, mm_

  interface operator(+)
     procedure madd_
  end interface operator(+)

  interface operator(-)
     procedure msub_
  end interface operator(-)
  
  interface operator(*)
     procedure mm_
  end interface operator(*)

  type matrix
     integer :: size
     real(8),allocatable,dimension(:, :) :: mat
   contains
     generic   :: init          => init0, init1
     procedure :: init0         => init0_, init1 => init1_
     procedure :: set_mat       => set_mat_
     procedure :: show_mat      => show_mat_
     procedure :: equal         => equal_
     generic   :: assignment(=) => equal
     final     :: fini
  end type matrix
  
contains
  subroutine gen_rand(size, a, val_min, val_max, seed)
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    integer,intent(in) :: size
    real(8),dimension(size, size),intent(out) :: a
    real(8),intent(in) :: val_min, val_max
    integer,intent(in) :: seed
    integer::ierr
    integer::brng, method
    type(vsl_stream_state)::stream

    brng   = vsl_brng_mt19937
    !     for older intel MKL
    !     method = vsl_method_duniform_std
    method = vsl_rng_method_uniform_std_accurate

    ierr = vslnewstream(stream, brng, seed)
    ierr = vdrnguniform(method, stream, size, a, val_min, val_max)
    ierr = vsldeletestream(stream)
  end subroutine gen_rand

  subroutine init0_(this)
    implicit none
    class(matrix),intent(inout) :: this
    integer :: size

    size = 8
    this%size = size
    allocate(this%mat(size, size))
    !$omp parallel
    !$omp workshare
    this%mat(:, :) = 0.0d0
    !$omp end workshare
    !$omp end parallel
  end subroutine init0_

  subroutine init1_(this, size)
    implicit none
    class(matrix),intent(inout) :: this
    integer,intent(in) :: size

    this%size = size
    allocate(this%mat(size, size))
    this%mat(:, :) = 0.0d0
  end subroutine init1_

  subroutine set_mat_(this, val_min, val_max, seed)
    implicit none
    class(matrix),intent(inout) :: this
    real(8),intent(in) :: val_min, val_max
    integer,intent(in) :: seed

    call gen_rand(this%size**2, this%mat, val_min, val_max, seed)
  end subroutine set_mat_
  
  subroutine show_mat_(this)
    class(matrix),intent(in) :: this
    integer :: i, j

    do i = 1, this%size
       do j = 1, this%size
          write(6, '(1pe14.5)', advance="no") this%mat(i, j)
       end do
       write(6, *)
    end do
  end subroutine show_mat_
  
  function madd_(a, b) result(c)
    implicit none
    type(matrix),intent(in)  :: a, b
    type(matrix) :: c
    integer :: i, j
    integer :: size

    size = a%size
    call c%init(a%size)
    !$omp parallel do
    do j = 1, size
       do i = 1, size
          c%mat(i, j) = a%mat(i, j) + b%mat(i, j)
       end do
    end do
  end function madd_

  function msub_(a, b) result(c)
    implicit none
    type(matrix),intent(in)  :: a, b
    type(matrix) :: c
    integer :: i, j
    integer :: size

    size = a%size
    call c%init(a%size)
    !$omp parallel do
    do j = 1, size
       do i = 1, size
          c%mat(i, j) = a%mat(i, j) - b%mat(i ,j)
       end do
    end do
  end function msub_
  
  function mm_(a, b) result(c)
    implicit none
    type(matrix),intent(in)  :: a, b
    type(matrix) :: c
    integer :: i, j, k
    integer :: size

    size = a%size
    call c%init(a%size)
    !$omp parallel do
    do j = 1, size
       do k = 1, size
          do i = 1, size
             c%mat(i, j) = c%mat(i, j) + a%mat(i, k)*b%mat(k ,j)
          end do
       end do
    end do
  end function mm_

  subroutine equal_(lhs, rhs)
    implicit none
    class(matrix),intent(out) :: lhs
    class(matrix),intent(in)  :: rhs

    call lhs%init1(rhs%size)
    
    !$omp parallel
    !$omp workshare
    lhs%mat(:, :) = rhs%mat(:, :)
    !$omp end workshare
    !$omp end parallel
  end subroutine equal_

  subroutine fini(this)
    implicit none
    type(matrix),intent(inout) :: this

    this%size = 0
    if (allocated(this%mat)) deallocate(this%mat)
  end subroutine fini
end module mod1

program main
  use mod1
  implicit none
  type(matrix) :: a, b, c, d
  integer :: size
  character(len=16) :: argv1
  integer :: i

  if (iargc() == 0) then
     size = 4
  else
     call getarg(1, argv1)
     read(argv1, *) size
  end if
  write(6, *) "size:", size

  call a%init(size)
  call b%init(size)
  call a%set_mat(-1.0d0, 1.0d0, 5555)
  call b%set_mat(-2.0d0, 2.0d0, 7777)

  c = a*b
  d = c+a-b

  write(6, *) "A:"
  call a%show_mat
  write(6, *) "B:"
  call b%show_mat
  write(6, *) "C = A*B:"
  call c%show_mat
  write(6, *) "D = C+A-B:"
  call d%show_mat
  
  stop
end program main

