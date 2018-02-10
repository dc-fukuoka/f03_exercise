module mod1
  implicit none
  type mm
     integer :: size
     real(8),allocatable,dimension(:, :) :: a, b, c
     contains
       generic   :: init  => init0, init1
       procedure :: init0 => init_, init1 => init__
       procedure :: set_mat => set_mat_
       procedure :: mm    => mm_
       procedure :: show_mat => show_mat_
       final   :: fini
    end type mm
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
    
  subroutine init_(this)
    implicit none
    class(mm),intent(inout) :: this
    integer :: size

    size = 8
    this%size = size
    allocate(this%a(size, size), this%b(size, size), this%c(size, size))
    this%a(:, :) = 0.0d0
    this%b(:, :) = 0.0d0
    this%c(:, :) = 0.0d0
  end subroutine init_
  
  subroutine init__(this, size)
    implicit none
    class(mm),intent(inout) :: this
    integer,intent(in) :: size

    this%size = size
    allocate(this%a(size, size), this%b(size, size), this%c(size, size))
    this%a(:, :) = 0.0d0
    this%b(:, :) = 0.0d0
    this%c(:, :) = 0.0d0
  end subroutine init__

  subroutine set_mat_(this, mat, val_min, val_max, seed)
    implicit none
    class(mm),intent(inout) :: this
    character(len=1),intent(in) :: mat
    real(8),intent(in) :: val_min, val_max
    integer,intent(in) :: seed

    if (mat == 'a') then
       call gen_rand(this%size**2, this%a, val_min, val_max, seed)
    else if (mat == 'b') then
       call gen_rand(this%size**2, this%b, val_min, val_max, seed)
    else if (mat == 'c') then
       call gen_rand(this%size**2, this%c, val_min, val_max, seed)
    end if
  end subroutine set_mat_
  subroutine show_mat_(this, mat)
    class(mm),intent(in) :: this
    character(len=1) :: mat
    integer :: i, j

    if (mat == 'a') then
       write(6, *) "A:"
       do i = 1, this%size
          do j = 1, this%size
             write(6, '(1pe14.5)', advance="no") this%a(i, j)
          end do
          write(6, *)
       end do
    else if (mat == 'b') then
       write(6, *) "B:"
       do i = 1, this%size
          do j = 1, this%size
             write(6, '(1pe14.5)', advance="no") this%b(i, j)
          end do
          write(6, *)
       end do
    else if (mat == 'c') then
       write(6, *) "C:"
       do i = 1, this%size
          do j = 1, this%size
             write(6, '(1pe14.5)', advance="no") this%c(i, j)
          end do
          write(6, *)
       end do
    end if
    
  end subroutine show_mat_

  subroutine mm_(this)
    implicit none
    class(mm),intent(inout) :: this
    integer :: i, j, k

    !$omp parallel do
    do j = 1, this%size
       do k = 1, this%size
          do i = 1, this%size
             this%c(i, j) = this%c(i, j) + this%a(i, k)*this%b(k ,j)
          end do
       end do
    end do
  end subroutine mm_

  subroutine fini(this)
    implicit none
    type(mm),intent(inout) :: this

    this%size = 0
    if (allocated(this%a)) deallocate(this%a)
    if (allocated(this%b)) deallocate(this%b)
    if (allocated(this%c)) deallocate(this%c)
  end subroutine fini
end module mod1

program main
  use mod1
  implicit none
  type(mm) :: mm1
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
  
  call mm1%init(size)
  call mm1%set_mat('a', -1.0d0, 1.0d0, 5555)
  call mm1%set_mat('b', -2.0d0, 2.0d0, 7777)
  call mm1%mm
  call mm1%show_mat('a')
  call mm1%show_mat('b')
  call mm1%show_mat('c')
  stop
end program main
