module matrix_calc
  implicit none

  integer,parameter  :: dp = kind(1.0d0)
  real(dp),parameter :: tol = 1.0d-10
  
  interface operator(+)
     procedure madd_
  end interface operator(+)

  interface operator(-)
     procedure msub_
  end interface operator(-)
  
  interface operator(*)
     procedure mm_
  end interface operator(*)

  interface operator(/)
     procedure mdiv_
  end interface operator(/)
#if 0
  interface assignment(=)
     procedure equal_
  end interface assignment(=)
#endif
  interface matrix
     module procedure init0_
     module procedure init1_
  end interface matrix

  private :: lu_decomp, inverse

  type matrix
     integer :: size
     real(dp),allocatable,dimension(:, :) :: mat
   contains
     generic          :: init          => init0, init1
     procedure,nopass :: init0         => init0_
     procedure,nopass :: init1         => init1_
     procedure        :: set_mat       => set_mat_
     procedure        :: show_mat      => show_mat_
#if 0
     procedure        :: equal         => equal_
     generic          :: assignment(=) => equal
#endif
     final            :: fini
  end type matrix
  
contains
  subroutine gen_rand(size, a, val_min, val_max, seed)
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(out) :: a
    real(dp),intent(in) :: val_min, val_max
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

  function init0_() result(this)
    implicit none
    type(matrix) :: this
    integer :: size

    size = 8
    this%size = size
    allocate(this%mat(size, size))
    !$omp parallel
    !$omp workshare
    this%mat(:, :) = 0.0d0
    !$omp end workshare
    !$omp end parallel
  end function init0_

  function init1_(size) result(this)
    implicit none
    type(matrix) :: this
    integer,intent(in) :: size

    this%size = size
    allocate(this%mat(size, size))
    !$omp parallel
    !$omp workshare
    this%mat(:, :) = 0.0d0
    !$omp end workshare
    !$omp end parallel
  end function init1_

  subroutine set_mat_(this, val_min, val_max, seed)
    implicit none
    class(matrix),intent(inout) :: this
    real(dp),intent(in) :: val_min, val_max
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
    c = matrix(size)
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
    c = matrix(size)
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
    c = matrix(size)
    !$omp parallel do
    do j = 1, size
       do k = 1, size
          do i = 1, size
             c%mat(i, j) = c%mat(i, j) + a%mat(i, k)*b%mat(k ,j)
          end do
       end do
    end do
  end function mm_
#if 0
  subroutine equal_(lhs, rhs)
    implicit none
    class(matrix),intent(out) :: lhs
    class(matrix),intent(in)  :: rhs

!    lhs = init1_(rhs%size)
    
    !$omp parallel
    !$omp workshare
    lhs%mat(:, :) = rhs%mat(:, :)
    !$omp end workshare
    !$omp end parallel
  end subroutine equal_
#endif

  ! ref: http://workspacememory.hatenablog.com/entry/2017/03/01/173753
  subroutine lu_decomp(size, a, ipivot, lu)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    integer,dimension(size),intent(out) :: ipivot
    real(dp),dimension(size, size),intent(out) :: lu
    integer :: i, j, k
    integer :: ip, tmp_ip
    real(dp) :: tmp, max0, w

    !$omp parallel
    !$omp workshare
    lu = a
    !$omp end workshare
    !$omp do
    do i = 1, size
       ipivot(i) = i
    end do
    !$omp end do
    !$omp end parallel
    do k = 1, size-1
       max0 = abs(lu(k, k))
       ip = k
       ! this loop is impossible to parallelize for me
       do i = k+1, size
          tmp = abs(lu(i, k))
          if (tmp > max0) then
             max0 = tmp
             ip  = i
          end if
       end do       

       if (max0 <= tol) then
          write(6, *) "one of diagonal component is smaller than", tol
          stop
       end if

       if (ip .ne. k) then
          !$omp parallel do private(tmp)
          do j = k, size
             tmp       = lu(ip, j)
             lu(ip, j) = lu(k,  j)
             lu(k,  j) = tmp
          end do
          tmp_ip     = ipivot(ip)
          ipivot(ip) = ipivot(k)
          ipivot(k)  = tmp_ip
          !$omp parallel do private(tmp)
          do j = 1, k-1
             tmp       = lu(k, j)
             lu(k,  j) = lu(ip, j)
             lu(ip, j) = tmp
          end do
       end if
       !$omp parallel do private(w)
       do i = k+1, size
          w        = lu(i, k)/lu(k, k)
          lu(i, k) = w

          do j = k+1, size
             lu(i, j) = lu(i, j) - w*lu(k, j)
          end do
       end do
    end do

    ! write(100, *) ipivot ! for a debug
    
  end subroutine lu_decomp

  subroutine inverse(size, a, a_inv)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in)  :: a
    real(dp),dimension(size, size),intent(out) :: a_inv
    real(dp),dimension(size, size) :: lu
    integer,dimension(size) :: ipivot
    integer :: i, j, k
    real(dp),dimension(size) :: unit_vec, y
    real(dp) :: tmp

    !$omp parallel
    !$omp workshare
    lu     = 0.0d0
    ipivot = 0
    a_inv  = 0.0d0
    !$omp end workshare
    !$omp end parallel
    call lu_decomp(size, a, ipivot, lu)
    
    do k = 1, size
       !$omp parallel
       !$omp workshare
       unit_vec = 0.0d0
       !$omp end workshare
       !$omp end parallel
       unit_vec(k) = 1.0d0
       
       ! forward substitution
       y(1) = unit_vec(ipivot(1))
       do i = 2, size
          tmp = 0.0d0
          !$omp parallel do reduction(+:tmp)
          do j = 1, i-1
             tmp = tmp + lu(i, j)*y(j)
          end do
          y(i) = unit_vec(ipivot(i)) - tmp
       end do
       
       ! backward substitution
       a_inv(size, k) = y(size)/lu(size, size)
       do i = size-1, 1, -1
          tmp = 0.0d0
       !$omp parallel do reduction(+:tmp)
          do j = i+1, size
             tmp = tmp + lu(i, j)*a_inv(j, k)
          end do
          a_inv(i, k) = (y(i) - tmp)/lu(i, i)
       end do
    end do
  end subroutine inverse

  function mdiv_(a, b) result(c)
    type(matrix),intent(in)  :: a, b
    type(matrix) :: c
    type(matrix) :: b_inv
    integer :: size

    size  = a%size
    b_inv = matrix(size)
    c     = matrix(size)
    call inverse(size, b%mat, b_inv%mat)
    c =  mm_(a, b_inv)
  end function mdiv_
  
  subroutine fini(this)
    implicit none
    type(matrix),intent(inout) :: this

    this%size = 0
    if (allocated(this%mat)) deallocate(this%mat)
  end subroutine fini
end module matrix_calc

program main
  use matrix_calc
  implicit none
  type(matrix) :: a, b, c, d, e, f
  integer :: size
  character(len=16) :: argv1

  if (iargc() == 0) then
     size = 4
  else
     call getarg(1, argv1)
     read(argv1, *) size
  end if
  write(6, *) "size:", size

  a = matrix(size)
  b = matrix(size)
  call a%set_mat(-1.0d0, 1.0d0, 5555)
  call b%set_mat(-2.0d0, 2.0d0, 7777)

  c = a*b
  d = c+a-b
  e = d
  f = c/b

  write(6, *) "A:"
  call a%show_mat
  write(6, *) "B:"
  call b%show_mat
  write(6, *) "C = A*B:"
  call c%show_mat
  write(6, *) "D = C+A-B:"
  call d%show_mat
  write(6, *) "E = D"
  call e%show_mat
  write(6, *) "F = C/B (=A)"
  call f%show_mat
  
  stop
end program main

