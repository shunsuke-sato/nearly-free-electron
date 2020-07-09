! Cvec: cartesian 
! Lvec: Lattice
! Rvec: Reciprocal lattice

module global_variables
  implicit none

! Mathematical  parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! Physical constants
  real(8),parameter :: AA2au=1d0/0.529177d0

! Grids
  integer :: NK1,NK2,NK3,NK
  real(8) :: a_Cvec(3,3),b_Cvec(3,3), Vcell

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call init_lattice

end program main
!-------------------------------------------------------------------------------
subroutine init_lattice
  use global_variables
  implicit none
  real(8) :: lattice_const
  real(8) :: tvec(3)

  lattice_const = 4.05d0*AA2au
  a_Cvec = 0d0
  a_Cvec(1:3,1) = (/ 0.0d0, 0.5d0, 0.5d0 /)
  a_Cvec(1:3,2) = (/ 0.5d0, 0.0d0, 0.5d0 /)
  a_Cvec(1:3,3) = (/ 0.5d0, 0.5d0, 0.0d0 /)
  a_Cvec = lattice_const*a_Cvec

  tvec(1)=a_Cvec(2,2)*a_Cvec(3,3)-a_Cvec(3,2)*a_Cvec(2,3)
  tvec(2)=a_Cvec(3,2)*a_Cvec(1,3)-a_Cvec(1,2)*a_Cvec(3,3)
  tvec(3)=a_Cvec(1,2)*a_Cvec(2,3)-a_Cvec(2,2)*a_Cvec(1,3)

  Vcell = abs(sum(a_Cvec(:,1)*tvec(:)))

  b_Cvec(1,1)=2d0*pi*(a_Cvec(2,2)*a_Cvec(3,3)-a_Cvec(3,2)*a_Cvec(2,3))
  b_Cvec(2,1)=2d0*pi*(a_Cvec(3,2)*a_Cvec(1,3)-a_Cvec(1,2)*a_Cvec(3,3))
  b_Cvec(3,1)=2d0*pi*(a_Cvec(1,2)*a_Cvec(2,3)-a_Cvec(2,2)*a_Cvec(1,3))

  b_Cvec(1,2)=2d0*pi*(a_Cvec(2,3)*a_Cvec(3,1)-a_Cvec(3,3)*a_Cvec(2,1))
  b_Cvec(2,2)=2d0*pi*(a_Cvec(3,3)*a_Cvec(1,1)-a_Cvec(1,3)*a_Cvec(3,1))
  b_Cvec(3,2)=2d0*pi*(a_Cvec(1,3)*a_Cvec(2,1)-a_Cvec(2,3)*a_Cvec(1,1))

  b_Cvec(1,3)=2d0*pi*(a_Cvec(2,1)*a_Cvec(3,2)-a_Cvec(3,1)*a_Cvec(2,2))
  b_Cvec(2,3)=2d0*pi*(a_Cvec(3,1)*a_Cvec(1,2)-a_Cvec(1,1)*a_Cvec(3,2))
  b_Cvec(3,3)=2d0*pi*(a_Cvec(1,1)*a_Cvec(2,2)-a_Cvec(2,1)*a_Cvec(1,2))

  b_Cvec=b_Cvec/sum(a_Cvec(:,1)*tvec(:))

!sanity check
  write(*,"(A)")"inter products of a_Cvec and b_Cvec (devided by 2*pi)"
  write(*,"(999e16.6e3)")sum(a_Cvec(:,1)*b_Cvec(:,1))/(2d0*pi) &
                        ,sum(a_Cvec(:,1)*b_Cvec(:,2))/(2d0*pi) &
                        ,sum(a_Cvec(:,1)*b_Cvec(:,3))/(2d0*pi)
  write(*,"(999e16.6e3)")sum(a_Cvec(:,2)*b_Cvec(:,1))/(2d0*pi) &
                        ,sum(a_Cvec(:,2)*b_Cvec(:,2))/(2d0*pi) &
                        ,sum(a_Cvec(:,2)*b_Cvec(:,3))/(2d0*pi)
  write(*,"(999e16.6e3)")sum(a_Cvec(:,3)*b_Cvec(:,1))/(2d0*pi) &
                        ,sum(a_Cvec(:,3)*b_Cvec(:,2))/(2d0*pi) &
                        ,sum(a_Cvec(:,3)*b_Cvec(:,3))/(2d0*pi)


end subroutine init_lattice
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
