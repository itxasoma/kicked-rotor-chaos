module constants
  implicit none
  real(8), parameter :: pi    = 4d0*atan(1d0)
  real(8), parameter :: twopi = 2d0*pi
contains
  pure real(8) function wrap_2pi(x) result(y)
    real(8), intent(in) :: x
    y = modulo(x, twopi)
  end function wrap_2pi
end module constants
