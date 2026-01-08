module standard_map
  use constants, only: wrap_2pi
  implicit none
contains
  subroutine step_standard_map(theta, p, K, periodic_p)
    real(8), intent(inout) :: theta, p
    real(8), intent(in)    :: K
    logical, intent(in)    :: periodic_p

    ! Standard map:
    ! p_{n+1} = p_n + K*sin(theta_n)
    ! theta_{n+1} = theta_n + p_{n+1}
    p     = p + K*sin(theta)
    theta = theta + p

    theta = wrap_2pi(theta)
    if (periodic_p) p = wrap_2pi(p)
  end subroutine step_standard_map
end module standard_map
