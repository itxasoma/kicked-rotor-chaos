module diagnostics
  use constants,    only: pi, twopi, wrap_2pi
  use standard_map, only: step_standard_map
  implicit none
contains

  subroutine lyapunov_standard_map(theta0, p0, K, n_steps, lambda)
    real(8), intent(in)  :: theta0, p0, K
    integer, intent(in)  :: n_steps
    real(8), intent(out) :: lambda

    real(8) :: theta, p
    real(8) :: u, v, u_new, v_new, normv, sumlog
    integer :: i
    logical, parameter :: periodic_p = .true.

    theta = wrap_2pi(theta0)
    p     = wrap_2pi(p0)

    ! Tangent vector
    u = 1d0
    v = 0d0
    sumlog = 0d0

    do i = 1, n_steps
      ! Tangent map for standard map
      v_new = v + K*cos(theta)*u
      u_new = u + v_new

      u = u_new
      v = v_new

      normv = sqrt(u*u + v*v)
      sumlog = sumlog + log(normv)

      u = u/normv
      v = v/normv

      call step_standard_map(theta, p, K, periodic_p)
    enddo

    lambda = sumlog / real(n_steps,8)
  end subroutine lyapunov_standard_map


  subroutine diffusion_msd(K, n_steps, n_ens, msd)
    real(8), intent(in)  :: K
    integer, intent(in)  :: n_steps, n_ens
    real(8), intent(out) :: msd(0:n_steps)

    real(8) :: theta, p, p0, r
    integer :: i, e
    logical, parameter :: periodic_p = .false.  ! cylinder: p unwrapped

    msd = 0d0

    do e = 1, n_ens
      call random_number(r); theta = twopi*r
      call random_number(r); p     = (r - 0.5d0)*twopi
      p0 = p

      msd(0) = msd(0) + (p - p0)**2
      do i = 1, n_steps
        call step_standard_map(theta, p, K, periodic_p)
        msd(i) = msd(i) + (p - p0)**2
      enddo
    enddo

    msd = msd / real(n_ens,8)
  end subroutine diffusion_msd

end module diagnostics
