program main_standardmap
  use constants,    only: pi, twopi
  use standard_map, only: step_standard_map
  use diagnostics,  only: lyapunov_standard_map, diffusion_msd
  implicit none

  character(len=*), parameter :: outdir = "out_standardmap"

  ! ---- K scan: explicit list ----
  integer, parameter :: nK = 33
  real(8), dimension(nK), parameter :: KLIST = (/ &
    0.00d0, 0.10d0, 0.20d0, 0.30d0, 0.40d0, 0.50d0, 0.60d0, 0.70d0, 0.80d0, 0.85d0, 0.90d0, &
    0.93d0, 0.95d0, 0.96d0, 0.97d0, 0.971d0, 0.972d0, 0.98d0, 0.99d0, 1.00d0, 1.02d0, 1.05d0, &
    1.10d0, 1.20d0, 1.40d0, 1.60d0, 1.80d0, 2.00d0, 2.50d0, 3.00d0, 4.00d0, 5.00d0, 5.50d0 /)

  ! ---- Poincare (torus) ----
  integer, parameter :: n_burn  = 100
  integer, parameter :: n_plot  = 200
  integer, parameter :: N_theta = 40
  integer, parameter :: N_p     = 40

  ! ---- Lyapunov vs K (ensemble) ----
  integer, parameter :: n_ens_lyap = 40
  integer, parameter :: n_lyap = 20000

  ! ---- Diffusion (cylinder) ----
  integer, parameter :: n_diff_steps = 5000
  integer, parameter :: n_ens        = 200

  ! ---- FTLE grid (torus) ----
  integer, parameter :: N_theta_ftle = 200
  integer, parameter :: N_p_ftle     = 200
  integer, parameter :: n_ftle_steps = 800

  integer :: j, m, n, i, e, u
  real(8) :: K, theta, p, lambda, theta0, p0
  real(8) :: lam_sum, lam2_sum, lam_mean, lam_std, r
  real(8), allocatable :: msd(:)

  character(len=256) :: path
  character(len=64)  :: kstr

  call random_seed()

  ! Create output directory (re-runnable)
  call execute_command_line("mkdir -p " // trim(outdir))

  ! -----------------------------
  ! 1) Lyapunov vs K (ensemble mean + std)
  ! -----------------------------
  path = trim(outdir)//"/lyapunov_vs_K.csv"
  open(newunit=u, file=trim(path), status="replace", action="write")
  write(u,'(A)') "K,lambda_mean,lambda_std"

  do j = 1, nK
    K = KLIST(j)
    call k_to_str(K, kstr)

    lam_sum  = 0.0d0
    lam2_sum = 0.0d0

    do e = 1, n_ens_lyap
      call random_number(r); theta0 = twopi*r
      call random_number(r); p0     = twopi*r
      call lyapunov_standard_map(theta0=theta0, p0=p0, K=K, n_steps=n_lyap, lambda=lambda)
      lam_sum  = lam_sum  + lambda
      lam2_sum = lam2_sum + lambda*lambda
    enddo

    lam_mean = lam_sum / real(n_ens_lyap,8)
    lam_std  = sqrt(max(0.0d0, lam2_sum/real(n_ens_lyap,8) - lam_mean*lam_mean))

    write(u,'(A,A,ES16.8,A,ES16.8)') trim(kstr), ",", lam_mean, ",", lam_std
  enddo
  close(u)

  ! -----------------------------
  ! 2) Poincaré sections (torus)
  ! -----------------------------
  do j = 1, nK
    K = KLIST(j)
    call k_to_str(K, kstr)

    path = trim(outdir)//"/poincare_K_"//trim(kstr)//".dat"
    open(newunit=u, file=trim(path), status="replace", action="write")
    write(u,'(A)') "# theta  p   (both wrapped to [0,2pi) )"

    do m = 0, N_theta-1
      theta = twopi*real(m,8)/real(N_theta,8)
      do n = 0, N_p-1
        p = twopi*real(n,8)/real(N_p,8)

        do i = 1, n_burn
          call step_standard_map(theta, p, K, periodic_p=.true.)
        enddo

        do i = 1, n_plot
          call step_standard_map(theta, p, K, periodic_p=.true.)
          write(u,'(F18.10,1X,F18.10)') theta, p
        enddo

        write(u,*)
      enddo
    enddo
    close(u)
  enddo

  ! -----------------------------
  ! 3) Diffusion MSD (cylinder)
  ! -----------------------------
  allocate(msd(0:n_diff_steps))

  do j = 1, nK
    K = KLIST(j)
    call k_to_str(K, kstr)

    call diffusion_msd(K=K, n_steps=n_diff_steps, n_ens=n_ens, msd=msd)

    path = trim(outdir)//"/msd_K_"//trim(kstr)//".csv"
    open(newunit=u, file=trim(path), status="replace", action="write")
    write(u,'(A)') "n,msd"
    do i = 0, n_diff_steps
      write(u,'(I0,A,ES18.10)') i, ",", msd(i)
    enddo
    close(u)
  enddo

  deallocate(msd)

  ! -----------------------------
  ! 4) FTLE grid (Lyapunov color map)
  ! -----------------------------
  do j = 1, nK
    K = KLIST(j)
    call k_to_str(K, kstr)

    path = trim(outdir)//"/ftle_grid_K_"//trim(kstr)//".csv"
    open(newunit=u, file=trim(path), status="replace", action="write")
    write(u,'(A)') "theta0,p0,ftle"

    do m = 0, N_theta_ftle-1
      theta = twopi*real(m,8)/real(N_theta_ftle,8)
      do n = 0, N_p_ftle-1
        p = twopi*real(n,8)/real(N_p_ftle,8)
        call lyapunov_standard_map(theta0=theta, p0=p, K=K, n_steps=n_ftle_steps, lambda=lambda)
        write(u,'(F18.10,A,F18.10,A,ES18.10)') theta, ",", p, ",", lambda
      enddo
    enddo
    close(u)
  enddo

contains

  subroutine k_to_str(K, s)
    real(8), intent(in) :: K
    character(len=*), intent(out) :: s
    character(len=32) :: tmp
    integer :: i, nout

    ! Fixed width to keep leading zero; then strip spaces
    write(tmp,'(F6.3)') K

    nout = 0
    s = ""
    do i = 1, len_trim(tmp)
      if (tmp(i:i) /= " ") then
        nout = nout + 1
        if (nout <= len(s)) s(nout:nout) = tmp(i:i)
      end if
    enddo

    ! Safety: convert "-0.200" to "m0.200" for filenames if needed
    if (len_trim(s) >= 1) then
      if (s(1:1) == "-") s(1:1) = "m"
    end if
  end subroutine k_to_str

end program main_standardmap
