!+ Monte Carlo simulation of a spaceborne lidar signal from the ocean
!
Program MC_water

! Description:
!   Forward Monte Carlo simulation of the return of a spaceborne lidar
!   from a homogeneous water column under a rough sea surface (532 nm
!   reference). Photons are launched from the apex of the receiver cone,
!   cross the Cox-Munk surface, random-walk in the water and contribute
!   to the signal through local estimates. The output gives the signal
!   per incident photon versus real depth and versus time-based depth,
!   by scattering order, plus the surface bin (glint) and
!   the number of collisions inside and outside the field of view.
!
! Method:
!   1. Input, seed and arrays: read_input, set_random_seed,
!      write_input_used (mc_module).
!   2. Geometry computed once: footprint ellipse (Shangguan et al. 2024,
!      Eq. 11-16), refracted optical axis (Snell), receiver directions.
!   3. Tables computed once: Fresnel LUTs, interface
!      transmittance on the refracted axis with shadowing.
!   4. Photon loop. Each photon: generation from the cone apex
!      (generate_photon); air->water surface interaction; glint local
!      estimate if the entry point is inside the footprint ellipse;
!      random walk with move-then-scatter ordering (free path sampled
!      with c, Beer-Lambert; weight multiplied by omega at each
!      collision; termination at ww <= wmin). At each subsurface
!      collision inside the FOV: local estimate (Shangguan et al. 2024,
!      Eq. 23), binned by real depth (Shangguan et al. 2025, Eq. 25) and
!      by time-based depth z_t = s_p cos(alpha_r) (Eq. 24 and 27). When
!      the photon reaches the surface from below: water->air interaction;
!      the transmitted fraction leaves without being scored (the local
!      estimate of the last collision already holds its contribution),
!      the reflected fraction continues.
!   5. Output: write_output.
!
! Input files:
!   Namelist file (first command-line argument, default input.nml),
!   groups &geometry, &atmosphere, &surface, &water, &numerics,
!   &receiver; read in
!   read_input. Every variable is required. Examples in input/.
!
! Output files:
!   SIGNAL_PROFILE_Z.DAT, SIGNAL_PROFILE_S.DAT, SCATTERING_COUNTS.DAT
!   (written in write_output) and input_used.nml (written in
!   write_input_used), in the current directory.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-09-22  Original code (time-based profile patch of
!                      2026-09-21). Valerio Ciattaglia
! v0.1     2026-11-05  Time-based bin on z_t = s_p cos(alpha_r);
!                      configurable output grid (zmax, dz).
!                      Valerio Ciattaglia
! v0.2     2026-11-05  Parameters from a namelist file; collision
!                      counters; phase-function arrays (0:nangle).
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Input and output moved to mc_module; omega only
!                      at collisions; s_tot cut at the exit point for
!                      reflected photons; FOV shift at the photon depth;
!                      analytic HG sampling. Valerio Ciattaglia
! v0.4     2026-11-07  Beam divergence beam_div from the input; header
!                      and layout to the European Standards.
!                      Valerio Ciattaglia
! v0.5     2026-11-08  Analog upwelling score removed: it counted again
!                      the event of the local estimate (B6).
!                      Valerio Ciattaglia
! v0.6     2026-11-08  Sea-surface slopes by analytic inversion; slope
!                      tables removed. Valerio Ciattaglia
! v0.7     2026-11-11  Glint call without n_water (B5a).
!                      Valerio Ciattaglia
! v0.7     2026-11-11  Glint toward the telescope centre x_sat instead
!                      of x_A (B5c). Valerio Ciattaglia
! v0.8     2026-11-11  Output grid from the receiver sampling interval
!                      (&receiver); binning unchanged. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran), C preprocessor for the
!                       version macro GIT_VERSION (set by the Makefile).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:

Use mc_module, Only :                                                  &
! Imported Scalar Variables with intent (in):
    theta_v, fov, beam_div, hlidar, r_tel, sampling_mode, tau0,        &
    n_water, internal_reflection, kext, omega, g, nphotons, wmin,      &
    zmax, dz, N_LUT, maxsca, D_l, sigma_g, sigma2, nz,         &
! Imported Routines:
    read_input, set_random_seed, write_input_used, write_output,       &
    compute_ellipse_params, generate_photon, photon_in_fov_simplified, &
    depth_bin, generate_fresnel_lut, get_lut_nearest,                  &
    lambda_saunders, surface_interaction, surface_interaction_wa,      &
    find_surface_crossing, surface_glint_signal_le,                    &
    local_estimate_signal, sample_hg_cos, scatter_direction

Implicit None

#ifndef GIT_VERSION
#define GIT_VERSION 'unknown'
#endif

! Local parameters:
Character(len=*), Parameter :: version = GIT_VERSION
                                    ! Code version (git describe).

! Local scalars:
!   Photon state.
Real(8) :: x                        ! Photon x [m].
Real(8) :: y                        ! Photon y [m].
Real(8) :: z                        ! Photon z [m], z < 0 in water (z up).
Real(8) :: s                        ! Free path of the current step [m].
Real(8) :: s_tot                    ! In-water path since entry [m].
Real(8) :: ww                       ! Photon weight.
Integer :: ns                       ! Scattering events so far.
Real(8) :: x_prev                   ! x before the step [m].
Real(8) :: y_prev                   ! y before the step [m].
Real(8) :: z_prev                   ! z before the step [m].
Real(8) :: x_E                      ! Surface crossing x [m].
Real(8) :: y_E                      ! Surface crossing y [m].
!   Scattering.
Real(8) :: cos_theta                ! Cosine of the polar angle.
Real(8) :: phi                      ! Azimuthal angle [rad].
Real(8) :: ranval                   ! Uniform random number.
!   Geometry.
Real(8) :: pi                       ! Pi.
Real(8) :: alpha                    ! Off-nadir angle [rad].
Real(8) :: Area                     ! Receiver area pi r_tel^2 [m^2].
Real(8) :: H_prime                  ! Altitude of the cone apex [m].
Real(8) :: m_axis                   ! Footprint semi-major axis [m].
Real(8) :: n_axis                   ! Footprint semi-minor axis [m].
Real(8) :: x_A                      ! Foot of the cone apex [m].
Real(8) :: x_C                      ! Aim point [m].
Real(8) :: x_sat                    ! Telescope centre x, at hlidar [m].
Real(8) :: geom_shift               ! Centre-to-aim-point shift [m].
Real(8) :: sin_alpha_r_axis         ! sin of the refracted axis angle.
Real(8) :: alpha_r_axis             ! Refracted axis angle [rad].
Real(8) :: tan_alpha_r_axis         ! tan of the refracted axis angle.
Real(8) :: cos_alpha_r_axis         ! cos of the refracted axis angle.
!   Sea surface.
Real(8) :: refl                     ! Reflectance at the sampled facet.
Real(8) :: trans                    ! Transmittance at the sampled facet.
Real(8) :: T_t                      ! T_wa on the refracted axis.
Real(8) :: Lambda_i_axis            ! Shadowing on the axis, water side.
Real(8) :: Lambda_t_axis            ! Shadowing on the axis, air side.
!   Signal and binning.
Real(8) :: signal                   ! Contribution of one event.
Real(8) :: signal_glint             ! Glint contribution of one photon.
Real(8) :: d_ret                    ! Return path in water d' [m].
Real(8) :: s_p                      ! Time-based position, Eq. 24 [m].
Real(8) :: z_t                      ! Time-based depth, Eq. 27 [m].
Real(8) :: s_cross                  ! In-water path to the exit [m].
Integer :: iphot                    ! Photon index.
Integer :: iz                       ! Real-depth bin.
Integer :: is                       ! Time-based bin.
Integer :: alloc_stat               ! Allocation status.

! Local arrays:
Real(8) :: u(3)                     ! Photon direction versor.
Real(8) :: k_ref(3)                 ! Reflected versor.
Real(8) :: k_tr(3)                  ! Refracted versor.
Real(8) :: u_wat_recv(3)            ! Toward the receiver, in water.
Real(8), Allocatable :: LUT_cos(:)       ! cos(incidence) grid.
Real(8), Allocatable :: LUT_R_aa(:)      ! Reflectance, air side.
Real(8), Allocatable :: LUT_T_aw(:)      ! Transmittance, air->water.
Real(8), Allocatable :: LUT_R_ww(:)      ! Reflectance, water side.
Real(8), Allocatable :: LUT_T_wa(:)      ! Transmittance, water->air.
Real(8), Allocatable :: photons(:,:)     ! Signal (order, real depth);
                                         ! order 0 = total.
Real(8), Allocatable :: photons_s(:,:)   ! Signal (order, z_t).
Integer(8), Allocatable :: n_in(:,:)     ! Collisions inside the FOV.
Integer(8), Allocatable :: n_out(:,:)    ! Collisions outside the FOV.

!- End of header ------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Input, random seed and arrays
!-----------------------------------------------------------------------
Call read_input(version)

Allocate(photons(0:maxsca,0:nz), photons_s(0:maxsca,0:nz),             &
         n_in(0:maxsca,0:nz), n_out(0:maxsca,0:nz), stat=alloc_stat)
If (alloc_stat /= 0) Then
  Write(*,'(A)') 'ERROR: allocation of the output arrays failed.'
  Stop 1
End If
Allocate(LUT_cos(N_LUT), LUT_R_aa(N_LUT), LUT_T_aw(N_LUT),             &
         LUT_R_ww(N_LUT), LUT_T_wa(N_LUT), stat=alloc_stat)
If (alloc_stat /= 0) Then
  Write(*,'(A)') 'ERROR: allocation of the tables failed.'
  Stop 1
End If

Call set_random_seed()
Call write_input_used()

!-----------------------------------------------------------------------
! 2. Geometry (computed once)
!-----------------------------------------------------------------------
pi    = acos(-1.0d0)
alpha = theta_v * pi / 180.0d0   ! Off-nadir angle [rad].
Area  = pi*r_tel*r_tel           ! Effective telescope aperture area.

! Lidar footprint ellipse, Shangguan et al. (2024), Eq. 11, 12, 14, 16.
Call compute_ellipse_params(hlidar, alpha, fov, D_l,                   &
                            H_prime, m_axis, n_axis, x_A, x_C,         &
                            geom_shift)

! Telescope centre: point of the optical axis at altitude hlidar, i.e.
! (H_prime - hlidar) below the virtual apex F, whose foot is x_A.
! Equivalent to x_C - hlidar tan(alpha).
x_sat = x_A + (H_prime - hlidar) * tan(alpha)

! FOV: refraction angle of the optical axis (Snell, computed once).
! The underwater FOV ellipse at depth |z| is the surface ellipse
! translated by |z| tan(alpha_r_axis) (photon_in_fov_simplified).
sin_alpha_r_axis = sin(alpha) / n_water
alpha_r_axis     = asin(sin_alpha_r_axis)
tan_alpha_r_axis = tan(alpha_r_axis)
cos_alpha_r_axis = cos(alpha_r_axis)

! Receiver directions (constant for all photons), z positive upward:
! u_wat_recv: receiver direction in water (refracted optical axis,
!             flat-surface approximation), from the scatter point
!             toward the receiver;
u_wat_recv(1) = -sin(alpha_r_axis)
u_wat_recv(2) =  0.0d0
u_wat_recv(3) =  cos(alpha_r_axis)

!-----------------------------------------------------------------------
! 3. Tables (computed once)
!-----------------------------------------------------------------------
! Surface optical properties: Fresnel LUTs and Cox-Munk slopes.
Call generate_fresnel_lut(N_LUT, n_water, LUT_cos, LUT_R_aa, LUT_T_aw, &
                          LUT_R_ww, LUT_T_wa)

! Local estimate: Fresnel T_wa at the refracted optical-axis angle,
! corrected by bidirectional shadowing, evaluated at the exit geometry
! (theta_i = alpha_r_axis in water, theta_t = alpha in air).
Call get_lut_nearest(cos(alpha_r_axis), LUT_cos, LUT_T_wa, T_t)
Lambda_i_axis = lambda_saunders(cos(alpha_r_axis), sigma2)
Lambda_t_axis = lambda_saunders(cos(alpha),        sigma2)
T_t = T_t / (1.0d0 + Lambda_i_axis + Lambda_t_axis)

photons_s(0:maxsca,0:nz) = 0.0d0
photons(0:maxsca,0:nz)   = 0.0d0
n_in(0:maxsca,0:nz)      = 0_8
n_out(0:maxsca,0:nz)     = 0_8

!-----------------------------------------------------------------------
! 4. Photon loop
!-----------------------------------------------------------------------
Do iphot = 1, nphotons

  ! Progress report every 10%.
  If (mod(iphot, max(1, nphotons/10)) == 0) Then
    Write(*,'(A,I3,A,I10,A,I10,A)')                                    &
      '  Progress: ', int(int(iphot, 8)*100_8/int(nphotons, 8)),       &
      '%  (photon ', iphot, ' / ', nphotons, ')'
  End If

  !---------------------------------------------------------------------
  ! 4.1 Photon generation
  !---------------------------------------------------------------------
  ! Sample entry point (x, y) and incident direction u.
  Call generate_photon(sigma_g, beam_div, sampling_mode, alpha,        &
                       H_prime, x_A, x, y, u)
  z = 0.0d0                             ! Photon enters at the surface.
  If (x <= -huge(1.0d0)*0.5d0) Cycle    ! Grazing ray, discarded.

  !---------------------------------------------------------------------
  ! 4.2 Surface interaction (air -> water) and glint
  !---------------------------------------------------------------------
  Call surface_interaction(u, N_LUT,                                   &
                           LUT_cos, LUT_R_aa, LUT_T_aw,                &
                           sigma2, n_water,                            &
                           refl, trans, k_ref, k_tr)

  ! Glint local estimate, if the surface hit point is visible to the
  ! receiver; u is still the incident direction from generate_photon.
  If ((x / m_axis)**2 + (y / n_axis)**2 <= 1.0d0) Then
    signal_glint = surface_glint_signal_le(u, x, y, x_sat, hlidar,     &
                     sigma2, N_LUT, LUT_cos, LUT_R_aa,                 &
                     tau0, alpha, Area)
    photons(0, 0)   = photons(0, 0)   + signal_glint   ! Real, iz = 0.
    photons_s(0, 0) = photons_s(0, 0) + signal_glint   ! Time, is = 0.
  End If

  ! Direction and weight after the surface interaction.
  u  = k_tr
  ww = 1.0d0 * trans
  ns = 0                                ! Scattering event counter.
  s_tot = 0.0d0                         ! Cumulative path length.

  !---------------------------------------------------------------------
  ! 4.3 Random walk
  !---------------------------------------------------------------------
  ! Move-then-scatter ordering: the free path is travelled along the
  ! CURRENT direction (k_tr on the first step), and the direction is
  ! updated only AT the collision, for the next flight.
  Do While (ww > wmin)

    ! Free-path sampling (Beer-Lambert).
    Call random_number(ranval)
    If (ranval < tiny(1.0d0)) ranval = tiny(1.0d0)   ! Avoid log(0).
    s = -log(ranval)/kext
    s_tot = s_tot + s

    ! Ballistic propagation along u; state saved for the surface
    ! crossing.
    x_prev = x
    y_prev = y
    z_prev = z
    x = x + u(1)*s
    y = y + u(2)*s
    z = z + u(3)*s

    If (z >= 0.0d0) Then
      !-----------------------------------------------------------------
      ! 4.3.1 Surface crossing from below
      !-----------------------------------------------------------------
      ! No collision: the step reaches the surface first, so no omega.
      ! Exact surface exit point.
      Call find_surface_crossing(x_prev, y_prev, z_prev,               &
                                 x, y, z, x_E, y_E)
      ! s_tot includes the full last step s, which overshoots z = 0:
      ! replace it with the segment from the previous position to the
      ! exit point E.
      s_cross = s_tot - s + sqrt((x_E - x_prev)**2                     &
                                 + (y_E - y_prev)**2                   &
                                 + z_prev**2)
      ! Reflected/transmitted split (including TIR), local Fresnel
      ! coefficients and bidirectional shadowing; u is the direction
      ! that carried the photon across the interface.
      Call surface_interaction_wa(u, N_LUT,                            &
                                  LUT_cos, LUT_R_ww, LUT_T_wa,         &
                                  sigma2, n_water,                     &
                                  refl, trans, k_ref, k_tr)

      ! The transmitted fraction leaves the system: its contribution to
      ! the receiver is already in the local estimate of the last
      ! collision.

      ! Reflected fraction: continues the random walk from E.
      If (internal_reflection) Then
        ww = ww * refl
      Else
        ww = 0.0d0                      ! The photon leaves the water.
      End If
      If (ww <= wmin) ww = 0.0d0        ! Terminate.
      u  = k_ref
      x  = x_E
      y  = y_E
      z  = -1.0d-6                      ! Epsilon below the surface [m].
      s_tot = s_cross

    Else
      !-----------------------------------------------------------------
      ! 4.3.2 Subsurface collision and local estimate
      !-----------------------------------------------------------------
      ! Collision at (x, y, z); u is the incoming (pre-scatter)
      ! direction.
      ns = ns + 1                       ! Scattering order counter.
      ww = ww*omega                     ! Absorption at the collision.
      iz = depth_bin(abs(z), dz, zmax, nz)   ! Real-depth bin, nz+1 if
                                             ! beyond zmax.
      If (iz <= nz) Then
        If (photon_in_fov_simplified(x, y, abs(z), m_axis, n_axis,     &
                                     tan_alpha_r_axis)) Then
          ! Local estimate signal (Shangguan et al. 2024, Eq. 23).
          signal = local_estimate_signal(ww, u, u_wat_recv,            &
                     abs(z), hlidar, alpha, alpha_r_axis, n_water,     &
                     kext, Area, T_t, tau0, g)
          ! Time-based bin (Shangguan et al. 2025, Eq. 24): forward path
          ! s_tot plus the return leg to the surface along the receiver
          ! axis; projected on the vertical (Eq. 27). iz holds the real
          ! depth z_p (Eq. 25).
          d_ret = abs(z) / cos_alpha_r_axis
          s_p   = 0.5d0 * (s_tot + d_ret)
          z_t   = s_p * cos_alpha_r_axis
          is    = depth_bin(z_t, dz, zmax, nz)
          n_in(0, iz) = n_in(0, iz) + 1_8
          If (ns <= maxsca) n_in(ns, iz) = n_in(ns, iz) + 1_8
          photons(0, iz) = photons(0, iz) + signal
          If (ns <= maxsca) photons(ns, iz) = photons(ns, iz) + signal
          If (is <= nz) Then
            photons_s(0, is) = photons_s(0, is) + signal
            If (ns <= maxsca) photons_s(ns, is) = photons_s(ns, is)    &
                                                  + signal
          End If
        Else
          ! Out-of-FOV collision: not scored, the photon continues.
          n_out(0, iz) = n_out(0, iz) + 1_8
          If (ns <= maxsca) n_out(ns, iz) = n_out(ns, iz) + 1_8
        End If
      End If

      ! Scattering: direction update for the NEXT free path.
      Call sample_hg_cos(g, cos_theta)
      Call random_number(ranval)
      phi = 2.0d0*pi*ranval             ! Azimuth angle.
      Call scatter_direction(u, cos_theta, phi)
    End If
  End Do

End Do

!-----------------------------------------------------------------------
! 5. Output
!-----------------------------------------------------------------------
Call write_output(photons, photons_s, n_in, n_out)

Deallocate(photons, photons_s, n_in, n_out)
Deallocate(LUT_cos, LUT_R_aa, LUT_T_aw, LUT_R_ww, LUT_T_wa)

End Program MC_water
