!+ Input parameters and procedures of the Monte Carlo lidar code MC_water
!
Module mc_module

!
! Description:
!   Holds the input parameters of MC_water (read from a namelist file),
!   the quantities derived from them, and every procedure of the code:
!   phase-function sampling, beam and footprint geometry, rough sea
!   surface (Fresnel, Cox-Munk slopes, Smith shadowing), signal estimators
!   (glint, local estimate) and input/output.
!   The main program MC_water only declares the photon state, prepares
!   the tables and runs the photon loop.
!
!   Procedures, by section:
!     1. Phase function and scattering : sample_hg_cos, scatter_direction
!     2. Beam and geometry             : compute_ellipse_params,
!                                        generate_photon,
!                                        photon_in_fov_simplified, depth_bin
!     3. Sea surface                   : generate_fresnel_lut, DCMP,
!                                        Fresnel_Trans_Scalar,
!                                        Fresnel_Refl_Scalar, Compute_CosT,
!                                        get_lut_nearest, erfinv_win,
!                                        sample_slope, lambda_saunders,
!                                        surface_interaction,
!                                        surface_interaction_wa,
!                                        find_surface_crossing
!     4. Signal estimators             : surface_glint_signal_le,
!                                        local_estimate_signal
!     5. Input and output              : read_input, check_input_complete,
!                                        compact_lower, is_name_char,
!                                        check_group, require,
!                                        set_random_seed, write_input_used,
!                                        write_output, write_header,
!                                        order_columns
!
!   Main references:
!     Shangguan, Guo, Liao and Lee (2024), Remote Sensing of Environment:
!       footprint geometry and local estimate (Eq. 11-23).
!     Shangguan et al. (2025), Remote Sensing of Environment: real and
!       time-based depth profiles (Eq. 24-27).
!     Nakajima and Tanaka (1983), J. Quant. Spectrosc. Radiat. Transfer
!       29(6): rough-surface slopes and shadowing (Eq. 14-20).
!     Cox and Munk (1954): mean-square slope of the sea surface.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
!
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.1     2026-11-05  depth_bin added (configurable output grid).
!                      Valerio Ciattaglia
! v0.2     2026-11-05  check_input_complete added (namelist input).
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Input parameters as module variables; input and
!                      output procedures moved here from MC.f90; unused
!                      procedures removed; analytic HG sampling.
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Laser divergence beam_div as input; headers and
!                      layout to the European Standards.
!                      Valerio Ciattaglia
! v0.5     2026-11-08  upwelling_surface_signal removed (B6).
!                      Valerio Ciattaglia
! v0.6     2026-11-08  Sea-surface slopes by analytic inversion
!                      (erfinv_win, sample_slope); slope tables removed.
!                      Valerio Ciattaglia
! v0.8     2026-11-11  Output grid from the receiver sampling interval:
!                      group &receiver (dt_bin, n_bins); zmax and dz
!                      derived; c_light. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!

Implicit None

! Global (i.e. public) Declarations:

! Global Parameters:
Real(8), Parameter :: c_light = 299792458.0d0   ! Speed of light in
                                                ! vacuum [m/s].

! Global Scalars:
!   Input parameters, namelist group &geometry.
Real(8) :: theta_v                  ! Off-nadir angle in air [deg].
Real(8) :: fov                      ! Full-angle receiver FOV [rad].
Real(8) :: beam_div                 ! Full-angle laser divergence at
                                    ! 1/e^2 of the intensity [rad].
Real(8) :: hlidar                   ! Orbital altitude H [m].
Real(8) :: r_tel                    ! Telescope radius [m].
Character(len=12) :: sampling_mode  ! 'gaussian' or 'uniform_cone'.
!   Input parameters, namelist group &atmosphere.
Real(8) :: tau0                     ! One-way optical depth at nadir [-].
!   Input parameters, namelist group &surface.
Real(8) :: n_water                  ! Refractive index of seawater [-].
Real(8) :: u10                      ! Wind speed at 10 m [m/s].
Logical :: flat_sea                 ! .true.: sigma2 = 0.
Logical :: internal_reflection      ! .true.: Fresnel reflection of
                                    ! upwelling photons into the water.
!   Input parameters, namelist group &water.
Real(8) :: kext                     ! Beam attenuation coefficient c [1/m].
Real(8) :: omega                    ! Single-scattering albedo b/c [-].
Real(8) :: g                        ! HG asymmetry factor [-].
!   Input parameters, namelist group &numerics.
Integer :: nphotons                 ! Number of photons.
Integer :: seed                     ! Random seed (0 = from the clock).
Real(8) :: wmin                     ! Weight cut-off [-].
Integer :: N_LUT                    ! Points of the Fresnel tables.
Integer :: maxsca                   ! Highest order recorded separately.
!   Input parameters, namelist group &receiver.
Real(8) :: dt_bin                   ! Receiver sampling interval [ns].
Integer :: n_bins                   ! Samples of the subsurface window.
!   Quantities derived from the input (set by read_input).
Real(8) :: dz                       ! Thickness of the output bins [m],
                                    ! from dt_bin.
Real(8) :: zmax                     ! Depth of the output grid,
                                    ! n_bins dz [m].
Real(8) :: D_l                      ! Telescope diameter 2 r_tel [m].
Real(8) :: sigma_g                  ! Beam angular std, beam_div/4 [rad].
Real(8) :: ksca                     ! Scattering coefficient b [1/m].
Real(8) :: sigma2                   ! Total mean-square slope [-].
Integer :: nz                       ! Number of subsurface bins, n_bins.
!   Run bookkeeping.
Character(len=512) :: input_file    ! Namelist input file.
Character(len=64)  :: code_version  ! Code version (git describe).

! Namelist groups of the input file.
Namelist /geometry/   theta_v, fov, beam_div, hlidar, r_tel, sampling_mode
Namelist /atmosphere/ tau0
Namelist /surface/    n_water, u10, flat_sea, internal_reflection
Namelist /water/      kext, omega, g
Namelist /numerics/   nphotons, seed, wmin, N_LUT, maxsca
Namelist /receiver/   dt_bin, n_bins

Contains

!-----------------------------------------------------------------------
! 1. Phase function and scattering
!-----------------------------------------------------------------------

!+ Samples cos(theta) from the Henyey-Greenstein phase function
!
Subroutine sample_hg_cos &
  ( g,  & ! in
    mu )  ! out

! Description:
!   Draws the cosine of the polar scattering angle from the
!   Henyey-Greenstein (HG) phase function with asymmetry factor g.
!
! Method:
!   Exact inversion of the HG cumulative distribution in mu = cos(theta).
!   With xi uniform in [0, 1):
!     mu = ( 1 + g^2 - ( (1 - g^2) / (1 - g + 2 g xi) )^2 ) / (2 g)
!   (Henyey and Greenstein 1941; standard Monte Carlo form, e.g. Wang,
!   Jacques and Zheng 1995, MCML; cited from memory, to be verified).
!   Limits: xi = 0 gives mu = -1, xi = 1 gives mu = +1.
!   For |g| < 1e-6 the function is isotropic: mu = 2 xi - 1.
!   The denominator 1 - g + 2 g xi is >= min(1 - g, 1 + g) > 0.
!   mu is clipped to [-1, 1] against rounding.
!   Test (runs/2026-11-07_revisione/test_hg): <mu> = 0.92001 +- 0.00007
!   for g = 0.92 with 1e7 samples.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.3     2026-11-07  Original code (replaces the tabulated phase
!                      function integra_phase). Valerio Ciattaglia
! v0.4     2026-11-07  Function turned into a subroutine (random number
!                      drawn: side effect); header to the standard.
!                      Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)  :: g      ! Asymmetry factor, -1 < g < 1.

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: mu     ! Sampled cos(theta).

! Local scalars:
Real(8) :: xi                  ! Uniform random number.
Real(8) :: frac                ! (1 - g^2)/(1 - g + 2 g xi).

!- End of header ------------------------------------------------------

Call random_number(xi)
If (abs(g) < 1.0d-6) Then
  mu = 2.0d0*xi - 1.0d0
Else
  frac = (1.0d0 - g*g) / (1.0d0 - g + 2.0d0*g*xi)
  mu   = (1.0d0 + g*g - frac*frac) / (2.0d0*g)
End If
mu = max(-1.0d0, min(1.0d0, mu))

End Subroutine sample_hg_cos

!+ Rotation of the photon direction at a scattering event
!
Subroutine scatter_direction &
  ( u,         & ! inout
    cos_theta, & ! in
    phi )        ! in

! Description:
!   Single scattering rotation of a photon direction (no free-flight
!   advance). Rotation part of the former propagate_step, extracted so
!   that the caller can advance the photon BEFORE scattering
!   (move-then-scatter ordering): the free path is travelled along the
!   incoming direction, and the direction is updated only at the
!   collision, for the next flight.
!
! Method:
!   Standard Monte Carlo rotation by the polar angle theta and the
!   azimuth phi around the incoming direction u:
!     u'x = sin(th)(ux uz cos(phi) - uy sin(phi))/sqrt(1-uz^2) + ux cos(th)
!     u'y = sin(th)(uy uz cos(phi) + ux sin(phi))/sqrt(1-uz^2) + uy cos(th)
!     u'z = -sin(th) cos(phi) sqrt(1-uz^2) + uz cos(th)
!   For |uz| > 0.9999 the direction is taken as vertical:
!     u' = (sin(th) cos(phi), sin(th) sin(phi), sign(uz) cos(th)).
!   The |u(3)| > 0.9999 guard keeps the argument of sqrt(1 - u(3)^2)
!   strictly positive and avoids division by zero near the vertical.
!   sin(theta) = sqrt(max(0, 1 - cos^2(theta))).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.3     2026-11-07  Takes cos(theta) instead of theta (analytic HG).
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)    :: cos_theta   ! Cosine of the polar angle.
Real(8), Intent(In)    :: phi         ! Azimuthal angle [rad].

! Array  arguments with intent(inout):
Real(8), Intent(InOut) :: u(3)        ! Direction cosines (unit vector).

! Local scalars:
Real(8) :: sin_theta                  ! Sine of the polar angle.
Real(8) :: uxx                        ! New x component.
Real(8) :: uyy                        ! New y component.
Real(8) :: uzz                        ! New z component.
Real(8) :: b                          ! 1/sqrt(1 - u(3)^2).

!- End of header ------------------------------------------------------

sin_theta = sqrt(max(0.0d0, 1.0d0 - cos_theta*cos_theta))
If (abs(u(3)) > 0.9999d0) Then
  uxx = sin_theta * cos(phi)
  uyy = sin_theta * sin(phi)
  uzz = sign(1.0d0, u(3)) * cos_theta
Else
  b   = 1.0d0 / sqrt(1.0d0 - u(3) * u(3))
  uxx = b * sin_theta * (u(1) * u(3) * cos(phi) - u(2) * sin(phi))     &
        + u(1) * cos_theta
  uyy = b * sin_theta * (u(2) * u(3) * cos(phi) + u(1) * sin(phi))     &
        + u(2) * cos_theta
  uzz =           -sin_theta * cos(phi) / b                            &
        + u(3) * cos_theta
End If
u(1) = uxx
u(2) = uyy
u(3) = uzz

End Subroutine scatter_direction

!-----------------------------------------------------------------------
! 2. Beam and geometry
!-----------------------------------------------------------------------

!+ Footprint ellipse of the receiver cone on the sea surface
!
Subroutine compute_ellipse_params &
  ( H,          & ! in
    alpha,      & ! in
    FOV,        & ! in
    D_l,        & ! in
    H_prime,    & ! out
    m_axis,     & ! out
    n_axis,     & ! out
    x_A,        & ! out
    x_C,        & ! out
    geom_shift )  ! out

! Description:
!   Geometric parameters of the slant-view lidar footprint ellipse on the
!   sea surface (z = 0).
!
!   Inputs:
!     H       : lidar altitude above the surface [m]
!     alpha   : off-nadir tilt angle [rad]
!     FOV     : full receiver field of view [rad]
!     D_l     : telescope diameter [m]
!   Outputs:
!     H_prime    : altitude of the receiver-cone apex F         - Eq. 12
!     m_axis     : ellipse semi-major axis                      - Eq. 11
!     n_axis     : ellipse semi-minor axis (via point M)        - Eq. 14-16
!     x_A        : horizontal position of foot A in the O frame (< 0)
!     x_C        : aim point (cone axis) in the O frame
!     geom_shift : O->C shift (= x_C, since x_O = 0 by definition)
!
! Method:
!   Geometry of Shangguan, Guo, Liao and Lee (2024), "Discrepancies
!   between time-based and real depth profiles in ocean lidar due to
!   multiple scattering", Remote Sensing of Environment:
!     H'  = H + D_l cos(alpha) / (2 tan(FOV/2))                  (Eq. 12)
!     m   = H' (tan(alpha + FOV/2) - tan(alpha - FOV/2)) / 2     (Eq. 11)
!     y_M = H' tan(FOV/2) / cos(alpha),  x_M = x_C               (Eq. 14)
!     n   = y_M / sqrt(1 - (x_M/m)^2)                            (Eq. 16)
!   Convention: origin O at the ellipse centre on the surface; the cone
!   apex F sits at altitude H' above foot A, offset horizontally by x_A.
!   At nadir (alpha -> 0) the footprint degenerates to a circle (m = n).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)  :: H            ! Lidar altitude [m].
Real(8), Intent(In)  :: alpha        ! Off-nadir angle [rad].
Real(8), Intent(In)  :: FOV          ! Full field of view [rad].
Real(8), Intent(In)  :: D_l          ! Telescope diameter [m].

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: H_prime      ! Altitude of the apex F [m].
Real(8), Intent(Out) :: m_axis       ! Semi-major axis [m].
Real(8), Intent(Out) :: n_axis       ! Semi-minor axis [m].
Real(8), Intent(Out) :: x_A          ! Foot A, O frame [m].
Real(8), Intent(Out) :: x_C          ! Aim point C, O frame [m].
Real(8), Intent(Out) :: geom_shift   ! O->C shift [m].

! Local scalars:
Real(8) :: tan_diff                  ! tan(alpha - FOV/2), vertex B.
Real(8) :: tan_sum                   ! tan(alpha + FOV/2), vertex D.
Real(8) :: half_fov                  ! FOV/2 [rad].
Real(8) :: x_O_inA                   ! Centre O in the A frame [m].
Real(8) :: x_C_inA                   ! Aim point C in the A frame [m].
Real(8) :: x_M                       ! Abscissa of point M [m].
Real(8) :: y_M                       ! Ordinate of point M [m].

!- End of header ------------------------------------------------------

half_fov = 0.5d0 * FOV
tan_diff = tan(alpha - half_fov)
tan_sum  = tan(alpha + half_fov)

!-----------------------------------------------------------------------
! 2.1 Apex altitude (Eq. 12) and semi-major axis (Eq. 11)
!-----------------------------------------------------------------------
H_prime = H + (D_l * cos(alpha)) / (2.0d0 * tan(half_fov))
m_axis  = 0.5d0 * H_prime * (tan_sum - tan_diff)

!-----------------------------------------------------------------------
! 2.2 Positions of O, A and C
!-----------------------------------------------------------------------
! Centre O in the frame with A at the origin.
x_O_inA = 0.5d0 * H_prime * (tan_diff + tan_sum)
! Positions in the frame with O at the origin.
x_A     = -x_O_inA
x_C_inA = H_prime * tan(alpha)
x_C     = x_C_inA - x_O_inA
geom_shift = x_C

!-----------------------------------------------------------------------
! 2.3 Semi-minor axis through point M (Eq. 14-16)
!-----------------------------------------------------------------------
x_M = x_C
y_M = H_prime * tan(half_fov) / cos(alpha)
n_axis = y_M / sqrt(1.0d0 - (x_M / m_axis)**2)

End Subroutine compute_ellipse_params

!+ Launch of one photon from the cone apex to the sea surface
!
Subroutine generate_photon &
  ( sigma_g,  & ! in
    beam_div, & ! in
    mode,     & ! in
    alpha,    & ! in
    H_apex,   & ! in
    x_A,      & ! in
    x_surf,   & ! out
    y_surf,   & ! out
    u_inc )     ! out

! Description:
!   Forward Monte Carlo launch of a single photon: samples the beam
!   angular divergence, propagates from the cone apex, and intersects the
!   sea surface (z = 0) deterministically.
!
!   Inputs:
!     sigma_g  : angular std of the beam per component [rad] ('gaussian')
!     beam_div : full-angle beam divergence [rad] ('uniform_cone')
!     mode     : 'gaussian' | 'uniform_cone'
!     alpha    : off-nadir tilt angle [rad]
!     H_apex   : altitude of the cone apex above the surface [m] (= H')
!     x_A      : horizontal offset of the apex in the O frame [m] (< 0)
!   Outputs:
!     x_surf, y_surf : surface impact coordinates in the O frame [m]
!     u_inc(3)       : incident propagation versor in the MC frame
!                      (z-up: into water is u_inc(3) < 0)
!
! Method:
!   Photons start from the cone apex at (x_A, 0, -H_apex) in an AUXILIARY
!   z-down frame; the apex altitude in the global MC z-up frame is
!   +H_apex. The propagation versor is computed and rotated in the
!   auxiliary z-down frame (sections 2-4), then explicitly flipped to the
!   global MC z-up convention in section 5 (u_inc(3) = -kz).
!   'gaussian'     : Rayleigh radial divergence, exact CDF inversion
!                    theta' = sigma_g sqrt(-2 ln R1) (Box-Muller polar
!                    core). Equivalent to a 2-D Gaussian in the angles with
!                    std sigma_g per component, intensity
!                    I(theta) ~ exp(-theta^2/(2 sigma_g^2)); I = I0/e^2 at
!                    theta = 2 sigma_g, hence sigma_g = beam_div/4 for a
!                    full-angle divergence measured at 1/e^2.
!   'uniform_cone' : cos(theta') uniform in [cos(beam_div/2), 1] (sharp
!                    edge at beam_div/2).
!   On a near-grazing ray (kz <= 0) the impact is flagged with -huge().
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  'uniform_cone' half-angle beam_div/2 (was FOV/2);
!                      header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8),          Intent(In)  :: sigma_g    ! Beam angular std [rad].
Real(8),          Intent(In)  :: beam_div   ! Full-angle divergence [rad].
Character(len=*), Intent(In)  :: mode       ! Sampling mode.
Real(8),          Intent(In)  :: alpha      ! Off-nadir angle [rad].
Real(8),          Intent(In)  :: H_apex     ! Apex altitude [m].
Real(8),          Intent(In)  :: x_A        ! Apex horizontal position [m].

! Scalar arguments with intent(out):
Real(8),          Intent(Out) :: x_surf     ! Surface impact x [m].
Real(8),          Intent(Out) :: y_surf     ! Surface impact y [m].

! Array  arguments with intent(out):
Real(8),          Intent(Out) :: u_inc(3)   ! Incident versor, z-up frame.

! Local parameters:
Real(8), Parameter :: PI = acos(-1.0d0)     ! Pi.

! Local scalars:
Real(8) :: R1          ! Uniform random number (radial).
Real(8) :: R2          ! Uniform random number (azimuth).
Real(8) :: theta_p     ! Angle from the cone axis [rad].
Real(8) :: phi_p       ! Azimuth around the cone axis [rad].
Real(8) :: vx          ! Versor x, cone-axis frame.
Real(8) :: vy          ! Versor y, cone-axis frame.
Real(8) :: vz          ! Versor z, cone-axis frame.
Real(8) :: kx          ! Versor x, tilted frame.
Real(8) :: ky          ! Versor y, tilted frame.
Real(8) :: kz          ! Versor z, tilted frame (z down).
Real(8) :: t_path      ! Path length to z = 0 [m].
Real(8) :: cos_min     ! cos(beam_div/2).

!- End of header ------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Random draw in the cone-axis frame
!-----------------------------------------------------------------------
Call random_number(R1)
Call random_number(R2)
If (R1 < tiny(1.0d0)) R1 = tiny(1.0d0)   ! Guard against log(0).

! Azimuth uniform on [0, 2 pi): PDF = 1/(2 pi) (isotropy).
phi_p = 2.0d0 * PI * R2

! Radial divergence angle, depending on the mode.
Select Case (trim(mode))
Case ('gaussian')
  ! Exact inversion of the Rayleigh CDF:
  !   CDF = 1 - exp(-theta'^2/(2 sigma_g^2))
  !   => theta' = sigma_g sqrt(-2 ln R1).
  theta_p = sigma_g * sqrt(-2.0d0 * log(R1))
Case ('uniform_cone')
  ! Uniform in solid angle within the cone of half-aperture beam_div/2.
  cos_min = cos(0.5d0 * beam_div)
  theta_p = acos(1.0d0 - R1 * (1.0d0 - cos_min))
Case Default
  theta_p = sigma_g * sqrt(-2.0d0 * log(R1))   ! Fallback: gaussian.
End Select

!-----------------------------------------------------------------------
! 2. Propagation versor in the auxiliary z-down frame
!-----------------------------------------------------------------------
! Cone axis = local z, physically pointing down: vz > 0 means downward
! in this auxiliary frame only.
vx = sin(theta_p) * cos(phi_p)
vy = sin(theta_p) * sin(phi_p)
vz = cos(theta_p)

!-----------------------------------------------------------------------
! 3. Rigid rotation R_y(alpha) in the auxiliary z-down frame
!-----------------------------------------------------------------------
!   R_y(alpha) = [  cos a   0   sin a ]
!                [    0     1     0   ]
!                [ -sin a   0   cos a ]
kx =  cos(alpha) * vx + sin(alpha) * vz
ky =  vy
kz = -sin(alpha) * vx + cos(alpha) * vz

!-----------------------------------------------------------------------
! 4. Intersection with the surface z = 0 (auxiliary frame)
!-----------------------------------------------------------------------
! Apex at (x_A, 0, -H_apex) in the auxiliary frame: t = H_apex/kz.
If (kz <= 0.0d0) Then
  x_surf   = -huge(1.0d0)
  y_surf   = -huge(1.0d0)
  u_inc(1) = 0.0d0
  u_inc(2) = 0.0d0
  u_inc(3) = 0.0d0
  Return
End If
t_path = H_apex / kz
x_surf = x_A + kx * t_path
y_surf =       ky * t_path

!-----------------------------------------------------------------------
! 5. Explicit flip to the global MC z-up frame
!-----------------------------------------------------------------------
! Only the z component changes: same horizontal axes in both frames.
u_inc(1) =  kx
u_inc(2) =  ky
u_inc(3) = -kz   ! Into water: u_inc(3) < 0 in the z-up MC frame.

End Subroutine generate_photon

!+ FOV test at depth, rigid-translation approximation
!
Pure Function photon_in_fov_simplified &
  ( x_ph,        & ! in
    y_ph,        & ! in
    z_depth,     & ! in
    m_axis,      & ! in
    n_axis,      & ! in
    tan_alpha_r) & ! in
  Result (in_fov)

! Description:
!   Checks whether a photon at depth z_depth falls within the FOV of the
!   receiver using the simplified rigid-translation approximation.
!
!   Inputs:
!     x_ph, y_ph  : horizontal coordinates of the photon in system O [m]
!     z_depth     : depth of the photon below the surface [m] (> 0)
!     m_axis      : semi-major axis of the FOV ellipse [m]
!     n_axis      : semi-minor axis of the FOV ellipse [m]
!     tan_alpha_r : tan of the refracted optical-axis angle
!   Output:
!     .true. if the photon is inside the FOV, .false. otherwise
!
! Method:
!   The underwater FOV at depth z_depth is the surface footprint ellipse
!   (semi-axes m, n) translated rigidly along x by the refraction
!   displacement of the optical axis:
!     Delta_x(z) = z tan(alpha_r),  sin(alpha_r) = sin(alpha)/n_water,
!     in FOV if ((x - Delta_x)/m)^2 + (y/n)^2 <= 1.
!   Approximations (all O(1e-4) or smaller):
!     - small beam divergence (sigma ~ 50 urad): incidence angles vary
!       little over the footprint;
!     - m/|x_A| ~ 1e-3: azimuthal direction omega ~ const over the
!       ellipse;
!     - geometric broadening DeltaM(z) ~ 0.04 m at 100 m depth << M.
!   No iteration required: O(1) evaluation time.
!   Convention: z-axis points UPWARD (global MC frame), origin at centre
!   O; z_depth is a POSITIVE depth below the surface (z_depth = -z).
!   Ref.: Shangguan et al. (2024), Eq. 17-22, simplified first order.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.1     2026-11-05  Shift tabulated at the bin centre.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Shift computed at the photon depth (argument
!                      tan_alpha_r instead of delta_x_fov).
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In) :: x_ph          ! Photon x [m].
Real(8), Intent(In) :: y_ph          ! Photon y [m].
Real(8), Intent(In) :: z_depth       ! Photon depth, positive [m].
Real(8), Intent(In) :: m_axis        ! Semi-major axis [m].
Real(8), Intent(In) :: n_axis        ! Semi-minor axis [m].
Real(8), Intent(In) :: tan_alpha_r   ! tan of the refracted axis angle.

! Function result:
Logical :: in_fov                    ! .true. inside the FOV.

! Local scalars:
Real(8) :: x_E_displaced             ! x relative to the shifted ellipse.
Real(8) :: ellipse_val               ! Ellipse equation value.

!- End of header ------------------------------------------------------

! Translate the point by the refraction displacement of the axis.
x_E_displaced = x_ph - z_depth * tan_alpha_r
! Test membership in the surface ellipse.
ellipse_val = (x_E_displaced / m_axis)**2 + (y_ph / n_axis)**2
in_fov = (ellipse_val <= 1.0d0)

End Function photon_in_fov_simplified

!+ Index of a depth on the output grid
!
Pure Function depth_bin &
  ( depth, & ! in
    dz,    & ! in
    zmax,  & ! in
    nz )   & ! in
  Result (ibin)

! Description:
!   Index of a depth on the output grid.
!   Grid: bin iz (1..nz) covers [(iz-1)*dz, iz*dz); bin 0 is reserved
!   for the surface and is never returned.
!
!   Inputs:
!     depth : positive depth (real or time-based) [m]
!     dz    : bin thickness [m]
!     zmax  : maximum depth of the grid [m], zmax = nz*dz
!     nz    : number of subsurface bins
!   Output:
!     bin index in 1..nz, or nz+1 if depth >= zmax (not scored)
!
! Method:
!   ibin = min(int(depth/dz) + 1, nz) if depth < zmax, else nz + 1.
!   The comparison is made in real arithmetic before int() to avoid
!   integer overflow for very large depths.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.1     2026-11-05  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In) :: depth   ! Depth to bin [m].
Real(8), Intent(In) :: dz      ! Bin thickness [m].
Real(8), Intent(In) :: zmax    ! Maximum grid depth [m].
Integer, Intent(In) :: nz      ! Number of subsurface bins.

! Function result:
Integer :: ibin                ! Bin index, nz+1 beyond zmax.

!- End of header ------------------------------------------------------

If (depth < zmax) Then
  ibin = min(int(depth / dz) + 1, nz)
Else
  ibin = nz + 1
End If

End Function depth_bin

!-----------------------------------------------------------------------
! 3. Sea surface
!-----------------------------------------------------------------------

!+ Look-up tables of Fresnel reflectance and transmittance
!
Subroutine generate_fresnel_lut &
  ( N_LUT,    & ! in
    n_water,  & ! in
    LUT_cos,  & ! out
    LUT_R_aa, & ! out
    LUT_T_aw, & ! out
    LUT_R_ww, & ! out
    LUT_T_wa )  ! out

! Description:
!   Builds look-up tables of scalar (unpolarised) Fresnel reflectance and
!   transmittance for an air <-> water interface as a function of
!   cos(theta). Gives fast access to the Fresnel coefficients for the two
!   directions: air->water and water->air.
!
!   Inputs:
!     N_LUT   : number of LUT entries
!     n_water : refractive index of seawater
!   Outputs:
!     LUT_cos  : cos(theta) nodes associated with the coefficients
!     LUT_R_aa : reflectance,   air-air     (reflection in air)
!     LUT_T_aw : transmittance, air->water  (transmission in water)
!     LUT_R_ww : reflectance,   water-water (reflection in water)
!     LUT_T_wa : transmittance, water->air  (transmission in air)
!
! Method:
!   Nodes uniform in cos(theta) from 0.02 (~87.7 deg) to 1.0 (0 deg):
!     cos_i = 0.02 + 0.98 (i - 1)/(N_LUT - 1).
!   Coefficients from Fresnel_Refl_Scalar and Fresnel_Trans_Scalar with
!   (n1, n2) = (1, n_water) and (n_water, 1).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.3     2026-11-07  n_water as argument (was 1.33).
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Integer, Intent(In)  :: N_LUT             ! Number of entries.
Real(8), Intent(In)  :: n_water           ! Refractive index of seawater.

! Array  arguments with intent(out):
Real(8), Intent(Out) :: LUT_cos(N_LUT)    ! cos(incidence) nodes.
Real(8), Intent(Out) :: LUT_R_aa(N_LUT)   ! Reflectance, air side.
Real(8), Intent(Out) :: LUT_T_aw(N_LUT)   ! Transmittance, air->water.
Real(8), Intent(Out) :: LUT_R_ww(N_LUT)   ! Reflectance, water side.
Real(8), Intent(Out) :: LUT_T_wa(N_LUT)   ! Transmittance, water->air.

! Local scalars:
Real(8)    :: cos_min                     ! First node.
Real(8)    :: cos_max                     ! Last node.
Real(8)    :: cos_theta                   ! Current node.
Complex(8) :: CN1                         ! Index of the incident medium.
Complex(8) :: CN2                         ! Index of the other medium.
Integer    :: i                           ! Node index.

!- End of header ------------------------------------------------------

cos_min = 0.02d0
cos_max = 1.0d0
Do i = 1, N_LUT
  cos_theta = cos_min + (cos_max - cos_min) * real(i - 1, 8)           &
                        / real(N_LUT - 1, 8)
  LUT_cos(i) = cos_theta
  ! First case: air -> water.
  CN1 = DCMP(1.0d0)
  CN2 = DCMP(n_water)
  LUT_R_aa(i) = Fresnel_Refl_Scalar(CN1, CN2, cos_theta)
  LUT_T_aw(i) = Fresnel_Trans_Scalar(CN1, CN2, cos_theta)
  ! Second case: water -> air.
  CN1 = DCMP(n_water)
  CN2 = DCMP(1.0d0)
  LUT_R_ww(i) = Fresnel_Refl_Scalar(CN1, CN2, cos_theta)
  LUT_T_wa(i) = Fresnel_Trans_Scalar(CN1, CN2, cos_theta)
End Do

End Subroutine generate_fresnel_lut

!+ Double-precision complex number from real and imaginary parts
!
Function DCMP &
  ( Re,  & ! in
    Im ) & ! in, optional
  Result (z)

! Description:
!   Builds a double-precision complex number from a mandatory real part
!   and an optional imaginary part.
!     DCMP(1.0d0)         -> (1.0d0,  0.0d0)
!     DCMP(1.0d0, -2.0d0) -> (1.0d0, -2.0d0)
!
! Method:
!   cmplx(Re, Im, kind=8); Im = 0 if absent.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)           :: Re   ! Real part.
Real(8), Intent(In), Optional :: Im   ! Imaginary part (0 if absent).

! Function result:
Complex(8) :: z                       ! Complex number.

!- End of header ------------------------------------------------------

If (present(Im)) Then
  z = cmplx(Re, Im, kind=8)
Else
  z = cmplx(Re, 0.0d0, kind=8)
End If

End Function DCMP

!+ Unpolarised Fresnel transmittance
!
Function Fresnel_Trans_Scalar &
  ( CN1,  & ! in
    CN2,  & ! in
    Coi ) & ! in
  Result (TF)

! Description:
!   Scalar (unpolarised) fraction of energy transmitted across the
!   interface between two media.
!
!   Inputs:
!     CN1 : refractive index of the incident medium    (complex(8))
!     CN2 : refractive index of the transmitted medium (complex(8))
!     Coi : cosine of the incidence angle
!   Output:
!     TF  : transmission coefficient (between 0 and 1)
!
! Method:
!   Averages the parallel (L) and perpendicular (P) polarisations:
!     t_L = 2 sqrt(n1 n2 cos_i cos_t) / (n2 cos_i + n1 cos_t)
!     t_P = 2 sqrt(n1 n2 cos_i cos_t) / (n1 cos_i + n2 cos_t)
!     TF  = ( |t_L|^2 + |t_P|^2 ) / 2
!   with cos_t from Compute_CosT (Snell).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Complex(8), Intent(In) :: CN1   ! Index of the incident medium.
Complex(8), Intent(In) :: CN2   ! Index of the transmitting medium.
Real(8),    Intent(In) :: Coi   ! cos(incidence).

! Function result:
Real(8) :: TF                   ! Transmittance.

! Local scalars:
Complex(8) :: Cot               ! cos(transmission).
Complex(8) :: Tc_L              ! Amplitude coefficient, parallel.
Complex(8) :: Tc_P              ! Amplitude coefficient, perpendicular.

!- End of header ------------------------------------------------------

! Transmitted-angle cosine via Snell's law.
Cot  = Compute_CosT(CN1, CN2, Coi)
! Amplitude transmission coefficients.
Tc_L = (2.0d0 * sqrt(CN1*CN2*Coi*Cot)) / (CN2*Coi + CN1*Cot)
Tc_P = (2.0d0 * sqrt(CN1*CN2*Coi*Cot)) / (CN1*Coi + CN2*Cot)
! Transmitted energy: average of the two polarisations.
TF = 0.5d0 * (abs(Tc_L)**2 + abs(Tc_P)**2)

End Function Fresnel_Trans_Scalar

!+ Unpolarised Fresnel reflectance
!
Function Fresnel_Refl_Scalar &
  ( CN1,  & ! in
    CN2,  & ! in
    Coi ) & ! in
  Result (RF)

! Description:
!   Scalar (unpolarised) fraction of energy reflected at the interface
!   between two media.
!
!   Inputs:
!     CN1 : refractive index of the incident medium    (complex(8))
!     CN2 : refractive index of the transmitted medium (complex(8))
!     Coi : cosine of the incidence angle
!   Output:
!     RF  : reflection coefficient (between 0 and 1)
!
! Method:
!   Averages the parallel (L) and perpendicular (P) polarisations:
!     r_L = (n2 cos_i - n1 cos_t) / (n2 cos_i + n1 cos_t)
!     r_P = (n1 cos_i - n2 cos_t) / (n1 cos_i + n2 cos_t)
!     RF  = ( |r_L|^2 + |r_P|^2 ) / 2
!   with cos_t from Compute_CosT (complex beyond the critical angle,
!   where RF = 1).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Complex(8), Intent(In) :: CN1   ! Index of the incident medium.
Complex(8), Intent(In) :: CN2   ! Index of the transmitting medium.
Real(8),    Intent(In) :: Coi   ! cos(incidence).

! Function result:
Real(8) :: RF                   ! Reflectance.

! Local scalars:
Complex(8) :: Cos_t             ! cos(transmission).
Complex(8) :: Tc_L              ! Amplitude coefficient, parallel.
Complex(8) :: Tc_P              ! Amplitude coefficient, perpendicular.

!- End of header ------------------------------------------------------

! Transmitted-angle cosine (used only for the amplitude).
Cos_t = Compute_CosT(CN1, CN2, Coi)
! Amplitude reflection coefficients.
Tc_L = (CN2*Coi - CN1*Cos_t) / (CN2*Coi + CN1*Cos_t)
Tc_P = (CN1*Coi - CN2*Cos_t) / (CN1*Coi + CN2*Cos_t)
! Reflected energy: average of the two polarisations.
RF = 0.5d0 * (abs(Tc_L)**2 + abs(Tc_P)**2)

End Function Fresnel_Refl_Scalar

!+ Cosine of the transmission angle from Snell's law
!
Function Compute_CosT &
  ( CN1,  & ! in
    CN2,  & ! in
    Coi ) & ! in
  Result (Cot)

! Description:
!   Cosine of the transmitted angle (Cot) from the generalised Snell's
!   law for complex refractive indices. Supports both denser and less
!   dense transmitted media.
!
!   Inputs:
!     CN1 : refractive index of the incident medium    (complex(8))
!     CN2 : refractive index of the transmitted medium (complex(8))
!     Coi : cosine of the incidence angle
!   Output:
!     Cot : cosine of the transmitted angle (complex(8))
!
! Method:
!   Snell: sin^2(theta_t) = (n1/n2)^2 (1 - cos^2(theta_i)),
!          cos(theta_t) = sqrt(1 - sin^2(theta_t)),
!   imaginary beyond the critical angle (total internal reflection).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Complex(8), Intent(In) :: CN1   ! Index of the incident medium.
Complex(8), Intent(In) :: CN2   ! Index of the transmitting medium.
Real(8),    Intent(In) :: Coi   ! cos(incidence).

! Function result:
Complex(8) :: Cot               ! cos(transmission).

! Local scalars:
Complex(8) :: sin2_t            ! sin^2(transmission).

!- End of header ------------------------------------------------------

sin2_t = (CN1 / CN2)**2 * (1.0d0 - Coi**2)
Cot = sqrt(1.0d0 - sin2_t)

End Function Compute_CosT

!+ LUT value at the node nearest to a given cosine
!
Subroutine get_lut_nearest &
  ( cos_theta, & ! in
    LUTcos_in, & ! in
    LUTval_in, & ! in
    coeff )      ! out

! Description:
!   Returns the LUT coefficient that best matches a given incidence-angle
!   cosine (nearest neighbour, no interpolation).
!
!   Inputs:
!     cos_theta    : cosine of the incidence angle [-1, 1]
!     LUTcos_in(:) : reference cosines (uniform grid)
!     LUTval_in(:) : coefficients associated with each cosine
!   Output:
!     coeff        : LUT coefficient (nearest neighbour)
!
! Method:
!   cos_theta is clipped to [-1, 1] (safety clipping). The grid step is
!   (LUTcos(n) - LUTcos(1))/(n - 1); the nearest node is guessed as
!   nint((cos - LUTcos(1))/step) + 1 and checked against its two
!   neighbours; on ties the lowest index wins. This gives the same node
!   as the linear scan of v0.0-v0.2.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code (linear scan). Valerio Ciattaglia
! v0.3     2026-11-07  Direct index on the uniform grid; commented
!                      linear-interpolation block removed.
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)  :: cos_theta      ! cos(incidence).

! Array  arguments with intent(in):
Real(8), Intent(In)  :: LUTcos_in(:)   ! Uniform cos grid.
Real(8), Intent(In)  :: LUTval_in(:)   ! Values on the grid.

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: coeff          ! Value at the nearest node.

! Local scalars:
Integer :: n_lut                       ! Number of nodes.
Integer :: i                           ! Node index.
Integer :: i_guess                     ! Nearest node from the step.
Integer :: idx_closest                 ! Nearest node.
Real(8) :: cos_clip                    ! cos_theta clipped to [-1, 1].
Real(8) :: step                        ! Grid step.
Real(8) :: diff                        ! Distance from a node.
Real(8) :: min_diff                    ! Smallest distance so far.

!- End of header ------------------------------------------------------

n_lut = size(LUTcos_in)
cos_clip = max(-1.0d0, min(1.0d0, cos_theta))

step = (LUTcos_in(n_lut) - LUTcos_in(1)) / real(n_lut - 1, 8)
i_guess = nint((cos_clip - LUTcos_in(1)) / step) + 1
i_guess = max(1, min(n_lut, i_guess))

! Neighbours of the guess, lowest index first.
idx_closest = max(1, i_guess - 1)
min_diff    = abs(LUTcos_in(idx_closest) - cos_clip)
Do i = idx_closest + 1, min(n_lut, i_guess + 1)
  diff = abs(LUTcos_in(i) - cos_clip)
  If (diff < min_diff) Then
    min_diff    = diff
    idx_closest = i
  End If
End Do

coeff = LUTval_in(idx_closest)

End Subroutine get_lut_nearest

!+ Approximate inverse error function
!
Pure Function erfinv_win &
  ( x ) & ! in
  Result (y)

! Description:
!   Inverse error function erf^-1(x), in Winitzki's closed approximation.
!
!   Inputs:
!     x : argument, -1 < x < 1
!   Output:
!     y : approximation of erf^-1(x)
!
! Method:
!   With a = 0.147 and L = ln(1 - x^2):
!     t = 2 / (pi a) + L / 2,
!     y = sign(x) sqrt( sqrt(t^2 - L / a) - t ).
!   Winitzki's approximation (cited from memory, to be verified).
!   Measured against the exact inverse (runs/2026-11-08_slope_analitica):
!   largest relative error 0.19%, near y = 2.5.
!   The inner difference can round to a small negative number for x near
!   0: it is clipped at 0 before the square root.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.6     2026-11-08  Original code. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In) :: x        ! Argument, -1 < x < 1.

! Function result:
Real(8) :: y                    ! Approximation of erf^-1(x).

! Local parameters:
Real(8), Parameter :: pi_loc = 3.141592653589793d0   ! Pi.
Real(8), Parameter :: a_win  = 0.147d0               ! Winitzki constant.

! Local scalars:
Real(8) :: ln_term              ! ln(1 - x^2).
Real(8) :: t                    ! 2 / (pi a) + ln_term / 2.

!- End of header ------------------------------------------------------

ln_term = log(1.0d0 - x*x)
t       = 2.0d0 / (pi_loc * a_win) + 0.5d0 * ln_term
y       = sign(sqrt(max(0.0d0,                                         &
                        sqrt(t*t - ln_term / a_win) - t)), x)

End Function erfinv_win

!+ Samples one component of the sea-surface slope
!
Subroutine sample_slope &
  ( sigma, & ! in
    s )      ! out

! Description:
!   Draws one slope component from the Gaussian distribution of the sea
!   surface, by analytic inversion of its cumulative distribution.
!
!   Inputs:
!     sigma : root of the total mean-square slope sigma2 (0 = flat sea)
!   Output:
!     s     : slope component
!
! Method:
!   1-D slope distribution (Nakajima and Tanaka 1983, Eq. 14):
!     p(s) = exp(-s^2/sigma^2) / (sqrt(pi) sigma).
!   sigma^2 is the total mean-square slope: each component has variance
!   sigma^2/2, so that <s_x^2 + s_y^2> = sigma^2, consistent with
!   lambda_saunders and with the Cox-Munk pdf of the glint.
!   Cumulative distribution and its inverse, with U uniform in (0, 1):
!     F(s) = (1 + erf(s/sigma)) / 2,   s = sigma erf^-1(2 U - 1).
!   Flat sea (sigma = 0): s = 0. The random number is drawn anyway, so
!   that flat-sea runs keep their random sequence.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.6     2026-11-08  Original code, replaces integra_slope and
!                      sample_from_cdf. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)  :: sigma   ! Root of sigma2.

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: s       ! Slope component.

! Local scalars:
Real(8) :: u_rnd                ! Uniform random number.

!- End of header ------------------------------------------------------

Call random_number(u_rnd)
If (sigma <= 0.0d0) Then
  s = 0.0d0                     ! Flat sea.
  Return
End If
If (u_rnd <= 0.0d0) u_rnd = epsilon(1.0d0)   ! Avoid log(0) in erfinv_win.
s = sigma * erfinv_win(2.0d0 * u_rnd - 1.0d0)

End Subroutine sample_slope

!+ Smith shadowing function for an isotropic Gaussian surface
!
Pure Function lambda_saunders &
  ( mu_zenith, & ! in
    sigma2 )   & ! in
  Result (lambda)

! Description:
!   Computes the Smith/Saunders scalar shadowing function for an
!   isotropic Gaussian sea surface.
!
!   Inputs:
!     mu_zenith : cosine of the zenith angle (any sign; |cos| used)
!     sigma2    : total mean-square slope of the surface; in the code
!                 sigma2 = 0.003 + 0.00512 u10 (Cox and Munk 1954)
!   Output:
!     Lambda(v) : shadowing risk parameter
!
! Method:
!   Ref.: Nakajima and Tanaka (1983), J. Quant. Spectrosc. Radiat.
!   Transfer 29(6).
!     v = |cos(theta)| / sqrt(sigma2 sin^2(theta)) = cot(theta)/sigma
!                                                                (Eq. 19)
!     Lambda = ( exp(-v^2)/(v sqrt(pi)) - erfc(v) ) / 2          (Eq. 20)
!   The paper uses mu > 0 implicitly (ray above the horizon). For
!   ascending photons (mu < 0) surface symmetry gives the same Lambda,
!   hence abs() here.
!   Uses the Fortran 2008 intrinsic erfc(v): full double precision.
!   Lambda -> 0 for v -> +inf (vertical ray), set to 0 for v >= 5;
!   Lambda -> +inf for v -> 0 (grazing ray), set to 1e30 for v < 1e-10.
!   Negative values from rounding near v = 5 are clamped to 0.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In) :: mu_zenith    ! Cosine of the zenith angle.
Real(8), Intent(In) :: sigma2       ! Total mean-square slope.

! Function result:
Real(8) :: lambda                   ! Shadowing function Lambda(v).

! Local parameters:
Real(8), Parameter :: INV_SQPI = 0.56418958354775628d0   ! 1/sqrt(pi).
Real(8), Parameter :: V_MAX    = 5.0d0      ! Lambda = 0 above this v.

! Local scalars:
Real(8) :: v                        ! cot(theta)/sigma.
Real(8) :: v2                       ! v^2.
Real(8) :: denom2                   ! sigma2 sin^2(theta).
Real(8) :: erfc_v                   ! erfc(v).

!- End of header ------------------------------------------------------

! Stability for the vertical ray (v -> +inf, Lambda -> 0).
denom2 = sigma2 * (1.0d0 - mu_zenith**2)
If (denom2 < 1.0d-14) Then
  lambda = 0.0d0
  Return
End If
v = abs(mu_zenith) / sqrt(denom2)
! Asymptotic tail (Lambda -> 0).
If (v >= V_MAX) Then
  lambda = 0.0d0
  Return
End If
! Grazing incidence (Lambda -> +inf).
If (v < 1.0d-10) Then
  lambda = 1.0d30
  Return
End If
v2     = v * v
erfc_v = erfc(v)
lambda = 0.5d0 * (INV_SQPI * exp(-v2) / v - erfc_v)
! Clamp numerical noise near V_MAX.
If (lambda < 0.0d0) lambda = 0.0d0

End Function lambda_saunders

!+ Air-to-water crossing of the rough sea surface
!
Subroutine surface_interaction &
  ( u,            & ! in
    N_LUT,        & ! in
    LUT_cos,      & ! in
    LUT_R_aa,     & ! in
    LUT_T_aw,     & ! in
    sigma2,       & ! in
    n_water,      & ! in
    refl,         & ! out
    trans,        & ! out
    k_ref,        & ! out
    k_tr )          ! out

! Description:
!   Computes reflection / transmission at a rough air-water interface for
!   a photon arriving from the air.
!
!   Inputs:
!     u               : incident versor of the photon
!     N_LUT           : length of the Fresnel LUT tables
!     LUT_cos(:)      : cos(theta) nodes of the Fresnel LUTs
!     LUT_R_aa(:)     : reflectance  air->air   (reflection in air)
!     LUT_T_aw(:)     : transmittance air->water (transmission in water)
!     sigma2          : total mean-square slope (Cox and Munk 1954)
!     n_water         : refractive index of seawater
!   Outputs:
!     refl     : Fresnel reflectance corrected for bidirectional shadowing
!     trans    : Fresnel transmittance corrected for bidirectional
!                shadowing
!     k_ref(3) : reflected unit vector
!     k_tr(3)  : transmitted unit vector (into water)
!
! Method:
!   Slope distribution sampled by analytic inversion, sample_slope
!   (Nakajima and Tanaka 1983, Eq. 14); facet normal
!   n = (-s_x, -s_y, 1)/|.|, toward the air.
!   Rejection sampling: facets with cos(theta_i) = -k_in.n <= 0 are
!   discarded (Heaviside). k_in and n are normalised with the same
!   sqrt(sum(v**2)) idiom.
!   Fresnel R and T from the LUTs at cos(theta_i).
!   Directions, eta = 1/n_water:
!     k_ref = k_in + 2 cos(theta_i) n
!     k_tr  = eta k_in + (eta cos(theta_i) - cos(theta_t)) n
!   Bidirectional shadowing (Nakajima and Tanaka 1983, Eq. 17; Smith
!   1967):
!     refl  / (1 + Lambda_i + Lambda_r)
!     trans / (1 + Lambda_i + Lambda_t)
!   Lambda computed from abs(k(3)) to account for the azimuthal rotation
!   of k_ref and k_tr relative to the incident plane. Calls
!   lambda_saunders (Fortran 2008 intrinsic erfc).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.3     2026-11-07  n_water as argument (was 1.33).
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
! v0.6     2026-11-08  Slopes from sample_slope, tables removed.
!                      Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Integer, Intent(In)  :: N_LUT                 ! Number of LUT entries.
Real(8), Intent(In)  :: sigma2                ! Total mean-square slope.
Real(8), Intent(In)  :: n_water               ! Index of seawater.

! Array  arguments with intent(in):
Real(8), Intent(In)  :: u(3)                  ! Incident versor.
Real(8), Intent(In)  :: LUT_cos(N_LUT)        ! LUT cos grid.
Real(8), Intent(In)  :: LUT_R_aa(N_LUT)       ! Reflectance, air side.
Real(8), Intent(In)  :: LUT_T_aw(N_LUT)       ! Transmittance, a->w.

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: refl                  ! Reflectance, shadowed.
Real(8), Intent(Out) :: trans                 ! Transmittance, shadowed.

! Array  arguments with intent(out):
Real(8), Intent(Out) :: k_ref(3)              ! Reflected versor.
Real(8), Intent(Out) :: k_tr(3)               ! Refracted versor.

! Local scalars:
Real(8)    :: slp_x                           ! Facet slope along x.
Real(8)    :: slp_y                           ! Facet slope along y.
Real(8)    :: cos_theta_i                     ! cos(local incidence).
Real(8)    :: cos_theta_t                     ! cos(local transmission).
Real(8)    :: eta                             ! n_air/n_water.
Complex(8) :: CN1                             ! Index of air.
Complex(8) :: CN2                             ! Index of water.
Real(8)    :: mag                             ! Vector norm.
Real(8)    :: Lambda_i                        ! Shadowing, incident leg.
Real(8)    :: Lambda_r                        ! Shadowing, reflected leg.
Real(8)    :: Lambda_t                        ! Shadowing, transmitted.
Logical    :: accepted                        ! Facet accepted.

! Local arrays:
Real(8)    :: n_raw(3)                        ! Unnormalised normal.
Real(8)    :: n_hat(3)                        ! Facet normal, toward air.
Real(8)    :: k_in(3)                         ! Normalised incident versor.

!- End of header ------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Facet sampling (rejection of back-facing facets)
!-----------------------------------------------------------------------
k_in = u / sqrt(sum(u**2))
accepted = .false.
Do While (.not. accepted)
  Call sample_slope(sqrt(sigma2), slp_x)
  Call sample_slope(sqrt(sigma2), slp_y)
  ! Local surface normal, oriented toward the air (+z).
  n_raw = [ -slp_x, -slp_y, 1.0d0 ]
  n_hat = n_raw / sqrt(sum(n_raw**2))
  ! Cosine of the incidence angle on the facet.
  cos_theta_i = -dot_product(k_in, n_hat)
  If (cos_theta_i > 0.0d0) accepted = .true.
End Do

!-----------------------------------------------------------------------
! 2. Fresnel coefficients from the LUTs
!-----------------------------------------------------------------------
Call get_lut_nearest(cos_theta_i, LUT_cos, LUT_R_aa, refl)
Call get_lut_nearest(cos_theta_i, LUT_cos, LUT_T_aw, trans)

!-----------------------------------------------------------------------
! 3. Reflected and refracted directions
!-----------------------------------------------------------------------
CN1         = DCMP(1.0d0)
CN2         = DCMP(n_water)
eta         = 1.0d0 / n_water
cos_theta_t = real(Compute_CosT(CN1, CN2, cos_theta_i), 8)
k_ref = k_in + 2.0d0 * cos_theta_i * n_hat
k_tr  = eta * k_in + (eta * cos_theta_i - cos_theta_t) * n_hat
mag = sqrt(sum(k_ref**2))
If (mag > 0.0d0) k_ref = k_ref / mag
mag = sqrt(sum(k_tr **2))
If (mag > 0.0d0) k_tr  = k_tr  / mag

!-----------------------------------------------------------------------
! 4. Bidirectional shadowing (Smith 1967)
!-----------------------------------------------------------------------
! abs(k(3)) for all three legs: k_ref and k_tr may carry an azimuthal
! component absent in the incident plane.
Lambda_i = lambda_saunders(abs(k_in(3)),  sigma2)
Lambda_r = lambda_saunders(abs(k_ref(3)), sigma2)
Lambda_t = lambda_saunders(abs(k_tr(3)),  sigma2)
refl  = refl  / (1.0d0 + Lambda_i + Lambda_r)
trans = trans / (1.0d0 + Lambda_i + Lambda_t)

End Subroutine surface_interaction

!+ Water-to-air crossing of the rough sea surface
!
Subroutine surface_interaction_wa &
  ( u,            & ! in
    N_LUT,        & ! in
    LUT_cos,      & ! in
    LUT_R_ww,     & ! in
    LUT_T_wa,     & ! in
    sigma2,       & ! in
    n_water,      & ! in
    refl,         & ! out
    trans,        & ! out
    k_ref,        & ! out
    k_tr )          ! out

! Description:
!   Reflection / transmission for a photon travelling WATER -> AIR across
!   a rough interface, with explicit handling of Total Internal
!   Reflection (TIR).
!
!   Inputs:
!     u               : incident versor of the photon (water side)
!     N_LUT           : length of the Fresnel LUT tables
!     LUT_cos(:)      : cos(theta) nodes of the Fresnel LUTs
!     LUT_R_ww(:)     : reflectance  water->water (reflection in water)
!     LUT_T_wa(:)     : transmittance water->air  (transmission in air)
!     sigma2          : total mean-square slope (Cox and Munk 1954)
!     n_water         : refractive index of seawater
!   Outputs:
!     refl     : Fresnel reflectance corrected for bidirectional shadowing
!     trans    : Fresnel transmittance corrected for bidirectional
!                shadowing
!     k_ref(3) : reflected unit vector (water side)
!     k_tr(3)  : transmitted unit vector (air side); zero vector if TIR
!
! Method:
!   Slope distribution sampled by analytic inversion, sample_slope
!   (Nakajima and Tanaka 1983, Eq. 14). Rejection sampling: facets with
!   cos(theta_i) <= 0 are
!   discarded (Heaviside).
!   Convention: n oriented toward the air (+z); cos(theta_i) = k_in.n > 0
!   for an upward-going photon hitting the illuminated side of the facet.
!   TIR if cos(theta_i) < cos_crit = sqrt(1 - (1/n_water)^2): refl = 1,
!   trans = 0, shadowing applied to refl only. Otherwise, eta = n_water:
!     k_ref = k_in - 2 cos(theta_i) n
!     k_tr  = eta k_in - (eta cos(theta_i) - cos(theta_t)) n
!   Bidirectional shadowing (Smith 1967):
!     refl  / (1 + Lambda_i + Lambda_r)
!     trans / (1 + Lambda_i + Lambda_t)   [non-TIR only]
!   so refl + trans < 1 when Lambda > 0.
!   Lambda computed from abs(k(3)) to account for the azimuthal rotation
!   of k_ref and k_tr relative to the incident plane.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.3     2026-11-07  n_water as argument (was 1.33).
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
! v0.6     2026-11-08  Slopes from sample_slope, tables removed.
!                      Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Integer, Intent(In)  :: N_LUT                 ! Number of LUT entries.
Real(8), Intent(In)  :: sigma2                ! Total mean-square slope.
Real(8), Intent(In)  :: n_water               ! Index of seawater.

! Array  arguments with intent(in):
Real(8), Intent(In)  :: u(3)                  ! Incident versor (water).
Real(8), Intent(In)  :: LUT_cos(N_LUT)        ! LUT cos grid.
Real(8), Intent(In)  :: LUT_R_ww(N_LUT)       ! Reflectance, water side.
Real(8), Intent(In)  :: LUT_T_wa(N_LUT)       ! Transmittance, w->a.

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: refl                  ! Reflectance, shadowed.
Real(8), Intent(Out) :: trans                 ! Transmittance, shadowed.

! Array  arguments with intent(out):
Real(8), Intent(Out) :: k_ref(3)              ! Reflected versor.
Real(8), Intent(Out) :: k_tr(3)               ! Refracted versor, 0 if TIR.

! Local scalars:
Real(8)    :: slp_x                           ! Facet slope along x.
Real(8)    :: slp_y                           ! Facet slope along y.
Real(8)    :: cos_theta_i                     ! cos(local incidence).
Real(8)    :: eta                             ! n_water/n_air.
Real(8)    :: mag                             ! Vector norm.
Real(8)    :: cos_crit                        ! cos(critical angle).
Complex(8) :: cos_theta_t                     ! cos(local transmission).
Complex(8) :: CN1                             ! Index of water.
Complex(8) :: CN2                             ! Index of air.
Real(8)    :: Lambda_i                        ! Shadowing, incident leg.
Real(8)    :: Lambda_r                        ! Shadowing, reflected leg.
Real(8)    :: Lambda_t                        ! Shadowing, transmitted.
Logical    :: tir                             ! Total internal reflection.
Logical    :: accepted                        ! Facet accepted.

! Local arrays:
Real(8)    :: n_raw(3)                        ! Unnormalised normal.
Real(8)    :: n_hat(3)                        ! Facet normal, toward air.
Real(8)    :: k_in(3)                         ! Normalised incident versor.

!- End of header ------------------------------------------------------

! Refraction constants; max() avoids a NaN in sqrt.
eta      = n_water
cos_crit = sqrt(max(0.0d0, 1.0d0 - (1.0d0/eta)**2))
k_in = u / sqrt(sum(u**2))

!-----------------------------------------------------------------------
! 1. Facet sampling (rejection of back-facing facets)
!-----------------------------------------------------------------------
accepted = .false.
Do While (.not. accepted)
  Call sample_slope(sqrt(sigma2), slp_x)
  Call sample_slope(sqrt(sigma2), slp_y)
  ! Local surface normal, oriented toward the air (+z).
  n_raw = [ -slp_x, -slp_y, 1.0d0 ]
  n_hat = n_raw / sqrt(sum(n_raw**2))
  ! No minus sign: k_in (upward) and n_hat (toward air) are both
  ! positive along z for an illuminated facet.
  cos_theta_i = dot_product(k_in, n_hat)
  If (cos_theta_i > 0.0d0) accepted = .true.
End Do

!-----------------------------------------------------------------------
! 2. Total internal reflection or Fresnel split
!-----------------------------------------------------------------------
tir = cos_theta_i < cos_crit
If (tir) Then
  refl  = 1.0d0
  trans = 0.0d0
  k_ref = k_in - 2.0d0 * cos_theta_i * n_hat
  k_tr  = [ 0.0d0, 0.0d0, 0.0d0 ]
Else
  Call get_lut_nearest(cos_theta_i, LUT_cos, LUT_R_ww, refl)
  Call get_lut_nearest(cos_theta_i, LUT_cos, LUT_T_wa, trans)
  ! Snell: n1 = n_water (water) -> n2 = 1 (air).
  CN1 = DCMP(n_water)
  CN2 = DCMP(1.0d0)
  cos_theta_t = Compute_CosT(CN1, CN2, cos_theta_i)
  k_ref = k_in - 2.0d0 * cos_theta_i * n_hat
  k_tr  = eta * k_in - (eta * cos_theta_i - real(cos_theta_t, 8)) * n_hat
End If
mag = sqrt(sum(k_ref**2))
If (mag > 0.0d0) k_ref = k_ref / mag
mag = sqrt(sum(k_tr **2))
If (mag > 0.0d0) k_tr  = k_tr  / mag

!-----------------------------------------------------------------------
! 3. Bidirectional shadowing (Smith 1967)
!-----------------------------------------------------------------------
Lambda_i = lambda_saunders(abs(k_in(3)),  sigma2)
Lambda_r = lambda_saunders(abs(k_ref(3)), sigma2)
refl = refl / (1.0d0 + Lambda_i + Lambda_r)
If (.not. tir) Then
  Lambda_t = lambda_saunders(abs(k_tr(3)), sigma2)
  trans = trans / (1.0d0 + Lambda_i + Lambda_t)
End If

End Subroutine surface_interaction_wa

!-----------------------------------------------------------------------
! Surface exit from below.
! When the random walk brings a photon to z >= 0, it has reached the
! water-air interface, modelled as a Cox-Munk rough surface (the same
! facet model used for downwelling). surface_interaction_wa returns the
! reflected and transmitted weight fractions (refl, trans, each divided
! by its shadowing factor, so refl + trans <= 1), together with the
! reflected and transmitted direction versors (k_ref, k_tr) for the
! sampled facet. It accounts for:
!   - the local facet normal (Cox-Munk slope sampling);
!   - Total Internal Reflection (trans = 0, k_ref = mirror of k_in);
!   - Fresnel transmittance/reflectance at the local incidence angle;
!   - Smith (1967) bidirectional shadowing.
! The main loop applies weight splitting at the interface:
!   - ww * refl continues the random walk in water with direction k_ref
!     (reflected facet direction, NOT a flat-surface mirror reflection);
!   - ww * trans leaves the system and is not scored: the probability
!     of reaching the receiver after the last collision is already in
!     the local estimate of that collision (Eq. 23). Scoring it again
!     would count the same event twice (B6, removed in v0.5).
! The value of this interface treatment is that the reflected fraction
! is not discarded: it continues the random walk and keeps contributing
! via the in-water local estimate.
!-----------------------------------------------------------------------

!+ Point where a photon step crosses the surface z = 0
!
Pure Subroutine find_surface_crossing &
  ( x0,  & ! in
    y0,  & ! in
    z0,  & ! in
    x1,  & ! in
    y1,  & ! in
    z1,  & ! in
    x_E, & ! out
    y_E )  ! out

! Description:
!   Finds the intersection of the photon path with z = 0. The photon
!   moved from (x0, y0, z0) with z0 < 0 to (x1, y1, z1) with z1 >= 0.
!
!   Inputs:
!     x0, y0, z0 : position before the step (z0 < 0)
!     x1, y1, z1 : position after the step  (z1 >= 0)
!   Outputs:
!     x_E, y_E   : surface crossing point coordinates at z = 0
!
! Method:
!   Linear interpolation along the straight step:
!     frac = -z0/(z1 - z0) in [0, 1],
!     x_E = x0 + frac (x1 - x0),  y_E = y0 + frac (y1 - y0).
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.3     2026-11-07  real(8) instead of kind dp. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Real(8), Intent(In)  :: x0    ! Start x [m].
Real(8), Intent(In)  :: y0    ! Start y [m].
Real(8), Intent(In)  :: z0    ! Start z [m], < 0.
Real(8), Intent(In)  :: x1    ! End x [m].
Real(8), Intent(In)  :: y1    ! End y [m].
Real(8), Intent(In)  :: z1    ! End z [m], >= 0.

! Scalar arguments with intent(out):
Real(8), Intent(Out) :: x_E   ! Crossing x [m].
Real(8), Intent(Out) :: y_E   ! Crossing y [m].

! Local scalars:
Real(8) :: frac               ! Fraction of the step up to z = 0.

!- End of header ------------------------------------------------------

frac = (-z0) / (z1 - z0)
x_E = x0 + frac * (x1 - x0)
y_E = y0 + frac * (y1 - y0)

End Subroutine find_surface_crossing

!-----------------------------------------------------------------------
! 4. Signal estimators
!-----------------------------------------------------------------------

!+ Local estimate of the surface glint (Nakajima and Tanaka 1983)
!
Function surface_glint_signal_le &
  ( u_in,     & ! in
    x_E,      & ! in
    y_E,      & ! in
    x_sat,    & ! in
    H_sat,    & ! in
    sigma2,   & ! in
    N_LUT,    & ! in
    LUT_cos,  & ! in
    LUT_R_aa, & ! in
    tau0,     & ! in
    alpha,    & ! in
    Area )    & ! in
  Result (signal_glint)

! Description:
!   Semi-analytic (local estimate) glint signal. Instead of checking
!   whether the randomly sampled k_ref falls inside the FOV cone (analog
!   approach, statistically useless for FOV <= 1 mrad), this function
!   computes the EXPECTED glint contribution at every surface hit inside
!   the FOV ellipse.
!
!   Inputs:
!     u_in(3)     : incident photon direction in air (u_in(3) < 0)
!     x_E, y_E    : surface hit point [m]
!     x_sat       : x of the telescope centre [m], at altitude H_sat:
!                   x_A + (H' - H) tan(alpha). NOT x_A, which is the
!                   foot of the virtual apex F at altitude H'
!                   (Shangguan et al. 2024, Eq. 12)
!     H_sat       : satellite altitude [m]     (= hlidar in MC_water)
!     sigma2      : total mean-square slope (Cox and Munk 1954)
!     N_LUT       : LUT size
!     LUT_cos(:)  : LUT cosine grid
!     LUT_R_aa(:) : LUT Fresnel air-air reflectance
!     tau0        : zenith one-way atmospheric optical depth
!     alpha       : off-nadir angle [rad] (used for T_atm)
!     Area        : receiver aperture area [m^2] = pi r_tel^2
!   Output:
!     glint contribution of the photon [same units as the weight]
!
! Method:
!   1. Receiver direction u_recv from the hit point to the satellite.
!   2. Required facet normal that specularly reflects u_in toward u_recv:
!        k_ref = u_in - 2 (u_in.n) n = u_recv  =>  n ~ u_recv - u_in
!      (n(3) > 0; u_in(3) < 0 and u_recv(3) > 0, so
!      |u_recv - u_in| >= 2 cos(alpha) >> 0: never degenerate).
!   3. Cox-Munk density of the SLOPES at s = (-n_x/n_z, -n_y/n_z):
!        P = exp(-(s_x^2 + s_y^2)/sigma2) / (pi sigma2)
!      (Nakajima and Tanaka 1983, p. 524, formula above Eq. 14).
!   4. Fresnel R at the local incidence angle (from the LUT).
!   5. Bidirectional shadowing G = 1/(1 + Lambda_i + Lambda_r).
!   Formula:
!     E = R_F P / (4 cos(theta_0) cos^4(theta_n)) G (A/R_sat^2) T_atm,
!     cos(theta_0) = -u_in(3), T_atm = exp(-2 tau0/cos(alpha)).
!   Derivation; theta_i and theta_n are LOCAL angles on the required
!   facet, not the off-nadir alpha:
!     a) probability that the photon hits a facet with slope in d2s,
!        facet area seen by the beam over the beam section:
!          P cos(theta_i) / (cos(theta_n) cos(theta_0)) d2s;
!     b) slopes to normal solid angle: d2s = dOmega_n / cos^3(theta_n),
!        i.e. p(mu_n) = P / mu_n^3 (Nakajima and Tanaka 1983, Eq. 14);
!     c) normal to reflected direction (the reflected ray turns twice
!        as fast as the normal): dOmega_n = dOmega_r / (4 cos(theta_i));
!     d) dOmega_r = A / R_sat^2; cos(theta_i) cancels.
!   Same as Nakajima and Tanaka (1983), Eq. 20, with p = P / mu_n^3,
!   turned from reflection function into fraction per photon.
!   R_sat = |r_sat - r_E|, the geometric distance from the hit point to
!   the satellite: the glint is reflected in air and never crosses the
!   interface, so no refraction and d(Omega_r) = A / R_sat^2 (unlike
!   R_eff = (n_water H + z)/cos(alpha) of local_estimate_signal, where
!   the ray refracts at the surface).
!   Flat sea (sigma2 = 0): the pdf diverges; the glint is non-zero ONLY
!   if the required slope is (numerically) zero, i.e. exactly specular
!   geometry; then E = R_F (A/R_sat^2) T_atm.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.3     2026-11-07  n_water as argument (was a local 1.33).
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
! v0.7     2026-11-11  R_sat = geometric distance to the satellite,
!                      without n_water (B5a); n_water argument removed.
!                      Valerio Ciattaglia
! v0.7     2026-11-11  Denominator 4 cos(theta_0) cos^4(theta_n) (B5b).
!                      Valerio Ciattaglia
! v0.7     2026-11-11  x_sat is the telescope centre, no longer x_A
!                      (B5c). Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In) :: x_E              ! Hit point x [m].
Real(8), Intent(In) :: y_E              ! Hit point y [m].
Real(8), Intent(In) :: x_sat            ! Telescope centre x [m].
Real(8), Intent(In) :: H_sat            ! Satellite altitude [m].
Real(8), Intent(In) :: sigma2           ! Total mean-square slope.
Integer, Intent(In) :: N_LUT            ! Number of LUT entries.
Real(8), Intent(In) :: tau0             ! One-way optical depth.
Real(8), Intent(In) :: alpha            ! Off-nadir angle [rad].
Real(8), Intent(In) :: Area             ! Receiver area [m^2].

! Array  arguments with intent(in):
Real(8), Intent(In) :: u_in(3)          ! Incident versor.
Real(8), Intent(In) :: LUT_cos(N_LUT)   ! LUT cos grid.
Real(8), Intent(In) :: LUT_R_aa(N_LUT)  ! Reflectance, air side.

! Function result:
Real(8) :: signal_glint                 ! Glint contribution.

! Local scalars:
Real(8) :: u_recv_norm                  ! Norm before normalisation.
Real(8) :: R_sat                        ! Slant range [m].
Real(8) :: n_req_norm                   ! Norm before normalisation.
Real(8) :: cos_theta_i                  ! cos(local incidence).
Real(8) :: cos_theta_n                  ! cos(facet tilt).
Real(8) :: cos_theta_0                  ! cos(incident zenith).
Real(8) :: slope_x                      ! Required slope along x.
Real(8) :: slope_y                      ! Required slope along y.
Real(8) :: pi                           ! Pi.
Real(8) :: R_fresnel                    ! Fresnel R (LUT).
Real(8) :: p_cox                        ! Cox-Munk slope pdf.
Real(8) :: lambda_i                     ! Shadowing, incident leg.
Real(8) :: lambda_r                     ! Shadowing, reflected leg.
Real(8) :: G_shadow                     ! Shadowing factor.
Real(8) :: cos_alpha                    ! cos(alpha), bounded from 0.

! Local arrays:
Real(8) :: u_recv(3)                    ! Versor toward the satellite.
Real(8) :: n_req(3)                     ! Required facet normal.

!- End of header ------------------------------------------------------

signal_glint = 0.0d0
cos_alpha    = max(cos(alpha), 1.0d-6)
pi = acos(-1.0d0)

!-----------------------------------------------------------------------
! 1. Exact receiver direction from the hit point to the satellite
!-----------------------------------------------------------------------
u_recv(1) = x_sat - x_E
u_recv(2) =       - y_E
u_recv(3) = H_sat
u_recv_norm = sqrt(u_recv(1)**2 + u_recv(2)**2 + u_recv(3)**2)
If (u_recv_norm < 1.0d-12) Return
R_sat  = u_recv_norm              ! Geometric distance E-satellite.
u_recv = u_recv / u_recv_norm

!-----------------------------------------------------------------------
! 2. Required facet normal n_req = normalize(u_recv - u_in)
!-----------------------------------------------------------------------
n_req(1) = u_recv(1) - u_in(1)
n_req(2) = u_recv(2) - u_in(2)
n_req(3) = u_recv(3) - u_in(3)
n_req_norm = sqrt(n_req(1)**2 + n_req(2)**2 + n_req(3)**2)
If (n_req_norm < 1.0d-12) Return
n_req = n_req / n_req_norm
If (n_req(3) <= 0.0d0) Return   ! Facet must point upward.

!-----------------------------------------------------------------------
! 3. Local incidence angle, facet tilt and required slope
!-----------------------------------------------------------------------
cos_theta_i = -(u_in(1)*n_req(1) + u_in(2)*n_req(2) + u_in(3)*n_req(3))
If (cos_theta_i <= 0.0d0) Return
cos_theta_i = min(1.0d0, cos_theta_i)
cos_theta_n = n_req(3)
If (cos_theta_n < 1.0d-6) Return
cos_theta_0 = -u_in(3)                  ! Zenith cosine of the incident ray.
If (cos_theta_0 <= 0.0d0) Return
slope_x = -n_req(1) / cos_theta_n
slope_y = -n_req(2) / cos_theta_n

!-----------------------------------------------------------------------
! 4. Fresnel R at the local incidence angle
!-----------------------------------------------------------------------
Call get_lut_nearest(cos_theta_i, LUT_cos, LUT_R_aa, R_fresnel)

!-----------------------------------------------------------------------
! 5. Flat sea (sigma2 -> 0): mirror-like, deterministic reflection
!-----------------------------------------------------------------------
If (sigma2 < 1.0d-12) Then
  If (abs(slope_x) > 1.0d-6 .or. abs(slope_y) > 1.0d-6) Return
  signal_glint = R_fresnel * (Area / (R_sat * R_sat))                  &
                 * exp(-2.0d0 * tau0 / cos_alpha)
  Return
End If

!-----------------------------------------------------------------------
! 6. Rough sea: Cox-Munk pdf, shadowing, local estimate
!-----------------------------------------------------------------------
p_cox = exp(-(slope_x**2 + slope_y**2) / sigma2) / (pi * sigma2)
! Zenith cosines of the incident (-u_in(3) > 0) and receiver rays.
lambda_i = lambda_saunders(-u_in(3),   sigma2)
lambda_r = lambda_saunders( u_recv(3), sigma2)
G_shadow = 1.0d0 / (1.0d0 + lambda_i + lambda_r)
signal_glint = R_fresnel * p_cox                                       &
               / (4.0d0 * cos_theta_0 * cos_theta_n**4)                &
               * G_shadow                                              &
               * (Area / (R_sat * R_sat))                              &
               * exp(-2.0d0 * tau0 / cos_alpha)

End Function surface_glint_signal_le

!+ Local estimate of the signal at a subsurface collision
!
Pure Function local_estimate_signal &
  ( ww,           & ! in
    u,            & ! in
    u_wat_recv,   & ! in
    z_pos,        & ! in
    hlidar,       & ! in
    alpha,        & ! in
    alpha_r_axis, & ! in
    n_water,      & ! in
    kext,         & ! in
    Area,         & ! in
    T_t,          & ! in
    tau0,         & ! in
    g )           & ! in
  Result (signal)

! Description:
!   Local estimate of the lidar signal from a collision at depth z_pos,
!   Shangguan et al. (2024), Eq. 23. The photon weight ww is NOT
!   modified: weight management stays in the main loop.
!
!   Inputs:
!     ww            : photon weight W_N (after the absorption update)
!     u(3)          : photon direction BEFORE this scatter event (unit)
!     u_wat_recv(3) : direction from the scatter point toward the
!                     receiver in water (unit), precomputed once:
!                     (-sin(alpha_r), 0, cos(alpha_r))
!     z_pos         : depth below the surface (positive, = abs(z)) [m]
!     hlidar        : orbital altitude [m]
!     alpha         : off-nadir viewing angle in air [rad]
!     alpha_r_axis  : refracted optical-axis angle in water [rad]
!     n_water       : water refractive index
!     kext          : beam attenuation coefficient c [1/m]
!     Area          : receiver collecting area A_r [m^2]
!     T_t           : interface transmittance T_wa at alpha_r (LUT)
!     tau0          : zenith one-way atmospheric optical depth
!     g             : Henyey-Greenstein asymmetry factor
!   Output:
!     E(z) signal contribution [same units as W_N]
!
! Method:
!   E(z) = W_N (A_r / R_eff^2) (beta(theta)/4 pi) exp(-c d') T_a T_t
!   R_eff = (n_w H + z) / cos(alpha)   effective slant range [m]
!   d'    = z / cos(alpha_r)           water path to the surface [m]
!   T_a   = exp(-2 tau0 / cos(alpha))  round-trip atmospheric
!           transmittance (down * up); ww does NOT pre-carry T_a_down.
!   T_t   = Fresnel T_wa at alpha_r    interface transmittance
!           (flat-surface approximation)
!   theta = angle between u and u_wat_recv (scattering angle for HG):
!   beta/4 pi = (1 - g^2) / (4 pi (1 + g^2 - 2 g cos(theta))^(3/2)).
!   Coordinate convention: z positive upward; photons in water have
!   z < 0; z_pos = abs(z) is always positive.
!   Ref.: Shangguan, Guo, Liao and Lee (2024), RSE, Eq. 23.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-07-15  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Pure; header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Real(8), Intent(In) :: ww              ! Photon weight.
Real(8), Intent(In) :: z_pos           ! Depth, positive [m].
Real(8), Intent(In) :: hlidar          ! Orbital altitude [m].
Real(8), Intent(In) :: alpha           ! Off-nadir angle in air [rad].
Real(8), Intent(In) :: alpha_r_axis    ! Refracted axis angle [rad].
Real(8), Intent(In) :: n_water         ! Index of seawater.
Real(8), Intent(In) :: kext            ! Attenuation coefficient [1/m].
Real(8), Intent(In) :: Area            ! Receiver area [m^2].
Real(8), Intent(In) :: T_t             ! Interface transmittance.
Real(8), Intent(In) :: tau0            ! One-way optical depth.
Real(8), Intent(In) :: g               ! Asymmetry factor.

! Array  arguments with intent(in):
Real(8), Intent(In) :: u(3)            ! Direction before the collision.
Real(8), Intent(In) :: u_wat_recv(3)   ! Direction toward the receiver.

! Function result:
Real(8) :: signal                      ! Local estimate E(z).

! Local parameters:
Real(8), Parameter :: pi_loc = 3.141592653589793d0   ! Pi.

! Local scalars:
Real(8) :: cos_theta                   ! cos(scattering angle).
Real(8) :: R_eff                       ! Effective slant range [m].
Real(8) :: beta_over_4pi               ! HG phase function/4 pi [1/sr].
Real(8) :: d_prime                     ! Return path in water [m].

!- End of header ------------------------------------------------------

! Scattering angle theta between u and u_wat_recv.
cos_theta = u(1)*u_wat_recv(1) + u(2)*u_wat_recv(2) + u(3)*u_wat_recv(3)
If (cos_theta >  1.0d0) cos_theta =  1.0d0
If (cos_theta < -1.0d0) cos_theta = -1.0d0

! beta(theta)/4 pi, HG [1/sr].
beta_over_4pi = ((1.0d0 - g*g) / (4.0d0 * pi_loc))                     &
                * (1.0d0 + g*g - 2.0d0*g*cos_theta)**(-1.5d0)
! Effective slant range and water path to the surface.
R_eff   = (n_water * hlidar + z_pos) / cos(alpha)
d_prime = z_pos / cos(alpha_r_axis)

! Shangguan et al. (2024), Eq. 23.
signal = ww * (Area / (R_eff * R_eff))                                 &
         * beta_over_4pi                                               &
         * exp(-kext * d_prime)                                        &
         * exp(-2.0d0 * tau0 / cos(alpha))                             &
         * T_t

End Function local_estimate_signal

!-----------------------------------------------------------------------
! 5. Input and output
!-----------------------------------------------------------------------

!+ Reads, checks and completes the namelist input
!
Subroutine read_input &
  ( version )   ! in

! Description:
!   Reads the namelist input file (first command-line argument, or
!   input.nml), checks that every variable is assigned and in range, and
!   sets the derived quantities. Stops with a message on any error.
!
! Method:
!   1. File name from get_command_argument; existence check.
!   2. Completeness: check_input_complete (no default values).
!   3. The six groups are read with a rewind before each one, so their
!      order in the file is free; check_group stops on a missing or
!      malformed group.
!   4. Domain checks with require; all violations are listed.
!   5. Derived quantities:
!        D_l = 2 r_tel,  sigma_g = beam_div/4 (intensity at 1/e^2 for
!        theta = 2 sigma_g, see generate_photon),  ksca = omega kext,
!        sigma2 = 0.003 + 0.00512 u10 (Cox and Munk 1954, clean surface,
!        total mean-square slope; 0 if flat_sea),
!        output grid from the receiver sampling interval dt_bin: a time
!        t after the surface echo is a one-way path in water
!        s_p = c t / (2 n_w), projected on the vertical as
!        z_t = s_p cos(alpha_r) (Shangguan et al. 2025, Eq. 24 and 27),
!        with cos(alpha_r) = sqrt(1 - (sin(theta_v)/n_w)^2) (Snell), so
!          dz = c dt_bin cos(alpha_r) / (2 n_w),
!          nz = n_bins,  zmax = nz dz.
!        The real-depth profile uses the same grid, to compare the two
!        profiles bin by bin; it is not a grid of the medium.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, input section of MC.f90.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module; derived quantities here;
!                      nangle removed. Valerio Ciattaglia
! v0.4     2026-11-07  beam_div; sigma_g = beam_div/4 (was fov/4);
!                      header to the standard. Valerio Ciattaglia
! v0.6     2026-11-08  nslope removed. Valerio Ciattaglia
! v0.8     2026-11-11  Group &receiver (dt_bin, n_bins); dz and zmax
!                      derived from dt_bin. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Character(len=*), Intent(In) :: version   ! Code version (git describe).

! Local scalars:
Integer :: unit_in                        ! Unit of the input file.
Integer :: io_stat                        ! I/O status.
Character(len=256) :: io_msg              ! I/O error message.
Logical :: file_exists                    ! .true. if the file exists.
Logical :: input_ok                       ! .false. if a check fails.
Real(8) :: pi                             ! Pi.
Real(8) :: cos_ar                         ! cos of the refracted axis.

!- End of header ------------------------------------------------------

code_version = version

!-----------------------------------------------------------------------
! 1. File name and existence
!-----------------------------------------------------------------------
If (command_argument_count() >= 1) Then
  Call get_command_argument(1, input_file)
Else
  input_file = 'input.nml'
End If
Inquire(file=trim(input_file), exist=file_exists)
If (.not. file_exists) Then
  Write(*,'(A)') 'ERROR: input file not found: '//trim(input_file)
  Write(*,'(A)') 'Usage: mc_water [input_file]   (default: input.nml)'
  Stop 1
End If
Write(*,'(A)') 'Input file: '//trim(input_file)

!-----------------------------------------------------------------------
! 2. Every variable must be assigned (no defaults)
!-----------------------------------------------------------------------
Call check_input_complete(trim(input_file),                            &
     [Character(len=24) ::                                             &
     'theta_v', 'fov', 'beam_div', 'hlidar', 'r_tel', 'sampling_mode', &
     'tau0',                                                           &
     'n_water', 'u10', 'flat_sea', 'internal_reflection',              &
     'kext', 'omega', 'g',                                             &
     'nphotons', 'seed', 'wmin', 'N_LUT', 'maxsca',                    &
     'dt_bin', 'n_bins'],                                              &
     input_ok)
If (.not. input_ok) Stop 1

!-----------------------------------------------------------------------
! 3. Namelist groups
!-----------------------------------------------------------------------
Open(newunit=unit_in, file=trim(input_file), status='old',             &
     action='read')
Rewind(unit_in)
Read(unit_in, nml=geometry, iostat=io_stat, iomsg=io_msg)
Call check_group('geometry', io_stat, io_msg)
Rewind(unit_in)
Read(unit_in, nml=atmosphere, iostat=io_stat, iomsg=io_msg)
Call check_group('atmosphere', io_stat, io_msg)
Rewind(unit_in)
Read(unit_in, nml=surface, iostat=io_stat, iomsg=io_msg)
Call check_group('surface', io_stat, io_msg)
Rewind(unit_in)
Read(unit_in, nml=water, iostat=io_stat, iomsg=io_msg)
Call check_group('water', io_stat, io_msg)
Rewind(unit_in)
Read(unit_in, nml=numerics, iostat=io_stat, iomsg=io_msg)
Call check_group('numerics', io_stat, io_msg)
Rewind(unit_in)
Read(unit_in, nml=receiver, iostat=io_stat, iomsg=io_msg)
Call check_group('receiver', io_stat, io_msg)
Close(unit_in)

!-----------------------------------------------------------------------
! 4. Domain checks
!-----------------------------------------------------------------------
input_ok = .true.
Call require(nphotons > 0, 'nphotons > 0', input_ok)
Call require(kext > 0.0d0, 'kext > 0', input_ok)
Call require(omega > 0.0d0 .and. omega < 1.0d0, '0 < omega < 1',       &
             input_ok)
Call require(g > -1.0d0 .and. g < 1.0d0, '-1 < g < 1', input_ok)
Call require(theta_v >= 0.0d0 .and. theta_v < 90.0d0,                  &
             '0 <= theta_v < 90', input_ok)
Call require(fov > 0.0d0, 'fov > 0', input_ok)
Call require(beam_div > 0.0d0, 'beam_div > 0', input_ok)
Call require(hlidar > 0.0d0, 'hlidar > 0', input_ok)
Call require(r_tel > 0.0d0, 'r_tel > 0', input_ok)
Call require(trim(sampling_mode) == 'gaussian' .or.                    &
             trim(sampling_mode) == 'uniform_cone',                    &
             "sampling_mode = 'gaussian' or 'uniform_cone'", input_ok)
Call require(n_water > 1.0d0, 'n_water > 1', input_ok)
Call require(u10 >= 0.0d0, 'u10 >= 0', input_ok)
Call require(tau0 >= 0.0d0, 'tau0 >= 0', input_ok)
Call require(wmin > 0.0d0 .and. wmin < 1.0d0, '0 < wmin < 1', input_ok)
Call require(dt_bin > 0.0d0, 'dt_bin > 0', input_ok)
Call require(n_bins >= 1, 'n_bins >= 1', input_ok)
Call require(N_LUT >= 2, 'N_LUT >= 2', input_ok)
Call require(maxsca >= 1, 'maxsca >= 1', input_ok)
Call require(seed >= 0, 'seed >= 0', input_ok)
If (.not. input_ok) Stop 1

!-----------------------------------------------------------------------
! 5. Derived quantities
!-----------------------------------------------------------------------
D_l     = 2.0d0*r_tel
sigma_g = beam_div/4.0d0
ksca    = omega*kext
sigma2  = 0.003d0 + 0.00512d0*u10
If (flat_sea) sigma2 = 0.0d0
! Output grid from the receiver sampling interval; the argument of sqrt
! is positive because n_water > 1 is checked above.
pi     = acos(-1.0d0)
cos_ar = sqrt(1.0d0 - (sin(theta_v*pi/180.0d0)/n_water)**2)
dz     = c_light*dt_bin*1.0d-9*cos_ar/(2.0d0*n_water)
nz     = n_bins
zmax   = real(nz, 8)*dz

End Subroutine read_input

!+ Checks that every required variable is assigned in the input file
!
Subroutine check_input_complete &
  ( fname, & ! in
    names, & ! in
    ok )     ! out

! Description:
!   Verifies that every required variable is assigned in the namelist
!   input file. Missing variables are listed on standard output.
!
!   Inputs:
!     fname : input file name
!     names : required variable names
!   Output:
!     ok    : .true. if all variables are assigned
!
! Method:
!   The file is read as plain text: comments (from '!') are removed, the
!   text is converted to lower case and blanks are dropped
!   (compact_lower); a variable is assigned if 'name=' appears preceded
!   by a non-name character (is_name_char) or at the start of the text.
!   Works for every type, logicals included, which have no "unassigned"
!   value.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code. Valerio Ciattaglia
! v0.4     2026-11-07  Internal functions moved to the module; header to
!                      the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Character(len=*), Intent(In)  :: fname      ! Input file name.

! Array  arguments with intent(in):
Character(len=*), Intent(In)  :: names(:)   ! Required variable names.

! Scalar arguments with intent(out):
Logical,          Intent(Out) :: ok         ! .true. if all assigned.

! Local scalars:
Character(len=:), Allocatable :: text       ! Whole file, cleaned.
Character(len=1024) :: line                 ! One line of the file.
Character(len=64)   :: key                  ! Lower-case 'name='.
Integer :: unit_f                           ! File unit.
Integer :: ios                              ! I/O status.
Integer :: k                                ! Name index.
Integer :: pos                              ! Position of a match.
Integer :: start                            ! Search start.
Integer :: icom                             ! Position of '!'.
Integer :: alloc_stat                       ! Allocation status.

! Local arrays:
Logical, Allocatable :: found(:)            ! Assignment found.

!- End of header ------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Read and clean the file
!-----------------------------------------------------------------------
text = ' '
Open(newunit=unit_f, file=fname, status='old', action='read')
Do
  Read(unit_f, '(A)', iostat=ios) line
  If (ios /= 0) Exit
  icom = index(line, '!')
  If (icom > 0) line(icom:) = ' '
  text = text // compact_lower(line) // ' '
End Do
Close(unit_f)

!-----------------------------------------------------------------------
! 2. Search each name
!-----------------------------------------------------------------------
Allocate(found(size(names)), stat=alloc_stat)
If (alloc_stat /= 0) Then
  Write(*,'(A)') 'ERROR: allocation failed in check_input_complete.'
  Stop 1
End If
found = .false.
Do k = 1, size(names)
  key = compact_lower(names(k)) // '='
  start = 1
  Do
    pos = index(text(start:), trim(key))
    If (pos == 0) Exit
    pos = start + pos - 1
    If (pos == 1) Then
      found(k) = .true.
    Else If (.not. is_name_char(text(pos-1:pos-1))) Then
      found(k) = .true.
    End If
    If (found(k)) Exit
    start = pos + 1
  End Do
End Do

!-----------------------------------------------------------------------
! 3. Report
!-----------------------------------------------------------------------
ok = all(found)
Do k = 1, size(names)
  If (.not. found(k)) Write(*,'(A)') 'ERROR: variable not assigned in ' &
                                     //trim(fname)//': '//trim(names(k))
End Do
Deallocate(found)

End Subroutine check_input_complete

!+ Lower-case copy of a string without blanks
!
Pure Function compact_lower &
  ( str )  & ! in
  Result (res)

! Description:
!   Returns str in lower case with blanks and tabs removed.
!
! Method:
!   Character codes A-Z are shifted by 32; blanks and tabs (code 9) are
!   skipped.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, internal to check_input_complete.
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Module procedure; header to the standard.
!                      Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Character(len=*), Intent(In)  :: str   ! Input string.

! Function result:
Character(len=:), Allocatable :: res   ! Compacted lower-case string.

! Local scalars:
Integer :: j                           ! Character index.
Integer :: ic                          ! Character code.

!- End of header ------------------------------------------------------

res = ''
Do j = 1, len_trim(str)
  ic = iachar(str(j:j))
  If (str(j:j) == ' ' .or. ic == 9) Cycle
  If (ic >= iachar('A') .and. ic <= iachar('Z')) ic = ic + 32
  res = res // achar(ic)
End Do

End Function compact_lower

!+ True for characters allowed in a Fortran name
!
Pure Function is_name_char &
  ( c )  & ! in
  Result (is_name)

! Description:
!   .true. for lower-case letters, digits and underscore.
!
! Method:
!   Character comparison; the text is already in lower case.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, internal to check_input_complete.
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Module procedure; header to the standard.
!                      Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Character(len=1), Intent(In) :: c   ! Character to test.

! Function result:
Logical :: is_name                  ! .true. for a name character.

!- End of header ------------------------------------------------------

is_name = (c >= 'a' .and. c <= 'z') .or. (c >= '0' .and. c <= '9')    &
          .or. c == '_'

End Function is_name_char

!+ Stops if a namelist group is missing or malformed
!
Subroutine check_group &
  ( group,   & ! in
    io_stat, & ! in
    io_msg )   ! in

! Description:
!   Stops the program with a message if the namelist group just read is
!   missing (io_stat < 0) or malformed (io_stat > 0).
!
! Method:
!   Test of the iostat of the namelist read; the message includes the
!   group name and the iomsg of the compiler.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, internal procedure of MC.f90.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Character(len=*), Intent(In) :: group     ! Name of the group read.
Integer,          Intent(In) :: io_stat   ! I/O status of the read.
Character(len=*), Intent(In) :: io_msg    ! I/O error message.

!- End of header ------------------------------------------------------

If (io_stat < 0) Then
  Write(*,'(A)') 'ERROR: namelist group &'//group//' not found in '     &
                 //trim(input_file)
  Stop 1
Else If (io_stat > 0) Then
  Write(*,'(A)') 'ERROR: cannot read namelist group &'//group//': '     &
                 //trim(io_msg)
  Stop 1
End If

End Subroutine check_group

!+ Reports an out-of-range input value
!
Subroutine require &
  ( condition, & ! in
    rule,      & ! in
    ok )         ! inout

! Description:
!   Prints the admissible range of an input value if the condition is
!   false, and sets ok to .false.; the caller stops after all checks, so
!   every violation is listed.
!
! Method:
!   Logical test of the condition passed by the caller.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, internal procedure of MC.f90.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module; ok as argument.
!                      Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Logical,          Intent(In)    :: condition   ! .true. if admissible.
Character(len=*), Intent(In)    :: rule        ! Admissible range.

! Scalar arguments with intent(inout):
Logical,          Intent(InOut) :: ok          ! .false. if violated.

!- End of header ------------------------------------------------------

If (.not. condition) Then
  Write(*,'(A)') 'ERROR: input value out of range, required: '//rule
  ok = .false.
End If

End Subroutine require

!+ Seeds the random number generator
!
Subroutine set_random_seed()

! Description:
!   Seeds the intrinsic generator from the input seed. seed = 0 is
!   replaced by a value drawn from the system clock, so that every run
!   is repeatable: the seed actually used is printed, written in the
!   headers of the output files and in input_used.nml.
!
! Method:
!   seed_arr(i) = seed + 37 (i - 1), i = 1..n, with n from
!   random_seed(size=n); random_seed(put=seed_arr).
!   seed = 0: seed = mod(clock count, 1e9) + 1.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, in MC.f90. Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Local scalars:
Integer    :: n_seed                   ! Size of the seed array.
Integer    :: i_seed                   ! Seed-array index.
Integer(8) :: clock_count              ! System clock count.
Integer    :: alloc_stat               ! Allocation status.

! Local arrays:
Integer, Allocatable :: seed_arr(:)    ! Seed array.

!- End of header ------------------------------------------------------

If (seed == 0) Then
  Call system_clock(clock_count)
  seed = int(mod(clock_count, 1000000000_8), kind=4) + 1
End If
Call random_seed(size=n_seed)
Allocate(seed_arr(n_seed), stat=alloc_stat)
If (alloc_stat /= 0) Then
  Write(*,'(A)') 'ERROR: allocation failed in set_random_seed.'
  Stop 1
End If
Do i_seed = 1, n_seed
  seed_arr(i_seed) = seed + 37*(i_seed - 1)
End Do
Call random_seed(put=seed_arr)
Deallocate(seed_arr)
Write(*,'(A,I0)') 'Random seed: ', seed

End Subroutine set_random_seed

!+ Writes the input actually used to input_used.nml
!
Subroutine write_input_used()

! Description:
!   Copy of the input actually used, with the actual random seed, in
!   input_used.nml in the run directory; the first line gives the code
!   version and the input file.
!
! Method:
!   Namelist write of the six groups.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, in MC.f90. Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
! v0.8     2026-11-11  Group &receiver. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Local scalars:
Integer :: unit_used   ! Unit of input_used.nml.

!- End of header ------------------------------------------------------

Open(newunit=unit_used, file='input_used.nml', status='replace',       &
     action='write')
Write(unit_used,'(A)') '! Input actually used by mc_water '            &
                       //trim(code_version)//', read from '            &
                       //trim(input_file)//'.'
Write(unit_used, nml=geometry)
Write(unit_used, nml=atmosphere)
Write(unit_used, nml=surface)
Write(unit_used, nml=water)
Write(unit_used, nml=numerics)
Write(unit_used, nml=receiver)
Close(unit_used)

End Subroutine write_input_used

!+ Writes the signal profiles and the collision counts
!
Subroutine write_output &
  ( photons,   & ! in
    photons_s, & ! in
    n_in,      & ! in
    n_out )      ! in

! Description:
!   Writes SIGNAL_PROFILE_Z.DAT (real depth, Shangguan et al. 2025,
!   Eq. 25), SIGNAL_PROFILE_S.DAT (time-based depth
!   z_t = s_p cos(alpha_r), Eq. 24 and 27) and SCATTERING_COUNTS.DAT
!   (collisions inside and outside the FOV), each with the run
!   parameters as '# key = value' header lines (write_header).
!
! Method:
!   Each signal value is divided by nphotons to give the mean
!   contribution per incident photon: the result is invariant under
!   changes to nphotons. Subsurface bins (iz >= 1) are also divided by
!   dz: signal per metre of depth, invariant under changes to dz. Bin 0
!   (surface: glint) is not.
!   Column layout: col 1 = bin index, col 2 = bin-centre depth
!   (iz - 0.5) dz [m], then the signal for scattering orders ns = 0
!   (all orders) through ns = maxsca. Counts: in(ns = 0..maxsca), then
!   out(ns = 0..maxsca), raw numbers of collisions.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.0     2026-09-22  Original code, output section of MC.f90.
!                      Valerio Ciattaglia
! v0.1     2026-11-05  Signal per metre; depth column.
!                      Valerio Ciattaglia
! v0.2     2026-11-05  Format from maxsca; SCATTERING_COUNTS.DAT.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
! v0.5     2026-11-08  Bin 0 holds the glint only. Valerio Ciattaglia
! v0.8     2026-11-11  dz derived from dt_bin; dt_bin and n_bins in the
!                      headers (write_header). Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Array  arguments with intent(in):
Real(8),    Intent(In) :: photons(0:maxsca, 0:nz)     ! Signal, real z.
Real(8),    Intent(In) :: photons_s(0:maxsca, 0:nz)   ! Signal, z_t.
Integer(8), Intent(In) :: n_in(0:maxsca, 0:nz)        ! Inside the FOV.
Integer(8), Intent(In) :: n_out(0:maxsca, 0:nz)       ! Outside the FOV.

! Local parameters:
Character(len=*), Parameter :: col_fmt = '(i6,f10.3,*(1pe11.3))'
                                         ! Format of the signal rows.
Character(len=*), Parameter :: cnt_fmt = '(i6,f10.3,*(i12))'
                                         ! Format of the count rows.

! Local scalars:
Integer :: unit_z     ! Unit of SIGNAL_PROFILE_Z.DAT.
Integer :: unit_s     ! Unit of SIGNAL_PROFILE_S.DAT.
Integer :: unit_c     ! Unit of SCATTERING_COUNTS.DAT.
Integer :: iz         ! Bin index.
Integer :: ns         ! Scattering order.

!- End of header ------------------------------------------------------

Open(newunit=unit_z, file='SIGNAL_PROFILE_Z.DAT',  status='replace')
Open(newunit=unit_s, file='SIGNAL_PROFILE_S.DAT',  status='replace')
Open(newunit=unit_c, file='SCATTERING_COUNTS.DAT', status='replace')

!-----------------------------------------------------------------------
! 1. Headers
!-----------------------------------------------------------------------
Call write_header(unit_z)
Write(unit_z,'(A,I0,A,F0.4,A)')                                        &
  '# iz=0: surface (glint, ns=0), not divided by dz;'                  &
  //' iz=1..', nz, ': subsurface bins of ', dz, ' m, signal per metre'
Write(unit_z,'(A)')                                                    &
  '# Depth axis: real depth z (Shangguan et al. 2025, Eq. 25)'
Write(unit_z,'(A)') '# col1=iz col2=z_centre[m] col3=ns0(total) '      &
                    //order_columns(4, 'ns')

Call write_header(unit_s)
Write(unit_s,'(A,I0,A,F0.4,A)')                                        &
  '# is=0: surface (glint, ns=0), not divided by dz;'                  &
  //' is=1..', nz, ': subsurface bins of ', dz, ' m, signal per metre'
Write(unit_s,'(A)')                                                    &
  '# Depth axis: time-based depth z_t = s_p*cos(alpha_r)'              &
  //' (Shangguan et al. 2025, Eq. 24 and 27)'
Write(unit_s,'(A)') '# col1=is col2=z_t_centre[m] col3=ns0(total) '    &
                    //order_columns(4, 'ns')

Call write_header(unit_c)
Write(unit_c,'(A,I0,A,F0.4,A)')                                        &
  '# Number of collisions (each counts 1, not divided by nphotons)'    &
  //' by real-depth bin; iz=1..', nz, ', bins of ', dz,                &
  ' m; iz=0 unused'
Write(unit_c,'(A)')                                                    &
  '# in = collision inside the FOV (scored), out = outside the FOV'    &
  //' (not scored, photon continues)'
Write(unit_c,'(A,I0,A,I0,A)')                                          &
  '# col1=iz col2=z_centre[m] col3..col', 3+maxsca,                    &
  '=in(ns=0..maxsca) col', 4+maxsca,                                   &
  '..=out(ns=0..maxsca); ns=0 is the total'

!-----------------------------------------------------------------------
! 2. Surface bin (iz = 0), not divided by dz
!-----------------------------------------------------------------------
Write(unit_z,col_fmt) 0, 0.0d0,                                        &
  (photons(ns,0)/real(nphotons, 8), ns=0,maxsca)
Write(unit_s,col_fmt) 0, 0.0d0,                                        &
  (photons_s(ns,0)/real(nphotons, 8), ns=0,maxsca)
Write(unit_c,cnt_fmt) 0, 0.0d0,                                        &
  (n_in(ns,0), ns=0,maxsca), (n_out(ns,0), ns=0,maxsca)

!-----------------------------------------------------------------------
! 3. Subsurface bins (iz >= 1), signal per metre
!-----------------------------------------------------------------------
Do iz = 1, nz
  Write(unit_z,col_fmt) iz, (real(iz, 8) - 0.5d0)*dz,                  &
    (photons(ns,iz)/(real(nphotons, 8)*dz), ns=0,maxsca)
  Write(unit_s,col_fmt) iz, (real(iz, 8) - 0.5d0)*dz,                  &
    (photons_s(ns,iz)/(real(nphotons, 8)*dz), ns=0,maxsca)
  Write(unit_c,cnt_fmt) iz, (real(iz, 8) - 0.5d0)*dz,                  &
    (n_in(ns,iz), ns=0,maxsca), (n_out(ns,iz), ns=0,maxsca)
End Do

Close(unit_z)
Close(unit_s)
Close(unit_c)

End Subroutine write_output

!+ Writes the run parameters as header lines of an output file
!
Subroutine write_header &
  ( unit )   ! in

! Description:
!   Writes the code version, the input file and every run parameter as
!   '# key = value' lines, read back by the plot scripts.
!
! Method:
!   One formatted line per parameter; alpha in degrees is written as
!   (theta_v pi/180)(180/pi), as in v0.2, so that the text is unchanged.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, internal procedure of MC.f90.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module; nangle removed.
!                      Valerio Ciattaglia
! v0.4     2026-11-07  beam_div line; header to the standard.
!                      Valerio Ciattaglia
! v0.6     2026-11-08  nslope removed. Valerio Ciattaglia
! v0.8     2026-11-11  dt_bin and n_bins lines; dz with 7 decimals.
!                      Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Subroutine arguments
! Scalar arguments with intent(in):
Integer, Intent(In) :: unit   ! Output unit.

! Local scalars:
Real(8) :: pi                 ! Pi.

!- End of header ------------------------------------------------------

pi = acos(-1.0d0)
Write(unit,'(A,A)')      '# code     = ', trim(code_version)
Write(unit,'(A,A)')      '# input    = ', trim(input_file)
Write(unit,'(A,I12)')    '# nphotons = ', nphotons
Write(unit,'(A,I12)')    '# seed     = ', seed
Write(unit,'(A,F8.3)')   '# alpha    = ',                              &
                         (theta_v * pi / 180.0d0) * (180.0d0 / pi)
Write(unit,'(A,ES12.5)') '# FOV      = ', fov
Write(unit,'(A,ES12.5)') '# beam_div = ', beam_div
Write(unit,'(A,ES12.5)') '# hlidar   = ', hlidar
Write(unit,'(A,F8.4)')   '# r_tel    = ', r_tel
Write(unit,'(A,A)')      '# sampling = ', trim(sampling_mode)
Write(unit,'(A,F8.4)')   '# tau0     = ', tau0
Write(unit,'(A,F8.4)')   '# n_water  = ', n_water
Write(unit,'(A,F8.3)')   '# U10      = ', u10
Write(unit,'(A,L1)')     '# flat_sea = ', flat_sea
Write(unit,'(A,F10.6)')  '# sigma2   = ', sigma2
Write(unit,'(A,L1)')     '# int_refl = ', internal_reflection
Write(unit,'(A,F8.5)')   '# kext     = ', kext
Write(unit,'(A,F8.4)')   '# omega    = ', omega
Write(unit,'(A,F8.4)')   '# g        = ', g
Write(unit,'(A,ES12.5)') '# wmin     = ', wmin
Write(unit,'(A,F10.4)')  '# dt_bin   = ', dt_bin
Write(unit,'(A,I8)')     '# n_bins   = ', n_bins
Write(unit,'(A,F10.3)')  '# zmax     = ', zmax
Write(unit,'(A,F12.7)')  '# dz       = ', dz
Write(unit,'(A,I8)')     '# N_LUT    = ', N_LUT
Write(unit,'(A,I8)')     '# maxsca   = ', maxsca

End Subroutine write_header

!+ Header labels of the output columns of orders 1..maxsca
!
Function order_columns &
  ( first_col, & ! in
    prefix )   & ! in
  Result (labels)

! Description:
!   Builds the header labels 'colN=prefixK' of the columns of scattering
!   orders K = 1..maxsca, starting from column first_col.
!
! Method:
!   Internal write of each label; labels separated by one blank.
!
! Current Code Owner: Valerio Ciattaglia
!
! History:
! Version  Date        Comment
! -------  ----        -------
! v0.2     2026-11-05  Original code, internal procedure of MC.f90.
!                      Valerio Ciattaglia
! v0.3     2026-11-07  Moved to the module. Valerio Ciattaglia
! v0.4     2026-11-07  Header to the standard. Valerio Ciattaglia
!
! Code Description:
!   Language:           Fortran 2008 (gfortran).
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Implicit None

! Function arguments
! Scalar arguments with intent(in):
Integer,          Intent(In)  :: first_col   ! Column of order 1.
Character(len=*), Intent(In)  :: prefix      ! Label prefix.

! Function result:
Character(len=:), Allocatable :: labels      ! Labels, blank-separated.

! Local scalars:
Character(len=32) :: one                     ! Label of one column.
Integer :: k                                 ! Scattering order.

!- End of header ------------------------------------------------------

labels = ''
Do k = 1, maxsca
  Write(one,'(A,I0,A,A,I0)') 'col', first_col + k - 1, '=', prefix, k
  labels = labels//trim(one)//' '
End Do

End Function order_columns

End Module mc_module

!- End of module header
