!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
      SUBROUTINE icarus(                      &
!cms reduced number of continuation lines
     &     debug, debugcol, npoints, sunlit,  &
     &     nlev, ncol, pfull, phalf,qv,         & 
     &     cc, dtau_s, top_height,  &
     &     top_height_direction, overlap,   &
     &     frac_out, skt, emsfc_lw,   &
     &     at,  &
     &     dem_s,  &
     &     fq_isccp, &
     &     totalcldarea, &
     &     meanptop, &
     &     meantaucld, &
     &     meanalbedocld, &
     &     meantb, &
     &     meantbclr, &
     &     boxtau, &
     &     boxptop, &
! BMD 07/26 - start
     &     dtau_s_liq, dtau_s_ice, reffl, reffi, &
     &     fq_modis_lwp_reffl, fq_modis_iwp_reffi, &
     &     cltmodis, clwmodis, climodis, &
     &     fq_modis_tau_ctp, fq_modis_tau_ctp_liq, fq_modis_tau_ctp_ice, &
     &     fq_modis_tau_reffl, fq_modis_tau_reffi, &
! BMD 07/26 - end
! BMD 09/26 - start
     &     modis_diag_liq_attempt, modis_diag_liq_fail, &
     &     modis_diag_liq_obshi, modis_diag_liq_obslo, &
     &     modis_diag_ice_attempt, modis_diag_ice_fail, &
     &     modis_diag_ice_obshi, modis_diag_ice_obslo &
! BMD 09/26 - end
     &)

     USE mo_kind,  ONLY: wp

!$Id: icarus.f,v 4.1 2010/05/27 16:30:18 hadmw Exp $

! *****************************COPYRIGHT****************************
! (c) 2009, Lawrence Livermore National Security Limited Liability 
! Corporation.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Lawrence Livermore National Security 
!       Limited Liability Corporation nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************

      IMPLICIT NONE

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER :: ncolprint

! BMD 07/26 - start
!     Bin-table parameters for two new MODIS-style joint histograms:
!     liquid water path (LWP) vs. liquid effective radius (ReffLiq), and
!     ice water path (IWP) vs. ice effective radius (ReffIce). These use
!     the Stephens (1978) bulk approximation CWP = (2/3)*rho*re*tau, the
!     same relation COSP's own modis_simulator.F90 uses, rather than the
!     brightness-temperature/cloud-top-pressure machinery below that is
!     specific to the ISCCP tau/CTP histogram. Bin edges mirror the ones
!     already adopted for the ICON MODIS-simulator port (cosp_config.F90,
!     BMD 07/2026) so the two implementations stay directly comparable.
!
!     Wired end-to-end: mo_cosp_modis_simulator.f90 is the dedicated MODIS
!     call site (passing this whole trailing argument group), while
!     mo_cosp_isccp_simulator.f90's pre-existing ISCCP-only call site omits
!     it -- safe because these 9 dummy arguments (plus the tau/CTP group
!     added in a second pass, further below) are all OPTIONAL and every
!     code block that touches them is gated on PRESENT(...).

      INTEGER, PARAMETER :: n_modis_lwp = 7, n_modis_iwp = 7 ! number of MODIS LWP, IWP bins
      INTEGER, PARAMETER :: n_modis_reffl = 8, n_modis_reffi = 8 ! number of MODIS effective radius (liquid,ice) bins

      REAL(wp), PARAMETER :: modis_lwp_bin_edges(n_modis_lwp+1) = &
     &  (/ 0._wp, 0.01_wp, 0.03_wp, 0.06_wp, 0.10_wp, 0.15_wp, 0.25_wp, 20._wp /)     ! kg m-2
      REAL(wp), PARAMETER :: modis_iwp_bin_edges(n_modis_iwp+1) = &
     &  (/ 0._wp, 0.02_wp, 0.05_wp, 0.10_wp, 0.20_wp, 0.40_wp, 1.00_wp, 20._wp /)     ! kg m-2
      REAL(wp), PARAMETER :: modis_reffl_bin_edges(n_modis_reffl+1) = &
     &  (/ 0._wp, 4.e-6_wp, 8.e-6_wp, 1.0e-5_wp, 1.25e-5_wp, 1.5e-5_wp, &
     &     2.0e-5_wp, 3.0e-5_wp, 1.0e-2_wp /)                                        ! m
      REAL(wp), PARAMETER :: modis_reffi_bin_edges(n_modis_reffi+1) = &
     &  (/ 0._wp, 5.e-6_wp, 1.0e-5_wp, 2.0e-5_wp, 3.0e-5_wp, 4.0e-5_wp, &
     &     5.0e-5_wp, 6.0e-5_wp, 1.0e-2_wp /)                                        ! m

      REAL(wp), PARAMETER :: modis_lwp_conversion = (2._wp/3._wp)*1000._wp   ! rho_liquid [kg/m3]
      REAL(wp), PARAMETER :: modis_ice_density    = 0.93_wp                  ! relative to liquid

!     Phase is decided once per subcolumn (not independently per phase): a subcolumn is liquid-dominated,
!     ice-dominated, or excluded from both
!     histograms.
      REAL(wp), PARAMETER :: modis_phase_threshold = 0.7_wp

!     How deep into the cloud (in optical depth, from cloud top) the phase
!     determination looks: matches real COSP's phase_TauLimit
!     (cosp_modis_interface.F90's COSP_MODIS_INIT in ~/COSPv2.0). Passive
!     VIS/NIR retrievals only "see" the topmost cloud optical depth, so
!     liquid_fraction below is deliberately NOT a whole-subcolumn average --
!     see the tau_phase_cum/tauf_phase_cum accumulation.
      REAL(wp), PARAMETER :: modis_phase_taulimit = 1._wp

!     Minimum detectable optical thickness for the whole MODIS block below
!     (cloud-fraction outputs AND the LWP/IWP-vs-Reff histograms): matches
!     real COSP's min_OpticalThickness (cosp_modis_interface.F90's
!     COSP_MODIS_INIT in ~/COSPv2.0), not the tiny tauchk epsilon used by the
!     ISCCP tau/CTP histogram above. A subcolumn with tau(j,ibox) below this
!     is not MODIS-visible at all, regardless of phase.
      REAL(wp), PARAMETER :: modis_min_optical_thickness = 0.3_wp

!     MODIS tau/cloud-top-pressure joint histogram: bin counts and edges
!     are NOT redeclared here -- COSPv2's numMODISTauBins/modis_histTau and
!     numMODISPresBins/modis_histPres (cosp_config.F90, ~/COSPv2.0) are
!     IDENTICAL to this file's own ISCCP tau/CTP bin edges (isccp_taumin/
!     itau(j) and the 180/310/440/560/680/800 mb breakpoints/ipres(j)
!     below), confirmed against ~/COSPv2.0/src/cosp.F90 (numMODISTauBins=
!     ntau, modis_histTau=tau_binBounds) -- so itau(j) is reused directly
!     as the MODIS tau bin index (see the fq_modis_tau_ctp block below);
!     only the pressure axis needs its own binning (imodis_pres), since
!     MODIS's own CTP retrieval differs from ISCCP's ptop (see below).
      INTEGER, PARAMETER :: n_modis_tau = 7, n_modis_pres = 7

!     CO2-slicing cloud-top-pressure retrieval parameters, matching real
!     COSP's CO2Slicing_TauLimit and CO2Slicing_PressureLimit
!     (cosp_modis_interface.F90's COSP_MODIS_INIT in ~/COSPv2.0).
!     CO2Slicing_PressureLimit is 70000 Pa upstream; kept here in mb (700)
!     to match this file's existing ptop/boxptop mb convention, so the
!     CO2-sliced retrieval and the ISCCP-fallback CTP can be compared and
!     binned without a unit conversion at every use site.
      REAL(wp), PARAMETER :: modis_co2slicing_taulimit = 1._wp
      REAL(wp), PARAMETER :: modis_co2slicing_preslimit_mb = 700._wp

!     Simulated-reflectance size retrieval, ported from COSPv2's
!     modis_subcolumn/get_g_nir/get_ssa_nir (modis_simulator.F90) via
!     ~/COSPv2.0 -- see the CONTAINS block at the end of this file for the
!     ported procedures themselves (fit_to_quadratic/fit_to_cubic,
!     get_g_nir/get_ssa_nir, two_stream/two_stream_reflectance,
!     adding_doubling, interpolate_to_min). Fully wired into the main body:
!     the "observed" reflectance is synthesized via two_stream/
!     adding_doubling over a per-level near-IR single-scattering albedo/
!     asymmetry parameter computed directly from the model's real reffl/
!     reffi through get_g_nir/get_ssa_nir (ports COSPv2's modis_optics,
!     cosp_optics.f90, tau-weighted liquid/ice combination -- confirmed
!     against CAM6's own COSP2 MODIS simulator, which sources its
!     equivalent sizeLIQ/sizeICE the same way), then inverted against
!     these trial-size libraries via interpolate_to_min, replacing the
!     earlier box_reffl/box_reffi model-truth shortcut. (Earlier passes of
!     this port instead used a near-IR piz/cg pair proxied from RRTMG-SW
!     band 16, which at 3.08-3.85 micron is much wider than the true
!     MODIS band 20 channel, 3.66-3.84 micron -- superseded now that the
!     narrowband get_g_nir/get_ssa_nir fits are evaluated directly.)
!
!     Trial-size grid bounds, in MICRONS (matching COSPv2's own unit
!     convention for these specific constants -- NOT the meters convention
!     the rest of this file's reffl/reffi/bin-edge code uses). Only the
!     final retrieved size, after interpolate_to_min, gets converted to
!     meters (multiply by 1.e-6_wp), exactly as COSPv2 does at its own
!     retrievedSize(i) = 1.0e-06_wp*interpolate_to_min(...) call site.
      INTEGER, PARAMETER :: num_trial_res = 15
      REAL(wp), PARAMETER :: modis_trial_re_water_min = 4._wp   ! micron
      REAL(wp), PARAMETER :: modis_trial_re_water_max = 30._wp  ! micron
      REAL(wp), PARAMETER :: modis_trial_re_ice_min = 5._wp     ! micron
      REAL(wp), PARAMETER :: modis_trial_re_ice_max = 90._wp    ! micron

!     Missing-value sentinel for a failed size retrieval (COSPv2's
!     re_fill), in microns -- matches the trial-grid unit convention above,
!     not the meters convention used once the final result is converted.
      REAL(wp), PARAMETER :: modis_re_fill = -999._wp           ! micron
! BMD 07/26 - end

!     -----
!     Input 
!     -----

      INTEGER :: npoints       !  number of model points in the horizontal
      INTEGER :: nlev          !  number of model levels in column
      INTEGER :: ncol          !  number of subcolumns

      INTEGER :: sunlit(npoints) !  1 for day points, 0 for night time

      REAL(wp) :: pfull(npoints,nlev)
                       !  pressure of full model levels (Pascals)
                  !  pfull(npoints,1) is top level of model
                  !  pfull(npoints,nlev) is bot of model

      REAL(wp) :: phalf(npoints,nlev+1)
                  !  pressure of half model levels (Pascals)
                  !  phalf(npoints,1) is top of model
                  !  phalf(npoints,nlev+1) is the surface pressure

      REAL(wp) :: qv(npoints,nlev)
                  !  water vapor specific humidity (kg vapor/ kg air)
                  !         on full model levels

      REAL(wp) :: cc(npoints,nlev)   
                  !  input cloud cover in each model level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds


      REAL(wp) :: dtau_s(npoints,nlev) 
                  !  mean 0.67 micron optical depth of stratiform
                !  clouds in each model level
                  !  NOTE:  this the cloud optical depth of only the
                  !  cloudy part of the grid box, it is not weighted
                  !  with the 0 cloud optical depth of the clear
                  !         part of the grid box

! BMD 07/26 - start
!     OPTIONAL: callers that don't want the MODIS-style histograms below
!     (e.g. mo_cosp_isccp_simulator.f90's ISCCP-only call site) can omit
!     this whole trailing group of 4 arguments. Every code block below that
!     touches any of them is gated on PRESENT(...), so it's illegal/unsafe
!     to reference one of these without checking PRESENT first.
      REAL(wp), OPTIONAL :: dtau_s_liq(npoints,nlev)
                  !  liquid-only component of dtau_s (0.67 micron optical
                  !  depth of the cloudy part of the grid box), used to
                  !  build the MODIS-style LWP-vs-ReffLiq joint histogram
      REAL(wp), OPTIONAL :: dtau_s_ice(npoints,nlev)
                  !  ice-only component of dtau_s, used to build the
                  !  MODIS-style IWP-vs-ReffIce joint histogram
      REAL(wp), OPTIONAL :: reffl(npoints,nlev)
                  !  liquid cloud droplet effective radius (m) on full
                  !  model levels
      REAL(wp), OPTIONAL :: reffi(npoints,nlev)
                  !  ice crystal effective radius (m) on full model levels
! BMD 07/26 - end

! BMD 09/26 - start
!     Diagnostic counters (fraction of sunlit subcolumns, boxarea-weighted,
!     same convention as cltmodis/clwmodis/climodis) added to test whether
!     clwmodis's near-collapse is caused by the near-IR size retrieval
!     chronically failing to bracket the observed reflectance for
!     liquid-classified subcolumns specifically (as opposed to, e.g., phase
!     misclassification biasing subcolumns toward the ice branch instead).
!     *_attempt: subcolumns where a liquid/ice-library retrieval was
!     actually attempted (obs_refl_active .and. PRESENT(dtau_s_liq)) --
!     the denominator for a fail-rate diagnostic.
!     *_fail: subset of *_attempt where interpolate_to_min could not
!     bracket a root (box_reffl/box_reffi came back modis_re_fill).
!     *_obshi/*_obslo: subset of *_fail split by which side of the trial
!     curve's achievable reflectance range the "observed" reflectance
!     (obs_refl_run, from the real per-level near-IR optics computed from
!     reffl/reffi via get_g_nir/get_ssa_nir) fell outside of: obshi =
!     exceeded the max reflectance any trial size can produce, obslo =
!     fell below the min. This distinguishes "retrieval saturates because
!     the observed near-IR optics don't match the phase-specific Platnick
!     trial library" from other possible causes.
!     OPTIONAL like the rest of the MODIS group: omitted by
!     mo_cosp_isccp_simulator.f90's ISCCP-only call site.
      REAL(wp), OPTIONAL :: modis_diag_liq_attempt(npoints)
      REAL(wp), OPTIONAL :: modis_diag_liq_fail(npoints)
      REAL(wp), OPTIONAL :: modis_diag_liq_obshi(npoints)
      REAL(wp), OPTIONAL :: modis_diag_liq_obslo(npoints)
      REAL(wp), OPTIONAL :: modis_diag_ice_attempt(npoints)
      REAL(wp), OPTIONAL :: modis_diag_ice_fail(npoints)
      REAL(wp), OPTIONAL :: modis_diag_ice_obshi(npoints)
      REAL(wp), OPTIONAL :: modis_diag_ice_obslo(npoints)
! BMD 09/26 - end

      INTEGER :: overlap                   !  overlap type
                              !  1=max
                              !  2=rand
                              !  3=max/rand

      INTEGER :: top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
                              !  optical depth to adjust cloud top pressure. Note
                              !  that this calculation is most appropriate to compare
                              !  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
                              !  3 = adjust top height using only the computed
                              !  infrared brightness temperature. Note that this
                              !  calculation is most appropriate to compare to ISCCP
                              !  IR only algortihm (i.e. you can compare to nighttime
                              !  ISCCP data with this option)

      INTEGER :: top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 !
                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 !
                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                 ! with interpolated temperature equal to the radiance 
                                 ! determined cloud-top temperature
                                 ! 
                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                 !                                 !
                                 ! 1 = old setting: matches all versions of 
                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                 !
                                 ! 2 = default setting: for version numbers 4.0 and higher
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL(wp) :: skt(npoints)                 !  skin Temperature (K)
      REAL(wp) :: emsfc_lw                     !  10.5 micron emissivity of surface (fraction)
      REAL(wp) :: at(npoints,nlev)             !  temperature in each model level (K)
      REAL(wp) :: dem_s(npoints,nlev)          !  10.5 micron longwave emissivity of stratiform
                                               !  clouds in each
                                               !  model level.  Same note applies as in dtau_s.

      REAL(wp) :: frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column



!     ------
!     Output
!     ------

      REAL(wp) :: fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL(wp) :: totalcldarea(npoints)        !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  NOTE: This diagnostic
                                        ! does not count model clouds with tau < isccp_taumin
                              ! Thus this diagnostic does not equal the sum over all entries of 
                              !  fq_isccp. However, this diagnostic does equal the sum over entries
                              !   of fq_isccp with itau = 2:7 (omitting itau = 1)
      
      
     ! The following three means are averages only over the cloudy areas with tau > isccp_taumin.  
     ! If no clouds with tau > isccp_taumin are in grid box all three quantities should equal zero.      
                              
      REAL(wp) :: meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
                              
      REAL(wp) :: meantaucld(npoints)          !  mean optical thickness 
                                        !  linear averaging in albedo performed.
      
      REAL(wp) :: meanalbedocld(npoints)        ! mean cloud albedo
                                        ! linear averaging in albedo performed
                                        
      REAL(wp) :: meantb(npoints)              ! mean all-sky 10.5 micron brightness temperature
      
      REAL(wp) :: meantbclr(npoints)           ! mean clear-sky 10.5 micron brightness temperature
      
      REAL(wp) :: boxtau(npoints,ncol)         !  optical thickness in each column
      
      REAL(wp) :: boxptop(npoints,ncol)        !  cloud top pressure (mb) in each column

! BMD 07/26 - start
!     OPTIONAL: see the note by dtau_s_liq above -- omit this whole group
!     together when the caller doesn't want MODIS output.
      REAL(wp), OPTIONAL :: fq_modis_lwp_reffl(npoints,n_modis_lwp,n_modis_reffl)
                  !  joint histogram of liquid water path vs. liquid
                  !  effective radius, fraction of model grid box
      REAL(wp), OPTIONAL :: fq_modis_iwp_reffi(npoints,n_modis_iwp,n_modis_reffi)
                  !  joint histogram of ice water path vs. ice
                  !  effective radius, fraction of model grid box
      REAL(wp), OPTIONAL :: cltmodis(npoints)
                  !  MODIS total cloud fraction: fraction of sunlit
                  !  subcolumns with tau >= modis_min_optical_thickness,
                  !  regardless of phase determination outcome (so this is
                  !  NOT simply clwmodis+climodis -- undetermined-phase
                  !  subcolumns count here but not in either of those)
      REAL(wp), OPTIONAL :: clwmodis(npoints)
                  !  MODIS liquid cloud fraction: fraction of sunlit
                  !  subcolumns classified liquid-dominated
      REAL(wp), OPTIONAL :: climodis(npoints)
                  !  MODIS ice cloud fraction: fraction of sunlit
                  !  subcolumns classified ice-dominated

!     OPTIONAL, own group: fq_modis_tau_ctp/_liq/_ice don't require
!     dtau_s_liq/dtau_s_ice/reffl/reffi at all (only dtau_s/frac_out/phalf,
!     all non-optional), but are gated the same OPTIONAL/PRESENT way as the
!     rest of the MODIS block since mo_cosp_modis_simulator.f90 -- the only
!     caller wanting MODIS output -- always supplies the whole group
!     together.
      REAL(wp), OPTIONAL :: fq_modis_tau_ctp(npoints,n_modis_tau,n_modis_pres)
                  !  joint histogram of optical thickness vs. cloud-top
                  !  pressure (all phases), fraction of model grid box
      REAL(wp), OPTIONAL :: fq_modis_tau_ctp_liq(npoints,n_modis_tau,n_modis_pres)
                  !  same, restricted to liquid-dominated subcolumns
      REAL(wp), OPTIONAL :: fq_modis_tau_ctp_ice(npoints,n_modis_tau,n_modis_pres)
                  !  same, restricted to ice-dominated subcolumns

!     OPTIONAL, own group: MODIS-style optical-thickness-vs-particle-size
!     joint histograms, ported from COSPv2's modis_column
!     (Optical_Thickness_vs_ReffLiq/_vs_ReffIce, modis_simulator.F90) via
!     ~/COSPv2.0 -- a different axis pairing than fq_modis_lwp_reffl/
!     fq_modis_iwp_reffi above (LWP/IWP vs Reff): here the FULL subcolumn
!     tau(j,ibox) (same itau(j) bin already used for fq_modis_tau_ctp) is
!     crossed with the retrieved size instead of the water/ice path.
!     Populated from the same box_reffl/box_reffi retrieval outcome as the
!     LWP/IWP-vs-Reff histograms -- no separate retrieval attempt needed.
      REAL(wp), OPTIONAL :: fq_modis_tau_reffl(npoints,n_modis_tau,n_modis_reffl)
                  !  joint histogram of optical thickness vs. liquid
                  !  effective radius, fraction of model grid box
      REAL(wp), OPTIONAL :: fq_modis_tau_reffi(npoints,n_modis_tau,n_modis_reffi)
                  !  joint histogram of optical thickness vs. ice
                  !  effective radius, fraction of model grid box
! BMD 07/26 - end


!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL(wp) :: dem(npoints,ncol),bb(npoints)     !  working variables for 10.5 micron longwave 
                              !  emissivity in part of
                              !  gridbox under consideration

      REAL(wp) :: ptrop(npoints)
      REAL(wp) :: attrop(npoints)
      REAL(wp) :: attropmin (npoints)
      REAL(wp) :: atmax(npoints)
      REAL(wp) :: btcmin(npoints)
      REAL(wp) :: transmax(npoints)

      INTEGER :: i,j,ilev,ibox,itrop(npoints)
      INTEGER :: ipres(npoints)
      INTEGER :: itau(npoints),ilev2
      INTEGER :: acc(nlev,ncol)
      INTEGER :: match(npoints,nlev-1)
      INTEGER :: nmatch(npoints)
      INTEGER :: levmatch(npoints,ncol)
      
      !variables needed for water vapor continuum absorption
      REAL(wp) :: fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      REAL(wp) :: taumin(npoints)
      REAL(wp) :: dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
      REAL(wp) :: press(npoints), dpress(npoints), atmden(npoints)
      REAL(wp) :: rvh20(npoints), wk(npoints), rhoave(npoints)
      REAL(wp) :: rh20s(npoints), rfrgn(npoints)
      REAL(wp) :: tmpexp(npoints),tauwv(npoints)
      
      CHARACTER*1 :: cchar(6),cchar_realtops(6)
      INTEGER :: icycle
      REAL(wp) :: tau(npoints,ncol)
      LOGICAL :: box_cloudy(npoints,ncol)
      REAL(wp) :: tb(npoints,ncol)
      REAL(wp) :: ptop(npoints,ncol)
      REAL(wp) :: emcld(npoints,ncol)
      REAL(wp) :: fluxtop(npoints,ncol)
      REAL(wp) :: trans_layers_above(npoints,ncol)
      REAL(wp) :: isccp_taumin,fluxtopinit(npoints),tauir(npoints)
      REAL(wp) :: albedocld(npoints,ncol)
      REAL(wp) :: boxarea

! BMD 07/26 - start
      REAL(wp) :: box_reffl(npoints), box_reffi(npoints)
                  !  retrieved (not model-truth) liquid/ice effective
                  !  radius (m) for the current subcolumn -- see the
                  !  simulated-reflectance retrieval block below. Negative
                  !  (modis_re_fill*1.e-6_wp) means the retrieval failed to
                  !  bracket a trial size; guarded before binning.
      REAL(wp) :: box_lwp(npoints), box_iwp(npoints)
      REAL(wp) :: liquid_fraction(npoints)
                  !  tauf_phase_cum/tau_phase_cum for the current subcolumn --
                  !  the single per-subcolumn phase-dominance fraction used
                  !  to decide liquid- vs. ice-dominated (see
                  !  modis_phase_threshold above and tau_phase_cum below)
      INTEGER  :: imodis_lwp(npoints), imodis_reffl(npoints)
      INTEGER  :: imodis_iwp(npoints), imodis_reffi(npoints)
      LOGICAL  :: retrieval_valid
                  !  BMD 08/26: whether THIS subcolumn's single size-
                  !  retrieval attempt (liquid library for liquid-dominated
                  !  and undetermined-phase subcolumns, ice library for
                  !  ice-dominated ones) succeeded -- gates cltmodis/
                  !  fq_modis_tau_ctp/fq_modis_tau_ctp_liq/_ice the same way
                  !  COSPv2's validRetrievalMask gates its "total" outputs.
                  !  See the detailed note at the retrieval-attempt block
                  !  below.

!     Simulated-reflectance MODIS size retrieval (Stage 4/5 of the port,
!     BMD 07/26 fourth pass, near-IR-optics sourcing updated BMD 09/26):
!     synthesizes the "observed" near-IR reflectance per subcolumn from the
!     model's real per-level optical properties (near-IR single-scattering
!     albedo/asymmetry computed from reffl/reffi via get_g_nir/get_ssa_nir,
!     see nir_w0_combined/nir_g_combined below, via two_stream/
!     adding_doubling), then inverts it against an independent trial-size
!     library (get_g_nir/get_ssa_nir/two_stream_reflectance/
!     interpolate_to_min, all ported in the CONTAINS block) -- replaces the
!     previous retau_liq/tau_liq model-truth shortcut for box_reffl/
!     box_reffi.
!
!     Reflectance-synthesis running state, one running (Refl,Tran) pair
!     per subcolumn, updated top-down over cloudy levels via the same
!     adding-layers recursion as adding_doubling (called with npts=2 each
!     time: the running state as "layer 1", this level as "layer 2") --
!     obs_refl_active tracks whether any cloudy level has been folded in
!     yet, since the very first contributing level initializes rather
!     than combines (mirrors adding_doubling's own Refl_cumulative(1) =
!     Refl(1) special case).
      LOGICAL  :: obs_refl_active(npoints,ncol)
      REAL(wp) :: obs_refl_run(npoints,ncol), obs_tran_run(npoints,ncol)
      REAL(wp) :: refl_layer, tran_layer, refl_combined, tran_combined
                  !  per-(j,ibox,ilev) scratch scalars for the reflectance
                  !  synthesis above

!     Per-level near-IR (true MODIS band 20, 3.66-3.84 micron) single-
!     scattering albedo/asymmetry parameter, computed from the model's
!     real reffl/reffi via get_g_nir/get_ssa_nir and tau-weighted between
!     liquid and ice -- ports COSPv2's modis_optics (cosp_optics.f90) via
!     ~/COSPv2.0, confirmed against CAM6's own COSP2 MODIS simulator
!     (cospsimulator_intr.f90's MODIS_waterSize/MODIS_iceSize -> the same
!     modis_optics), which sources its equivalent near-IR optics the same
!     way. Supersedes the RRTMG-SW-band-16-sourced dpiz_s_nir/dcg_s_nir
!     (BMD 09/26) -- see the PRESENT(dtau_s_liq) block below.
      REAL(wp) :: nir_w0_liq(npoints,nlev), nir_g_liq(npoints,nlev)
      REAL(wp) :: nir_w0_ice(npoints,nlev), nir_g_ice(npoints,nlev)
      REAL(wp) :: nir_w0_combined(npoints,nlev), nir_g_combined(npoints,nlev)

!     Trial-size grid and its near-IR optical-property library (Platnick
!     Collection-6 fits, get_g_nir/get_ssa_nir), computed once per icarus
!     call (not per subcolumn -- these don't depend on j/ibox) and reused
!     for every liquid/ice-dominated subcolumn's retrieval below. Units
!     are microns, matching modis_trial_re_water/ice_min/max's convention.
      REAL(wp) :: trial_re_w(num_trial_res), trial_re_i(num_trial_res)
      REAL(wp) :: trial_g_w(num_trial_res), trial_w0_w(num_trial_res)
      REAL(wp) :: trial_g_i(num_trial_res), trial_w0_i(num_trial_res)
      REAL(wp) :: predicted_refl_w(num_trial_res), predicted_refl_i(num_trial_res)
                  !  per-subcolumn predicted-reflectance-vs-trial-size
                  !  curve, recomputed each liquid-/ice-dominated
                  !  subcolumn since it depends on that subcolumn's own
                  !  total tau(j,ibox)

      REAL(wp) :: tau_phase_cum(npoints,ncol), tauf_phase_cum(npoints,ncol)
                  !  running (cloud-top-down) cumulative total optical depth
                  !  and optical-depth-weighted liquid-fraction product used
                  !  to build liquid_fraction above, ports COSPv2's
                  !  weight_by_extinction (modis_simulator.F90) restricted to
                  !  the top modis_phase_taulimit of cloud optical depth --
                  !  a passive VIS/NIR sensor only sees the cloud-top phase,
                  !  not a full-column average
      LOGICAL  :: phase_active(npoints,ncol)
                  !  .true. until a subcolumn's cumulative tau reaches
                  !  modis_phase_taulimit; freezes tau_phase_cum/
                  !  tauf_phase_cum thereafter (the vectorized-over-npoints
                  !  equivalent of weight_by_extinction's scalar "exit")
      REAL(wp) :: phase_tau_incr, phase_f_liq, phase_deltatau
                  !  per-(j,ibox,ilev) scratch scalars for the
                  !  tau_phase_cum/tauf_phase_cum accumulation: this level's
                  !  total (liquid+ice) tau increment, this level's liquid
                  !  fraction, and (on the level where the cumulative tau
                  !  crosses modis_phase_taulimit) the partial increment
                  !  actually used, mirroring weight_by_extinction's deltaX

      REAL(wp) :: ctp_tau_cum(npoints,ncol), ctp_pres_cum(npoints,ncol)
                  !  running (cloud-top-down) cumulative total optical depth
                  !  and optical-depth-weighted pressure product, ports
                  !  COSPv2's cloud_top_pressure (modis_simulator.F90)
                  !  restricted to the top modis_co2slicing_taulimit of
                  !  cloud optical depth -- the CO2-slicing MODIS
                  !  cloud-top-pressure retrieval, distinct from ISCCP's
                  !  brightness-temperature-based ptop(j,ibox)
      LOGICAL  :: ctp_active(npoints,ncol)
                  !  .true. until a subcolumn's cumulative tau reaches
                  !  modis_co2slicing_taulimit; freezes ctp_tau_cum/
                  !  ctp_pres_cum thereafter, same role as phase_active above
      REAL(wp) :: ctp_deltatau
                  !  per-(j,ibox,ilev) scratch scalar: on the level where
                  !  the cumulative tau crosses modis_co2slicing_taulimit,
                  !  the partial tau increment actually used (mirrors
                  !  phase_deltatau/weight_by_extinction's deltaX, but for
                  !  cloud_top_pressure's trapezoidal-rule pressure integral)
      REAL(wp) :: modis_ctp_j
                  !  current subcolumn's blended MODIS cloud-top pressure
                  !  (mb): the CO2-sliced ctp_pres_cum/ctp_tau_cum result,
                  !  or ptop(j,ibox) (ISCCP's CTP) when the CO2-sliced value
                  !  falls below modis_co2slicing_preslimit_mb (low cloud) --
                  !  ports COSPv2's modis_subcolumn low-cloud fallback
      INTEGER  :: imodis_pres(npoints)
                  !  MODIS CTP bin index for the current subcolumn, binned
                  !  from modis_ctp_j using the same 180/310/440/560/680/800
                  !  mb edges as ipres(j) above (see n_modis_pres note) --
                  !  NOT the same value as ipres(j), since modis_ctp_j is a
                  !  different retrieval than ISCCP's ptop(j,ibox)
! BMD 07/26 - end

      INTEGER :: debug       ! set to non-zero value to print out inputs
                    ! with step debug
      INTEGER :: debugcol    ! set to non-zero value to print out column
                    ! decomposition with step debugcol
      INTEGER :: rangevec(npoints),rangeerror

      INTEGER :: k1,k2
      REAL(wp) :: rec2p13,tauchk,logp,logp1,logp2,atd
      REAL(wp) :: output_missing_value

      CHARACTER*10 :: ftn09
      
      DATA isccp_taumin / 0.3_wp /
      DATA output_missing_value / -1.E+30 /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

!     ------ End duplicate definitions common to wrapper routine

      tauchk = -1._wp*log(0.9999999_wp)
      rec2p13=1._wp/2.13_wp

      ncolprint=0

      IF ( debug.ne.0 ) THEN
          j=1
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'debug='
          write(6,'(8I10)') debug
          write(6,'(a10)') 'debugcol='
          write(6,'(8I10)') debugcol
          write(6,'(a10)') 'npoints='
          write(6,'(8I10)') npoints
          write(6,'(a10)') 'nlev='
          write(6,'(8I10)') nlev
          write(6,'(a10)') 'ncol='
          write(6,'(8I10)') ncol
          write(6,'(a11)') 'top_height='
          write(6,'(8I10)') top_height
          write(6,'(a21)') 'top_height_direction='
          write(6,'(8I10)') top_height_direction
          write(6,'(a10)') 'overlap='
          write(6,'(8I10)') overlap
          write(6,'(a10)') 'emsfc_lw='
          write(6,'(8f10.2)') emsfc_lw
        DO j=1,npoints,debug
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'sunlit='
          write(6,'(8I10)') sunlit(j)
          write(6,'(a10)') 'pfull='
          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
          write(6,'(a10)') 'phalf='
          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
          write(6,'(a10)') 'qv='
          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
          write(6,'(a10)') 'cc='
          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
          write(6,'(a10)') 'dtau_s='
          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
          write(6,'(a10)') 'skt='
          write(6,'(8f10.2)') skt(j)
          write(6,'(a10)') 'at='
          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
          write(6,'(a10)') 'dem_s='
          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
        ENDDO
      ENDIF

!     ---------------------------------------------------!

      IF (ncolprint.ne.0) THEN
      DO j=1,npoints,1000
        write(6,'(a10)') 'j='
        write(6,'(8I10)') j
      ENDDO
      ENDIF

      IF (top_height .eq. 1 .or. top_height .eq. 3) THEN 

      DO j=1,npoints 
          ptrop(j)=5000._wp
          attropmin(j) = 400._wp
          atmax(j) = 0._wp
          attrop(j) = 120._wp
          itrop(j) = 1
      ENDDO 

      DO 12 ilev=1,nlev
        DO j=1,npoints 
         IF (pfull(j,ilev) .lt. 40000._wp .and.   &
     &          pfull(j,ilev) .gt.  5000._wp .and. &
     &          at(j,ilev) .lt. attropmin(j)) THEN
                ptrop(j) = pfull(j,ilev)
                attropmin(j) = at(j,ilev)
                attrop(j) = attropmin(j)
                itrop(j)=ilev
           ENDIF
        ENDDO
12    CONTINUE

      DO 13 ilev=1,nlev
        DO j=1,npoints 
           IF (at(j,ilev) .gt. atmax(j) .and. &
     &              ilev  .ge. itrop(j)) atmax(j)=at(j,ilev)
        ENDDO
13    CONTINUE

      ENDIF


      IF (top_height .eq. 1 .or. top_height .eq. 3) THEN
          DO j=1,npoints
              meantb(j) = 0._wp
              meantbclr(j) = 0._wp 
          ENDDO
      ELSE
          DO j=1,npoints
              meantb(j) = output_missing_value
              meantbclr(j) = output_missing_value
          ENDDO
      ENDIF
      
!     -----------------------------------------------------!

!     ---------------------------------------------------!

      DO ilev=1,nlev
        DO j=1,npoints

          rangevec(j)=0

          IF (cc(j,ilev) .lt. 0._wp .or. cc(j,ilev) .gt. 1._wp) THEN
!           error = cloud fraction less than zero
!           error = cloud fraction greater than 1
            rangevec(j)=rangevec(j)+1
          ENDIF


          IF (dtau_s(j,ilev) .lt. 0._wp) THEN
!           ' error = stratiform cloud opt. depth less than zero'
            rangevec(j)=rangevec(j)+4
          ENDIF


          IF (dem_s(j,ilev) .lt. 0._wp .or. dem_s(j,ilev) .gt. 1._wp) THEN
!             ' error = stratiform cloud emissivity less than zero'
!             ' error = stratiform cloud emissivity greater than 1'
            rangevec(j)=rangevec(j)+16
          ENDIF
        ENDDO 

        rangeerror=0
        DO j=1,npoints
            rangeerror=rangeerror+rangevec(j)
        ENDDO

        IF (rangeerror.ne.0) THEN 
              write (6,*) 'Input variable out of range'
              write (6,*) 'rangevec:'
              write (6,*) rangevec
!              call flush(6)
              STOP
        ENDIF
      ENDDO
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      DO 15 ibox=1,ncol
        DO j=1,npoints 
            tau(j,ibox)=0._wp
          albedocld(j,ibox)=0._wp
          boxtau(j,ibox)=output_missing_value
          boxptop(j,ibox)=output_missing_value
          box_cloudy(j,ibox)=.false.
! BMD 07/26 - start
          tau_phase_cum(j,ibox)=0._wp
          tauf_phase_cum(j,ibox)=0._wp
          phase_active(j,ibox)=.true.
          ctp_tau_cum(j,ibox)=0._wp
          ctp_pres_cum(j,ibox)=0._wp
          ctp_active(j,ibox)=.true.
          obs_refl_active(j,ibox)=.false.
          obs_refl_run(j,ibox)=0._wp
          obs_tran_run(j,ibox)=1._wp
! BMD 07/26 - end
        ENDDO
15    CONTINUE

      !compute total cloud optical depth for each column     
      DO ilev=1,nlev
            !increment tau for each of the boxes
            DO ibox=1,ncol
              DO j=1,npoints 
                 IF (frac_out(j,ibox,ilev).eq.1._wp) THEN
                        tau(j,ibox)=tau(j,ibox) &
     &                     + dtau_s(j,ilev)
                 ENDIF
              ENDDO
            ENDDO ! ibox
      ENDDO ! ilev

! BMD 07/26 - start
      !Gated: dtau_s_liq/dtau_s_ice/reffl/reffi are OPTIONAL, and it is
      !illegal to reference an absent optional dummy argument -- callers
      !that don't pass this group (e.g. mo_cosp_isccp_simulator.f90) must
      !skip this block entirely.
      IF (PRESENT(dtau_s_liq)) THEN

      !Phase discrimination, ported from COSPv2's weight_by_extinction
      !(modis_simulator.F90) via ~/COSPv2.0: the optical-depth-weighted mean
      !liquid fraction over only the top modis_phase_taulimit (1.0) of cloud
      !optical depth from cloud top -- a passive
      !VIS/NIR sensor only sees the cloud-top phase. ilev=1..nlev is already
      !TOA-to-surface (see pfull/phalf convention above), matching
      !weight_by_extinction's implicit top-down ordering, so no extra flip
      !is needed. phase_active freezes each (j,ibox) once its cumulative tau
      !reaches modis_phase_taulimit, taking only the partial contribution
      !from the crossing level -- the vectorized-over-npoints equivalent of
      !weight_by_extinction's scalar "exit".
      DO ilev=1,nlev
            DO ibox=1,ncol
              DO j=1,npoints
                 IF (frac_out(j,ibox,ilev).eq.1._wp .and. &
     &               phase_active(j,ibox)) THEN
                    phase_tau_incr = dtau_s_liq(j,ilev) + dtau_s_ice(j,ilev)
                    IF (phase_tau_incr .gt. 0._wp) THEN
                       phase_f_liq = dtau_s_liq(j,ilev) / phase_tau_incr
                       IF (tau_phase_cum(j,ibox)+phase_tau_incr .gt. &
     &                     modis_phase_taulimit) THEN
                          phase_deltatau = modis_phase_taulimit &
     &                       - tau_phase_cum(j,ibox)
                          tau_phase_cum(j,ibox) = tau_phase_cum(j,ibox) &
     &                       + phase_deltatau
                          tauf_phase_cum(j,ibox) = tauf_phase_cum(j,ibox) &
     &                       + phase_deltatau*phase_f_liq
                          phase_active(j,ibox) = .false.
                       ELSE
                          tau_phase_cum(j,ibox) = tau_phase_cum(j,ibox) &
     &                       + phase_tau_incr
                          tauf_phase_cum(j,ibox) = tauf_phase_cum(j,ibox) &
     &                       + phase_tau_incr*phase_f_liq
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
            ENDDO ! ibox
      ENDDO ! ilev
      ENDIF
! BMD 07/26 - end

! BMD 07/26 - start
      !CO2-slicing cloud-top-pressure retrieval, ported from COSPv2's
      !cloud_top_pressure (modis_simulator.F90) via ~/COSPv2.0: an
      !optical-depth-weighted trapezoidal-rule pressure integral (pressure
      !assumed to vary linearly across each layer, unlike the phase
      !discrimination above which assumes a constant value per layer),
      !restricted to the top modis_co2slicing_taulimit (1.0) of cloud
      !optical depth from cloud top. Layer-edge pressures come from phalf
      !(already Pa, TOA-to-surface, same convention as pfull above);
      !divided by 100 to match this file's mb convention (ptop/boxptop),
      !so the result can be blended with ptop(j,ibox) (ISCCP's CTP, the
      !low-cloud fallback) and binned with the same edges as ipres(j)
      !below without further unit conversion. Only needs dtau_s/frac_out/
      !phalf (all non-optional) -- gated on PRESENT(fq_modis_tau_ctp)
      !specifically, independent of the dtau_s_liq/reffl group above,
      !since the (phase-independent) tau/CTP histogram doesn't need the
      !liquid/ice split at all.
      IF (PRESENT(fq_modis_tau_ctp)) THEN
      DO ilev=1,nlev
            DO ibox=1,ncol
              DO j=1,npoints
                 IF (frac_out(j,ibox,ilev).eq.1._wp .and. &
     &               ctp_active(j,ibox)) THEN
                    IF (dtau_s(j,ilev) .gt. 0._wp) THEN
                       IF (ctp_tau_cum(j,ibox)+dtau_s(j,ilev) .gt. &
     &                     modis_co2slicing_taulimit) THEN
                          ctp_deltatau = modis_co2slicing_taulimit &
     &                       - ctp_tau_cum(j,ibox)
                          ctp_tau_cum(j,ibox) = ctp_tau_cum(j,ibox) &
     &                       + ctp_deltatau
                          ctp_pres_cum(j,ibox) = ctp_pres_cum(j,ibox) &
     &                       + (phalf(j,ilev)/100._wp)*ctp_deltatau &
     &                       + ((phalf(j,ilev+1)-phalf(j,ilev))/100._wp) &
     &                         *ctp_deltatau**2/(2._wp*dtau_s(j,ilev))
                          ctp_active(j,ibox) = .false.
                       ELSE
                          ctp_tau_cum(j,ibox) = ctp_tau_cum(j,ibox) &
     &                       + dtau_s(j,ilev)
                          ctp_pres_cum(j,ibox) = ctp_pres_cum(j,ibox) &
     &                       + dtau_s(j,ilev) &
     &                         *(phalf(j,ilev)+phalf(j,ilev+1)) &
     &                         /2._wp/100._wp
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
            ENDDO ! ibox
      ENDDO ! ilev
      ENDIF
! BMD 07/26 - end

! BMD 09/26 - start
      !Per-level near-IR (true MODIS band 20) single-scattering albedo/
      !asymmetry parameter, ported from COSPv2's modis_optics (cosp_optics.f90)
      !via ~/COSPv2.0: get_g_nir/get_ssa_nir evaluated at the model's real
      !reffl/reffi (meters -> microns via *1.e6_wp, matching CAM6's own
      !COSP2 MODIS_waterSize/MODIS_iceSize sourcing in cospsimulator_intr.F90),
      !then tau-weighted between liquid and ice. Supersedes dpiz_s_nir/
      !dcg_s_nir (RRTMG-SW band 16, 3.08-3.85 micron -- much wider than
      !MODIS band 20's actual 3.66-3.84 micron). Gated on PRESENT(dtau_s_liq)
      !like the phase-discrimination block above, since it needs the same
      !dtau_s_liq/dtau_s_ice/reffl/reffi group; masked to levels with
      !nonzero liquid+ice tau to avoid a 0/0 in the w0-weighted g average.
      IF (PRESENT(dtau_s_liq)) THEN
      WHERE (dtau_s_liq+dtau_s_ice .gt. 0._wp)
         nir_w0_liq      = get_ssa_nir(.true.,  reffl*1.e6_wp)
         nir_g_liq       = get_g_nir(.true.,  reffl*1.e6_wp)
         nir_w0_ice      = get_ssa_nir(.false., reffi*1.e6_wp)
         nir_g_ice       = get_g_nir(.false., reffi*1.e6_wp)
         nir_w0_combined = (dtau_s_liq*nir_w0_liq + dtau_s_ice*nir_w0_ice) &
     &                    / (dtau_s_liq+dtau_s_ice)
         nir_g_combined  = (dtau_s_liq*nir_w0_liq*nir_g_liq &
     &                    + dtau_s_ice*nir_w0_ice*nir_g_ice) &
     &                    / (nir_w0_combined*(dtau_s_liq+dtau_s_ice))
      ELSEWHERE
         nir_w0_combined = 0._wp
         nir_g_combined  = 0._wp
      ENDWHERE
      ENDIF
! BMD 09/26 - end

! BMD 07/26 - start
      !Simulated-reflectance synthesis (Stage 4 of the MODIS size-retrieval
      !port, BMD 07/26 fourth pass), ported from COSPv2's
      !compute_toa_reflectace/two_stream/adding_doubling (modis_simulator.F90)
      !via ~/COSPv2.0: synthesizes the "observed" near-IR reflectance per
      !subcolumn from the model's real per-level tau/single-scattering
      !albedo/asymmetry (nir_w0_combined/nir_g_combined above), combining
      !cloudy levels top-down via the adding-layers recursion (same math as
      !adding_doubling in the CONTAINS block, called here 2 layers at a
      !time: the running state as "layer 1", the current level as
      !"layer 2" -- mathematically identical to running the full N-layer
      !recursion at once, since it's already strictly sequential). No tau
      !limit here -- unlike the phase/CTP retrievals above,
      !compute_toa_reflectace processes the whole cloudy profile. Gated on
      !PRESENT(dtau_s_liq) specifically, independent of the
      !fq_modis_tau_ctp group above, since this only needs
      !dtau_s/nir_w0_combined/nir_g_combined/frac_out.
      IF (PRESENT(dtau_s_liq)) THEN
      DO ilev=1,nlev
            DO ibox=1,ncol
              DO j=1,npoints
                 IF (frac_out(j,ibox,ilev).eq.1._wp .and. &
     &               dtau_s(j,ilev).gt.0._wp) THEN
                    CALL two_stream(dtau_s(j,ilev), nir_g_combined(j,ilev), &
     &                 nir_w0_combined(j,ilev), refl_layer, tran_layer)
                    IF (.not. obs_refl_active(j,ibox)) THEN
                       obs_refl_run(j,ibox) = refl_layer
                       obs_tran_run(j,ibox) = tran_layer
                       obs_refl_active(j,ibox) = .true.
                    ELSE
                       CALL adding_doubling(2, &
     &                    (/ obs_refl_run(j,ibox), refl_layer /), &
     &                    (/ obs_tran_run(j,ibox), tran_layer /), &
     &                    refl_combined, tran_combined)
                       obs_refl_run(j,ibox) = refl_combined
                       obs_tran_run(j,ibox) = tran_combined
                    ENDIF
                 ENDIF
              ENDDO
            ENDDO ! ibox
      ENDDO ! ilev
      ENDIF
! BMD 07/26 - end

          IF (ncolprint.ne.0) THEN

              DO j=1,npoints ,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write(6,'(i2,1X,8(f7.2,1X))')  &
     &          ilev,                          &
     &          (tau(j,ibox),ibox=1,ncolprint)
              ENDDO
          ENDIF 
!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done IF top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      IF (top_height .eq. 1 .or. top_height .eq. 3) THEN
        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644_wp
        wtmh20 = 18.01534_wp
        Navo = 6.023E+23_wp
        grav = 9.806650E+02_wp
        pstd = 1.013250E+06_wp
        t0 = 296._wp
        IF (ncolprint .ne. 0)  &
     &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
        DO 125 ilev=1,nlev
          DO j=1,npoints 
               !press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10._wp
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10._wp
               !atmden = g/cm2 = kg/m2 / 10 
               atmden(j) = dpress(j)/grav
               rvh20(j) = qv(j,ilev)*wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
               tmpexp(j) = exp(-0.02_wp*(at(j,ilev)-t0))
               tauwv(j) = wk(j)*1.e-20_wp*(  &
     &           (0.0224697_wp*rh20s(j)*tmpexp(j)) +  &
     &                (3.41817e-7_wp*rfrgn(j)) )*0.98_wp
               dem_wv(j,ilev) = 1._wp - exp( -1._wp * tauwv(j))
          ENDDO
               IF (ncolprint .ne. 0) THEN
               DO j=1,npoints ,1000
               write(6,'(a10)') 'j='
               write(6,'(8I10)') j
               write(6,'(i2,1X,3(f8.3,3X))') ilev, &
     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100._wp), &
     &           tauwv(j),dem_wv(j,ilev)
               ENDDO
             ENDIF
125     CONTINUE

        !initialize variables
        DO j=1,npoints 
          fluxtop_clrsky(j) = 0._wp
          trans_layers_above_clrsky(j)=1._wp
        ENDDO

        DO ilev=1,nlev
          DO j=1,npoints 
 
            ! Black body emission at temperature of the layer

              bb(j)=1._wp / ( exp(1307.27_wp/at(j,ilev)) - 1._wp )
              !bb(j)= 5.67e-8*at(j,ilev)**4

              ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

                fluxtop_clrsky(j) = fluxtop_clrsky(j)  &
     &            + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j) 
            
                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

                trans_layers_above_clrsky(j)= &
     &            trans_layers_above_clrsky(j)*(1._wp-dem_wv(j,ilev))
          ENDDO   
            IF (ncolprint.ne.0) THEN
             DO j=1,npoints ,1000
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
              write (6,'(a)') &
     &        'emiss_layer,100._wp*bb(j),100._wp*f,total_trans:'
              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100._wp*bb(j), &
     &             100._wp*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
             ENDDO   
            ENDIF
        ENDDO   !loop over level
        
        DO j=1,npoints 
          !add in surface emission
          bb(j)=1._wp/( exp(1307.27_wp/skt(j)) - 1._wp )
          !bb(j)=5.67e-8*skt(j)**4

          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j)  &
     &     * trans_layers_above_clrsky(j)
     
          !clear sky brightness temperature
          meantbclr(j) = 1307.27_wp/(log(1._wp+(1._wp/fluxtop_clrsky(j))))
        ENDDO

        IF (ncolprint.ne.0) THEN
        DO j=1,npoints ,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'
          write (6,'(a)') 'emsfc,100._wp*bb(j),100._wp*f,total_trans:'
          write (6,'(5(f7.2,1X))') emsfc_lw,100._wp*bb(j), &
     &      100._wp*fluxtop_clrsky(j), &
     &       trans_layers_above_clrsky(j), meantbclr(j)
        ENDDO
      ENDIF

        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------

        IF (ncolprint.ne.0) THEN

        DO j=1,npoints ,1000
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'ts:'
            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
            write (6,'(a)') 'ta_rev:'
            write (6,'(8f7.2)') &
     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
        ENDDO
        ENDIF 
        !loop over columns 
        DO ibox=1,ncol
          DO j=1,npoints
            fluxtop(j,ibox)=0._wp
            trans_layers_above(j,ibox)=1._wp
          ENDDO
        ENDDO

        DO ilev=1,nlev
              DO j=1,npoints 
                ! Black body emission at temperature of the layer

              bb(j)=1._wp / ( exp(1307.27_wp/at(j,ilev)) - 1._wp )
              !bb(j)= 5.67e-8*at(j,ilev)**4
              ENDDO

            DO ibox=1,ncol
              DO j=1,npoints 

              ! emissivity for point in this layer
                IF (frac_out(j,ibox,ilev).eq.1._wp) THEN
                dem(j,ibox)= 1._wp - &
     &             ( (1._wp - dem_wv(j,ilev)) * (1._wp -  dem_s(j,ilev)) )
                ELSE
                dem(j,ibox)=  dem_wv(j,ilev)
                ENDIF

                ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

                fluxtop(j,ibox) = fluxtop(j,ibox) &
     &            + dem(j,ibox) * bb(j) &
     &            * trans_layers_above(j,ibox) 
            
                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

                trans_layers_above(j,ibox)= &
     &            trans_layers_above(j,ibox)*(1._wp-dem(j,ibox))

              ENDDO ! j
            ENDDO ! ibox

            IF (ncolprint.ne.0) THEN
              DO j=1,npoints,1000
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'emiss_layer:'
              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
              write (6,'(a)') '100.*bb(j):'
              write (6,'(8f7.2)') (100._wp*bb(j),ibox=1,ncolprint)
              write (6,'(a)') '100.*f:'
              write (6,'(8f7.2)')  &
     &         (100._wp*fluxtop(j,ibox),ibox=1,ncolprint)
              write (6,'(a)') 'total_trans:'
              write (6,'(8f7.2)') &
     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
            ENDDO
          ENDIF

        ENDDO ! ilev

          DO j=1,npoints 
            !add in surface emission
            bb(j)=1._wp/( exp(1307.27_wp/skt(j)) - 1._wp )
            !bb(j)=5.67e-8*skt(j)**4
          ENDDO

        DO ibox=1,ncol
          DO j=1,npoints 

            !add in surface emission

            fluxtop(j,ibox) = fluxtop(j,ibox) &
     &         + emsfc_lw * bb(j) &
     &         * trans_layers_above(j,ibox) 
            
          ENDDO
        ENDDO

        !calculate mean infrared brightness temperature
        DO ibox=1,ncol
          DO j=1,npoints 
            meantb(j) = meantb(j)+1307.27_wp/(log(1._wp+(1._wp/fluxtop(j,ibox))))
          ENDDO
        ENDDO
          DO j=1, npoints
            meantb(j) = meantb(j) / REAL(ncol,wp)
          ENDDO        

        IF (ncolprint.ne.0) THEN

          DO j=1,npoints ,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'
          write (6,'(a)') 'emiss_layer:'
          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
          write (6,'(a)') '100.*bb(j):'
          write (6,'(8f7.2)') (100._wp*bb(j),ibox=1,ncolprint)
          write (6,'(a)') '100.*f:'
          write (6,'(8f7.2)') (100._wp*fluxtop(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'meantb(j):'
          write (6,'(8f7.2)') (meantb(j),ibox=1,ncolprint)
      
          ENDDO
      ENDIF
    
        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature

        !account for partially transmitting cloud recompute flux 
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the 
        !adjustment performed in this section
        !
      !If it turns out that the cloud brightness temperature
      !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
      !
        !Note that this is discussed on pages 85-87 of 
        !the ISCCP D level documentation (Rossow et al. 1996)
           
          DO j=1,npoints  
            !compute minimum brightness temperature and optical depth
            btcmin(j) = 1._wp /  ( exp(1307.27_wp/(attrop(j)-5._wp)) - 1._wp ) 
          ENDDO 
        DO ibox=1,ncol
          DO j=1,npoints  
            transmax(j) = (fluxtop(j,ibox)-btcmin(j)) &
     &                /(fluxtop_clrsky(j)-btcmin(j))
          !note that the initial setting of tauir(j) is needed so that
          !tauir(j) has a realistic value should the next if block be
          !bypassed
            tauir(j) = tau(j,ibox) * rec2p13
            taumin(j) = -1._wp * log(max(min(transmax(j),0.9999999_wp),0.001_wp))

          ENDDO 

          IF (top_height .eq. 1) THEN
            DO j=1,npoints  
              IF (transmax(j) .gt. 0.001_wp .and. &
     &          transmax(j) .le. 0.9999999_wp) THEN
                fluxtopinit(j) = fluxtop(j,ibox)
              tauir(j) = tau(j,ibox) *rec2p13
              ENDIF
            ENDDO
            DO icycle=1,2
              DO j=1,npoints  
                IF (tau(j,ibox) .gt. (tauchk            )) THEN 
                IF (transmax(j) .gt. 0.001_wp .and. &
     &            transmax(j) .le. 0.9999999_wp) THEN
                  emcld(j,ibox) = 1._wp - exp(-1._wp * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -  &  
     &              ((1._wp-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.E-06_wp, &
     &              (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= 1307.27_wp &
     &              / (log(1._wp + (1._wp/fluxtop(j,ibox))))
                  IF (tb(j,ibox) .gt. 260._wp) THEN
                  tauir(j) = tau(j,ibox) / 2.56_wp
                  ENDIF                   
                 ENDIF
                ENDIF
              ENDDO
            ENDDO
                
          ENDIF
        
          DO j=1,npoints
            IF (tau(j,ibox) .gt. (tauchk            )) THEN 
                !cloudy box 
                !NOTE: tb is the cloud-top temperature not infrared brightness temperature 
                !at this point in the code
                tb(j,ibox)= 1307.27_wp/ (log(1._wp + (1._wp/fluxtop(j,ibox))))
                IF (top_height.eq.1.and.tauir(j).lt.taumin(j)) THEN
                         tb(j,ibox) = attrop(j) - 5._wp
                   tau(j,ibox) = 2.13_wp*taumin(j)
                ENDIF
            ELSE
                !clear sky brightness temperature
                tb(j,ibox) = meantbclr(j)
            ENDIF
          ENDDO ! j
        ENDDO ! ibox

        IF (ncolprint.ne.0) THEN

          DO j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'attrop:'
          write (6,'(8f7.2)') (attrop(j))
          write (6,'(a)') 'btcmin:'
          write (6,'(8f7.2)') (btcmin(j))
          write (6,'(a)') 'fluxtop_clrsky*100:'
          write (6,'(8f7.2)') &
     &      (100._wp*fluxtop_clrsky(j))
          write (6,'(a)') '100.*f_adj:'
          write (6,'(8f7.2)') (100._wp*fluxtop(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'transmax:'
          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
          write (6,'(a)') 'tau:'
          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'emcld:'
          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') &
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'total_emiss:'
          write (6,'(8f7.2)') &
     &        (1.0_wp-trans_layers_above(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') &
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'ppout:'
          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
          ENDDO ! j
      ENDIF

      ENDIF
!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      DO 30 ibox=1,ncol
        !segregate according to optical thickness
        IF (top_height .eq. 1 .or. top_height .eq. 3) THEN  
          !find level whose temperature
          !most closely matches brightness temperature
          DO j=1,npoints 
            nmatch(j)=0
          ENDDO
          DO 29 k1=1,nlev-1
            IF (top_height_direction .eq. 2) THEN
              ilev = nlev - k1 
            ELSE
              ilev = k1
            ENDIF
            !cdir nodep
            DO j=1,npoints 
             IF (ilev .ge. itrop(j)) THEN
              IF ((at(j,ilev)   .ge. tb(j,ibox) .and. &
     &          at(j,ilev+1) .le. tb(j,ibox)) .or. &
     &          (at(j,ilev) .le. tb(j,ibox) .and.  &
     &          at(j,ilev+1) .ge. tb(j,ibox))) THEN 
                nmatch(j)=nmatch(j)+1
                match(j,nmatch(j))=ilev
              ENDIF  
             ENDIF                         
            ENDDO
29        CONTINUE

          DO j=1,npoints 
            IF (nmatch(j) .ge. 1) THEN
              k1 = match(j,nmatch(j))
              k2 = k1 + 1
              logp1 = log(pfull(j,k1))
              logp2 = log(pfull(j,k2))
              atd = max(tauchk,abs(at(j,k2) - at(j,k1)))
              logp=logp1+(logp2-logp1)*abs(tb(j,ibox)-at(j,k1))/atd
              ptop(j,ibox) = exp(logp)
              IF(abs(pfull(j,k1)-ptop(j,ibox)) .lt. &
     &            abs(pfull(j,k2)-ptop(j,ibox))) THEN
                 levmatch(j,ibox)=k1
              ELSE
                 levmatch(j,ibox)=k2
              ENDIF   
            ELSE
              IF (tb(j,ibox) .le. attrop(j)) THEN
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              ENDIF
              IF (tb(j,ibox) .ge. atmax(j)) THEN
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              ENDIF                                
            ENDIF
          ENDDO ! j

        ELSE ! if (top_height .eq. 1 .or. top_height .eq. 3) 
 
          DO j=1,npoints     
            ptop(j,ibox)=0._wp
          ENDDO
          DO ilev=1,nlev
            DO j=1,npoints     
              IF ((ptop(j,ibox) .eq. 0._wp ) &
     &           .and.(frac_out(j,ibox,ilev) .ne. 0._wp)) THEN
                ptop(j,ibox)=phalf(j,ilev)
              levmatch(j,ibox)=ilev
              ENDIF
            ENDDO
          ENDDO
        ENDIF                            
          
        DO j=1,npoints
          IF (tau(j,ibox) .le. (tauchk            )) THEN
            ptop(j,ibox)=0._wp
            levmatch(j,ibox)=0      
          ENDIF 
        ENDDO

30    CONTINUE
              
!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy 
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.  
!
      !compute isccp frequencies

      !reset frequencies
      DO 37 ilev=1,7
      DO 38 ilev2=1,7
        DO j=1,npoints ! 
             IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN 
                fq_isccp(j,ilev,ilev2)= 0._wp
             ELSE
                fq_isccp(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
38    CONTINUE
37    CONTINUE

! BMD 07/26 - start
      !reset MODIS-style LWP-vs-ReffLiq and IWP-vs-ReffIce joint-histogram
      !frequencies. Unlike the ISCCP tau/CTP histogram, these bins rely on
      !solar reflectance retrievals only, so they are gated on sunlit alone
      !(no top_height==3 IR-only override, which is ISCCP-specific).
      !Gated on PRESENT: fq_modis_lwp_reffl etc. are OPTIONAL, omitted by
      !callers that don't want MODIS output (e.g. the ISCCP-only call site
      !in mo_cosp_isccp_simulator.f90) -- referencing them unconditionally
      !would dereference an absent argument.
      IF (PRESENT(fq_modis_lwp_reffl)) THEN
      DO ilev=1,n_modis_lwp
      DO ilev2=1,n_modis_reffl
        DO j=1,npoints
             IF (sunlit(j).eq.1) THEN
                fq_modis_lwp_reffl(j,ilev,ilev2)= 0._wp
             ELSE
                fq_modis_lwp_reffl(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
      ENDDO
      ENDDO

      DO ilev=1,n_modis_iwp
      DO ilev2=1,n_modis_reffi
        DO j=1,npoints
             IF (sunlit(j).eq.1) THEN
                fq_modis_iwp_reffi(j,ilev,ilev2)= 0._wp
             ELSE
                fq_modis_iwp_reffi(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
      ENDDO
      ENDDO

      !reset MODIS-style tau-vs-ReffLiq and tau-vs-ReffIce joint-histogram
      !frequencies -- same sunlit-only gate, same OPTIONAL group (populated
      !alongside fq_modis_lwp_reffl/fq_modis_iwp_reffi from the same
      !box_reffl/box_reffi retrieval outcome, just crossed with the tau bin
      !instead of the LWP/IWP bin).
      DO ilev=1,n_modis_tau
      DO ilev2=1,n_modis_reffl
        DO j=1,npoints
             IF (sunlit(j).eq.1) THEN
                fq_modis_tau_reffl(j,ilev,ilev2)= 0._wp
             ELSE
                fq_modis_tau_reffl(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
      ENDDO
      ENDDO

      DO ilev=1,n_modis_tau
      DO ilev2=1,n_modis_reffi
        DO j=1,npoints
             IF (sunlit(j).eq.1) THEN
                fq_modis_tau_reffi(j,ilev,ilev2)= 0._wp
             ELSE
                fq_modis_tau_reffi(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
      ENDDO
      ENDDO

      !reset MODIS total/liquid/ice cloud fraction
      DO j=1,npoints
           IF (sunlit(j).eq.1) THEN
              cltmodis(j)= 0._wp
              clwmodis(j)= 0._wp
              climodis(j)= 0._wp
           ELSE
              cltmodis(j)= output_missing_value
              clwmodis(j)= output_missing_value
              climodis(j)= output_missing_value
           ENDIF
      ENDDO
      ENDIF
! BMD 07/26 - end

! BMD 09/26 - start
      !reset MODIS retrieval bracket-failure diagnostic counters (see the
      !declaration comment above for what these measure).
      IF (PRESENT(modis_diag_liq_attempt)) THEN
      DO j=1,npoints
           IF (sunlit(j).eq.1) THEN
              modis_diag_liq_attempt(j)= 0._wp
              modis_diag_liq_fail(j)= 0._wp
              modis_diag_liq_obshi(j)= 0._wp
              modis_diag_liq_obslo(j)= 0._wp
              modis_diag_ice_attempt(j)= 0._wp
              modis_diag_ice_fail(j)= 0._wp
              modis_diag_ice_obshi(j)= 0._wp
              modis_diag_ice_obslo(j)= 0._wp
           ELSE
              modis_diag_liq_attempt(j)= output_missing_value
              modis_diag_liq_fail(j)= output_missing_value
              modis_diag_liq_obshi(j)= output_missing_value
              modis_diag_liq_obslo(j)= output_missing_value
              modis_diag_ice_attempt(j)= output_missing_value
              modis_diag_ice_fail(j)= output_missing_value
              modis_diag_ice_obshi(j)= output_missing_value
              modis_diag_ice_obslo(j)= output_missing_value
           ENDIF
      ENDDO
      ENDIF
! BMD 09/26 - end

! BMD 07/26 - start
      !reset MODIS-style tau/CTP joint-histogram frequencies (all phases,
      !liquid-only, ice-only). Same sunlit-only gate as the LWP/IWP reset
      !above (no top_height==3 IR-only override), and its own PRESENT
      !check since fq_modis_tau_ctp/_liq/_ice are their own OPTIONAL group.
      IF (PRESENT(fq_modis_tau_ctp)) THEN
      DO ilev=1,n_modis_tau
      DO ilev2=1,n_modis_pres
        DO j=1,npoints
             IF (sunlit(j).eq.1) THEN
                fq_modis_tau_ctp(j,ilev,ilev2)= 0._wp
                fq_modis_tau_ctp_liq(j,ilev,ilev2)= 0._wp
                fq_modis_tau_ctp_ice(j,ilev,ilev2)= 0._wp
             ELSE
                fq_modis_tau_ctp(j,ilev,ilev2)= output_missing_value
                fq_modis_tau_ctp_liq(j,ilev,ilev2)= output_missing_value
                fq_modis_tau_ctp_ice(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
      ENDDO
      ENDDO
      ENDIF
! BMD 07/26 - end

      !reset variables need for averaging cloud properties
      DO j=1,npoints 
        IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN 
             totalcldarea(j) = 0._wp
             meanalbedocld(j) = 0._wp
             meanptop(j) = 0._wp
             meantaucld(j) = 0._wp
        ELSE
             totalcldarea(j) = output_missing_value
             meanalbedocld(j) = output_missing_value
             meanptop(j) = output_missing_value
             meantaucld(j) = output_missing_value
        ENDIF
      ENDDO ! j

      boxarea = 1._wp/REAL(ncol,wp)

! BMD 07/26 - start
      !Trial-size grid and its near-IR optical-property library (Stage 5 of
      !the MODIS size-retrieval port), ported from COSPv2's COSP_MODIS_INIT
      !(cosp_modis_interface.F90) via ~/COSPv2.0: computed once per icarus
      !call, not per subcolumn, since it doesn't depend on j/ibox. Reused
      !below by every liquid-/ice-dominated subcolumn's retrieval.
      IF (PRESENT(dtau_s_liq)) THEN
         trial_re_w(1:num_trial_res) = modis_trial_re_water_min &
     &      + (modis_trial_re_water_max-modis_trial_re_water_min) &
     &        /REAL(num_trial_res-1,wp) * (/(i,i=0,num_trial_res-1)/)
         trial_re_i(1:num_trial_res) = modis_trial_re_ice_min &
     &      + (modis_trial_re_ice_max-modis_trial_re_ice_min) &
     &        /REAL(num_trial_res-1,wp) * (/(i,i=0,num_trial_res-1)/)
         trial_g_w(1:num_trial_res)  = get_g_nir(.true.,trial_re_w(1:num_trial_res))
         trial_w0_w(1:num_trial_res) = get_ssa_nir(.true.,trial_re_w(1:num_trial_res))
         trial_g_i(1:num_trial_res)  = get_g_nir(.false.,trial_re_i(1:num_trial_res))
         trial_w0_i(1:num_trial_res) = get_ssa_nir(.false.,trial_re_i(1:num_trial_res))
      ENDIF
! BMD 07/26 - end

      DO 39 ibox=1,ncol
        DO j=1,npoints 

          IF (tau(j,ibox) .gt. (tauchk            ) &
     &      .and. ptop(j,ibox) .gt. 0._wp) THEN
              box_cloudy(j,ibox)=.true.
          ENDIF

          IF (box_cloudy(j,ibox)) THEN

              IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN

                boxtau(j,ibox) = tau(j,ibox)

                IF (tau(j,ibox) .ge. isccp_taumin) THEN

                   totalcldarea(j) = totalcldarea(j) + boxarea
                
                   !convert optical thickness to albedo
                   albedocld(j,ibox) &
     &                   = (tau(j,ibox)**0.895_wp)/((tau(j,ibox)**0.895_wp)+6.82_wp)
         
                   !contribute to averaging
                   meanalbedocld(j) = meanalbedocld(j) &
     &                                +albedocld(j,ibox)*boxarea

                ENDIF
            ENDIF
          ENDIF

          IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN 

           IF (box_cloudy(j,ibox)) THEN
          
              !convert ptop to millibars
              ptop(j,ibox)=ptop(j,ibox) / 100._wp
            
              !save for output cloud top pressure and optical thickness
              boxptop(j,ibox) = ptop(j,ibox)
    
              IF (tau(j,ibox) .ge. isccp_taumin) THEN
               meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea
              ENDIF

              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0

              !determine optical depth category
              IF (tau(j,ibox) .lt. isccp_taumin) THEN
                  itau(j)=1
              ELSE IF (tau(j,ibox) .ge. isccp_taumin &
     &          .and. tau(j,ibox) .lt. 1.3_wp) THEN
                itau(j)=2
              ELSE IF (tau(j,ibox) .ge. 1.3_wp  &
     &          .and. tau(j,ibox) .lt. 3.6_wp) THEN
                itau(j)=3
              ELSE IF (tau(j,ibox) .ge. 3.6_wp &
     &          .and. tau(j,ibox) .lt. 9.4_wp) THEN
                  itau(j)=4
              ELSE IF (tau(j,ibox) .ge. 9.4_wp &
     &          .and. tau(j,ibox) .lt. 23._wp) THEN
                  itau(j)=5
              ELSE IF (tau(j,ibox) .ge. 23._wp &
     &          .and. tau(j,ibox) .lt. 60._wp) THEN
                  itau(j)=6
              ELSE IF (tau(j,ibox) .ge. 60._wp) THEN
                  itau(j)=7
              ENDIF

              !determine cloud top pressure category
              IF (    ptop(j,ibox) .gt. 0._wp  &
     &          .and.ptop(j,ibox) .lt. 180._wp) THEN
                  ipres(j)=1
              ELSE IF(ptop(j,ibox) .ge. 180._wp &
     &          .and.ptop(j,ibox) .lt. 310._wp) THEN
                  ipres(j)=2
              ELSE IF(ptop(j,ibox) .ge. 310._wp &
     &          .and.ptop(j,ibox) .lt. 440._wp) THEN
                  ipres(j)=3
              ELSE IF(ptop(j,ibox) .ge. 440._wp &
     &          .and.ptop(j,ibox) .lt. 560._wp) THEN
                  ipres(j)=4
              ELSE IF(ptop(j,ibox) .ge. 560._wp &
     &          .and.ptop(j,ibox) .lt. 680._wp) THEN
                  ipres(j)=5
              ELSE IF(ptop(j,ibox) .ge. 680._wp &
     &          .and.ptop(j,ibox) .lt. 800._wp) THEN
                  ipres(j)=6
              ELSE IF(ptop(j,ibox) .ge. 800._wp) THEN
                  ipres(j)=7
              ENDIF 

              !update frequencies
              IF(ipres(j) .gt. 0.and.itau(j) .gt. 0) THEN
              fq_isccp(j,itau(j),ipres(j))= &
     &          fq_isccp(j,itau(j),ipres(j))+ boxarea
              ENDIF
            ENDIF
          ENDIF

! BMD 07/26 - start
          !MODIS-style cltmodis/clwmodis/climodis, tau/CTP (total/liquid/
          !ice), LWP-vs-ReffLiq, IWP-vs-ReffIce, and (BMD 08/26) tau-vs-
          !ReffLiq/tau-vs-ReffIce joint histograms. Matches COSPv2's
          !modis_column logic (modis_simulator.F90, ~/COSPv2.0): each
          !subcolumn is assigned a SINGLE dominant phase (liquid, ice, or
          !undetermined -- see modis_phase_threshold above, standing in for
          !that file's phaseDiscrimination_Threshold/phaseIsUndetermined),
          !then the whole subcolumn's TOTAL optical depth tau(j,ibox) (not
          !just its liquid- or ice-only component) is combined with that
          !phase's characteristic effective radius via the Stephens (1978)
          !bulk approximation CWP = (2/3)*rho*re*tau -- the same relation
          !modis_simulator.F90's lwpWRK/iwpWRK apply to the full retrievedTau.
          !Gate is modis_min_optical_thickness (0.3), not the near-zero
          !tauchk the ISCCP machinery above uses -- matches real MODIS's
          !actual detection limit (see modis_min_optical_thickness above).
          !Outer PRESENT check: fq_modis_lwp_reffl/cltmodis/etc. are
          !OPTIONAL, omitted by callers that don't want MODIS output (e.g.
          !mo_cosp_isccp_simulator.f90) -- referencing them unconditionally
          !would dereference an absent argument.
          IF (PRESENT(fq_modis_lwp_reffl)) THEN
          IF (sunlit(j).eq.1 .and. &
     &        tau(j,ibox) .ge. modis_min_optical_thickness) THEN

            !cloud-top-weighted phase fraction (tau_phase_cum/tauf_phase_cum,
            !ported from weight_by_extinction above). Moved ahead of
            !cltmodis/fq_modis_tau_ctp (BMD 08/26): the size-retrieval
            !attempt just below needs the phase classification up front to
            !pick a trial library, and that retrieval's success now gates
            !the total-cloud outputs (see the next comment).
            IF ( tau_phase_cum(j,ibox) .gt. (tauchk) ) THEN
              liquid_fraction(j) = tauf_phase_cum(j,ibox) &
     &          / tau_phase_cum(j,ibox)
            ELSE
              liquid_fraction(j) = -1._wp
            ENDIF

            !BMD 08/26 -- size-retrieval attempt, run once per subcolumn
            !regardless of phase outcome; its success gates the total-cloud
            !outputs below. This replaces the earlier design (547ce56's
            !commit message) that left cltmodis/fq_modis_tau_ctp
            !unconditional on retrieval success, believing that "matched
            !the reference" -- direct inspection of ~/COSPv2.0's
            !modis_simulator.F90 shows the opposite: modis_subcolumn
            !(:159-165) sets retrievedTau/cloudMask from tau alone, but
            !modis_column (:311-330) builds
            !validRetrievalMask = particle_size > 0 and ANDs it into
            !cloudMask before computing Cloud_Fraction_Total_Mean and
            !masking tauWRK/ctpWRK for
            !Optical_Thickness_vs_Cloud_Top_Pressure (:402-414) -- i.e. a
            !subcolumn whose near-IR reflectance retrieval fails to bracket
            !a trial size is dropped from the TOTAL tau/CTP histogram too,
            !not just from the LWP/IWP-vs-Reff ones. modis_subcolumn
            !(:204-215) also shows which library each phase uses:
            !phaseIsLiquid AND phaseIsUndetermined both invert against the
            !LIQUID trial library; only phaseIsIce uses the ice library --
            !so an undetermined-phase subcolumn still gets a (liquid-
            !library) retrieval attempt purely to decide whether it counts
            !toward the total, even though it never contributes to
            !clwmodis or the liquid-specific histograms below (those
            !require phase==liquid specifically, not just "a liquid-
            !library retrieval was attempted"). clwmodis/climodis
            !themselves are ALSO gated on box_reffl(j)/box_reffi(j) >= 0
            !(see the liquid/ice phase branches below), matching COSPv2's
            !Cloud_Fraction_Water_Mean/Cloud_Fraction_Ice_Mean, whose
            !waterCloudMask/iceCloudMask (modis_simulator.F90:314-317) are
            !ANDed with the same validRetrievalMask as everything else in
            !this block -- no divergence left here.
            IF (liquid_fraction(j) .gt. (1._wp-modis_phase_threshold)) THEN
              !liquid-dominated (>=modis_phase_threshold) or undetermined
              !((1-modis_phase_threshold, modis_phase_threshold)): liquid
              !trial library, matching modis_subcolumn's phaseIsLiquid .OR.
              !phaseIsUndetermined branch
              IF (PRESENT(dtau_s_liq) .and. obs_refl_active(j,ibox)) THEN
                 predicted_refl_w(1:num_trial_res) = two_stream_reflectance( &
     &              tau(j,ibox), trial_g_w(1:num_trial_res), &
     &              trial_w0_w(1:num_trial_res))
                 box_reffl(j) = 1.e-6_wp * interpolate_to_min( &
     &              trial_re_w(1:num_trial_res), &
     &              predicted_refl_w(1:num_trial_res), obs_refl_run(j,ibox))

! BMD 09/26 - start
                 IF (PRESENT(modis_diag_liq_attempt)) THEN
                    modis_diag_liq_attempt(j) = modis_diag_liq_attempt(j) + boxarea
                    IF (box_reffl(j) .lt. 0._wp) THEN
                       modis_diag_liq_fail(j) = modis_diag_liq_fail(j) + boxarea
                       IF (obs_refl_run(j,ibox) .gt. &
     &                     maxval(predicted_refl_w(1:num_trial_res))) THEN
                          modis_diag_liq_obshi(j) = modis_diag_liq_obshi(j) + boxarea
                       ELSE IF (obs_refl_run(j,ibox) .lt. &
     &                     minval(predicted_refl_w(1:num_trial_res))) THEN
                          modis_diag_liq_obslo(j) = modis_diag_liq_obslo(j) + boxarea
                       ENDIF
                    ENDIF
                 ENDIF
! BMD 09/26 - end
              ELSE
                 box_reffl(j) = -1._wp
              ENDIF
              retrieval_valid = (box_reffl(j) .ge. 0._wp)
            ELSE IF (liquid_fraction(j) .ge. 0._wp .and. &
     &        liquid_fraction(j) .le. (1._wp-modis_phase_threshold)) THEN
              !ice-dominated: ice trial library
              IF (PRESENT(dtau_s_liq) .and. obs_refl_active(j,ibox)) THEN
                 predicted_refl_i(1:num_trial_res) = two_stream_reflectance( &
     &              tau(j,ibox), trial_g_i(1:num_trial_res), &
     &              trial_w0_i(1:num_trial_res))
                 box_reffi(j) = 1.e-6_wp * interpolate_to_min( &
     &              trial_re_i(1:num_trial_res), &
     &              predicted_refl_i(1:num_trial_res), obs_refl_run(j,ibox))

! BMD 09/26 - start
                 IF (PRESENT(modis_diag_ice_attempt)) THEN
                    modis_diag_ice_attempt(j) = modis_diag_ice_attempt(j) + boxarea
                    IF (box_reffi(j) .lt. 0._wp) THEN
                       modis_diag_ice_fail(j) = modis_diag_ice_fail(j) + boxarea
                       IF (obs_refl_run(j,ibox) .gt. &
     &                     maxval(predicted_refl_i(1:num_trial_res))) THEN
                          modis_diag_ice_obshi(j) = modis_diag_ice_obshi(j) + boxarea
                       ELSE IF (obs_refl_run(j,ibox) .lt. &
     &                     minval(predicted_refl_i(1:num_trial_res))) THEN
                          modis_diag_ice_obslo(j) = modis_diag_ice_obslo(j) + boxarea
                       ENDIF
                    ENDIF
                 ENDIF
! BMD 09/26 - end
              ELSE
                 box_reffi(j) = -1._wp
              ENDIF
              retrieval_valid = (box_reffi(j) .ge. 0._wp)
            ELSE
              !liquid_fraction(j) < 0: no extinction at all within the
              !top-of-cloud phase-weighting window (an ECHAM-specific edge
              !case -- COSPv2's own weight_by_extinction has no true
              !sentinel here, it falls back to a 0.0 liquid fraction, i.e.
              !phaseIsIce; this file keeps the -1 sentinel as a distinct
              !"no phase info" case instead of silently reclassifying it as
              !ice -- a pre-existing, deliberate divergence, unrelated to
              !this retrieval-gating change). No retrieval is attempted.
              box_reffl(j) = -1._wp
              box_reffi(j) = -1._wp
              retrieval_valid = .false.
            ENDIF

            !Total-cloud outputs require the retrieval above to have
            !succeeded (see the BMD 08/26 note above).
            IF (retrieval_valid) THEN

            cltmodis(j) = cltmodis(j) + boxarea

            !MODIS-style tau/CTP joint histogram (all phases). itau(j) was
            !already computed above for fq_isccp -- reused directly here,
            !not recomputed, since COSPv2's MODIS tau bins are identical to
            !this file's ISCCP tau bins (see n_modis_tau note above).
            !modis_ctp_j is a genuinely different retrieval than ptop(j,ibox)
            !though: the CO2-sliced ctp_pres_cum/ctp_tau_cum result, falling
            !back to ptop(j,ibox) (ISCCP's CTP) for low clouds -- ports
            !COSPv2's modis_subcolumn low-cloud fallback.
            IF (PRESENT(fq_modis_tau_ctp)) THEN
              IF ( ctp_tau_cum(j,ibox) .gt. 0._wp ) THEN
                modis_ctp_j = ctp_pres_cum(j,ibox) / ctp_tau_cum(j,ibox)
              ELSE
                modis_ctp_j = 0._wp
              ENDIF
              IF ( modis_ctp_j .gt. modis_co2slicing_preslimit_mb ) THEN
                modis_ctp_j = ptop(j,ibox)
              ENDIF

              imodis_pres(j) = 0
              IF (    modis_ctp_j .gt. 0._wp  &
     &          .and.modis_ctp_j .lt. 180._wp) THEN
                  imodis_pres(j)=1
              ELSE IF(modis_ctp_j .ge. 180._wp &
     &          .and.modis_ctp_j .lt. 310._wp) THEN
                  imodis_pres(j)=2
              ELSE IF(modis_ctp_j .ge. 310._wp &
     &          .and.modis_ctp_j .lt. 440._wp) THEN
                  imodis_pres(j)=3
              ELSE IF(modis_ctp_j .ge. 440._wp &
     &          .and.modis_ctp_j .lt. 560._wp) THEN
                  imodis_pres(j)=4
              ELSE IF(modis_ctp_j .ge. 560._wp &
     &          .and.modis_ctp_j .lt. 680._wp) THEN
                  imodis_pres(j)=5
              ELSE IF(modis_ctp_j .ge. 680._wp &
     &          .and.modis_ctp_j .lt. 800._wp) THEN
                  imodis_pres(j)=6
              ELSE IF(modis_ctp_j .ge. 800._wp) THEN
                  imodis_pres(j)=7
              ENDIF

              IF (imodis_pres(j) .gt. 0 .and. itau(j) .gt. 0) THEN
                fq_modis_tau_ctp(j,itau(j),imodis_pres(j)) = &
     &            fq_modis_tau_ctp(j,itau(j),imodis_pres(j)) + boxarea
              ENDIF
            ENDIF

            ENDIF
            !ENDIF above closes the retrieval_valid guard around
            !cltmodis/fq_modis_tau_ctp (BMD 08/26). imodis_pres(j) is only
            !computed inside that guard, but that's safe: the liquid/ice
            !tau/CTP writes below are themselves now gated on
            !box_reffl(j)/box_reffi(j) >= 0, which (given the phase
            !classification is the same test used above) can only be true
            !when retrieval_valid was already true and imodis_pres(j) was
            !freshly computed this same (j,ibox) iteration.

            IF (liquid_fraction(j) .ge. modis_phase_threshold) THEN
              !liquid-dominated subcolumn -> LWP-vs-ReffLiq and tau-vs-
              !ReffLiq. box_reffl was already computed above (single
              !retrieval attempt per subcolumn, shared with the total-cloud
              !gate) -- not recomputed here.

              IF (box_reffl(j) .ge. 0._wp) THEN
              !clwmodis now INSIDE the box_reffl(j).ge.0._wp guard (BMD
              !08/26, previously unconditional on phase alone): matches
              !COSPv2's Cloud_Fraction_Water_Mean, whose waterCloudMask
              !(modis_simulator.F90:314-315) is phase==phaseIsLiquid ANDed
              !with the same validRetrievalMask as everything else in this
              !block -- closing the divergence noted (and deliberately left
              !open) in the previous retrieval-gating change.
              clwmodis(j) = clwmodis(j) + boxarea
              box_lwp(j) = modis_lwp_conversion * box_reffl(j) * tau(j,ibox)

              imodis_lwp(j) = 0
              imodis_reffl(j) = 0

              IF (box_lwp(j) .ge. modis_lwp_bin_edges(1) .and. &
     &          box_lwp(j) .lt. modis_lwp_bin_edges(2)) THEN
                imodis_lwp(j)=1
              ELSE IF (box_lwp(j) .ge. modis_lwp_bin_edges(2) .and. &
     &          box_lwp(j) .lt. modis_lwp_bin_edges(3)) THEN
                imodis_lwp(j)=2
              ELSE IF (box_lwp(j) .ge. modis_lwp_bin_edges(3) .and. &
     &          box_lwp(j) .lt. modis_lwp_bin_edges(4)) THEN
                imodis_lwp(j)=3
              ELSE IF (box_lwp(j) .ge. modis_lwp_bin_edges(4) .and. &
     &          box_lwp(j) .lt. modis_lwp_bin_edges(5)) THEN
                imodis_lwp(j)=4
              ELSE IF (box_lwp(j) .ge. modis_lwp_bin_edges(5) .and. &
     &          box_lwp(j) .lt. modis_lwp_bin_edges(6)) THEN
                imodis_lwp(j)=5
              ELSE IF (box_lwp(j) .ge. modis_lwp_bin_edges(6) .and. &
     &          box_lwp(j) .lt. modis_lwp_bin_edges(7)) THEN
                imodis_lwp(j)=6
              ELSE IF (box_lwp(j) .ge. modis_lwp_bin_edges(7)) THEN
                imodis_lwp(j)=7
              ENDIF

              IF (box_reffl(j) .ge. modis_reffl_bin_edges(1) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(2)) THEN
                imodis_reffl(j)=1
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(2) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(3)) THEN
                imodis_reffl(j)=2
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(3) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(4)) THEN
                imodis_reffl(j)=3
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(4) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(5)) THEN
                imodis_reffl(j)=4
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(5) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(6)) THEN
                imodis_reffl(j)=5
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(6) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(7)) THEN
                imodis_reffl(j)=6
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(7) .and. &
     &          box_reffl(j) .lt. modis_reffl_bin_edges(8)) THEN
                imodis_reffl(j)=7
              ELSE IF (box_reffl(j) .ge. modis_reffl_bin_edges(8)) THEN
                imodis_reffl(j)=8
              ENDIF

              IF (imodis_lwp(j) .gt. 0 .and. imodis_reffl(j) .gt. 0) THEN
                fq_modis_lwp_reffl(j,imodis_lwp(j),imodis_reffl(j)) = &
     &            fq_modis_lwp_reffl(j,imodis_lwp(j),imodis_reffl(j)) + boxarea
              ENDIF

              !MODIS-style tau-vs-ReffLiq joint histogram (BMD 08/26),
              !ported from COSPv2's Optical_Thickness_vs_ReffLiq
              !(modis_simulator.F90:438-440,455): the FULL subcolumn tau
              !(itau(j), same bin as fq_modis_tau_ctp) crossed with the
              !just-computed imodis_reffl(j) -- reuses both indices, no new
              !computation needed.
              IF (PRESENT(fq_modis_tau_reffl) .and. &
     &            itau(j) .gt. 0 .and. imodis_reffl(j) .gt. 0) THEN
                fq_modis_tau_reffl(j,itau(j),imodis_reffl(j)) = &
     &            fq_modis_tau_reffl(j,itau(j),imodis_reffl(j)) + boxarea
              ENDIF

              !liquid-phase tau/CTP: reuses itau(j)/imodis_pres(j) already
              !computed above, same PRESENT gate as the total tau/CTP
              !write. Now INSIDE the box_reffl(j).ge.0._wp guard (BMD
              !08/26, previously outside it) -- see the detailed
              !retrieval-gating note above: fq_modis_tau_ctp_liq can no
              !longer contain area that fq_modis_tau_ctp (total) excludes.
              IF (PRESENT(fq_modis_tau_ctp) .and. &
     &            imodis_pres(j) .gt. 0 .and. itau(j) .gt. 0) THEN
                fq_modis_tau_ctp_liq(j,itau(j),imodis_pres(j)) = &
     &            fq_modis_tau_ctp_liq(j,itau(j),imodis_pres(j)) + boxarea
              ENDIF
              ENDIF
              !ENDIF above closes the box_reffl(j).ge.0._wp guard, now
              !wrapping LWP-vs-ReffLiq, tau-vs-ReffLiq, AND tau/CTP-Liq
              !together (BMD 08/26).

            ELSE IF (liquid_fraction(j) .ge. 0._wp .and. &
     &        liquid_fraction(j) .le. (1._wp-modis_phase_threshold)) THEN
              !ice-dominated subcolumn -> IWP-vs-ReffIce and tau-vs-ReffIce.
              !box_reffi was already computed above -- see the liquid-
              !branch note.

              IF (box_reffi(j) .ge. 0._wp) THEN
              !climodis now INSIDE the box_reffi(j).ge.0._wp guard (BMD
              !08/26) -- same fix as clwmodis above, matching COSPv2's
              !Cloud_Fraction_Ice_Mean/iceCloudMask coupling to
              !validRetrievalMask.
              climodis(j) = climodis(j) + boxarea
              box_iwp(j) = modis_lwp_conversion * modis_ice_density &
     &          * box_reffi(j) * tau(j,ibox)

              imodis_iwp(j) = 0
              imodis_reffi(j) = 0

              IF (box_iwp(j) .ge. modis_iwp_bin_edges(1) .and. &
     &          box_iwp(j) .lt. modis_iwp_bin_edges(2)) THEN
                imodis_iwp(j)=1
              ELSE IF (box_iwp(j) .ge. modis_iwp_bin_edges(2) .and. &
     &          box_iwp(j) .lt. modis_iwp_bin_edges(3)) THEN
                imodis_iwp(j)=2
              ELSE IF (box_iwp(j) .ge. modis_iwp_bin_edges(3) .and. &
     &          box_iwp(j) .lt. modis_iwp_bin_edges(4)) THEN
                imodis_iwp(j)=3
              ELSE IF (box_iwp(j) .ge. modis_iwp_bin_edges(4) .and. &
     &          box_iwp(j) .lt. modis_iwp_bin_edges(5)) THEN
                imodis_iwp(j)=4
              ELSE IF (box_iwp(j) .ge. modis_iwp_bin_edges(5) .and. &
     &          box_iwp(j) .lt. modis_iwp_bin_edges(6)) THEN
                imodis_iwp(j)=5
              ELSE IF (box_iwp(j) .ge. modis_iwp_bin_edges(6) .and. &
     &          box_iwp(j) .lt. modis_iwp_bin_edges(7)) THEN
                imodis_iwp(j)=6
              ELSE IF (box_iwp(j) .ge. modis_iwp_bin_edges(7)) THEN
                imodis_iwp(j)=7
              ENDIF

              IF (box_reffi(j) .ge. modis_reffi_bin_edges(1) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(2)) THEN
                imodis_reffi(j)=1
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(2) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(3)) THEN
                imodis_reffi(j)=2
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(3) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(4)) THEN
                imodis_reffi(j)=3
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(4) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(5)) THEN
                imodis_reffi(j)=4
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(5) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(6)) THEN
                imodis_reffi(j)=5
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(6) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(7)) THEN
                imodis_reffi(j)=6
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(7) .and. &
     &          box_reffi(j) .lt. modis_reffi_bin_edges(8)) THEN
                imodis_reffi(j)=7
              ELSE IF (box_reffi(j) .ge. modis_reffi_bin_edges(8)) THEN
                imodis_reffi(j)=8
              ENDIF

              IF (imodis_iwp(j) .gt. 0 .and. imodis_reffi(j) .gt. 0) THEN
                fq_modis_iwp_reffi(j,imodis_iwp(j),imodis_reffi(j)) = &
     &            fq_modis_iwp_reffi(j,imodis_iwp(j),imodis_reffi(j)) + boxarea
              ENDIF

              !MODIS-style tau-vs-ReffIce joint histogram (BMD 08/26) --
              !see the tau-vs-ReffLiq note in the liquid branch above.
              IF (PRESENT(fq_modis_tau_reffi) .and. &
     &            itau(j) .gt. 0 .and. imodis_reffi(j) .gt. 0) THEN
                fq_modis_tau_reffi(j,itau(j),imodis_reffi(j)) = &
     &            fq_modis_tau_reffi(j,itau(j),imodis_reffi(j)) + boxarea
              ENDIF

              !ice-phase tau/CTP: same reuse pattern as the liquid branch
              !above -- now also INSIDE the box_reffi(j).ge.0._wp guard
              !(BMD 08/26), matching the liquid branch's move.
              IF (PRESENT(fq_modis_tau_ctp) .and. &
     &            imodis_pres(j) .gt. 0 .and. itau(j) .gt. 0) THEN
                fq_modis_tau_ctp_ice(j,itau(j),imodis_pres(j)) = &
     &            fq_modis_tau_ctp_ice(j,itau(j),imodis_pres(j)) + boxarea
              ENDIF
              ENDIF
              !ENDIF above closes the box_reffi(j).ge.0._wp guard, now
              !wrapping IWP-vs-ReffIce, tau-vs-ReffIce, AND tau/CTP-Ice
              !together (BMD 08/26).
            ENDIF
            !ELSE: mixed-phase/undetermined subcolumn, or no phase-resolved
            !tau at all -- excluded from the liquid/ice-specific outputs,
            !matching ICON's phaseIsUndetermined/phaseIsNone handling.
            !(Undetermined subcolumns DO still count toward cltmodis/
            !fq_modis_tau_ctp above if their liquid-library retrieval
            !succeeded -- see the BMD 08/26 note.)

          ENDIF
          ENDIF
! BMD 07/26 - end
        ENDDO ! j
39    CONTINUE
      
      !compute mean cloud properties
      DO j=1,npoints 
        IF (totalcldarea(j) .gt. 0._wp) THEN
          ! code above guarantees that totalcldarea > 0 
          ! only if sunlit .eq. 1 .or. top_height = 3 
          ! and applies only to clouds with tau > isccp_taumin
          meanptop(j) = meanptop(j) / totalcldarea(j)
          meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
          meantaucld(j) = (6.82_wp/((1._wp/meanalbedocld(j))-1._wp))**(1._wp/0.895_wp)
        ELSE
          ! this code is necessary so that in the case that totalcldarea = 0.,
          ! that these variables, which are in-cloud averages, are set to missing
          ! note that totalcldarea will be 0. if all the clouds in the grid box have
          ! tau < isccp_taumin 
          meanptop(j) = output_missing_value
          meanalbedocld(j) = output_missing_value
          meantaucld(j) = output_missing_value
        ENDIF
      ENDDO ! j

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
      IF (debugcol.ne.0) THEN
!     
         DO j=1,npoints,debugcol

            !produce character output
            DO ilev=1,nlev
              DO ibox=1,ncol
                   acc(ilev,ibox)=0
              ENDDO
            ENDDO

            DO ilev=1,nlev
              DO ibox=1,ncol
                   acc(ilev,ibox)=INT(frac_out(j,ibox,ilev)*2._wp)
                   IF (levmatch(j,ibox) .eq. ilev)  &
     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
              ENDDO
            ENDDO

          write(ftn09,11) j
11        format('ftn09.',i4.4)
          open(9, FILE=ftn09, FORM='FORMATTED')
             write(9,'(a1)') ' '
             write(9,'(10i5)') &
     &                  (ilev,ilev=5,nlev,5)
             write(9,'(a1)') ' '
             DO ibox=1,ncol
               write(9,'(40(a1),1x,40(a1))') &
     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) &
     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
             ENDDO
             close(9)

             IF (ncolprint.ne.0) THEN
               write(6,'(a1)') ' '
                    write(6,'(a2,1X,5(a7,1X),a50)') &
     &                  'ilev',&
     &                  'pfull','at',&
     &                  'cc*100','dem_s','dtau_s',&
     &                  'cchar'

!               DO 4012 ilev=1,nlev
!                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
!                   write(6,'(i2,1X,5(f7.2,1X),50(a1))') 
!     &                  ilev,
!     &                  pfull(j,ilev)/100.,at(j,ilev),
!     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
!     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
!4012           continue
               write (6,'(a)') 'skt(j):'
               write (6,'(8f7.2)') skt(j)
               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
               write (6,'(a)') 'tau:'
               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
               write (6,'(a)') 'tb:'
               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
               write (6,'(a)') 'ptop:'
               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
             ENDIF 
        ENDDO
      ENDIF

! BMD 07/26 - start
!     Simulated-reflectance MODIS size-retrieval procedures, ported from
!     ~/COSPv2.0/src/simulator/MODIS_simulator/modis_simulator.F90.
!     Internal (CONTAINS'd) rather than separate external procedures: some
!     of these are ELEMENTAL and called with array arguments (the 15-point
!     trial-size grid) from the main body, which requires an explicit
!     interface visible at the call site -- host association gives that
!     automatically, without a separate INTERFACE block. IMPLICIT NONE is
!     inherited from the host subroutine, not repeated below. Some local
!     variable names below (tau, g, w0, f, beta, den, rk, th, e1, e2, ...)
!     are ported verbatim from modis_simulator.F90 and are NOT the same as
!     this subroutine's own host-scope variables of similar names (e.g.
!     the host's tau(npoints,ncol)) -- Fortran scoping shadows these
!     correctly within each internal procedure's own body, but readers
!     should not assume any connection to the host's tau array.
      CONTAINS

      PURE FUNCTION fit_to_cubic(x, coefficients) RESULT(y)
      REAL(wp), INTENT(IN) :: x
      REAL(wp), INTENT(IN) :: coefficients(4)
      REAL(wp) :: y
      y = coefficients(1) + x*(coefficients(2) &
     &  + x*(coefficients(3) + x*coefficients(4)))
      END FUNCTION fit_to_cubic

      PURE FUNCTION fit_to_quadratic(x, coefficients) RESULT(y)
      REAL(wp), INTENT(IN) :: x
      REAL(wp), INTENT(IN) :: coefficients(3)
      REAL(wp) :: y
      y = coefficients(1) + x*(coefficients(2) + x*coefficients(3))
      END FUNCTION fit_to_quadratic

!     Polynomial fits for near-IR (3.7 micron, Collection-6) asymmetry
!     parameter (get_g_nir) and single-scattering albedo (get_ssa_nir) as
!     functions of particle size, ported verbatim (fits from Steve
!     Platnick). re is in MICRONS. Adapted: takes a LOGICAL is_liquid
!     instead of an INTEGER phase enum -- this file already classifies
!     phase via liquid_fraction/modis_phase_threshold directly, so
!     phaseIsLiquid/phaseIsIce/etc. were not ported.
      ELEMENTAL FUNCTION get_g_nir(is_liquid, re) RESULT(g_nir)
      LOGICAL, INTENT(IN) :: is_liquid
      REAL(wp), INTENT(IN) :: re
      REAL(wp) :: g_nir
      REAL(wp), PARAMETER :: ice_coefficients(3) = &
     &  (/ 0.7490_wp, 6.5153e-3_wp, -5.4136e-5_wp /)
      REAL(wp), PARAMETER :: small_water_coefficients(3) = &
     &  (/ 1.0364_wp, -8.8800e-2_wp, 7.0000e-3_wp /)
      REAL(wp), PARAMETER :: big_water_coefficients(4) = &
     &  (/ 0.6035_wp, 2.8993e-2_wp, -1.1051e-3_wp, 1.5134e-5_wp /)

      IF (is_liquid) THEN
         IF (re .lt. 7._wp) THEN
            g_nir = fit_to_quadratic(re, small_water_coefficients)
            IF (re .lt. modis_trial_re_water_min) g_nir = &
     &         fit_to_quadratic(modis_trial_re_water_min, small_water_coefficients)
         ELSE
            g_nir = fit_to_cubic(re, big_water_coefficients)
            IF (re .gt. modis_trial_re_water_max) g_nir = &
     &         fit_to_cubic(modis_trial_re_water_max, big_water_coefficients)
         ENDIF
      ELSE
         g_nir = fit_to_quadratic(re, ice_coefficients)
         IF (re .lt. modis_trial_re_ice_min) g_nir = &
     &      fit_to_quadratic(modis_trial_re_ice_min, ice_coefficients)
         IF (re .gt. modis_trial_re_ice_max) g_nir = &
     &      fit_to_quadratic(modis_trial_re_ice_max, ice_coefficients)
      ENDIF
      END FUNCTION get_g_nir

      ELEMENTAL FUNCTION get_ssa_nir(is_liquid, re) RESULT(ssa_nir)
      LOGICAL, INTENT(IN) :: is_liquid
      REAL(wp), INTENT(IN) :: re
      REAL(wp) :: ssa_nir
      REAL(wp), PARAMETER :: ice_coefficients(4) = &
     &  (/ 0.9625_wp, -1.8069e-2_wp, 3.3281e-4_wp, -2.2865e-6_wp /)
      REAL(wp), PARAMETER :: water_coefficients(3) = &
     &  (/ 1.0044_wp, -1.1397e-2_wp, 1.3300e-4_wp /)

      IF (is_liquid) THEN
         ssa_nir = fit_to_quadratic(re, water_coefficients)
         IF (re .lt. modis_trial_re_water_min) ssa_nir = &
     &      fit_to_quadratic(modis_trial_re_water_min, water_coefficients)
         IF (re .gt. modis_trial_re_water_max) ssa_nir = &
     &      fit_to_quadratic(modis_trial_re_water_max, water_coefficients)
      ELSE
         ssa_nir = fit_to_cubic(re, ice_coefficients)
         IF (re .lt. modis_trial_re_ice_min) ssa_nir = &
     &      fit_to_cubic(modis_trial_re_ice_min, ice_coefficients)
         IF (re .gt. modis_trial_re_ice_max) ssa_nir = &
     &      fit_to_cubic(modis_trial_re_ice_max, ice_coefficients)
      ENDIF
      END FUNCTION get_ssa_nir

!     Single-aggregate-layer reflectance via the delta-Eddington two-stream
!     approximation (Joseph et al. 1976; code from Lazaros Oreopoulos via
!     Steve Platnick), ported verbatim. Used for the trial-size PREDICTED
!     reflectance curve (whole subcolumn treated as one layer with its
!     total tau) -- NOT the same code path as two_stream/adding_doubling
!     below, which build the multi-level OBSERVED reflectance from the
!     model's real per-level vertical structure. This split (simplified
!     single-layer retrieval assumption vs. real multi-layer truth) is
!     itself part of the genuine retrieval algorithm being simulated, not
!     something to reconcile.
      ELEMENTAL FUNCTION two_stream_reflectance(tauint, gint, w0int) RESULT(refl)
      REAL(wp), INTENT(IN) :: tauint, gint, w0int
      REAL(wp) :: refl
      INTEGER, PARAMETER :: beam = 2
      REAL(wp), PARAMETER :: xmu = 0.866_wp, minConservativeW0 = 0.9999999_wp
      REAL(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
     &  a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, &
     &  beta, e1, e2, ef1, ef2

      f   = gint**2
      tau = (1._wp - w0int*f) * tauint
      w0  = (1._wp - f) * w0int / (1._wp - w0int*f)
      g   = (gint - f) / (1._wp - f)

      gamma1 =  (7._wp - w0*(4._wp + 3._wp*g)) / 4._wp
      gamma2 = -(1._wp - w0*(4._wp - 3._wp*g)) / 4._wp
      gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
      gamma4 =   1._wp - gamma3

      IF (w0int .gt. minConservativeW0) THEN
         !Conservative scattering (beam=2 only, matches modis_simulator.F90)
         refl = gamma1*tau / (1._wp + gamma1*tau)
      ELSE
         !Non-conservative scattering
         a1 = gamma1*gamma4 + gamma2*gamma3
         a2 = gamma1*gamma3 + gamma2*gamma4
         rk = sqrt(gamma1**2 - gamma2**2)

         r1 = (1._wp - rk*xmu) * (a2 + rk*gamma3)
         r2 = (1._wp + rk*xmu) * (a2 - rk*gamma3)
         r3 = 2._wp*rk * (gamma3 - a2*xmu)
         r4 = (1._wp - (rk*xmu)**2) * (rk + gamma1)
         r5 = (1._wp - (rk*xmu)**2) * (rk - gamma1)

         t1 = (1._wp + rk*xmu) * (a1 + rk*gamma4)
         t2 = (1._wp - rk*xmu) * (a1 - rk*gamma4)
         t3 = 2._wp*rk * (gamma4 + a1*xmu)
         t4 = r4
         t5 = r5

         beta = -r5/r4

         e1 = min(rk*tau, 500._wp)
         e2 = min(tau/xmu, 500._wp)

         ef1 = exp(-e1)
         ef2 = exp(-2*e1)
         refl = (gamma2*(1._wp-ef2)) / ((rk+gamma1)*(1._wp-beta*ef2))
      ENDIF
      END FUNCTION two_stream_reflectance

!     Per-layer reflectance/transmittance, same delta-Eddington math as
!     two_stream_reflectance above but returning both ref and tra (needed
!     by adding_doubling below to combine multiple real layers) --
!     ported verbatim.
      PURE SUBROUTINE two_stream(tauint, gint, w0int, ref, tra)
      REAL(wp), INTENT(IN) :: tauint, gint, w0int
      REAL(wp), INTENT(OUT) :: ref, tra
      INTEGER, PARAMETER :: beam = 2
      REAL(wp), PARAMETER :: xmu = 0.866_wp, minConservativeW0 = 0.9999999_wp
      REAL(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
     &  a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, &
     &  beta, e1, e2, ef1, ef2

      f   = gint**2
      tau = (1._wp - w0int*f) * tauint
      w0  = (1._wp - f) * w0int / (1._wp - w0int*f)
      g   = (gint - f) / (1._wp - f)

      gamma1 =  (7._wp - w0*(4._wp + 3._wp*g)) / 4._wp
      gamma2 = -(1._wp - w0*(4._wp - 3._wp*g)) / 4._wp
      gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
      gamma4 =   1._wp - gamma3

      IF (w0int .gt. minConservativeW0) THEN
         !Conservative scattering (beam=2 only)
         ref = gamma1*tau / (1._wp + gamma1*tau)
         tra = 1._wp - ref
      ELSE
         !Non-conservative scattering
         a1 = gamma1*gamma4 + gamma2*gamma3
         a2 = gamma1*gamma3 + gamma2*gamma4
         rk = sqrt(gamma1**2 - gamma2**2)

         r1 = (1._wp - rk*xmu) * (a2 + rk*gamma3)
         r2 = (1._wp + rk*xmu) * (a2 - rk*gamma3)
         r3 = 2._wp*rk * (gamma3 - a2*xmu)
         r4 = (1._wp - (rk*xmu)**2) * (rk + gamma1)
         r5 = (1._wp - (rk*xmu)**2) * (rk - gamma1)

         t1 = (1._wp + rk*xmu) * (a1 + rk*gamma4)
         t2 = (1._wp - rk*xmu) * (a1 - rk*gamma4)
         t3 = 2._wp*rk * (gamma4 + a1*xmu)
         t4 = r4
         t5 = r5

         beta = -r5/r4

         e1 = min(rk*tau, 500._wp)
         e2 = min(tau/xmu, 500._wp)

         ef1 = exp(-e1)
         ef2 = exp(-2*e1)
         ref = (gamma2*(1._wp-ef2)) / ((rk+gamma1)*(1._wp-beta*ef2))
         tra = (2._wp*rk*ef1) / ((rk+gamma1)*(1._wp-beta*ef2))
      ENDIF
      END SUBROUTINE two_stream

!     Combine per-layer (Refl,Tran) pairs into a single stack (Refl_tot,
!     Tran_tot) via the adding-layers formula, sequentially from index 1
!     (TOA-most, matching this file's own pfull/phalf TOA-to-surface
!     convention) downward -- ported verbatim.
      PURE SUBROUTINE adding_doubling(npts, Refl, Tran, Refl_tot, Tran_tot)
      INTEGER, INTENT(IN) :: npts
      REAL(wp), INTENT(IN) :: Refl(npts), Tran(npts)
      REAL(wp), INTENT(OUT) :: Refl_tot, Tran_tot
      INTEGER :: i
      REAL(wp) :: Refl_cumulative(npts), Tran_cumulative(npts)

      Refl_cumulative(1) = Refl(1)
      Tran_cumulative(1) = Tran(1)

      DO i=2,npts
         Refl_cumulative(i) = Refl_cumulative(i-1) + Refl(i) &
     &     * (Tran_cumulative(i-1)**2) &
     &     / (1._wp - Refl_cumulative(i-1)*Refl(i))
         Tran_cumulative(i) = (Tran_cumulative(i-1)*Tran(i)) &
     &     / (1._wp - Refl_cumulative(i-1)*Refl(i))
      ENDDO

      Refl_tot = Refl_cumulative(npts)
      Tran_tot = Tran_cumulative(npts)
      END SUBROUTINE adding_doubling

!     Given a predicted-reflectance-vs-trial-size curve (y as a function
!     of x), find the trial size whose predicted reflectance brackets the
!     observed reflectance yobs, linearly interpolating between the two
!     bracketing trial points -- ported verbatim. Returns modis_re_fill
!     (micron) if no bracket is found (the observed reflectance falls
!     outside what any trial size in [min,max] can produce): this is the
!     actual mechanism real MODIS retrievals go missing by, NOT a direct
!     comparison of some "true" size against modis_trial_re_water_min/max
!     (see the note on this distinction from the design discussion this
!     was ported from).
      PURE FUNCTION interpolate_to_min(x, y, yobs) RESULT(x_min)
      REAL(wp), INTENT(IN) :: x(num_trial_res), y(num_trial_res)
      REAL(wp), INTENT(IN) :: yobs
      REAL(wp) :: x_min
      REAL(wp) :: diff(num_trial_res)
      INTEGER :: minDiffLoc, lowerBound, upperBound

      diff(1:num_trial_res) = y(1:num_trial_res) - yobs
      minDiffLoc = minloc(abs(diff), dim=1)

      IF (minDiffLoc .eq. 1) THEN
         lowerBound = minDiffLoc
         upperBound = minDiffLoc + 1
      ELSE IF (minDiffLoc .eq. num_trial_res) THEN
         lowerBound = minDiffLoc - 1
         upperBound = minDiffLoc
      ELSE
         IF (diff(minDiffLoc-1)*diff(minDiffLoc) .lt. 0._wp) THEN
            lowerBound = minDiffLoc - 1
            upperBound = minDiffLoc
         ELSE
            lowerBound = minDiffLoc
            upperBound = minDiffLoc + 1
         ENDIF
      ENDIF

      IF (diff(lowerBound)*diff(upperBound) .lt. 0._wp) THEN
         !Root is bracketed -- interpolate linearly
         x_min = x(upperBound) - diff(upperBound) &
     &     * (x(upperBound) - x(lowerBound)) &
     &     / (diff(upperBound) - diff(lowerBound))
      ELSE
         x_min = modis_re_fill
      ENDIF
      END FUNCTION interpolate_to_min
! BMD 07/26 - end

      END SUBROUTINE icarus 
