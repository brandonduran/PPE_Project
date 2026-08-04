!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are permitted
! provided that the following conditions are met:
!
!     * Redistributions of source code must retain the above copyright notice, this list
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used
!       to endorse or promote products derived from this software without specific prior written
!       permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! BMD 07/26 - start
! Cloned from mo_cosp_isccp_simulator.f90 as the starting point for a MODIS-style
! simulator wrapper: same ECHAM-grid <-> icarus level-flip pattern, extended with
! the liquid/ice-split optical depth, liquid/ice effective radius inputs and the
! LWP-vs-ReffLiq / IWP-vs-ReffIce / tau-vs-CTP (all phases, liquid, ice)
! joint-histogram outputs added to the icarus engine in cosp_icarus.f90. Wired
! end-to-end: the Lmodis_sim namelist switch and stream-field declarations are
! in mo_cosp_simulator.f90, dtau_s_liq/dtau_s_ice are sourced there from
! cisccp_cldtau3d_liq/_ice (in turn from mo_psrad_interface.f90's
! cld_tau_swl_vr/cld_tau_swi_vr), and mo_cosp_metoff_cosp.f90 passes
! everything through. Still needs end-to-end run/output validation.
!
! The underlying icarus engine is monolithic: it always computes the ISCCP
! tau/CTP histogram (fq_isccp) and its associated cloud-top-pressure/mean
! quantities as part of computing per-subcolumn tau, regardless of whether the
! caller wants them. Rather than exposing those ISCCP-shaped outputs as dummy
! arguments here -- which would mean this subroutine could be handed the
! caller's real fq_isccp/totalcldarea/etc. arrays and rewrite them a second
! time whenever both Lisccp_sim and Lmodis_sim are on, redundantly duplicating
! the ISCCP tau/CTP calculation already done by cosp_isccp_simulator -- this
! subroutine keeps them as local, discarded work arrays. cosp_modis_simulator's
! own interface therefore only exposes the MODIS-specific joint histograms
! (LWP-vs-ReffLiq, IWP-vs-ReffIce, and -- BMD 07/26, second pass -- tau/CTP,
! tau/CTP-Liq, tau/CTP-Ice, this last group its own genuinely new retrieval,
! CO2-slicing rather than ISCCP's brightness-temperature CTP, see
! cosp_icarus.f90); the ISCCP tau/CTP histogram itself is computed exactly
! once, in cosp_isccp_simulator, as it should be.
! BMD 07/26 - end

MODULE mo_cosp_modis_simulator

  USE mo_cosp_constants
  USE mo_kind, only:wp

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------- SUBROUTINE COSP_MODIS_SIMULATOR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE cosp_modis_simulator(                          &
             kbdim,  klev,  Ncolumns,                     &
!gbx
             isccp_top_height, isccp_top_height_direction, &
             isccp_overlap, sh, tca, dtau_s, T, dem_s,    &
             p, ph, sunlit, skt, isccp_emsfc_lw,          &
!sgx
             frac_out,                                    &
!modis
             dtau_s_liq, dtau_s_ice, reffl, reffi,        &
             fq_modis_lwp_reffl, fq_modis_iwp_reffi,      &
             cltmodis, clwmodis, climodis,                &
             fq_modis_tau_ctp, fq_modis_tau_ctp_liq,      &
             fq_modis_tau_ctp_ice                         &
             )

  ! Arguments
   INTEGER, INTENT(IN) :: kbdim, klev, Ncolumns

!gbx
   INTEGER, INTENT(IN) :: isccp_top_height
   INTEGER, INTENT(IN) :: isccp_top_height_direction
   INTEGER, INTENT(IN) :: isccp_overlap

   REAL(wp), INTENT(INOUT) :: p(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: ph(kbdim,klev)
   REAL(wp), INTENT(IN)    :: isccp_emsfc_lw
   REAL(wp), INTENT(INOUT) :: sh(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: tca(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: dem_s(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: dtau_s(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: sunlit(kbdim)
   REAL(wp), INTENT(INOUT) :: T(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: skt(kbdim)

! subgrid sgx
   REAL(wp), INTENT(IN) :: frac_out(kbdim,Ncolumns,klev)

! modis
!  liquid/ice-split counterpart of dtau_s, and liquid/ice effective radius (m).
!  dtau_s_liq/dtau_s_ice are the phase-partitioned equivalent of the
!  0.67 micron cld_tau_sw band used to build dtau_s (see mo_psrad_interface.f90,
!  cld_tau_swl_vr/cld_tau_swi_vr); reffl/reffi are the same liquid/ice effective
!  radii already threaded through mo_cosp_simulator.f90 as Reff(:,:,I_LSCLIQ/ICE)
!  for the lidar path (in turn derived from cosp_reffl/cosp_reffi), just not yet
!  passed down to the ISCCP/MODIS call site.
   REAL(wp), INTENT(INOUT) :: dtau_s_liq(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: dtau_s_ice(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: reffl(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: reffi(kbdim,klev)

!  joint-histogram outputs. Bin counts (7 LWP x 8 ReffLiq, 7 IWP x 8 ReffIce)
!  must match n_modis_lwp/n_modis_iwp/n_modis_reffl/n_modis_reffi declared in
!  cosp_icarus.f90 -- hardcoded here as literals rather than shared module
!  parameters, matching how fq_isccp(kbdim,7,7) in mo_cosp_isccp_simulator.f90
!  hardcodes the ISCCP 7x7 bin count rather than referencing a shared constant.
   REAL(wp), INTENT(INOUT) :: fq_modis_lwp_reffl(kbdim,7,8)
   REAL(wp), INTENT(INOUT) :: fq_modis_iwp_reffi(kbdim,7,8)

!  1D total/liquid/ice cloud fraction (fraction of sunlit subcolumns), the
!  simplest MODIS diagnostics to sanity-check against obs -- see the naming
!  convention (cltmodis/clwmodis/climodis) confirmed against the CESM/ICON
!  driver cospsimulator_intr.F90.
   REAL(wp), INTENT(INOUT) :: cltmodis(kbdim)
   REAL(wp), INTENT(INOUT) :: clwmodis(kbdim)
   REAL(wp), INTENT(INOUT) :: climodis(kbdim)

!  joint histograms of optical thickness vs. cloud-top pressure (all
!  phases, liquid-only, ice-only). Bin counts (7 tau x 7 CTP) must match
!  n_modis_tau/n_modis_pres declared in cosp_icarus.f90 -- same
!  hardcoded-literal convention as fq_modis_lwp_reffl/fq_modis_iwp_reffi
!  above. cloud-top pressure here is COSPv2's own CO2-slicing retrieval
!  (see cosp_icarus.f90), NOT ISCCP's brightness-temperature ptop, though
!  it falls back to it for low clouds.
   REAL(wp), INTENT(INOUT) :: fq_modis_tau_ctp(kbdim,7,7)
   REAL(wp), INTENT(INOUT) :: fq_modis_tau_ctp_liq(kbdim,7,7)
   REAL(wp), INTENT(INOUT) :: fq_modis_tau_ctp_ice(kbdim,7,7)

  ! Local variables
  REAL(wp) :: pfull(kbdim, klev)
  REAL(wp) :: phalf(kbdim, klev + 1)
  REAL(wp) :: qv(kbdim, klev)
  REAL(wp) :: cc(kbdim, klev)
  REAL(wp) :: dtau_s_l(kbdim, klev)
  REAL(wp) :: at(kbdim, klev)
  REAL(wp) :: dem_s_l(kbdim, klev)
  REAL(wp) :: frac_out_l(kbdim, Ncolumns, klev)
  INTEGER :: sunlit_l(kbdim)

  REAL(wp) :: dtau_s_liq_l(kbdim, klev)
  REAL(wp) :: dtau_s_ice_l(kbdim, klev)
  REAL(wp) :: reffl_l(kbdim, klev)
  REAL(wp) :: reffi_l(kbdim, klev)

! icarus always computes the ISCCP tau/CTP histogram and its associated
! cloud-top-pressure/mean quantities as a side effect of computing per-
! subcolumn tau -- these are discarded local work arrays, not returned to the
! caller, so the real ISCCP output (owned by cosp_isccp_simulator) is computed
! exactly once rather than twice. See the module-level comment above.
  REAL(wp) :: fq_isccp_unused(kbdim,7,7), totalcldarea_unused(kbdim)
  REAL(wp) :: meanptop_unused(kbdim), meantaucld_unused(kbdim)
  REAL(wp) :: meantb_unused(kbdim), meantbclr_unused(kbdim)
  REAL(wp) :: boxtau_unused(kbdim,Ncolumns), boxptop_unused(kbdim,Ncolumns)
  REAL(wp) :: meanalbedocld_unused(kbdim)

  ! Flip inputs. Levels from TOA to surface
  pfull  = p(:,klev:1:-1)
  phalf(:,1)         = 0.0_wp ! Top level
  phalf(:,2:klev+1) = ph(:,klev:1:-1)
  qv     = sh(:,klev:1:-1)
  cc     = 0.999999_wp*tca(:,klev:1:-1)
  dtau_s_l = dtau_s(:,klev:1:-1)
  at     = T(:,klev:1:-1)
  dem_s_l  = dem_s(:,klev:1:-1)
  frac_out_l(1:Kbdim,:,1:klev) = frac_out(1:Kbdim,:,klev:1:-1)
  sunlit_l = int(sunlit)

  dtau_s_liq_l = dtau_s_liq(:,klev:1:-1)
  dtau_s_ice_l = dtau_s_ice(:,klev:1:-1)
  reffl_l      = reffl(:,klev:1:-1)
  reffi_l      = reffi(:,klev:1:-1)

  CALL icarus(0,0,kbdim,sunlit_l,klev,ncolumns, &
            pfull,phalf,qv,cc,dtau_s_l, &
            isccp_top_height,isccp_top_height_direction, &
            isccp_overlap,frac_out_l, &
            skt,isccp_emsfc_lw,at,dem_s_l, &
            fq_isccp_unused,totalcldarea_unused, &
            meanptop_unused,meantaucld_unused,meanalbedocld_unused, &
            meantb_unused,meantbclr_unused,boxtau_unused,boxptop_unused, &
            dtau_s_liq_l,dtau_s_ice_l,reffl_l,reffi_l, &
            fq_modis_lwp_reffl,fq_modis_iwp_reffi, &
            cltmodis,clwmodis,climodis, &
            fq_modis_tau_ctp,fq_modis_tau_ctp_liq,fq_modis_tau_ctp_ice)

  ! No output flip needed for the LWP/IWP/ReffLiq/ReffIce histograms: unlike
  ! fq_isccp's pressure axis (see cosp_isccp_simulator, which owns that
  ! flip), those bin axes are physical-value thresholds with no dependence
  ! on model level ordering, so icarus's bin index already matches the
  ! intended output convention directly.
  !
  ! fq_modis_tau_ctp/_liq/_ice DO need the same pressure-axis flip as
  ! fq_isccp, though: imodis_pres (cosp_icarus.f90) is binned with the
  ! identical edges/index convention as ISCCP's own ipres, so it inherits
  ! the same "index 1 = low pressure" vs. the model's expected output
  ! ordering mismatch that fq_isccp's own flip below corrects. tau (dim 2)
  ! is untouched, matching fq_isccp's flip pattern exactly.
  fq_modis_tau_ctp(:,:,:)     = fq_modis_tau_ctp(:,:,7:1:-1)
  fq_modis_tau_ctp_liq(:,:,:) = fq_modis_tau_ctp_liq(:,:,7:1:-1)
  fq_modis_tau_ctp_ice(:,:,:) = fq_modis_tau_ctp_ice(:,:,7:1:-1)

END SUBROUTINE cosp_modis_simulator

END MODULE mo_cosp_modis_simulator
