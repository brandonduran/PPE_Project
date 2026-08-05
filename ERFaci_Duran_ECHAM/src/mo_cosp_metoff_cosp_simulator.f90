!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! (c) COPYRIGHT British Crown / Met Office 2008  
! Please refer to Met_Office_license.txt for details. 
MODULE mo_cosp_metoff_cosp_simulator
  USE mo_cosp_lidar
  USE mo_cosp_stats
  USE mo_cosp_isccp_simulator
! BMD 07/26 - start
  USE mo_cosp_modis_simulator
! BMD 07/26 - end
  USE mo_kind, only:wp

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_SIMULATOR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE cosp_simulator(                         &
       kbdim,Ncolumns,klev,                         &
       Lisccp_sim,  Llidar_sim,  Llidar_cfad, Lstats,            &
! BMD 07/26 - start
       Lmodis_sim,                                  &
! BMD 07/26 - end
!gbx
       lidar_ice_type,                               &
       isccp_top_height, isccp_top_height_direction, &
       isccp_overlap, isccp_emsfc_lw,               &
       ph, p, T, Reff, tca,                         &             
       sh, dem_s, dtau_s,                           &
       sunlit, skt, zlev, zlev_half, land,          &
!sgx
       frac_out,                                    & 
       sg_mr_hydro,                                 & 
!sglidar
       beta_mol, beta_tot, refl, tau_tot,           &
!stlidar
       srbval, cfad_sr, lidarcld,                   &
       cldlayer, parasolrefl,                       &
!isccp
       fq_isccp,  totalcldarea,                     &
       meanptop,  meantaucld,                       &
       meantb, meantbclr, boxtau,  boxptop,         &
       meanalbedocld,                                &
! BMD 07/26 - start
!modis
       dtau_s_liq, dtau_s_ice,                      &
       fq_modis_lwp_reffl, fq_modis_iwp_reffi,      &
       cltmodis, clwmodis, climodis,                &
       fq_modis_tau_ctp, fq_modis_tau_ctp_liq,      &
       fq_modis_tau_ctp_ice,                        &
       dpiz_s_nir, dcg_s_nir                        &
! BMD 07/26 - end
                                                     )

  ! Arguments
  INTEGER, INTENT(IN) :: kbdim,Ncolumns,klev
  LOGICAL, INTENT(IN) ::  Lisccp_sim
  LOGICAL, INTENT(IN) ::  Llidar_sim
  LOGICAL, INTENT(IN) ::  Llidar_cfad
  LOGICAL, INTENT(IN) ::  Lstats
! BMD 07/26 - start
  LOGICAL, INTENT(IN) ::  Lmodis_sim
! BMD 07/26 - end


!gbx

  INTEGER, INTENT(IN) :: lidar_ice_type
  INTEGER, INTENT(IN) :: isccp_top_height
  INTEGER, INTENT(IN) :: isccp_top_height_direction
  INTEGER, INTENT(IN) :: isccp_overlap

  REAL(wp), INTENT(IN) :: isccp_emsfc_lw 

  ! check intent statements !
  REAL(wp), INTENT(INOUT) :: p(kbdim,klev) 
  REAL(wp), INTENT(INOUT) :: ph(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: T(kbdim,klev) 
  REAL(wp), INTENT(INOUT) :: Reff(kbdim,klev,N_hydro)
  REAL(wp), INTENT(INOUT) :: tca(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: sh(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: dem_s(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: dtau_s(kbdim,klev)
! BMD 07/26 - start
!  liquid/ice-split counterpart of dtau_s, needed by cosp_modis_simulator's
!  LWP/IWP-vs-Reff and tau/CTP-Liq/Ice histograms. Sourced in
!  mo_cosp_simulator.f90's call_cospsimulator from cisccp_cldtau3d_liq/_ice
!  (in turn from cld_tau_swl_vr/cld_tau_swi_vr in mo_psrad_interface.f90),
!  passed straight through this pass-through layer.
  REAL(wp), INTENT(INOUT) :: dtau_s_liq(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: dtau_s_ice(kbdim,klev)
! BMD 07/26 - end
  REAL(wp), INTENT(INOUT) :: sunlit(kbdim)
  REAL(wp), INTENT(INOUT) :: skt(kbdim)
  REAL(wp), INTENT(INOUT) :: zlev(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: zlev_half(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: land(kbdim)

!sglidar
  REAL(wp) :: beta_mol(kbdim,klev)
  REAL(wp) :: beta_tot(kbdim,Ncolumns,klev)
  REAL(wp) :: refl(kbdim,Ncolumns,PARASOL_NREFL)
  REAL(wp) :: tau_tot(kbdim,Ncolumns,klev)

! stlidar lidar stats
  REAL(wp), INTENT(INOUT) :: srbval(SR_BINS)
  REAL(wp), INTENT(INOUT) :: cfad_sr(kbdim,SR_BINS,Nlr) 
  REAL(wp), INTENT(INOUT) :: lidarcld(kbdim,Nlr)
  REAL(wp), INTENT(INOUT) :: cldlayer(kbdim,LIDAR_NCAT)
  REAL(wp), INTENT(INOUT) :: parasolrefl(kbdim,PARASOL_NREFL)

! isccp
  REAL(wp) :: fq_isccp(kbdim,7,7), totalcldarea(kbdim)
  REAL(wp) :: meanptop(kbdim), meantaucld(kbdim)
  REAL(wp) :: meantb(kbdim), meantbclr(kbdim)
  REAL(wp) :: boxtau(kbdim,Ncolumns), boxptop(kbdim,Ncolumns)
  REAL(wp) :: meanalbedocld(kbdim)

! BMD 07/26 - start
! modis
!  bin counts (7,8 / 7,8) must match n_modis_lwp/n_modis_reffl/etc. declared
!  in cosp_icarus.f90 -- see the note in mo_cosp_modis_simulator.f90
  REAL(wp) :: fq_modis_lwp_reffl(kbdim,7,8)
  REAL(wp) :: fq_modis_iwp_reffi(kbdim,7,8)
!  total/liquid/ice cloud fraction (fraction of sunlit subcolumns)
  REAL(wp) :: cltmodis(kbdim)
  REAL(wp) :: clwmodis(kbdim)
  REAL(wp) :: climodis(kbdim)
!  joint histograms of optical thickness vs. cloud-top pressure (all
!  phases, liquid-only, ice-only) -- bin counts (7,7) must match
!  n_modis_tau/n_modis_pres declared in cosp_icarus.f90
  REAL(wp) :: fq_modis_tau_ctp(kbdim,7,7)
  REAL(wp) :: fq_modis_tau_ctp_liq(kbdim,7,7)
  REAL(wp) :: fq_modis_tau_ctp_ice(kbdim,7,7)
!  near-IR (3.7 micron) single-scattering albedo/asymmetry, passed straight
!  through to cosp_modis_simulator (ISCCP doesn't need these)
  REAL(wp) :: dpiz_s_nir(kbdim,klev)
  REAL(wp) :: dcg_s_nir(kbdim,klev)
! BMD 07/26 - end

! subgrid sgx
  REAL(wp) :: frac_out(kbdim,Ncolumns,klev)

! sghydro
  REAL(wp) :: sg_mr_hydro(kbdim,Ncolumns,klev,N_hydro)

  !+++++++++ Lidar model ++++++++++
  IF (Llidar_sim) THEN
    CALL cosp_lidar(                     &
         kbdim,Ncolumns,klev,            &
!gbx
         lidar_ice_type, Reff,           &
         ph,  sg_mr_hydro, p, T,        &
!sgx, sglidar   
!         frac_out,                       &
         beta_mol, beta_tot, tau_tot, refl)
  ENDIF

  !+++++++++ ISCCP simulator ++++++++++
  IF (Lisccp_sim) THEN
    CALL cosp_isccp_simulator(                     &
      kbdim,  klev, Ncolumns,                      & 
!gbx
      isccp_top_height, isccp_top_height_direction,&
      isccp_overlap, sh, tca, dtau_s, T, dem_s,    &
      p, ph , sunlit, skt, isccp_emsfc_lw,         &
!sgx
      frac_out,                                    &
!isccp
      fq_isccp,  totalcldarea,                     &
      meanptop,  meantaucld,                       &
      meantb, meantbclr, boxtau,  boxptop,         &
      meanalbedocld )
  ENDIF

! BMD 07/26 - start
  !+++++++++ MODIS simulator ++++++++++
  ! cosp_modis_simulator only returns the MODIS-specific joint histograms
  ! (fq_modis_lwp_reffl/fq_modis_iwp_reffi and, BMD 07/26 second pass,
  ! fq_modis_tau_ctp/_liq/_ice). It still runs the full icarus
  ! engine internally (tau/CTP is computed as part of the same subcolumn pass
  ! that produces the MODIS histograms), but the ISCCP-shaped outputs from
  ! that internal call are local/discarded there, not threaded back through
  ! to this subroutine's real fq_isccp/totalcldarea/etc. -- so the ISCCP
  ! tau/CTP histogram itself is only ever produced once, by
  ! cosp_isccp_simulator above, regardless of whether Lisccp_sim and
  ! Lmodis_sim are both on.
  !
  ! reffl/reffi are sliced directly from the Reff(:,:,N_hydro) array already
  ! threaded to this subroutine for the lidar path (I_LSCLIQ/I_LSCICE from
  ! mo_cosp_constants) -- no new dummy argument needed for those.
  IF (Lmodis_sim) THEN
    CALL cosp_modis_simulator(                     &
      kbdim,  klev, Ncolumns,                      &
!gbx
      isccp_top_height, isccp_top_height_direction,&
      isccp_overlap, sh, tca, dtau_s, T, dem_s,    &
      p, ph , sunlit, skt, isccp_emsfc_lw,         &
!sgx
      frac_out,                                    &
!modis
      dtau_s_liq, dtau_s_ice,                      &
      Reff(:,:,I_LSCLIQ), Reff(:,:,I_LSCICE),      &
      fq_modis_lwp_reffl, fq_modis_iwp_reffi,      &
      cltmodis, clwmodis, climodis,                &
      fq_modis_tau_ctp, fq_modis_tau_ctp_liq,      &
      fq_modis_tau_ctp_ice,                        &
      dpiz_s_nir, dcg_s_nir )
  ENDIF
! BMD 07/26 - end

  !+++++++++++ Summary statistics +++++++++++
  IF (Lstats) THEN
    CALL cosp_stats(kbdim,Ncolumns,klev,       &
      Llidar_sim, Llidar_cfad,                 &
!gbx
      ph, zlev, zlev_half, land,               &
!sglidar
      beta_mol, beta_tot, refl,                &
!stlidar
      srbval, cfad_sr, lidarcld,               &
      cldlayer, parasolrefl                    )
  ENDIF

  !+++++++++++ change of units after computation of statistics +++++++++++
!!  if (cfg%Llidar_sim) then
!!    where((sglidar%beta_tot > 0.0_wp) .and. (sglidar%beta_tot /= R_UNDEF)) 
!!        sglidar%beta_tot = log10(sglidar%beta_tot)
!!    elsewhere
!!        sglidar%beta_tot = R_UNDEF
!!    end where
!!  endif

END SUBROUTINE cosp_simulator

END MODULE mo_cosp_metoff_cosp_simulator
