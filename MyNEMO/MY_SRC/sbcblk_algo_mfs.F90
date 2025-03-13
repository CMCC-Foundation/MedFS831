MODULE sbcblk_algo_mfs
   !!======================================================================
   !!                       ***  MODULE  sbcblk_algo_mfs  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  3.3  !   2010-05 (P. Oddo) Original Code
   !!            4.0  !   2020-04 (S. Causio, S.Ciliberti, from E. Clementi)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_mfs  : bulk formulation as ocean surface boundary condition
   !!                   (forced mode, mfs bulk formulae)
   !!   blk_oce_mfs  : ocean: computes momentum, heat and freshwater fluxes
   !!   turb_cd_2z   : Computes iturb surface drag at 2m
   !!   psi_m_mfs           : universal profile stability function for momentum
   !!   psi_h_mfs           : universal profile stability function for temperature
   !and humidity
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE fldread   ! read input fields
   USE sbcwave, ONLY   :  cdn_wave ! wave module
#if defined key_si3 || defined key_cice
   USE sbc_ice         ! Surface boundary condition: ice fields
#endif
   !
   USE iom             ! I/O manager library
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lib_fortran     ! to use key_nosignedzero
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   IMPLICIT NONE
   PRIVATE

   PUBLIC ::   blk_oce_mfs   ! called by sbcblk.F90
   PUBLIC ::   turb_cd_2z    ! routine called in sbcblk_mfs module
   !!----------------------------------------------------------------------

#  include "do_loop_substitute.h90"

CONTAINS

     SUBROUTINE blk_oce_mfs( kt, pwndi, pwndj, pair, pprec, pprec_fqh, pslp, pclc, prh, pst)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE bulk_oce_mfs  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!      (momentum, heat, freshwater, runoff is added later in the code)
      !!
      !! ** Method  : (1) READ Atmospheric data from NetCDF files:
      !!      the 10m wind velocity (i-component) (m/s) at T-point
      !!      the 10m wind velocity (j-component) (m/s) at T-point
      !!      the 2m Dew point Temperature        (Kelvin)
      !!      the Cloud Cover                     (%)
      !!      the 2m air temperature              (Kelvin)
      !!      the Mean Sea Level Preesure         (hPa)
      !!      the Climatological Precipitation    (kg/m2/s)
      !!              (2) CALL blk_oce_mfs
      !!
      !!      Computes:
      !!      Solar Radiation using Reed formula (1975, 1977)
      !!      Net Long wave radiation using Brunt–Berliand as in Rosati and
      !!      Miyakoda (1988)
      !!      (replaced original MFS formulation of Bignami et al., 1995)
      !!      Latent and Sensible heat using Kondo (1975)
      !!      Drag coeff using Hllerman and Rosenstein (1983)
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the mesh
      !!                      referential
      !!                      i.e. the (i,j) referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point over free ocean,
      !!                            or leads in presence of sea-ice
      !!              - qns, qsr    non-slor and solar heat flux
      !!              - emp         evaporation minus precipitation
      !!----------------------------------------------------------------------
        REAL(wp), INTENT(inout   ), DIMENSION(jpi,jpj) ::   pwndi  ! atmospheric wind at T-point              [m/s]
        REAL(wp), INTENT(inout   ), DIMENSION(jpi,jpj) ::   pwndj  ! atmospheric wind at T-point              [m/s]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   prh     ! dew point at 2m at T-points              [K]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pair  ! potential temperature at 2m at T-points        [Kelvin]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pslp   ! sea-level pressure [Pa]
        REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   pprec  ! precipitation
        REAL(wp), INTENT(in   )                 ::   pprec_fqh  ! precipitation frequency
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pclc   ! cloud cover
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pst   ! sea surfacetemperature

        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sh_now   ! specific humidity at T-point
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  catm     ! Cover
        !  !---------------------------------------------------------------------
        !! Local fluxes variables
        !!---------------------------------------------------------------------
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  qbw     ! Net Long wave
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ha      ! Sesnible
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  elat    ! Latent
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  evap    ! evaporation rate
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  prec    ! precipitation rate
        INTEGER, INTENT( in  ) ::   kt
        !!
        INTEGER  :: ierror                          ! return error code
        INTEGER  :: ifpr     ! dummy loop indice
        INTEGER  :: jj,ji    ! dummy loop arguments
        INTEGER  ::   ios    ! Local integer output status for namelist read
        !!--------------------------------------------------------------------
        !!--------------------------------------------------------------------
        !! Variables for specific humidity computation
        !!--------------------------------------------------------------------
        REAL(wp) :: onsea,par1,par2
        DATA onsea,par1,par2 / 0.98, 640380., -5107.4 /
        !!                      par1 [Kg/m3], par2 [K]
        REAL(wp) :: rn_pfac, rn_cfac, rn_efac, rn_vfac
        !!---------------------------------------------------------------------
        !!---------------------------------------------------------------------
        IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
            !                                      ! ====================== !
            ALLOCATE( sh_now(jpi,jpj), catm(jpi,jpj),  &
                    &        qbw(jpi,jpj),    ha(jpi,jpj), elat(jpi,jpj),     &
                    &        evap(jpi,jpj), prec(jpi,jpj),    STAT=ierror )

            IF( ierror /= 0 )   CALL ctl_warn('sbc_blk_mfs: failed to allocate arrays')

            catm(:,:)   = 0.0_wp    ! initializze cloud cover variable
            sh_now(:,:) = 0.0_wp    ! initializze specifif humidity variable
        ENDIF

        DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

                ! Calculate Specific Humidity
                !-------------------------------------------------
                ! Modif Emanuela 2023-02-10 from Maicu v3.6 
                !sh_now(ji,jj) = (1_wp/1.22_wp) * onsea * par1 *EXP(par2/prh(ji,jj))
                sh_now(ji,jj) = (1_wp/1.22_wp) * par1*EXP(par2/prh(ji,jj))

                ! Normalize Clouds
                !-------------------------------------------------
                catm(ji,jj)   = max(0.0_wp,min(1.0_wp,pclc(ji,jj)*0.01_wp))

                ! Convert precipitations from Total prec [m] to
                ! rain fall rate in kg/m2/s
                ! rainfall rate = tot prec [m]* water dens [1000kg/m3] / time
                ! period
                ! [s]

                IF (pprec(ji,jj).lt.0) pprec(ji,jj)=0.0
                pprec(ji,jj) = pprec(ji,jj) *1000._wp/(pprec_fqh*3600)

        END_2D
        ! wind module at 10m
        !--------------------------------------------------
        wndm(:,:) = SQRT(  pwndi(:,:) * pwndi(:,:)   &
                &             + pwndj(:,:) * pwndj(:,:)  )

            ! Force to zero the output of fluxes
            !-------------------------------------------------
            qsr(:,:)  = 0.0 ; qbw(:,:)  = 0.0 ;
            ha(:,:)   = 0.0 ; elat(:,:) = 0.0 ;
            evap(:,:) = 0.0 ; utau(:,:) = 0.0 ;
            vtau(:,:) = 0.0 ; prec(:,:) = 0.0

            CALL lbc_lnk( 'sbcblk', pwndi(:,:), 'T', -1. , pwndj(:,:) , 'T', -1.)

            CALL turb_mfs(pst, pair , sh_now,     &     ! input dynamic
                    pwndi, pwndj,               &     ! input dynamic
                    pslp, catm ,               &     ! input dynamic
                    qsr,qbw,ha,elat,evap,utau,vtau)                          !output

            ! Shortwave radiation
            !--------------------------------------------------
            qsr(:,:) = qsr(:,:) * tmask(:,:,1)
            ! total non solar heat flux over water
            !--------------------------------------------------
            qns(:,:) = -1. * ( qbw(:,:) + ha(:,:) + elat(:,:) )
            qns(:,:) = qns(:,:)*tmask(:,:,1)

            ! mask the wind module at 10m
            !--------------------------------------------------
            wndm(:,:) = wndm(:,:) * tmask(:,:,1)

            !   wind stress module (taum) into T-grid
            !--------------------------------------------------
            taum(:,:) = SQRT( utau(:,:) * utau(:,:) + vtau(:,:) * vtau(:,:) ) *tmask(:,:,1)

            CALL lbc_lnk( 'sbcblk', taum, 'T', 1. )

            ! Interpolate utau, vtau into the grid_V and grid_V
            !-------------------------------------------------
            ! Note the use of 0.5*(2-umask) in order to unmask the stress along
            ! coastlines
            ! Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over
            ! ice shelves
            DO_2D( 1, 0, 1, 0 )
                    utau(ji,jj) = 0.5_wp * ( 2._wp - umask(ji,jj,1) ) * (utau(ji,jj) * tmask(ji,jj,1) &
                            &                                + utau(ji+1,jj) *tmask(ji+1,jj,1) )        &
                            &                 * MAX(tmask(ji,jj,1),tmask(ji+1,jj,1))
                    vtau(ji,jj) = 0.5_wp * ( 2._wp - vmask(ji,jj,1) ) * (vtau(ji,jj) * tmask(ji,jj,1) &
                            &                                + vtau(ji,jj+1) *tmask(ji,jj+1,1) )        &
                            &                 * MAX(tmask(ji,jj,1),tmask(ji,jj+1,1))
            END_2D

            CALL lbc_lnk( 'sbcblk', utau(:,:), 'U', -1.,vtau(:,:), 'V', -1. )

            ! for basin budget and coherence
            !--------------------------------------------------
            !CDIR COLLAPSE
            emp (:,:) = evap(:,:) - pprec(:,:) * tmask(:,:,1)
            prec(:,:) = pprec(:,:) * tmask(:,:,1)
            !CDIR COLLAPSE

            CALL iom_put( "qlw_oce"    ,  -qbw  )               ! output downward longwave heat over the ocean
            CALL iom_put( "qsb_oce"    ,  -ha   )               ! output downward sensible heat over the ocean
            CALL iom_put( "qla_oce"    ,  -elat )               ! output downward latent   heat over the ocean
            CALL iom_put( "qns_oce"    ,   qns  )               ! output downward non solar heat over the ocean
            CALL iom_put( "evap_oce"   ,   evap )               ! output water evaporation flux
            CALL iom_put( "precip  "   ,   prec )               ! output precipitation flux

  END SUBROUTINE blk_oce_mfs

  SUBROUTINE turb_mfs(sst,tnow,shnow,unow,vnow,mslnow,cldnow,qsw,qbw,ha,elat,&
             evap,taux,tauy)
    !!----------------------------------------------------------------------
    !!                    ***  ROUTINE fluxes_mfs  ***
    !!
    !! --- it provides SURFACE HEAT and MOMENTUM FLUXES in MKS :
    !!
    !!  1)   Water flux (WFLUX)                 [ watt/m*m ]
    !!  2)   Short wave flux (QSW)              [ watt/m*m ] Reed 1977
    !!  3)   Long wave flux backward (QBW)      [ watt/m*m ]
    !!  4)   Latent heat of evaporation (ELAT)  [ watt/m*m ]
    !!  5)   Sensible heat flux   (HA)          [ watt/m*m ]
    !!  6)   Wind stress x-component   (TAUX)   [ newton/m*m ]
    !!  7)   Wind stress y-component   (TAUY)   [ newton/m*m ]
    !!
    !!----------------------------------------------------------------------


  REAL(wp), INTENT(in   ), DIMENSION (jpi,jpj) :: sst, unow
  REAL(wp), INTENT(in   ), DIMENSION (jpi,jpj) :: vnow, cldnow, mslnow
  REAL(wp), INTENT(out  ), DIMENSION (jpi,jpj) :: qsw, qbw, ha, elat
  REAL(wp), INTENT(out  ), DIMENSION (jpi,jpj) :: evap, taux,tauy
  REAL(wp), INTENT(in), DIMENSION (jpi,jpj) :: tnow , shnow
  INTEGER :: ji,jj
  REAL(wp)  :: wair, vtnow, ea, deltemp, s, stp , fh , fe
  REAL(wp)  :: esre, cseep

  REAL(wp), DIMENSION(jpi,jpj) ::   rspeed, cdx, ce, shms
  REAL(wp), DIMENSION(jpi,jpj) ::   rhom, sstk, ch, rel_windu, rel_windv

    !!----------------------------------------------------------------------
    !!     coefficients ( in MKS )  :
    !!----------------------------------------------------------------------

  REAL(wp), PARAMETER ::  ps = 1013.    ! --- surface air pressure
  REAL(wp), PARAMETER ::  expsi=0.622   ! --- expsi
  REAL(wp), PARAMETER ::  rd=287.       ! --- dry air gas constant
  REAL(wp), PARAMETER ::  cp=1005.      ! --- specific heat capacity
  REAL(wp), PARAMETER ::  onsea=0.98    ! --- specific humidity factors
  REAL(wp), PARAMETER ::  par1=640380.  ! [Kg/m3]
  REAL(wp), PARAMETER ::  par2=-5107.4  ! [K]

    !---------------------------------------------------------------------
    !--- define Kondo parameters
    !---------------------------------------------------------------------

  REAL(wp), DIMENSION(5) :: a_h = (/0.0,0.927,1.15,1.17,1.652/)
  REAL(wp), DIMENSION(5) :: a_e = (/0.0,0.969,1.18,1.196,1.68/)
  REAL(wp), DIMENSION(5) :: b_h = (/1.185,0.0546,0.01,0.0075,-0.017/)
  REAL(wp), DIMENSION(5) :: b_e = (/1.23,0.0521,0.01,0.008,-0.016/)
  REAL(wp), DIMENSION(5) :: c_h = (/0.0,0.0,0.0,-0.00045,0.0/)
  REAL(wp), DIMENSION(5) :: c_e = (/0.0,0.0,0.0,-0.0004,0.0/)
  REAL(wp), DIMENSION(5) :: p_h = (/-0.157,1.0,1.0,1.0,1.0/)
  REAL(wp), DIMENSION(5) :: p_e = (/-0.16,1.0,1.0,1.0,1.0/)
  INTEGER :: kku                        !index varing with wind speed
    !!----------------------------------------------------------------------
    !! ------------------ (i)      short wave
    !!----------------------------------------------------------------------

  CALL qshort(cldnow,qsw)

  rel_windu(:,:) = 0.0_wp
  rel_windv(:,:) = 0.0_wp


  DO_2D( 0, 1, 0, 1 )
            rel_windu(ji,jj) = unow(ji,jj) - ( ssu_m(ji-1,jj) + ssu_m(ji,jj) ) &
                    &  /max(1.0, umask(ji-1,jj,1) + umask(ji,jj,1))
            rel_windv(ji,jj) = vnow(ji,jj) - ( ssv_m(ji,jj-1) + ssv_m(ji,jj) ) &
                    &  /max(1.0, vmask(ji,jj-1,1) + vmask(ji,jj,1))
  END_2D

  CALL lbc_lnk( 'sbcblk',rel_windu(:,:), 'T', -1. , rel_windv(:,:), 'T', -1.)
  rspeed(:,:)= SQRT(rel_windu(:,:)*rel_windu(:,:)   &
          &                   + rel_windv(:,:)*rel_windv(:,:))

  sstk(:,:) = sst(:,:) + rt0                          !- SST data in Kelvin degrees
  shms(:,:) = (1._wp/1.22_wp)*onsea*par1*EXP(par2/sstk(:,:)) !- Saturation Specific Humidity

    !!----------------------------------------------------------------------
    !! ------------------ (ii)    net long wave
    !!----------------------------------------------------------------------
  DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
          wair = shnow(ji,jj) / (1 - shnow(ji,jj))    ! mixing ratio of the air (Wallace and Hobbs)
          vtnow = (tnow(ji,jj)*(expsi+wair))/(expsi*(1.+wair))   ! virtual temperature of air
          rhom(ji,jj) = 100._wp*(ps/rd)/vtnow                       ! density of the moist air

            !Bignami et al. 1995 used in the Mediterranean Sea
            ea   = (wair  / (wair  + 0.622_wp )) * mslnow(ji,jj)*0.01_wp
            ! Emanuela correction 2023-02-10 Following Maicu
            !qbw(ji,jj) = emic*stefan*( sstk(ji,jj)**4. )                          &
            !           - ( stefan*( tnow(ji,jj)**4. ) * ( 0.653_wp + 0.00535_wp*ea ) )  &
            !           * ( 1._wp + 0.1762_wp*( cldnow(ji,jj)**2._wp ) )
            qbw(ji,jj) = stefan*( sstk(ji,jj)**4. ) &
                       - ( stefan*( tnow(ji,jj)**4. ) * ( 0.653_wp + 0.00535_wp*ea ) )  &
                       * ( 1._wp + 0.1762_wp*( cldnow(ji,jj)**2._wp ) )


            ! Rosati and Miyakoda (1988) used in the Black Sea
            !ea   = (wair  / (wair  + 0.622_wp )) * (mslnow(ji,jj)*0.01_wp) 
            !qbw(ji,jj) = emic*stefan*( sstk(ji,jj)**4. )*(0.39_wp -0.05_wp*sqrt(ea))    &
            !           * (1._wp-cldnow(ji,jj)*0.8_wp) + 4._wp*emic*stefan*(sstk(ji,jj)**3) &
            !           * (sstk(ji,jj)-tnow(ji,jj))
  END_2D

  DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !!----------------------------------------------------------------------
            !! ------------------ (iii)   sensible heat
            !!----------------------------------------------------------------------

            !! --- calculates the term :      ( Ts - Ta )
            !!----------------------------------------------------------------------

            deltemp = sstk(ji,jj) - tnow (ji,jj)

            !!----------------------------------------------------------------------
            !! --- variable turbulent exchange coefficients ( from Kondo 1975 )
            !! --- calculate the Neutral Transfer Coefficent using an empiric formula
            !! --- by Kondo et al. Then it applies the diabatic approximation.
            !!----------------------------------------------------------------------

            s = deltemp/(wndm(ji,jj)**2._wp)   !! --- calculate S
            stp = s*abs(s)/(abs(s)+0.01_wp)    !! --- calculate the Stability Parameter

            !!----------------------------------------------------------------------
            !! --- for stable condition (sst-t_air < 0):
            !!----------------------------------------------------------------------

            IF (s.lt.0. .and. ((stp.gt.-3.3).and.(stp.lt.0.))) THEN
                fh = 0.1_wp+0.03_wp*stp+0.9_wp*exp(4.8_wp*stp)
                fe = fh
            ELSE IF (s.lt.0. .and. stp.le.-3.3) THEN
                fh = 0._wp
                fe = fh
            ELSE                                       ! --- for unstable condition
                fh = 1.0_wp+0.63_wp*sqrt(stp)
                fe = fh
            ENDIF

            !!----------------------------------------------------------------------
            !! --- calculate the coefficient CH,CE,CD
            !!----------------------------------------------------------------------

            IF (wndm(ji,jj) >= 0. .AND. wndm(ji,jj) <= 2.2)       THEN
                kku=1
            ELSE IF (wndm(ji,jj) > 2.2 .AND. wndm(ji,jj) <= 5.0)  THEN
                kku=2
            ELSE IF (wndm(ji,jj) > 5.0 .AND. wndm(ji,jj) <= 8.0)  THEN
                kku=3
            ELSE IF (wndm(ji,jj) > 8.0 .AND. wndm(ji,jj) <= 25.0) THEN
                kku=4
            ELSE IF (wndm(ji,jj) > 25.0 )                         THEN
                kku=5
            ENDIF

            ch(ji,jj) = ( a_h(kku) + b_h(kku) * wndm(ji,jj) ** p_h(kku)      &
                    + c_h(kku) * (wndm(ji,jj)-8._wp ) **2) * fh

            ce(ji,jj) = ( a_e(kku) + b_e(kku) * wndm(ji,jj) ** p_e(kku)      &
                    + c_e(kku) * (wndm(ji,jj)-8._wp ) **2) * fe

            ch(ji,jj) = ch(ji,jj) / 1000.0_wp
            ce(ji,jj) = ce(ji,jj) / 1000.0_wp

            IF (wndm(ji,jj)<0.3) THEN
                ch(ji,jj) = 1.3e-03_wp * fh
                ce(ji,jj) = 1.5e-03_wp * fe
            ELSE IF(wndm(ji,jj)>50.0) THEN
                ch(ji,jj) = 1.25e-03_wp * fh
                ce(ji,jj) = 1.30e-03_wp * fe
            ENDIF

            !!----------------------------------------------------------------------
            !! --- calculates the SENSIBLE HEAT FLUX in MKS ( watt/m*m )
            !!----------------------------------------------------------------------

            ha(ji,jj) = rhom(ji,jj)*cp*ch(ji,jj)*wndm(ji,jj)*deltemp

            !!----------------------------------------------------------------------
            !! ------------------ (iv)  latent heat
            !! --- calculates the LATENT HEAT FLUX  ( watt/m*m )
            !! --- ELAT = L*rho*Ce*|V|*[qs(Ts)-qa(t2d)]
            !!----------------------------------------------------------------------

            esre  = shms(ji,jj) - shnow(ji,jj)   ! --- calculates the term : qs(Ta)-qa(t2d)

            cseep = ce(ji,jj) * wndm(ji,jj) * esre     ! --- calculates the term : Ce*|V|*[qs(Ts)-qa(t2d)]

            evap(ji,jj) = (cseep * rhom(ji,jj))  ! in [kg/m2/sec] !! --- calculates the EVAPORATION RATE [m/yr]
          
            elat(ji,jj) = rhom(ji,jj) * cseep * heatlat(sst(ji,jj))

            !!----------------------------------------------------------------------
            !! --- calculates the Drag Coefficient
            !!----------------------------------------------------------------------

            !!----------------------------------------------------------------------
            !! --- deltemp should be (Ts - Ta) in the formula estimating
            !! --- drag coefficient
            !!----------------------------------------------------------------------

            IF ( .NOT. ln_cdgw ) THEN
               ! Emanuela correction 2023-02-10 Following Maicu
               ! deltemp should be (Ta - Ts) in the formula estimating
               ! cdx(ji,jj) = cd_HR(wndm(ji,jj),deltemp)
               cdx(ji,jj) = cd_HR(wndm(ji,jj),-deltemp)
            ENDIF

  END_2D

  IF (ln_cdgw ) THEN
        CALL turb_cd_2z(2._wp,10._wp,sstk,tnow+2._wp*0.0098_wp,shms,shnow,rspeed,cdx)
  ENDIF
  
  IF( iom_use('Cd_oce') ) CALL iom_put( "Cd_oce",  cdx * tmask(:,:,1) )


    !!----------------------------------------------------------------------
    !! --- calculates the wind stresses in MKS ( newton/m*m )
    !! ---            taux= rho*Cd*|V|u      tauy= rho*Cd*|V|v
    !!----------------------------------------------------------------------

  taux(:,:)= rhom(:,:) * cdx(:,:) * rspeed(:,:) * rel_windu(:,:)
  tauy(:,:)= rhom(:,:) * cdx(:,:) * rspeed(:,:) * rel_windv(:,:)
    !
  END SUBROUTINE turb_mfs

  SUBROUTINE turb_cd_2z( zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd )
        !!----------------------------------------------------------------------
        !!               ***  modified from ROUTINE  turb_core  ***
        !!
        !! ** Purpose :   Computes turbulent transfert coefficients of surface
        !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
        !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
        !!
        !! ** Method : Monin Obukhov Similarity Theory
        !!             + Large & Yeager (2004,2008) closure: CD_n10 = f(U_n10)
        !!
        !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
        !!
        !! ** update: Laurent Brodeau, June 2014:
        !!    - handles both cases zt=zu and zt/=zu
        !!    - optimized: less 2D arrays allocated and less operations
        !!    - better first guess of stability by checking air-sea difference of virtual temperature
        !!       rather than temperature difference only...
        !!      Large & Yeager 2008. Drag-coefficient reduction for Cyclone conditions
        !!    - using code-wide physical constants defined into "phycst.mod" rather than redifining them
        !!      => 'vkarmn' and 'grav'
        !! ** Last update: Emanuela Clementi, Nov 2018:
        !!    - removed case not ln_cdgw
        !!    - removed function cd_neutral_10m (not needed)
        !!    - removed Ce, Ch, ,t_zu, sh_zu from output fields: not needed for mfs
        !!----------------------------------------------------------------------
        REAL(wp), INTENT(in   )                     ::   zt       ! height for T_zt and q_zt                  [m]
        REAL(wp), INTENT(in   )                     ::   zu       ! height for dU                             [m]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature              [Kelvin]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   T_zt     ! potential air temperature            [Kelvin]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_sat    ! sea surface specific humidity         [kg/kg]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                 [kg/kg]
        REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   dU       ! relative wind module at zu              [m/s]
        REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum       (tau)
        !
        INTEGER ::   j_itt
        INTEGER , PARAMETER ::   nb_itt = 5       ! number of itterations
        LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given a different height than U
        !
        REAL(wp), DIMENSION(jpi,jpj) ::   U_zu          ! relative wind at zu [m/s]
        REAL(wp), DIMENSION(jpi,jpj) ::   Ce_n10        ! 10m neutral latent coefficient
        REAL(wp), DIMENSION(jpi,jpj) ::   Ch_n10        ! 10m neutral sensible coefficient
        REAL(wp), DIMENSION(jpi,jpj) ::   Ch            ! transfer coefficient for sensible heat (Q_sens)
        REAL(wp), DIMENSION(jpi,jpj) ::   Ce            ! transfert coefficient for evaporation   (Q_lat)
        REAL(wp), DIMENSION(jpi,jpj) ::   T_zu          ! air temp. shifted at zu                     [K]
        REAL(wp), DIMENSION(jpi,jpj) ::   q_zu          ! spec. hum.shifted at zu               [kg/kg]
        REAL(wp), DIMENSION(jpi,jpj) ::   sqrt_Cd_n10   ! root square of Cd_n10
        REAL(wp), DIMENSION(jpi,jpj) ::   sqrt_Cd       ! root square of Cd
        REAL(wp), DIMENSION(jpi,jpj) ::   zeta_u        ! stability parameter at height zu
        REAL(wp), DIMENSION(jpi,jpj) ::   zeta_t        ! stability parameter at height zt
        REAL(wp), DIMENSION(jpi,jpj) ::   zpsi_h_u, zpsi_m_u
        REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
        REAL(wp), DIMENSION(jpi,jpj) ::   stab          ! 1st stability test integer
        REAL(wp),ALLOCATABLE, DIMENSION(:,:)  :: dT_AS
        !!----------------------------------------------------------------------
        l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp ) ! testing "zu == zt" is risky with double precision

        U_zu = MAX( 0.5_wp , dU )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

        ! emanuela needed for online wave coupling
        dT_AS=T_zt-sst

        !! First guess of stability:
        ztmp0 = T_zt*(1._wp + 0.608_wp*q_zt) - sst*(1._wp + 0.608_wp*q_sat) ! air-sea difference of virtual pot. temp. at zt
        stab  = 0.5_wp + sign(0.5,ztmp0)                           ! stab = 1 if dTv  > 0  => STABLE, 0 if unstable

        !! Neutral coefficients at 10m read from wave output:
        cdn_wave(:,:) = cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) )
        ztmp0   (:,:) = cdn_wave(:,:)

        sqrt_Cd_n10 = SQRT( ztmp0 )
        Ce_n10  = 1.e-3_wp*( 34.6_wp * sqrt_Cd_n10 )
        Ch_n10  = 1.e-3_wp*sqrt_Cd_n10*(18._wp*stab + 32.7_wp*(1._wp - stab))

        !! Initializing transf. coeff. with their first guess neutral equivalents:
        Cd = ztmp0   ;   Ce = Ce_n10   ;   Ch = Ch_n10   ;   sqrt_Cd = sqrt_Cd_n10

        !! Initializing values at z_u with z_t values:
        T_zu = T_zt   ;   q_zu = q_zt

        !!  * Now starting iteration loop
        DO j_itt=1, nb_itt
            !
            ztmp1 = T_zu - sst   ! Updating air/sea differences
            ztmp2 = q_zu - q_sat

            ! Updating turbulent scales :   (L&Y 2004 eq. (7))
            ztmp1  = Ch/sqrt_Cd*ztmp1    ! theta*
            ztmp2  = Ce/sqrt_Cd*ztmp2    ! q*

            ztmp0 = T_zu*(1._wp + 0.608_wp*q_zu) ! virtual potential temperature at zu

            ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
            ztmp0 =  (vkarmn*grav/ztmp0*(ztmp1*(1._wp+0.608_wp*q_zu) + 0.608_wp*T_zu*ztmp2)) / (Cd*U_zu*U_zu)
            !                                                           ( Cd*U_zu*U_zu is U*^2 at zu)
            !! Stability parameters :
            zeta_u   = zu*ztmp0   ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
            zpsi_h_u = psi_h_mfs( zeta_u )
            zpsi_m_u = psi_m_mfs( zeta_u )

            !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
            IF ( .NOT. l_zt_equal_zu ) THEN
                zeta_t = zt*ztmp0 ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
                stab = LOG(zu/zt) - zpsi_h_u + psi_h_mfs(zeta_t)  ! stab just used as temp array!!!
                T_zu = T_zt + ztmp1/vkarmn*stab    ! ztmp1 is still theta*
                q_zu = q_zt + ztmp2/vkarmn*stab    ! ztmp2 is still q*
                q_zu = max(0., q_zu)
            END IF

            sqrt_Cd = vkarmn / ( vkarmn / sqrt_Cd_n10 - zpsi_m_u )
            Cd      = sqrt_Cd * sqrt_Cd
            !
            ztmp0 = (LOG(zu/10._wp) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
            ztmp2 = sqrt_Cd / sqrt_Cd_n10
            ztmp1 = 1._wp + Ch_n10*ztmp0
            Ch  = Ch_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
            !
            ztmp1 = 1._wp + Ce_n10*ztmp0
            Ce  = Ce_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)
            !
        END DO
        !
   END SUBROUTINE turb_cd_2z


   FUNCTION psi_m_mfs(pzeta)   !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
        !-------------------------------------------------------------------------------
        ! universal profile stability function for momentum
        !-------------------------------------------------------------------------------
        REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
        REAL(wp), DIMENSION(jpi,jpj)             :: psi_m_mfs
       !
        INTEGER  ::   ji, jj    ! dummy loop indices
        REAL(wp) :: zta, zx2, zx, zpsi_unst,  zstab   ! local scalars
        !-------------------------------------------------------------------------------
        !
        DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !
            zta = pzeta(ji,jj) !MIN( pzeta(ji,jj) , 15._wp ) !! Very stable conditions (Lpositif and big!) AIMIE should we do this???
            !
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25

            zpsi_unst = 2._wp*LOG(ABS( (1._wp + zx )*0.5_wp ))   & ! AIMIE ABS was not in MFS
               &            + LOG(ABS( (1._wp + zx2)*0.5_wp ))   &
               &          - 2._wp*ATAN(zx) + rpi*0.5_wp
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_m_mfs(ji,jj) = -5._wp* zta  * zstab &   !Stable
               &              + (1._wp - zstab) * zpsi_unst !Unstable
            !
        END_2D
        !
   END FUNCTION psi_m_mfs


   FUNCTION psi_h_mfs( pzeta )    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
        !-------------------------------------------------------------------------------
        ! universal profile stability function for temperature and humidity
        !-------------------------------------------------------------------------------
        REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pzeta
        REAL(wp), DIMENSION(jpi,jpj)             ::   psi_h_mfs
        !
        INTEGER  ::   ji, jj     ! dummy loop indices
        REAL(wp) :: zta, zx2, zx, zpsi_unst, zstab  ! local scalars
        !-------------------------------------------------------------------------------
        !
        DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !
            zta = pzeta(ji,jj) ! MIN( pzeta(ji,jj) , 15._wp ) !! Very stable conditions (Lpositif and large!) AIMIE should we include this???!!!
            !
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT(zx2)
        
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_h_mfs(ji,jj) = -5._wp*zta*zstab   &  ! Stable
                             & + (1._wp - zstab) * zpsi_unst    ! Unstable
            !
        END_2D
        !
   END FUNCTION psi_h_mfs


   REAL(wp) FUNCTION cd_HR(speed,delt)
        !!----------------------------------------------------------------------
        !! --- calculates the Drag Coefficient as a function of the abs. value
        !! --- of the wind velocity ( Hellermann and Rosenstein )
        !!----------------------------------------------------------------------

        REAL(wp), INTENT(in) :: speed,delt
        REAL(wp), PARAMETER  :: a1=0.934e-3 , a2=0.788e-4, a3=0.868e-4
        REAL(wp), PARAMETER  :: a4=-0.616e-6, a5=-.120e-5, a6=-.214e-5

        cd_HR = a1 + a2*speed + a3*delt + a4*speed*speed        &
                + a5*delt*delt  + a6*speed*delt

   END FUNCTION cd_HR

   REAL(wp) function HEATLAT(t)
        !!----------------------------------------------------------------------
        !! --- calculates the Latent Heat of Vaporization ( J/kg ) as function of
        !! --- the temperature ( Celsius degrees )
        !! --- ( from A. Gill  pag. 607 )
        !!
        !! --- Constant Latent Heat of Vaporization ( Rosati,Miyakoda 1988 )
        !!     L = 2.501e+6  (MKS)
        !!----------------------------------------------------------------------

        REAL(wp) , intent(in) :: t

        heatlat = 2.5008e+6_wp -2.3e+3_wp*t

   END FUNCTION HEATLAT


   SUBROUTINE qshort(cldnow,qsw)
        !!----------------------------------------------------------------------
        !!                    ***  ROUTINE qshort  ***
        !!
        !! ** Purpose :   Compute Solar Radiation
        !!
        !! ** Method  :   Compute Solar Radiation according Astronomical
        !!                formulae
        !!
        !! References :   Reed RK (1975) and Reed RK (1977)
        !!
        !! Note: alat,alon - (lat, lon)  in radians
        !!----------------------------------------------------------------------
        REAL(wp) :: hour

        REAL(wp), DIMENSION(jpi,jpj) :: alat,alon
        REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: cldnow
        REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: qsw
        REAL(wp), DIMENSION(12) :: alpham

        REAL(wp), PARAMETER ::   eclips=23.439* (3.141592653589793_wp / 180._wp)
        REAL(wp), PARAMETER ::   solar = 1350.
        REAL(wp), PARAMETER ::   tau = 0.7
        REAL(wp), PARAMETER ::   aozone = 0.09
        REAL(wp), PARAMETER ::   yrdays = 360.
        REAL(wp) :: days, th0,th02,th03, sundec, thsun, coszen, qatten
        REAL(wp) :: qzer, qdir,qdiff,qtot,tjul,sunbet
        REAL(wp) :: albedo
        INTEGER :: jj, ji

        !!----------------------------------------------------------------------
        !! --- albedo monthly values from Payne (1972) as means of the values
        !! --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
        !! --- band of the Mediterranean Sea ) :
        !!----------------------------------------------------------------------

        data alpham /0.095,0.08,0.065,0.065,0.06,0.06,0.06,0.06,        &
                0.065,0.075,0.09,0.10/

        !!----------------------------------------------------------------------
        !!   days is the number of days elapsed until the day=nday_year
        !!----------------------------------------------------------------------
        days = nday_year -1._wp
        th0  = 2._wp*rpi*days/yrdays
        th02 = 2._wp*th0
        th03 = 3._wp*th0

        !! --- sun declination :
        !!----------------------------------------------------------------------
        sundec = 0.006918_wp - 0.399912_wp*cos(th0) + 0.070257_wp*sin(th0) -   &
                0.006758_wp*cos(th02) + 0.000907_wp*sin(th02) -   &
                0.002697_wp*cos(th03) + 0.001480_wp*sin(th03)

        hour = (( nsec_year / rday ) - INT (nsec_year / rday)) * rjjhh

        DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                alon(ji,jj) = glamt(ji,jj) * rad
                alat(ji,jj) = gphit(ji,jj) * rad        
                !! --- sun hour angle :
                !!----------------------------------------------------------------------
                thsun = (hour -12._wp)*15._wp*rad + alon(ji,jj)

                !! --- cosine of the solar zenith angle :
                !!----------------------------------------------------------------------
                coszen =sin(alat(ji,jj))*sin(sundec)                 &
                        +cos(alat(ji,jj))*cos(sundec)*cos(thsun)

                IF(coszen .LE. 5.035D-04) THEN
                    coszen = 0.0_wp
                    qatten = 0.0_wp
                ELSE
                    qatten = tau**(1./coszen)
                END IF

                qzer  = coszen * solar *                                 &
                        (1.0_wp+1.67E-2_wp*cos(rpi*2._wp*(days-3.0_wp)/365.0_wp))**2._wp
                qdir  = qzer * qatten
                qdiff = ((1._wp-aozone)*qzer - qdir) * 0.5_wp
                
                qtot  =  qdir + qdiff
                tjul = (days -81._wp)*rad

                !! --- sin of the solar noon altitude in radians :
                !!----------------------------------------------------------------------
                sunbet=sin(alat(ji,jj))*sin(eclips*sin(tjul)) +   &
                        cos(alat(ji,jj))*cos(eclips*sin(tjul))

                !! --- solar noon altitude in degrees :
                !!----------------------------------------------------------------------

                sunbet = asin(sunbet)/rad

                !!----------------------------------------------------------------------
                !! --- calculates the albedo according to Payne (1972)
                !!----------------------------------------------------------------------

                albedo = alpham(nmonth)

                !!----------------------------------------------------------------------
                !! --- ( radiation as from Reed(1977), Simpson and Paulson(1979) )
                !! --- calculates SHORT WAVE FLUX ( watt/m*m )
                !! --- ( Rosati,Miyakoda 1988 ; eq. 3.8 )
                !!----------------------------------------------------------------------
                IF(cldnow(ji,jj).LT.0.3) THEN
                    qsw(ji,jj) = qtot * (1._wp-albedo)
                ELSE
                    qsw(ji,jj) = qtot*(1._wp-0.62_wp*cldnow(ji,jj)              &
                            + .0019_wp*sunbet)*(1._wp-albedo)
                ENDIF
      END_2D
   END SUBROUTINE qshort


    !!======================================================================

END MODULE sbcblk_algo_mfs
