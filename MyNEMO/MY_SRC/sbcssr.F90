MODULE sbcssr
   !!======================================================================
   !!                       ***  MODULE  sbcssr  ***
   !! Surface module :  heat and fresh water fluxes a restoring term toward observed SST/SSS
   !!======================================================================
   !! History :  3.0  !  2006-06  (G. Madec)  Original code
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.6  !  2018-07  (R. Escudier) Relax can be done only around obs time
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssr       : add to sbc a restoring term toward SST/SSS climatology
   !!   sbc_ssr_init  : initialisation of surface restoring
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition
   USE phycst         ! physical constants
   USE sbcrnf         ! surface boundary condition : runoffs
   !
   USE fldread        ! read input fields
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssr        ! routine called in sbcmod
   PUBLIC   sbc_ssr_init   ! routine called in sbcmod
   PUBLIC   sbc_ssr_alloc  ! routine called in sbcmod

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   erp   !: evaporation damping   [kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qrp   !: heat flux damping        [w/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   coefice   !: under ice relaxation coefficient

   !                                   !!* Namelist namsbc_ssr *
   INTEGER, PUBLIC ::   nn_sstr         ! SST/SSS restoring indicator
   INTEGER, PUBLIC ::   nn_sssr         ! SST/SSS restoring indicator
   REAL(wp)        ::   rn_dqdt         ! restoring factor on SST and SSS
   REAL(wp)        ::   rn_deds         ! restoring factor on SST and SSS
   LOGICAL         ::   ln_sssr_bnd     ! flag to bound erp term 
   REAL(wp)        ::   rn_sssr_bnd     ! ABS(Max./Min.) value of erp term [mm/day]
   INTEGER         ::   nn_sssr_ice     ! Control of restoring under ice
   ! AIMIE modify namelist Relax can be done only around obs time
   INTEGER         ::   nn_hssr         ! Number of hours around observation time for SST/SSS restoration
   INTEGER         ::   nn_ssr_wgts     ! Use weights on the relaxation (stronger around obs time) (0=no weight,1=gaussian,2=cubic)
   LOGICAL         ::   ln_esst         ! Only relax SST if error (read in sn_esst) is less than rn_esst
   REAL(wp)        ::   rn_esst         ! Maximum of SST error for relaxation (if ln_esst)

   REAL(wp) , ALLOCATABLE, DIMENSION(:) ::   buffer   ! Temporary buffer for exchange
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sst   ! structure of input SST (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sss   ! structure of input SSS (file informations, fields read)
   ! AIMIE modify namelist Relax can be done only around obs time
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_esst  ! structure of input error SST (file informations, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcssr.F90 14834 2021-05-11 09:24:44Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssr( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssr  ***
      !!
      !! ** Purpose :   Add to heat and/or freshwater fluxes a damping term
      !!                toward observed SST and/or SSS.
      !!
      !! ** Method  : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS
      !!              - at each nscb time step
      !!                   add a retroaction term on qns    (nn_sstr = 1)
      !!                   add a damping term on sfx        (nn_sssr = 1)
      !!                   add a damping term on emp        (nn_sssr = 2)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zerp     ! local scalar for evaporation damping
      REAL(wp) ::   zqrp     ! local scalar for heat flux damping
      REAL(wp) ::   zsrp     ! local scalar for unit conversion of rn_deds factor
      REAL(wp) ::   zerp_bnd ! local scalar for unit conversion of rn_epr_max factor
      INTEGER  ::   ierror   ! return error code
      !!
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      ! AIMIE modify namelist Relax can be done only around obs time
      TYPE(FLD_N) ::   sn_sst, sn_sss, sn_esst        ! informations about the fields to be read
      INTEGER     ::   isec                     ! Seconds since Jan 1st nit000
      REAL(wp)    ::   dt                       ! Difference in seconds between obs time and current time
      REAL(wp)    ::   weight_sst               ! Weight for SST relaxation
      !!----------------------------------------------------------------------
      !
      IF( nn_sstr + nn_sssr /= 0 ) THEN
         !
         IF( nn_sstr == 1)   CALL fld_read( kt, nn_fsbc, sf_sst )   ! Read SST data and provides it at kt
      ! AIMIE modify namelist Relax can be done only around obs time
         IF( ln_esst )       CALL fld_read( kt, nn_fsbc, sf_esst )  !Read SST error and provides it at kt
         IF( nn_sssr >= 1)   CALL fld_read( kt, nn_fsbc, sf_sss )   ! Read SSS data and provides it at kt

      ! AIMIE modify to have Relax can be done only around obs time
         isec = nsec_year + nsec1jan000  ! Get current time in seconds since Jan 1st nit000 (to compare with naa/nbb)
         !
         !                                         ! ========================= !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN      !    Add restoring term     !
            !                                      ! ========================= !
            !
            qrp(:,:) = 0._wp ! necessary init
            erp(:,:) = 0._wp
            !
            IF( nn_sstr == 1 ) THEN                                   !* Temperature restoring term
! AIMIE modify calculation Relax can be done only around obs time               
               dt = MIN(ABS(isec - sf_sst(1)%nrec(2,sf_sst(1)%naa)),ABS(isec - sf_sst(1)%nrec(2,sf_sst(1)%nbb)))
               IF ( nn_ssr_wgts == 1 ) THEN                              ! Apply gaussian weighted relaxation (with respect to observation time)
                  weight_sst = EXP(-(dt/(2._wp*nn_hssr*3600._wp))*(dt/(nn_hssr*3600._wp)))
               ELSE IF ( nn_ssr_wgts == 2 ) THEN                         ! Apply cubic weighted relaxation (with respect to observation time)
                  weight_sst = EXP(-(dt/(2._wp*nn_hssr*3600._wp))*(dt/(nn_hssr*3600._wp))*(dt/(nn_hssr*3600._wp)))
               ELSE
                  weight_sst = 1._wp
               ENDIF
               !
               IF ( nn_ssr_wgts > 0 .OR. dt < nn_hssr*3600._wp ) THEN        ! If nn_sstr==0, only do the restoring around observation date (nn_hssr hours)
                  IF ( ln_esst ) THEN ! Only apply relaxation for error of SST < rn_esst
                     DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                           IF (sf_esst(1)%fnow(ji,jj,1) < rn_esst) THEN
                              zqrp = rn_dqdt * weight_sst * ( sst_m(ji,jj) - sf_sst(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1)
                           ELSE
                              zqrp = 0._wp
                           ENDIF
                           qns(ji,jj) = qns(ji,jj) + zqrp
                           qrp(ji,jj) = zqrp
                     END_2D
!Call iom not in NEMO 4.2
!                     CALL iom_put( "qrp", qrp )                             ! heat flux damping
                  ELSE
                     DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                           zqrp = rn_dqdt * weight_sst * ( sst_m(ji,jj) - sf_sst(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1)
                           qns(ji,jj) = qns(ji,jj) + zqrp
                           qrp(ji,jj) = zqrp
                     END_2D
!Call iom not in NEMO 4.2
!                    CALL iom_put( "qrp", qrp )                             !heat flux damping
                  ENDIF
               ENDIF
               !
!end romain
            ENDIF
            !
            IF( nn_sssr /= 0 .AND. nn_sssr_ice /= 1 ) THEN
              ! use fraction of ice ( fr_i ) to adjust relaxation under ice if nn_sssr_ice .ne. 1
              ! n.b. coefice is initialised and fixed to 1._wp if nn_sssr_ice = 1
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  SELECT CASE ( nn_sssr_ice )
                    CASE ( 0 )    ;  coefice(ji,jj) = 1._wp - fr_i(ji,jj)              ! no/reduced damping under ice
                    CASE  DEFAULT ;  coefice(ji,jj) = 1._wp + ( nn_sssr_ice - 1 ) * fr_i(ji,jj) ! reinforced damping (x nn_sssr_ice) under ice )
                  END SELECT
               END_2D
            ENDIF
            !
            IF( nn_sssr == 1 ) THEN                                   !* Salinity damping term (salt flux only (sfx))
               zsrp = rn_deds / rday                                  ! from [mm/day] to [kg/m2/s]
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  zerp = zsrp * ( 1. - 2.*rnfmsk(ji,jj) )   &      ! No damping in vicinity of river mouths
                     &        *   coefice(ji,jj)            &      ! Optional control of damping under sea-ice
                     &        * ( sss_m(ji,jj) - sf_sss(1)%fnow(ji,jj,1) ) * tmask(ji,jj,1)
                  sfx(ji,jj) = sfx(ji,jj) + zerp                 ! salt flux
                  erp(ji,jj) = zerp / MAX( sss_m(ji,jj), 1.e-20 ) ! converted into an equivalent volume flux (diagnostic only)
               END_2D
               !
            ELSEIF( nn_sssr == 2 ) THEN                               !* Salinity damping term (volume flux (emp) and associated heat flux (qns)
               zsrp = rn_deds / rday                                  ! from [mm/day] to [kg/m2/s]
               zerp_bnd = rn_sssr_bnd / rday                          !       -              -    
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  zerp = zsrp * ( 1. - 2.*rnfmsk(ji,jj) )   &      ! No damping in vicinity of river mouths
                     &        *   coefice(ji,jj)            &      ! Optional control of damping under sea-ice
                     &        * ( sss_m(ji,jj) - sf_sss(1)%fnow(ji,jj,1) )   &
                     &        / MAX(  sss_m(ji,jj), 1.e-20   ) * tmask(ji,jj,1)
                  IF( ln_sssr_bnd )   zerp = SIGN( 1.0_wp, zerp ) * MIN( zerp_bnd, ABS(zerp) )
                  emp(ji,jj) = emp (ji,jj) + zerp
                  qns(ji,jj) = qns(ji,jj) - zerp * rcp * sst_m(ji,jj)
                  erp(ji,jj) = zerp
                  qrp(ji,jj) = qrp(ji,jj) - zerp * rcp * sst_m(ji,jj)
               END_2D
            ENDIF
            ! outputs
            CALL iom_put( 'hflx_ssr_cea', qrp(:,:) )
            IF( nn_sssr == 1 )   CALL iom_put( 'sflx_ssr_cea',  erp(:,:) * sss_m(:,:) )
            IF( nn_sssr == 2 )   CALL iom_put( 'vflx_ssr_cea', -erp(:,:) )
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_ssr

 
   SUBROUTINE sbc_ssr_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssr_init  ***
      !!
      !! ** Purpose :   initialisation of surface damping term
      !!
      !! ** Method  : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS if required
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zerp     ! local scalar for evaporation damping
      REAL(wp) ::   zqrp     ! local scalar for heat flux damping
      REAL(wp) ::   zsrp     ! local scalar for unit conversion of rn_deds factor
      REAL(wp) ::   zerp_bnd ! local scalar for unit conversion of rn_epr_max factor
      INTEGER  ::   ierror   ! return error code
      !!
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      ! AIMIE restore to obs Romain
      TYPE(FLD_N) ::   sn_sst, sn_sss, sn_esst        ! informations about the fields to be read
      NAMELIST/namsbc_ssr/ cn_dir, nn_sstr, nn_sssr, rn_dqdt, rn_deds, sn_sst, sn_esst, &
              & sn_sss, ln_sssr_bnd, rn_sssr_bnd, nn_sssr_ice, nn_hssr,nn_ssr_wgts,     &
              & ln_esst, rn_esst
      INTEGER     ::  ios
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssr : SST and/or SSS damping term '
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      ! 
      READ  ( numnam_ref, namsbc_ssr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_ssr in reference namelist' )

      READ  ( numnam_cfg, namsbc_ssr, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_ssr in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_ssr )

      IF(lwp) THEN                 !* control print
         WRITE(numout,*) '   Namelist namsbc_ssr :'
         WRITE(numout,*) '      SST restoring term (Yes=1)             nn_sstr        = ', nn_sstr
         WRITE(numout,*) '         dQ/dT (restoring magnitude on SST)         rn_dqdt     = ', rn_dqdt, ' W/m2/K'
         WRITE(numout,*) '         Hours around obs time for SST/SSS restore  nn_hssr     = ', nn_hssr, ' hours'
         WRITE(numout,*) '         flag for relaxation weights                nn_ssr_wgts = ', nn_ssr_wgts
         WRITE(numout,*) '         flag for SST error limit                   ln_esst     = ', ln_esst
         WRITE(numout,*) '         Maximum of SST error                       rn_esst     = ', rn_esst, ' degC'
         WRITE(numout,*) '      SSS damping term (Yes=1, salt   flux)  nn_sssr        = ', nn_sssr
         WRITE(numout,*) '                       (Yes=2, volume flux) '
         WRITE(numout,*) '         dE/dS (restoring magnitude on SST)     rn_deds     = ', rn_deds, ' mm/day'
         WRITE(numout,*) '         flag to bound erp term                 ln_sssr_bnd = ', ln_sssr_bnd
         WRITE(numout,*) '         ABS(Max./Min.) erp threshold           rn_sssr_bnd = ', rn_sssr_bnd, ' mm/day'
         WRITE(numout,*) '      Cntrl of surface restoration under ice nn_sssr_ice    = ', nn_sssr_ice
         WRITE(numout,*) '          ( 0 = no restoration under ice)'
         WRITE(numout,*) '          ( 1 = restoration everywhere  )'
         WRITE(numout,*) '          (>1 = enhanced restoration under ice  )'
      ENDIF
      !
      IF( nn_sstr == 1 ) THEN      !* set sf_sst structure & allocate arrays
         !
         ALLOCATE( sf_sst(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst structure' )
         ALLOCATE( sf_sst(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst now array' )
         !
         ! fill sf_sst with sn_sst and control print
         CALL fld_fill( sf_sst, (/ sn_sst /), cn_dir, 'sbc_ssr', 'SST restoring term toward SST data', 'namsbc_ssr', no_print )
         IF( sf_sst(1)%ln_tint )   ALLOCATE( sf_sst(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst data array' )
         !
      ENDIF
      !
      IF( nn_sssr >= 1 ) THEN      !* set sf_sss structure & allocate arrays
         !
         ALLOCATE( sf_sss(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss structure' )
         ALLOCATE( sf_sss(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss now array' )
         !
         ! fill sf_sss with sn_sss and control print
         CALL fld_fill( sf_sss, (/ sn_sss /), cn_dir, 'sbc_ssr', 'SSS restoring term toward SSS data', 'namsbc_ssr', no_print )
         IF( sf_sss(1)%ln_tint )   ALLOCATE( sf_sss(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss data array' )
         !
      ENDIF
      !
      IF( ln_esst ) THEN      !* set sf_esst structure & allocate arrays
         !
         ALLOCATE( sf_esst(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_esst structure' )
         ALLOCATE( sf_esst(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_esst now array' )
         !
         ! fill sf_esst with sn_esst and control print
         CALL fld_fill( sf_esst, (/ sn_esst /), cn_dir, 'sbc_ssr', 'SST restoring term toward SST data', 'namsbc_ssr', no_print )
         IF( sf_esst(1)%ln_tint )   ALLOCATE( sf_esst(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_esst data array' )
         !
      ENDIF
      !
      coefice(:,:) = 1._wp         !  Initialise coefice to 1._wp ; will not need to be changed if nn_sssr_ice=1
      !                            !* Initialize qrp and erp if no restoring 
      IF( nn_sstr /= 1                   )   qrp(:,:) = 0._wp
      IF( nn_sssr /= 1 .OR. nn_sssr /= 2 )   erp(:,:) = 0._wp
      !
   END SUBROUTINE sbc_ssr_init
         
   INTEGER FUNCTION sbc_ssr_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION sbc_ssr_alloc  ***
      !!----------------------------------------------------------------------
      sbc_ssr_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( erp ) ) THEN
         ALLOCATE( qrp(jpi,jpj), erp(jpi,jpj), coefice(jpi,jpj), STAT= sbc_ssr_alloc )
         !
         IF( lk_mpp                  )   CALL mpp_sum ( 'sbcssr', sbc_ssr_alloc )
         IF( sbc_ssr_alloc /= 0 )   CALL ctl_warn('sbc_ssr_alloc: failed to allocate arrays.')
         !
      ENDIF
   END FUNCTION
      
   !!======================================================================
END MODULE sbcssr
