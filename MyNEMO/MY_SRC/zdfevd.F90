MODULE zdfevd
   !!======================================================================
   !!                       ***  MODULE  zdfevd  ***
   !! Ocean physics: parameterization of convection through an enhancement
   !!                of vertical eddy mixing coefficient
   !!======================================================================
   !! History :  OPA  !  1997-06  (G. Madec, A. Lazar)  Original code
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  !  2009-03  (M. Leclair, G. Madec, R. Benshila) test on both before & after
   !!            4.0  !  2017-04  (G. Madec)  evd applied on avm (at t-point) 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_evd       : increase the momentum and tracer Kz at the location of
   !!                   statically unstable portion of the water column (ln_zdfevd=T)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics variables
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! for iom_put
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_evd    ! called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdfevd.F90 15298 2021-09-28 10:06:42Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_evd( kt, Kmm, Krhs, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_evd  ***
      !!                   
      !! ** Purpose :   Local increased the vertical eddy viscosity and diffu-
      !!      sivity coefficients when a static instability is encountered.
      !!
      !! ** Method  :   tracer (and momentum if nn_evdm=1) vertical mixing 
      !!              coefficients are set to rn_evd (namelist parameter) 
      !!              if the water column is statically unstable.
      !!                The test of static instability is performed using
      !!              Brunt-Vaisala frequency (rn2 < -1.e-12) of to successive
      !!              time-step (Leap-Frog environnement): before and
      !!              now time-step.
      !!
      !! ** Action  :   avt, avm   enhanced where static instability occurs
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kt             ! ocean time-step indexocean time step
      INTEGER                    , INTENT(in   ) ::   Kmm, Krhs      ! time level indices
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avm, p_avt   !  momentum and tracer Kz (w-points)
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      ! NOTE: [tiling] use a SAVE array to store diagnostics, then send after all tiles are finished. This is necessary because p_avt/p_avm are modified on adjacent tiles when using nn_hls > 1. zavt_evd/zavm_evd are then zero on some points when subsequently calculated for these tiles.
      REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   zavt_evd, zavm_evd
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'zdf_evd : Enhanced Vertical Diffusion (evd)'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
            IF(lwp) WRITE(numout,*)
         ENDIF

         ALLOCATE( zavt_evd(jpi,jpj,jpk) )
         IF( nn_evdm == 1 ) ALLOCATE( zavm_evd(jpi,jpj,jpk) )
      ENDIF
      !
      !
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
         zavt_evd(ji,jj,jk) = p_avt(ji,jj,jk)         ! set avt prior to evd application
      END_3D
      !
      SELECT CASE ( nn_evdm )
      !
      CASE ( 1 )           !==  enhance tracer & momentum Kz  ==!   (if rn2<-1.e-12)
         !
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
            zavm_evd(ji,jj,jk) = p_avm(ji,jj,jk)      ! set avm prior to evd application
         END_3D
         !
!! change last digits results
!         WHERE( MAX( rn2(2:jpi,2:jpj,2:jpkm1), rn2b(2:jpi,2:jpj,2:jpkm1) )  <= -1.e-12 ) THEN
!            p_avt(2:jpi,2:jpj,2:jpkm1) = rn_evd * wmask(2:jpi,2:jpj,2:jpkm1)
!            p_avm(2:jpi,2:jpj,2:jpkm1) = rn_evd * wmask(2:jpi,2:jpj,2:jpkm1)
!         END WHERE
         !
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
            IF(  MIN( rn2(ji,jj,jk), rn2b(ji,jj,jk) ) <= -1.e-12 ) THEN
               p_avt(ji,jj,jk) = rn_evd * wmask(ji,jj,jk)
               p_avm(ji,jj,jk) = rn_evd * wmask(ji,jj,jk)
            ENDIF
         END_3D
         !
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
            zavm_evd(ji,jj,jk) = p_avm(ji,jj,jk) - zavm_evd(ji,jj,jk)   ! change in avm due to evd
         END_3D
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN                       ! Do only on the last tile
            CALL iom_put( "avm_evd", zavm_evd )                ! output this change
            DEALLOCATE( zavm_evd )
         ENDIF
         !
      CASE DEFAULT         !==  enhance tracer Kz  ==!   (if rn2<-1.e-12) 
!! change last digits results
!         WHERE( MAX( rn2(2:jpi,2:jpj,2:jpkm1), rn2b(2:jpi,2:jpj,2:jpkm1) )  <= -1.e-12 )
!            p_avt(2:jpi,2:jpj,2:jpkm1) = rn_evd * wmask(2:jpi,2:jpj,2:jpkm1)
!         END WHERE

         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
            IF(  MIN( rn2(ji,jj,jk), rn2b(ji,jj,jk) ) <= -1.e-12 )   &
               p_avt(ji,jj,jk) = rn_evd * wmask(ji,jj,jk)
         END_3D
         !
      END SELECT 

#if defined key_mfs
      !----------------------------------------------------------------------
      ! Enhanced vertical diffusivity in front of Dardanelles Strait
      !----------------------------------------------------------------------
      IF( TRIM( cn_cfg ) == "med24" .OR. TRIM( cn_cfg ) == "MED24" ) THEN
        IF( kt == nit000 ) THEN
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'Enhanced Vertical Diffusivity in front of Dardanelles Strait (MED24):'
          IF(lwp) WRITE(numout,*) 'from level 1 to level 15 (commented for DAR as river)'
          IF(lwp) WRITE(numout,*) '~~~~~~~ '
          IF(lwp) WRITE(numout,*)
        ENDIF

     ! Enhance vertical diffusivity at w-points
       DO jk=1,15
         DO jj=mj0(236+nn_hls), mj1(238+nn_hls)
           DO ji=mi0(1061+nn_hls-1), mi1(1063+nn_hls-1)
              p_avt (ji,jj,jk) = MIN(p_avt (ji,jj,jk)*5000._wp,rn_evd)
           END DO
         END DO
       END DO
      ELSE
       IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'NO Enhanced Vertical Diffusivity in front of Dardanelles Strait:'
         IF(lwp) WRITE(numout,*) 'cfg_name has to be defined'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
       ENDIF
      ENDIF
#endif
      !
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
         zavt_evd(ji,jj,jk) = p_avt(ji,jj,jk) - zavt_evd(ji,jj,jk)   ! change in avt due to evd
      END_3D
      !
      IF( l_trdtra ) CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_tem, jptra_evd, zavt_evd )
      !
      IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN                       ! Do only on the last tile
         CALL iom_put( "avt_evd", zavt_evd )              ! output this change
         DEALLOCATE( zavt_evd )
      ENDIF
      !
   END SUBROUTINE zdf_evd

   !!======================================================================
END MODULE zdfevd
