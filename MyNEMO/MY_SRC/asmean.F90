! ToDo: check why these variables are not exported in v42
!!!#define t_aml ts(:,:,1,jp_tem,:)
!!!#define q_aml ts(:,:,1,jp_sal,:)
MODULE asmean
   !!======================================================================
   !!                       ***  MODULE  sbcssm  ***
   !! Surface module :  provide time-mean ocean surface variables
   !!======================================================================
   !! History - Coded for OceanVar
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   as_mean        : calculate sea surface mean currents, temperature,
   !!                    and salinity over nn_fsbc time-step
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE par_oce         ! ocean space and time domain
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! ocean space and time domain
   USE in_out_manager  ! ocean space and time domain
   USE sbcapr

   IMPLICIT NONE
   PRIVATE

   PUBLIC   as_mean    ! routine called by step.F90

   REAL(wp), DIMENSION(:, :, :, : ), PUBLIC, ALLOCATABLE :: amuv_m, amts_m
   REAL(wp), DIMENSION(:, :       ), PUBLIC, ALLOCATABLE :: amh_m
   REAL(wp), DIMENSION(:, :, :    ), PUBLIC, ALLOCATABLE :: amtq_m

   !! * Substitutions
#  include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE as_mean( kt, KKnn, ktfreq, nv, ln_mv, ln_ib )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_oce  ***
      !!
      !! ** Purpose :   provide ocean variables to assimilation misfit
      !!                computation
      !!
      !! ** Method  :   compute mean fields for the period (kt - ktfreq ) to kt
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, ktfreq, nv, KKnn    ! ocean time step
      LOGICAL, INTENT(in) ::   ln_mv(nv), ln_ib
      !
      REAL(wp) ::   zcoef       ! temporary scalar
      !!---------------------------------------------------------------------
      !                                                   ! ---------------------------------------- !
      IF( ktfreq == 1 ) THEN                              !      Instantaneous surface fields        !
         !                                                ! ---------------------------------------- !
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'as_mean : assimilation mean fields, ktfreq=1 : instantaneous values'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
         !
         IF( ln_mv(1) ) amts_m(:,:,:,1) = ts(:,:,:,jp_tem,KKnn)
         IF( ln_mv(2) ) amts_m(:,:,:,2) = ts(:,:,:,jp_sal,KKnn)
         IF( ln_mv(3) ) amh_m = ssh(:,:,KKnn)
         IF( .NOT. ln_mv(3) ) THEN
                              IF (ln_ib .OR. .NOT. ln_apr_dyn) THEN
                              amh_m = ssh(:,:,KKnn)
                              ELSE
                              amh_m=ssh(:,:,KKnn) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) )
                              ENDIF
         ENDIF

         IF( ln_mv(4) ) amuv_m(:,:,:,1) = uu(:,:,:,KKnn)
         IF( ln_mv(5) ) amuv_m(:,:,:,2) = vv(:,:,:,KKnn)
         IF( ln_mv(6) ) amtq_m(:,:,  1) = ts(:,:,1,jp_tem,KKnn) !t_aml(:,:,KKnn)
         IF( ln_mv(7) ) amtq_m(:,:,  2) = ts(:,:,1,jp_sal,KKnn) !q_aml(:,:,KKnn)
         !
      ELSE
         !                                                ! ---------------------------------------- !
         IF( kt == nit000) THEN                           !       Initialisation: 1st time-step      !
            !                                             ! ---------------------------------------- !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'as_mean : sea surface mean fields'
            !
               IF(lwp) WRITE(numout,*) '~~~~~~~   mean fields initialised to instantaneous values'
               zcoef = REAL( ktfreq - 1, wp )
               IF( ln_mv(1) ) amts_m (:,:,:,1) = zcoef * ts(:,:,:,jp_tem,KKnn)
               IF( ln_mv(2) ) amts_m (:,:,:,2) = zcoef * ts(:,:,:,jp_sal,KKnn)
               IF( ln_mv(3) ) THEN
                              IF (ln_ib) THEN
                              amh_m = zcoef * ssh(:,:,KKnn)
                              ELSE
                              amh_m=zcoef *( ssh(:,:,KKnn) - 0.5 * (ssh_ib(:,:) + ssh_ibb(:,:) ))
                              ENDIF
               ENDIF
               IF( ln_mv(4) ) amuv_m (:,:,:,1) = zcoef * uu(:,:,:,KKnn)
               IF( ln_mv(5) ) amuv_m (:,:,:,2) = zcoef * vv(:,:,:,KKnn)
               IF( ln_mv(6) ) amtq_m (:,:,  1) = zcoef * ts(:,:,1,jp_tem,KKnn) !t_aml(:,:,KKnn)
               IF( ln_mv(7) ) amtq_m (:,:,  2) = zcoef * ts(:,:,1,jp_sal,KKnn) !q_aml(:,:,KKnn)
            !                                             ! ---------------------------------------- !
         ELSEIF( MOD( kt - 2 , ktfreq ) == 0 ) THEN       !   Initialisation: New mean computation   !
            !                                             ! ---------------------------------------- !
            IF( ln_mv(1) ) amts_m(:,:,:,1) = 0._wp
            IF( ln_mv(2) ) amts_m(:,:,:,2) = 0._wp
            IF( ln_mv(3) ) amh_m = 0._wp
            IF( ln_mv(4) ) amuv_m(:,:,:,1) = 0._wp
            IF( ln_mv(5) ) amuv_m(:,:,:,2) = 0._wp
            IF( ln_mv(6) ) amtq_m(:,:,  1) = 0._wp
            IF( ln_mv(7) ) amtq_m(:,:,  2) = 0._wp

         ENDIF
         !                                                ! ---------------------------------------- !
         !                                                !        Cumulate at each time step        !
         !                                                ! ---------------------------------------- !
         IF( ln_mv(1) ) amts_m(:,:,:,1) = amts_m(:,:,:,1) + ts(:,:,:,jp_tem,KKnn)
         IF( ln_mv(2) ) amts_m(:,:,:,2) = amts_m(:,:,:,2) + ts(:,:,:,jp_sal,KKnn)
         IF( ln_mv(3) ) THEN
                         IF( ln_ib ) THEN
                         amh_m = amh_m + ssh(:,:,KKnn)
                         ELSE
                         amh_m = amh_m + ( ssh(:,:,KKnn) - 0.5 * (ssh_ib(:,:) +ssh_ibb(:,:) ))
                         ENDIF
         ENDIF
         IF( ln_mv(4) ) amuv_m(:,:,:,1) = amuv_m(:,:,:,1) + uu(:,:,:,KKnn)
         IF( ln_mv(5) ) amuv_m(:,:,:,2) = amuv_m(:,:,:,2) + vv(:,:,:,KKnn)
         IF( ln_mv(6) ) amtq_m(:,:,  1) = amtq_m(:,:,  1) + ts(:,:,1,jp_tem,KKnn) ! t_aml(:,:,KKnn)
         IF( ln_mv(7) ) amtq_m(:,:,  2) = amtq_m(:,:,  2) + ts(:,:,1,jp_sal,KKnn) !q_aml(:,:,KKnn)


         !                                                ! ---------------------------------------- !
         IF( MOD( kt - 1 , ktfreq ) == 0 ) THEN           !
            !                                             ! ---------------------------------------- !
            zcoef = 1. / REAL( ktfreq, wp )
            IF( ln_mv(1) ) amts_m(:,:,:,1) = amts_m(:,:,:,1) * zcoef
            IF( ln_mv(2) ) amts_m(:,:,:,2) = amts_m(:,:,:,2) * zcoef
            IF( ln_mv(3) ) amh_m = amh_m * zcoef
            IF( ln_mv(4) ) amuv_m(:,:,:,1) = amuv_m(:,:,:,1) * zcoef
            IF( ln_mv(5) ) amuv_m(:,:,:,2) = amuv_m(:,:,:,2) * zcoef
            IF( ln_mv(6) ) amtq_m(:,:  ,1) = amtq_m(:,:,  1) * zcoef
            IF( ln_mv(7) ) amtq_m(:,:  ,2) = amtq_m(:,:,  2) * zcoef
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE as_mean

   !!======================================================================
END MODULE asmean
