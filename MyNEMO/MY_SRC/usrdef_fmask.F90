MODULE usrdef_fmask
   !!======================================================================
   !!                     ***  MODULE usrdef_fmask   ***
   !!
   !!                      ===  ORCA configuration  ===
   !!                            (2 and 1 degrees)
   !!
   !! User defined : alteration of land/sea f-point mask in some straits
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec, S. Flavoni)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_fmask  : alteration of f-point land/ocean mask in some straits
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! Massively Parallel Processing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_fmask    ! routine called by dommsk.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_fmask.F90 13435 2020-08-25 14:48:42Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_fmask( cd_cfg, kcfg, pfmsk )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   User defined alteration of the lateral boundary 
      !!              condition on velocity.
      !!
      !! ** Method  :   Local change of the value of fmask at lateral ocean/land 
      !!              boundary in straits in order to increase the viscous 
      !!              boundary layer and thus reduce the transport through the 
      !!              corresponding straits.
      !!                Here only alterations in ORCA R2 and R1 cases
      !!
      !! ** Action :   fmask : land/ocean mask at f-point with increased value 
      !!                       in some user defined straits
      !!----------------------------------------------------------------------
      CHARACTER(len=*)          , INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER                   , INTENT(in   ) ::   kcfg     ! configuration identifier
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pfmsk    ! Ocean/Land f-point mask including lateral boundary cond.
      !
      INTEGER  ::   iif, iil, ii0, ii1, ii   ! local integers
      INTEGER  ::   ijf, ijl, ij0, ij1       !   -       -
      INTEGER  ::   isrow                    ! index for ORCA1 starting row
      !!----------------------------------------------------------------------
      !
      IF( TRIM( cd_cfg ) == "orca" .OR. TRIM( cd_cfg ) == "ORCA" ) THEN      !==  ORCA Configurations  ==!
         !
         SELECT CASE ( kcfg )
         !
         CASE( 2 )                           !  R2 case
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R2: increase lateral friction near the following straits:'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
            !
            IF(lwp) WRITE(numout,*) '      Gibraltar '
            ij0 = 101 + nn_hls       ;   ij1 = 101 + nn_hls           ! Gibraltar strait  : partial slip (pfmsk=0.5)
            ii0 = 139 + nn_hls - 1   ;   ii1 = 140 + nn_hls - 1
            pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            ij0 = 102 + nn_hls       ;   ij1 = 102 + nn_hls
            ii0 = 139 + nn_hls - 1   ;   ii1 = 140 + nn_hls - 1
            pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            !
            IF(lwp) WRITE(numout,*) '      Bab el Mandeb '
            ij0 =  87 + nn_hls       ;   ij1 = 88  + nn_hls          ! Bab el Mandeb : partial slip (pfmsk=1)
            ii0 = 160 + nn_hls - 1   ;   ii1 = 160 + nn_hls - 1
            pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            ij0 =  88 + nn_hls       ;   ij1 =  88 + nn_hls
            ii0 = 159 + nn_hls - 1   ;   ii1 = 159 + nn_hls - 1
            pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            !
            ! We keep this as an example but it is instable in this case 
            !IF(lwp) WRITE(numout,*) '      Danish straits '
            !         ij0 = 115   ;   ij1 = 115 ! Danish straits  : strong slip (pfmsk > 2)
            !         ii0 = 145   ;   ii1 = 146   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
            !         ij0 = 116   ;   ij1 = 116
            !         ii0 = 145   ;   ii1 = 146   ;   pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
            !
         CASE( 1 )                           ! R1 case
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R1: increase lateral friction near the following straits:'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'   
!!gm    ! This dirty section will be suppressed by simplification process:
!!gm    ! all this will come back in input files
!!gm    ! Currently these hard-wired indices relate to configuration with extend grid (jpjglo=332)
            !
            isrow = 332 - (Nj0glo + 1)   ! was 332 - jpjglo -> jpjglo_old_version = Nj0glo + 1
            !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   orca_r1: increase friction near the following straits : '
            IF(lwp) WRITE(numout,*) '      Gibraltar '
            ii0 = 282 + nn_hls - 1       ;   ii1 = 283 + nn_hls - 1        ! Gibraltar Strait 
            ij0 = 241 + nn_hls - isrow   ;   ij1 = 241 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Bhosporus '
            ii0 = 314 + nn_hls - 1       ;   ii1 = 315 + nn_hls - 1        ! Bhosporus Strait 
            ij0 = 248 + nn_hls - isrow   ;   ij1 = 248 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Makassar (Top) '
            ii0 =  48 + nn_hls - 1       ;   ii1 =  48 + nn_hls - 1        ! Makassar Strait (Top) 
            ij0 = 189 + nn_hls - isrow   ;   ij1 = 190 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
            !
            IF(lwp) WRITE(numout,*) '      Lombok '
            ii0 =  44 + nn_hls - 1       ;   ii1 =  44 + nn_hls - 1        ! Lombok Strait 
            ij0 = 164 + nn_hls - isrow   ;   ij1 = 165 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Ombai '
            ii0 =  53 + nn_hls - 1       ;   ii1 =  53 + nn_hls - 1        ! Ombai Strait 
            ij0 = 164 + nn_hls - isrow   ;   ij1 = 165 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      Timor Passage '
            ii0 =  56 + nn_hls - 1       ;   ii1 =  56 + nn_hls - 1        ! Timor Passage 
            ij0 = 164 + nn_hls - isrow   ;   ij1 = 165 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  
            !
            IF(lwp) WRITE(numout,*) '      West Halmahera '
            ii0 =  58 + nn_hls - 1       ;   ii1 =  58 + nn_hls - 1        ! West Halmahera Strait 
            ij0 = 181 + nn_hls - isrow   ;   ij1 = 182 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
            !
            IF(lwp) WRITE(numout,*) '      East Halmahera '
            ii0 =  55 + nn_hls - 1       ;   ii1 =  55 + nn_hls - 1        ! East Halmahera Strait 
            ij0 = 181 + nn_hls - isrow   ;   ij1 = 182 + nn_hls - isrow
            pfmsk( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
            !
         CASE DEFAULT
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'usr_def_fmask : ORCA_R', kcfg,' : NO alteration of fmask in specific straits '
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'   
         END SELECT
#ifdef key_mfs

      ELSEIF( TRIM( cd_cfg ) == "med24" .OR. TRIM( cd_cfg ) == "MED24" ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'usr_def_fmask : increase lateral friction near the following strait (MED24) :'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'

! AIMIE issue with DO LOOP NEMO 4.2 due to mi0 to mi1 in order to avoid the MPI
! error put a function in the whole defined area in order to have 0 in the
! land, 1 in the sea and 4 in the coastline mask :
! function : (fmask-1)*fmask+fmask 

         IF(lwp) WRITE(numout,*) 'Messina strait' !
! EAS6 tidal: full area sea + coastline increase to 4 '
! EAS7 : reduced area + only coastline point increase to 8
! EAS8 : only coastline point increase to 4
         ii0 = 809    ;   ii1 = 811
         ij0 = 192    ;   ij1 = 194
         pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = pfmsk(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) * &
                                                              & (pfmsk(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) - 1.0_wp)+ &
                                                              & pfmsk(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk )

         IF(lwp) WRITE(numout,*) 'Gibraltar strait' !
! EAS6 : no increase in the lateral friction at Gibraltar
! EAS7 : only coastline point increase to 8
! EAS8 : only coastline point increase to 4
         ij0 = 135    ;   ij1 = 141
         ii0 = 297    ;   ii1 = 305
         pfmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = pfmsk(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) *  &
                                                               & (pfmsk(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) - 1.0_wp)+ &
                                                               & pfmsk(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) 
#endif
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'usr_def_fmask : NO alteration of fmask in specific straits '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      CALL lbc_lnk( 'usrdef_fmask', pfmsk, 'F', 1._wp )      ! Lateral boundary conditions on fmask
      !
   END SUBROUTINE usr_def_fmask
   
   !!======================================================================
END MODULE usrdef_fmask
