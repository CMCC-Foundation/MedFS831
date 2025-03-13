MODULE coastdist

   USE utils
   USE netcdf

   IMPLICIT NONE
   PUBLIC 

   CONTAINS

   SUBROUTINE cofdis( pdct, jk )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cofdis  ***
      !!
      !! ** Purpose :   Compute the distance between ocean T-points and the
      !!      ocean model coastlines.
      !!
      !! ** Method  :   For each model level, the distance-to-coast is 
      !!      computed as follows : 
      !!       - The coastline is defined as the serie of U-,V-,F-points
      !!      that are at the ocean-land bound.
      !!       - For each ocean T-point, the distance-to-coast is then 
      !!      computed as the smallest distance (on the sphere) between the 
      !!      T-point and all the coastline points.
      !!       - For land T-points, the distance-to-coast is set to zero.
      !!
      !! ** Action  : - pdct, distance to the coastline (argument)
      !!              - NetCDF file 'dist.coast.nc' 
      !!----------------------------------------------------------------------
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT( out ) ::   pdct   ! distance to the coastline
      !!
      INTEGER ::   ji, jj, jl,jk   ! dummy loop indices
      INTEGER ::   iju, ijt, icoast, itime, ierr, icot   ! local integers
      CHARACTER (len=32) ::   clname                     ! local name
      REAL(wp) ::   zdate0                               ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:) ::  zxt, zyt, zzt, zmask
      REAL(wp), POINTER, DIMENSION(:  ) ::  zxc, zyc, zzc, zdis    ! temporary workspace
      LOGICAL , ALLOCATABLE, DIMENSION(:,:) ::  llcotu, llcotv, llcotf   ! 2D logical workspace

      !!----------------------------------------------------------------------
      !
      ALLOCATE( zxt(jpi,jpj) , zyt(jpi,jpj) , zzt(jpi,jpj) , zmask(jpi,jpj)    )
      ALLOCATE(zxc(3*jpi*jpj), zyc(3*jpi*jpj), zzc(3*jpi*jpj), zdis(3*jpi*jpj)     )
      ALLOCATE( llcotu(jpi,jpj), llcotv(jpi,jpj), llcotf(jpi,jpj)  )
      ALLOCATE( gphiu(jpi,jpj), gphiv(jpi,jpj), gphif(jpi,jpj)  )
      ALLOCATE( glamu(jpi,jpj), glamv(jpi,jpj), glamf(jpi,jpj), glamt(jpi,jpj)  )
      ALLOCATE( umask(jpi,jpj), vmask(jpi,jpj), fmask(jpi,jpj)  )
      !

      CALL check_nf90( nf90_get_var( ncin, gphit_id, gphit, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, gphiu_id, gphiu, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, gphiv_id, gphiv, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, gphif_id, gphif, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, glamt_id, glamt, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, glamu_id, glamu, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, glamv_id, glamv, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, glamf_id, glamf, (/ 1,1 /), (/ jpi, jpj /) ) )
      CALL check_nf90( nf90_get_var( ncin, tmask_id, tmask, (/ 1,1,jk /), (/ jpi, jpj,1 /) ) )
      CALL check_nf90( nf90_get_var( ncin, umask_id, umask, (/ 1,1,jk /), (/ jpi, jpj,1 /) ) )
      CALL check_nf90( nf90_get_var( ncin, vmask_id, vmask, (/ 1,1,jk /), (/ jpi, jpj,1 /) ) )
      CALL check_nf90( nf90_get_var( ncin, fmask_id, fmask, (/ 1,1,jk /), (/ jpi, jpj,1 /) ) )

      pdct(:,:) = 0._wp
      zxt(:,:) = COS( rad * gphit(:,:) ) * COS( rad * glamt(:,:) )
      zyt(:,:) = COS( rad * gphit(:,:) ) * SIN( rad * glamt(:,:) )
      zzt(:,:) = SIN( rad * gphit(:,:) )


         ! Define the coastline points (U, V and F)
         DO jj = 2, jpj-1
            DO ji = 2, jpi-1
               zmask(ji,jj) =  ( tmask(ji,jj+1) + tmask(ji+1,jj+1) &
                   &           + tmask(ji,jj  ) + tmask(ji+1,jj  ) )
               llcotu(ji,jj) = ( tmask(ji,jj ) + tmask(ji+1,jj  ) == 1._wp ) 
               llcotv(ji,jj) = ( tmask(ji,jj  ) + tmask(ji  ,jj+1) == 1._wp ) 
               llcotf(ji,jj) = ( zmask(ji,jj) > 0._wp ) .AND. ( zmask(ji,jj) < 4._wp )
            END DO
         END DO

         ! Lateral boundaries conditions
         llcotu(:, 1 ) = umask(:,  2  ) == 1
         llcotu(:,jpj) = umask(:,jpj-1) == 1
         llcotv(:, 1 ) = vmask(:,  2  ) == 1
         llcotv(:,jpj) = vmask(:,jpj-1) == 1
         llcotf(:, 1 ) = fmask(:,  2  ) == 1
         llcotf(:,jpj) = fmask(:,jpj-1) == 1

         llcotu( 1 ,:) = umask(  2  ,:) == 1
         llcotu(jpi,:) = umask(jpi-1,:) == 1
         llcotv( 1 ,:) = vmask(  2  ,:) == 1
         llcotv(jpi,:) = vmask(jpi-1,:) == 1
         llcotf( 1 ,:) = fmask(  2  ,:) == 1
         llcotf(jpi,:) = fmask(jpi-1,:) == 1

         ! Compute cartesian coordinates of coastline points
         ! and the number of coastline points
         icoast = 0
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( llcotf(ji,jj) ) THEN
                  icoast = icoast + 1
                  zxc(icoast) = COS( rad*gphif(ji,jj) ) * COS( rad*glamf(ji,jj) )
                  zyc(icoast) = COS( rad*gphif(ji,jj) ) * SIN( rad*glamf(ji,jj) )
                  zzc(icoast) = SIN( rad*gphif(ji,jj) )
               ENDIF
               IF( llcotu(ji,jj) ) THEN
                  icoast = icoast+1
                  zxc(icoast) = COS( rad*gphiu(ji,jj) ) * COS( rad*glamu(ji,jj) )
                  zyc(icoast) = COS( rad*gphiu(ji,jj) ) * SIN( rad*glamu(ji,jj) )
                  zzc(icoast) = SIN( rad*gphiu(ji,jj) )
               ENDIF
               IF( llcotv(ji,jj) ) THEN
                  icoast = icoast+1
                  zxc(icoast) = COS( rad*gphiv(ji,jj) ) * COS( rad*glamv(ji,jj) )
                  zyc(icoast) = COS( rad*gphiv(ji,jj) ) * SIN( rad*glamv(ji,jj) )
                  zzc(icoast) = SIN( rad*gphiv(ji,jj) )
               ENDIF
            END DO
         END DO

         ! Distance for the T-points
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj) == 0._wp ) THEN
                  pdct(ji,jj) = 0._wp
               ELSE
                  DO jl = 1, icoast
                     zdis(jl) = ( zxt(ji,jj) - zxc(jl) )**2   &
                        &     + ( zyt(ji,jj) - zyc(jl) )**2   &
                        &     + ( zzt(ji,jj) - zzc(jl) )**2
                  END DO
                  pdct(ji,jj) = ra * SQRT( MINVAL( zdis(1:icoast) ) )
               ENDIF
            END DO
         END DO

      DEALLOCATE( zxt , zyt , zzt , zmask    )
      DEALLOCATE(zxc, zyc, zzc, zdis    )
      DEALLOCATE( llcotu, llcotv, llcotf  )
      DEALLOCATE( gphiu, gphiv, gphif  )
      DEALLOCATE( glamu, glamv, glamf, glamt  )
      DEALLOCATE( umask, vmask, fmask  )

   END SUBROUTINE cofdis

END MODULE coastdist
