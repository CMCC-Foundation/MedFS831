PROGRAM mk_cstdst_file
  !================================================================================
  !               *** PROGRAM mk_cstdst_file ****
  !================================================================================
  !
  !  Purpose: Create a file containing coastal distance
  !
  !  Method:  1) Read in tmask from mesh_mask file to use as a template
  !           2) Calculate 2d/3d distance from coast
  !           3) Write the array to output file
  !
  !  History: Original code: Mario Adani (Feb 2025) - it uses adapted feature of NEMO
  !-------------------------------------------------------------------------------

  ! Declare variables
  USE netcdf
  USE utils
  USE coastdist

  IMPLICIT NONE

  INTEGER  :: ji, jj, jk                         ! dummpy loop variables
  CHARACTER(LEN=200) :: meshfile = '../../medsea-nemo42/tools/DOMAINcfg/mesh_mask.nc'   ! mesh file
  CHARACTER(LEN=200) :: outfile = '../../DATA/coastal_distance.nc'                       ! output file


  CALL grid_info(meshfile)
  WRITE(numout, *) 'jpi = ',jpi
  WRITE(numout, *) 'jpj = ',jpj
  WRITE(numout, *) 'jpk = ',jpk

  ALLOCATE( coast2d(jpi, jpj) )
  ALLOCATE( coast3d(jpi, jpj, jpk) )

  !Create output file
  CALL make_outfile( outfile )

  CALL check_nf90( nf90_get_var( ncin, gphit_id, gphit, (/ 1,1 /), (/ jpi, jpj /) ) )

  !Loop through levels and read in tmask for each level as starting point for
  !coefficient array
  DO jk = 1, jpk
     print*,'Level = ',jk
     CALL cofdis(coast3d(:,:,jk),jk)
     ! Write out resto for this level
     CALL check_nf90( nf90_put_var( ncout, coast3d_id, coast3d(:,:,jk), (/ 1,1,jk /), (/ jpi, jpj,1 /) ) )
  ENDDO
  coast2d(:,:) = coast3d(:,:,1)
  CALL check_nf90( nf90_put_var( ncout, coast2d_id, coast2d, (/ 1,1 /), (/ jpi, jpj /) ) )

  ! Close the output file
  CALL check_nf90( nf90_close(ncout) )


END PROGRAM mk_cstdst_file
