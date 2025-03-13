#define nprd_sla 12
#define nprd_sst 12
#define nprd_sss 12
#define NPRD_SLA 12
#define NPRD_SST 12
#define NPRD_SSS 12


MODULE obsmisfits
   !!======================================================================
   !!                       ***  MODULE   obsmisfits  ***
   !! Compute observation misfits for use with OceanVar
   !!=====================================================================
   !! History :  1.0  !  2010-09  (A. Storto)  Original code
   !!----------------------------------------------------------------------

   !!- This module provides an interface to compute observation misifts online
   !!- for successive use within GLOceanVar.
   !!-
   !!- * Note that the observational vector is parallelized following
   !!- the domain decomposition in NEMO.

   !!----------------------------------------------------------------------
   !!   obsmis_init   : initialize misfits computation
   !!   obs_misfits   : performs computation
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE ioipsl, ONLY :   ymds2ju
   USE netcdf
!   USE prepbc
   USE sbc_oce
   USE sbcapr
   USE asmean
   USE mppini
! DISABLED FOR LAM CONFIG
!!!!!!!!!!!   USE diurnal_bulk    ! diurnal warm layer
!!!!!!!!!!!   USE cool_skin       ! Cool skin


   IMPLICIT NONE
   PUBLIC
   SAVE

   INTEGER :: nobsm_freq = 10, nomts=0, nn_sststart=-1, nn_sstend=-1
   INTEGER :: nn_save_background = -1
   LOGICAL :: ln_obsmisfits = .TRUE.
   LOGICAL, PARAMETER :: ln_prepbc = .FALSE.
   INTEGER , PARAMETER, PUBLIC ::   nprds = 12
!   INTEGER , PARAMETER, PUBLIC ::   nprd_sla = nprds, nprd_sst = nprds, nprd_sss = nprds
   INTEGER , PARAMETER:: numuni = 1031
   INTEGER ,ALLOCATABLE :: itstep(:),ntobs(:, :),ntstart(:, :),ntend(:, :)
   REAL(wp),ALLOCATABLE :: mdt(:,:)
   LOGICAL :: ln_noaa_sst = .TRUE., ln_addsstre = .TRUE.
   LOGICAL :: ln_sshcorr = .FALSE.
   LOGICAL :: ll_slaobs
   INTEGER, PARAMETER :: jpav = 7
   INTEGER, PARAMETER :: np_sip = 6
   LOGICAL :: ln_meanvars(jpav)
   LOGICAL :: ln_asmean

   REAL ( wp ), ALLOCATABLE :: dsst(:,:) , e1e2_i(:,:)
   INTEGER, PARAMETER :: nreyn_x=1440, nreyn_y=720, npiq=4
   REAL(wp), DIMENSION(nreyn_x, nreyn_y) :: osst, oerr, oice
   INTEGER , DIMENSION(nreyn_x, nreyn_y, npiq) :: orib, orjb
   REAL(wp), DIMENSION(nreyn_x, nreyn_y, npiq) :: orpb

   REAL(wp), ALLOCATABLE, DIMENSION(:) :: south, west, north, east
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: tsn_n, tsn_e
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: tq2_n, tq2_e
   REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: sla_n, sla_e
   INTEGER  :: imsn, imwe
   LOGICAL  :: ln_nf, ln_sf
   LOGICAL  :: ln_ef, ln_wf
   LOGICAL  :: l_mfc = .TRUE.
   LOGICAL  :: ln_omdebug = .FALSE.
   LOGICAL  :: ln_add_ib_to_ssh = .FALSE.

   LOGICAL  :: ln_smooth_noaa_sst = .TRUE.
   INTEGER  :: nn_sstts = 0
   REAL(wp), ALLOCATABLE :: sst_wgt(:)

   REAL(wp) :: atot
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: seaice


!... For restdiff
   LOGICAL :: ln_restdiff = .FALSE.
   INTEGER :: nrd_ts        = 5
   INTEGER :: nn_rd_rest    = 560
   INTEGER :: nn_rd_restd(99)

!... For mpp
   INTEGER :: nono,noso,noea,nowe
   INTEGER ,  DIMENSION(  :), ALLOCATABLE :: nleit, nlejt !first, last indoor index for each j-domain
   INTEGER :: nproc
   INTEGER, DIMENSION(  :), ALLOCATABLE ::  nlcit,nlcjt,nimppt,njmppt

#include "obs_params.h90"
#include "obs_str.h90"

CONTAINS

   SUBROUTINE obs_misfits( kt, Knn )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obs_misfits  ***
      !!
      !! ** Purpose :   performs misfits computation at each timestep
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT( in  ) ::   kt, Knn   ! ocean time step
      INTEGER :: nt , jobs , i1 , j1 , k1 , i2 , j2 , k2, kp, jk, knl
      INTEGER :: i1m , j1m , i2m , j2m , i_aux, j_aux, jpred, nn_tmp, npobs
      INTEGER, PARAMETER :: jp_tem = 1, jp_sal = 2

      REAL(wp) :: zsshcorr = 0._wp, rdiv, sal_tmp

      IF( kt == nit000 ) CALL obsmis_init
      IF( .not. ln_obsmisfits ) RETURN
      IF( kt == nn_save_background ) THEN
          WRITE(numuni,*) ' SAVE BACKGROUND AT TIMESTEP ', kt
          CALL save_background(Knn)
      ENDIF

      IF( ln_prepbc ) THEN
!          CALL prep_bc(kt, ln_noaa_sst, nn_sststart, nn_sstend)
      ELSE
          IF( sla%no .GT. 0 ) sla%bcp(:,:)=0._wp
          IF( sst%no .GT. 0 ) sst%bcp(:,:)=0._wp
      ENDIF

      IF( ln_asmean ) CALL as_mean( kt, Knn, nobsm_freq, jpav, ln_meanvars,ln_add_ib_to_ssh )

      IF( .NOT. ln_meanvars(1) ) amts_m(:,:,:,1) = ts(:,:,:,jp_tem,Knn)
      IF( .NOT. ln_meanvars(2) ) amts_m(:,:,:,2) = ts(:,:,:,jp_sal,Knn)
      IF( .NOT. ln_meanvars(3) ) THEN
                                 IF (ln_add_ib_to_ssh .OR. .NOT. ln_apr_dyn  ) THEN
                                        amh_m = ssh(:,:,Knn)
                                 ELSE
                                         amh_m =ssh(:,:,Knn)- 0.5 * (ssh_ib(:,:) + ssh_ibb(:,:) )
                                 ENDIF
      ENDIF
      IF( .NOT. ln_meanvars(4) ) amuv_m(:,:,:,1) = uu(:,:,:,Knn)
      IF( .NOT. ln_meanvars(5) ) amuv_m(:,:,:,2) = vv(:,:,:,Knn)
      IF( .NOT. ln_meanvars(6) ) amtq_m(:,:,1) = ts(:,:,1,jp_tem,Knn)
      IF( .NOT. ln_meanvars(7) ) amtq_m(:,:,2) = ts(:,:,1,jp_sal,Knn)

      IF( MOD(kt-1,nobsm_freq) == 0) THEN
      nomts=nomts+1

      WRITE(numuni,*) ' OBSMISFITS TIMESTEP = ',nomts,'  MOD TS = ',kt
      CALL flush(numuni)

      nt = (kt-nit000)/nobsm_freq + 1

      WRITE(numuni,*) ' OBSMISFITS TIMESTEP nt = ',nt
      WRITE(numuni,*) ' NOBS INS = ',ntobs(nt,NnINS)
      WRITE(numuni,*) ' AREA, ISTART, JSTART == ', &
                      narea, nimppt(narea), njmppt(narea)
      CALL flush(numuni)

      npobs = maxval( ntobs(nt,:) )
      IF( lk_mpp ) then
          CALL mpp_max('obs_misfits',npobs)
          IF(npobs .gt. 0 .and. jpnij .gt. 1) CALL prep_mppmisfits(Knn)
      ENDIF

      IF( ntobs(nt,NnINS) .GT. 0 ) THEN

        DO jobs=ntstart(nt,NnINS),ntend(nt,NnINS)

          IF( ins%kty(jobs) .LT. 5 ) THEN

             k1 = ins%kb(jobs)
             k2 = k1 + 1

             IF( ins%par(jobs) .EQ. kktemp .OR. ins%par(jobs) .EQ. kksst) kp = jp_tem
             IF( ins%par(jobs) .EQ. kksal ) kp = jp_sal

             if (ins%sd1(jobs) .eq. ins%sd2(jobs) ) then
               ins%bac(jobs)  = osum(ins%pq(jobs,1:npq),&
               & amts_m(ins%ib(jobs,1:npq)-nimppt(narea)+1,ins%jb(jobs,1:npq)-njmppt(narea)+1,k1,kp)) + &
                              osum(ins%pq(jobs,npq+1:2*npq),&
               & amts_m(ins%ib(jobs,1:npq)-nimppt(narea)+1,ins%jb(jobs,1:npq)-njmppt(narea)+1,k2,kp))
             else
               ins%bac(jobs) = mpp_misfits_ts(ins%sd1(jobs),ins%sd2(jobs),ins%ib(jobs,:), ins%jb(jobs,:), &
               & kp, k1, ins%pq(jobs,:), ins%moi(jobs,:), ins%moj(jobs,:),ins%nind(jobs) )
             endif
             if(ln_omdebug) then
                WRITE(numuni,*) '*' , jobs, ins%par(jobs), ins%sd1(jobs), ins%sd2(jobs), &
                & ins%bac(jobs), nimppt(narea), &
                & ins%ib(jobs,1:npq) - nimppt(narea) + 1, ins%jb(jobs,1:npq) - njmppt(narea) + 1, &
                & k1,kp, amts_m(ins%ib(jobs,1)-nimppt(narea)+1,ins%jb(jobs,1)-njmppt(narea)+1,k1,kp)
             endif
             ! Compute pot. temp. from model salinity (nearest neighbour) if needed
             if ( kp == jp_tem .AND. ins%val(jobs) .gt. 150._wp ) then
               sal_tmp = osum(ins%pq(jobs,1:npq),&
               & amts_m(ins%ib(jobs,1:npq)-nimppt(narea)+1,ins%jb(jobs,1:npq)-njmppt(narea)+1,k1,jp_sal)) + &
                         osum(ins%pq(jobs,npq+1:2*npq),&
               & amts_m(ins%ib(jobs,1:npq)-nimppt(narea)+1,ins%jb(jobs,1:npq)-njmppt(narea)+1,k2,jp_sal))
               ins%val(jobs) = potemp( sal_tmp, ins%val(jobs) - 200._wp, dep_to_p(ins%dpt(jobs),ins%lat(jobs)), 0.0_wp )
             endif

          ELSEIF ( ins%kty(jobs) .EQ. 5 ) THEN

             IF( ins%par(jobs) .EQ. kkt2m ) kp = 1
             IF( ins%par(jobs) .EQ. kkq2m ) kp = 2

             if (ins%sd1(jobs) .eq. ins%sd2(jobs) ) then
               ins%bac(jobs)  = osum(ins%pq(jobs,1:npq),&
               & amtq_m(ins%ib(jobs,1:npq)-nimppt(narea)+1,ins%jb(jobs,1:npq)-njmppt(narea)+1,kp))
             else
               ins%bac(jobs) = mpp_misfits_tq(ins%sd1(jobs),ins%sd2(jobs),ins%ib(jobs,:), ins%jb(jobs,:), &
               & kp, ins%pq(jobs,:), ins%moi(jobs,:), ins%moj(jobs,:),ins%nind(jobs) )
             endif

          ELSE
            WRITE(numuni,*) 'Unsupported kty : ',ins%kty(jobs)
            CALL ABORT()
          ENDIF
          ins%res(jobs)  = ins%val(jobs) - ins%bac(jobs)
        ENDDO

      ENDIF

      WRITE(numuni,*) ' NOBS SLA = ',ntobs(nt,NnSLA)
      CALL flush(numuni)

      IF( ln_sshcorr .and. ll_slaobs ) THEN

          ! Firstly, Unbias the SSH !
          zsshcorr = SUM( e1e2_i(:,:) * amh_m(:,:) ) / atot
          IF (lk_mpp) CALL mpp_sum('obs_misfits', zsshcorr )

          WRITE(numuni,*) ' ssh correction = ', zsshcorr
          CALL flush(numuni)

      ENDIF

      IF( ntobs(nt,NnSLA) .GT. 0 ) THEN

        DO jobs=ntstart(nt,NnSLA),ntend(nt,NnSLA)
            if (sla%sd1(jobs) .eq. sla%sd2(jobs) ) then
            sla%bac(jobs)  = &
              & osum(sla%pq(jobs,:), amh_m(sla%ib(jobs,:) - nimppt(narea) + 1, sla%jb(jobs,:) - njmppt(narea) + 1)) - &
              & osum(sla%pq(jobs,:), mdt(sla%ib(jobs,:) - nimppt(narea) + 1, sla%jb(jobs,:) - njmppt(narea) + 1))
          else
            sla%bac(jobs) = mpp_misfits_sla(sla%sd1(jobs),sla%sd2(jobs),sla%ib(jobs,:), sla%jb(jobs,:), &
              & sla%pq(jobs,:), sla%moi(jobs,:), sla%moj(jobs,:), sla%nind(jobs), Knn )

          endif

          if (sla%sd1(jobs) .eq. sla%sd2(jobs) ) then
           DO jk=1,jpk
            sla%tb(jobs,jk)  = osum(sla%pq(jobs,:),&
            & amts_m(sla%ib(jobs,:)-nimppt(narea)+1,sla%jb(jobs,:)-njmppt(narea)+1,jk,1))
            sla%sb(jobs,jk)  = osum(sla%pq(jobs,:),&
            & amts_m(sla%ib(jobs,:)-nimppt(narea)+1,sla%jb(jobs,:)-njmppt(narea)+1,jk,2))
           ENDDO
          else
           DO jk=1,jpk
            sla%tb(jobs,jk)  = mpp_misfits_ts(sla%sd1(jobs),sla%sd2(jobs),sla%ib(jobs,:), sla%jb(jobs,:), &
            & 1, jk, sla%pq(jobs,:), sla%moi(jobs,:), sla%moj(jobs,:),sla%nind(jobs) )
            sla%sb(jobs,jk)  = mpp_misfits_ts(sla%sd1(jobs),sla%sd2(jobs),sla%ib(jobs,:), sla%jb(jobs,:), &
            & 2, jk, sla%pq(jobs,:), sla%moi(jobs,:), sla%moj(jobs,:),sla%nind(jobs) )
           ENDDO
          endif

          sla%bcp(jobs,:) = -9.E+12

          sla%res(jobs)  = sla%val(jobs) - ( sla%bac(jobs) - zsshcorr )

        ENDDO

      ENDIF

      WRITE(numuni,*) ' NOBS SST = ',ntobs(nt,NnSST)
      CALL flush(numuni)

      IF( .NOT. ln_noaa_sst .AND. ntobs(nt,NnSST) .GT. 0 ) THEN

        DO jobs=ntstart(nt,NnSST),ntend(nt,NnSST)

          if (sst%sd1(jobs) .eq. sst%sd2(jobs) ) then
            sst%bac(jobs)  = osum(sst%pq(jobs,1:npq),&
            & amts_m(sst%ib(jobs,1:npq)-nimppt(narea)+1,sst%jb(jobs,1:npq)-njmppt(narea)+1,1,1))
          else
            sst%bac(jobs) = mpp_misfits_ts(sst%sd1(jobs),sst%sd2(jobs),sst%ib(jobs,:), sst%jb(jobs,:), &
            & 1, 1, sst%pq(jobs,:), sst%moi(jobs,:), sst%moj(jobs,:),sst%nind(jobs) )
          endif

           sst%res(jobs) = sst%val(jobs) - sst%bac(jobs)

        ENDDO
      ENDIF

      WRITE(numuni,*) ' NOBS SSS = ',ntobs(nt,NnSSS)
      CALL flush(numuni)

      IF( ntobs(nt,NnSSS) .GT. 0 ) THEN

        DO jobs=ntstart(nt,NnSSS),ntend(nt,NnSSS)
          if (sss%sd1(jobs) .eq. sss%sd2(jobs) ) then
            sss%bac(jobs)  = osum(sss%pq(jobs,1:npq),&
            & amts_m(sss%ib(jobs,1:npq)-nimppt(narea)+1,sss%jb(jobs,1:npq)-njmppt(narea)+1,1,2))
          else
            sss%bac(jobs) = mpp_misfits_ts(sss%sd1(jobs),sss%sd2(jobs),sss%ib(jobs,:), sss%jb(jobs,:), &
            & 2, 1, sss%pq(jobs,:), sss%moi(jobs,:), sss%moj(jobs,:),sss%nind(jobs) )
            sss%bcp(jobs,4) = mpp_misfits_ts(sss%sd1(jobs),sss%sd2(jobs),sss%ib(jobs,:), sss%jb(jobs,:), &
            & 1, 1, sss%pq(jobs,:), sss%moi(jobs,:), sss%moj(jobs,:),sss%nind(jobs) )
            sss%bcp(jobs,5:7) = -1.e+20
          endif
          sss%res(jobs) = sss%val(jobs) - sss%bac(jobs)
        ENDDO

      ENDIF

      WRITE(numuni,*) ' Done nt ',nt
      CALL flush(numuni)

      ENDIF

      ! For NOAA/SST average using all timesteps
      IF( ln_noaa_sst) THEN
        nn_tmp = kt - nit000 + 1
        IF( nn_tmp .ge. nn_sststart .and. nn_tmp .le. nn_sstend ) THEN
          IF( .NOT. ln_smooth_noaa_sst ) THEN
            dsst = dsst + ts(:,:,1,jp_tem,Knn)/REAL( nn_sstend-nn_sststart+1 , wp )
          ELSE
            nn_sstts = nn_sstts+1
            dsst = dsst + ts(:,:,1,jp_tem,Knn)*sst_wgt(nn_sstts)
          ENDIF
        ENDIF
      ENDIF


      IF( kt == nitend ) CALL obsmis_term

   END SUBROUTINE obs_misfits

   SUBROUTINE obsmis_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obsmis_init  ***
      !!
      !! ** Purpose :   initialize misfits computation
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      !MPP
      INTEGER ::    iimax, ijmax,jj,ji,jl,ii,ij,jn
      INTEGER :: nreci,nrecj,iresti,irestj
      INTEGER, DIMENSION(:,:), ALLOCATABLE ::  ilcit,ilcjt
      INTEGER, DIMENSION(:,:), ALLOCATABLE ::  iimppt, ijmppt, ijpi,ijpj, iproc
      INTEGER, DIMENSION(  :), ALLOCATABLE ::  iipos,  ijpos
      LOGICAL, DIMENSION(:,:), ALLOCATABLE ::  llisoce

      INTEGER :: tot_sla, it, nt, ios
      CHARACTER(LEN=60) :: cfdiam
      LOGICAL :: ln_meant, ln_meanh, ln_meanuv, ln_mean2m
      NAMELIST/namobs_misfits/ ln_obsmisfits, nobsm_freq, ln_noaa_sst, nn_save_background, &
      & nn_sststart, nn_sstend, ln_addsstre, ln_sshcorr, &
      & ln_meant, ln_meanh, ln_meanuv, ln_mean2m, ln_smooth_noaa_sst, &
      & ln_restdiff, nrd_ts, nn_rd_rest, nn_rd_restd, &
      & ln_omdebug,ln_add_ib_to_ssh
        !nn_thick_cons,ln_seaice,ln_only_sic, ln_icevar, ,ll_iv_iau,nn_iv_iau,nn_iv_verb, nn_istop_iv, nn_ivtime, &
      ln_meant = .FALSE.
      ln_meanh = .FALSE.
      ln_meanuv= .FALSE.
!      ln_seaice= .FALSE.

      nono=mpinei(jpno)
      noso=mpinei(jpso)
      nowe=mpinei(jpwe)
      noea=mpinei(jpea)

      ALLOCATE( nleit(jpnij), nlejt(jpnij) )
      ALLOCATE(nlcit(jpnij),nlcjt(jpnij),nimppt(jpnij),njmppt(jpnij))
      ! find decomposition
      ALLOCATE(iimppt(jpni,jpnj),ijmppt(jpni,jpnj),ijpi(jpni,jpnj), &
       & ijpj(jpni,jpnj),llisoce(jpni,jpnj),  iproc(jpni,jpnj),&
       & iipos(jpni*jpnj), ijpos(jpni*jpnj) )

       ALLOCATE(ilcit(jpni,jpnj),ilcjt(jpni,jpnj))

       CALL mpp_basesplit( jpiglo, jpjglo, nn_hls, jpni, jpnj, iimax,ijmax, iimppt, ijmppt, ijpi, ijpj )
       CALL mpp_is_ocean( llisoce )
       CALL mpp_getnum( llisoce, iproc, iipos, ijpos )

       DO jl =1,jpnij
             ii = iipos(jl)
             ij = ijpos(jl)
             nimppt(jl)= iimppt(ii,ij)
             njmppt(jl)= ijmppt(ii,ij)
             nleit(jl) = ijmppt(ii,ij) - 1 +      1      + nn_hls
             nlejt(jl) = ijmppt(ii,ij) - 1 + ijpj(ii,ij) - nn_hls
             nlcit (jl) = ijpi(ii,ij)  - 2 * nn_hls
             nlcjt (jl) = ijpj(ii,ij)  - 2 * nn_hls
       END DO

      nproc = narea - 1

      WRITE(cfdiam,'(A,I4.4,A)') 'misfits_',narea-1,'.out'
      OPEN( numuni, FILE=cfdiam)

      WRITE(numuni,* ) 'nproc = ',nproc
      WRITE(numuni,* ) 'iproc = ',iproc
      WRITE(numuni,* ) 'nimppt = ',nimppt
      WRITE(numuni,* ) 'njmppt = ',njmppt
      WRITE(numuni,* ) 'nlcit = ',nlcit
      WRITE(numuni,* ) 'nlcjt = ',nlcjt


   DEALLOCATE( iimppt, ijmppt, ijpi, ijpj, llisoce, iproc, iipos, ijpos)
   DEALLOCATE(ilcjt,ilcit)

!
      READ  ( numnam_ref, namobs_misfits, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namobs_misfits in reference namelist' )
      !
      READ  ( numnam_cfg, namobs_misfits, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namobs_misfits in configuration namelist' )
      IF(lwm) WRITE ( numond, namobs_misfits )


      IF( .not. ln_obsmisfits ) RETURN

      ALLOCATE( dsst(jpi,jpj), e1e2_i(jpi,jpj) )

      ALLOCATE( amts_m(0:jpi-1,0:jpj-1,jpk,2) )
      ALLOCATE( amuv_m(0:jpi-1,0:jpj-1,jpk,2) )
      ALLOCATE( amh_m (0:jpi-1,0:jpj      ) )
      ALLOCATE( amtq_m(0:jpi-1,0:jpj-1,    2) )

      imsn = 3*jpi*jpk*2+3*jpi*3
      imwe = 3*jpj*jpk*2+3*jpj*3
      ALLOCATE ( south(imsn), north(imsn) )
      ALLOCATE ( west (imwe), east (imwe) )
      ALLOCATE ( tsn_n(jpi,3,jpk,2), tsn_e(3,jpj,jpk,2) )
      ALLOCATE ( sla_n(jpi,3), sla_e(3,jpj) )
      ALLOCATE ( tq2_n(jpi,3,2), tq2_e(3,jpj,2) )

      ln_sf = ( njmpp .eq. 1 .OR. noso .lt. 0  )
      ln_nf = ( nono  .lt. 0 .OR. ( njmpp + Nj_0 ) .ge. jpjglo .OR.  nono.ge.jpnij )
      ln_wf = ( nowe  .lt. 0 .OR. nowe .ge. jpnij )
      ln_ef = ( noea  .lt. 0 .OR. noea .ge. jpnij )

      ln_meanvars(1) = ln_meant
      ln_meanvars(2) = ln_meant
      ln_meanvars(3) = ln_meanh
      ln_meanvars(4) = ln_meanuv
      ln_meanvars(5) = ln_meanuv
      ln_meanvars(6) = ln_mean2m
      ln_meanvars(7) = ln_mean2m

      ln_asmean = ( COUNT( (/ln_meant,ln_meanh,ln_meanuv,ln_mean2m/) ) .GT. 0 )

      WRITE(numuni,* ) '*** OBSMISFITS ***'
      WRITE(numuni,* ) 'ln_obsmisfits = ',ln_obsmisfits
      WRITE(numuni,* ) 'ln_noaa_sst   = ',ln_noaa_sst
      WRITE(numuni,* ) 'nn_sststart   = ',nn_sststart
      WRITE(numuni,* ) 'nn_sstend     = ',nn_sstend
      WRITE(numuni,* ) 'ln_addsstre   = ',ln_addsstre
      WRITE(numuni,* ) 'nobsm_freq    = ',nobsm_freq
      WRITE(numuni,* ) 'processor     = ',narea-1
      WRITE(numuni,* ) 'domain        = ',narea
      WRITE(numuni,* ) 'ln_meanvars   = ',ln_meanvars
!      WRITE(numuni,* ) 'ln_seaice     = ',ln_seaice
      CALL flush(numuni)

     IF( ( ln_noaa_sst ).AND. (nn_sststart.EQ.-1 .OR. nn_sstend.EQ.-1) ) THEN
         WRITE(*,*) 'sststart AND sstend to be properly setup!'
          STOP 'obsmis_init E R R O R'
      ENDIF

      IF( ln_smooth_noaa_sst ) THEN
        WRITE(numuni,* ) ' Smoothing SST parameters'
        nt = nn_sstend - nn_sststart + 1
        WRITE(numuni,* ) ' SST Timestep :', nt
        ALLOCATE( sst_wgt(nt) )
        nt = nt / 3
        WRITE(numuni,* ) ' Third SST Timestep :', nt
        sst_wgt(1:nt) = 1._wp
        sst_wgt((nt+1):(2*nt)) = 2._wp
        sst_wgt((2*nt+1):(3*nt)) = 1._wp
        sst_wgt = sst_wgt / SUM(sst_wgt)
        WRITE(numuni,* ) ' WEIGTHS'
        DO it=1,nt*3
            WRITE(numuni,* ) it,sst_wgt(it)
        ENDDO
      ENDIF

      ! Read observations
      WRITE(numuni,*) ' About to read obs...'
      CALL flush(numuni)
      CALL readobs
      WRITE(numuni,*) ' Done'
      WRITE(numuni,*)
      CALL flush(numuni)

      ! Timestep arranging
      WRITE(numuni,*) ' About to prepare obs...'
      CALL flush(numuni)
      CALL prepobs
      WRITE(numuni,*) ' Done'
      WRITE(numuni,*)
      CALL flush(numuni)

      ! Read MDT
      tot_sla = sla%no
      WRITE(numuni,*) ' Local SLA obs = ',tot_sla
      CALL flush(numuni)
      IF ( lk_mpp ) CALL mpp_sum('obs_misfits', tot_sla)
      WRITE(numuni,*) ' Total SLA obs = ',tot_sla
      CALL flush(numuni)
      WRITE(numuni,*) ' About to read MDT...'
      CALL flush(numuni)
      CALL readmdt(tot_sla)
      WRITE(numuni,*) ' Done'
      WRITE(numuni,*)
      CALL flush(numuni)

      ! Read SST
      IF( ln_noaa_sst ) Then
          Call Read_ReynSST
          dsst = 0._wp
      ENDIF

      ! Prepare for SSH correction
      if (ln_sshcorr ) THEN
        e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
        atot = SUM( e1e2_i(:,:) )
        IF( lk_mpp )   CALL  mpp_sum('obs_misfits', atot )  ! sum over the global domain
        WRITE(numuni,*) ' Global Ocean Area : ',atot
      Endif

      WRITE(numuni,*) ' Leaving obsmis_init...'
      CALL flush(numuni)

   END SUBROUTINE obsmis_init

   SUBROUTINE Arrange_ReynSST
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE Arrange_ReynSST  ***
      !!
      !! ** Purpose :   Arrange Gridded SST data for 3DVAR
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER  :: JX, JY, nr, jpred
      INTEGER  :: myfi, myfj, myli, mylj, nn_ice
      REAL(wp) :: zice, zmice
      INTEGER  :: i1,i2,j1,j2
      INTEGER  :: i1l,i2l,j1l,j2l
      CHARACTER(LEN=60) :: CFILE

      sst%no= nreyn_x*nreyn_y
      nr = 0
      nn_ice = 0

      allocate(sst%lon(sst%no),&
      &sst%lat(sst%no),&
      &sst%val(sst%no),&
      &sst%tdist(sst%no),&
      &sst%nind(sst%no),&
      &sst%flg(sst%no),&
      &sst%flc(sst%no),&
      &sst%ksat(sst%no),&
      &sst%eve(sst%no),&
      &sst%inc(sst%no),&
      &sst%bac(sst%no),&
      &sst%bia(sst%no),&
      &sst%err(sst%no),&
      &sst%res(sst%no),&
      &sst%b_a(sst%no),&
      &sst%ib (sst%no,npq),&
      &sst%jb (sst%no,npq),&
      &sst%pq(sst%no,npq),&
      &sst%bcp(sst%no,1:nprd_sst),&
      &sst%track(sst%no),&
      &sst%tim(sst%no))

      myfi = nimppt(narea)
      myfj = njmppt(narea)
      myli = nimppt(narea) + nleit(narea) - 1
      mylj = njmppt(narea) + nlejt(narea) - 1

      WRITE(numuni,*) ' SST REPORT FROM PROC',narea
      WRITE(numuni,*)
      WRITE(numuni,*) ' My x limits are ', myfi, myli
      WRITE(numuni,*) ' My y limits are ', myfj, mylj
      WRITE(numuni,*)
      WRITE(numuni,*) ' Ice model MIN ', MINVAL(fr_i)
      WRITE(numuni,*) ' Ice model MAX ', MAXVAL(fr_i)
      WRITE(numuni,*)

      DO JY=1, nreyn_y
        DO JX=1, nreyn_x

	   i1 = minval( orib(JX,JY,:) )
	   i2 = maxval( orib(JX,JY,:) )
	   j1 = minval( orjb(JX,JY,:) )
	   j2 = maxval( orjb(JX,JY,:) )

           IF( i1 .GE. myfi .AND. i2 .LE. myli .AND. &
             & j1 .GE. myfj .AND. j2 .LE. mylj .AND. &
             & OSST(JX, JY) .GT. -2. ) THEN

             zice =  osum(orpb(JX,JY,:),oice(orib(JX,JY,:),orjb(JX,JY,:) ) )
             zmice =  osum(orpb(JX,JY,:),&
	   & oice(orib(JX,JY,:)-nimppt(narea)+1,orjb(JX,JY,:)-njmppt(narea)+1))

             ! Reject obs within observed or modelled ice-covered areas
             IF( zice .lt. 0.01_wp .AND. zmice .lt. 0.01_wp) THEN

                 nr = nr + 1
                 sst%nind(nr) = nr
                 sst%lon( nr ) = 0.125_wp + (JX-1)*0.25_wp

		 if ( sst%lon( nr ) .GT. 180._wp ) &
                 & sst%lon( nr )=sst%lon( nr )-360._wp

                 sst%lat( nr ) = -89.875_wp + (JY-1)*0.25_wp
                 sst%val( nr ) = OSST(JX, JY)
                 sst%tdist( nr ) = 0._wp
                 sst%flg( nr ) = 1
                 sst%flc( nr ) = 1
                 sst%ksat( nr ) = 6
                 sst%eve( nr ) = 0

                 sst%ib ( nr ,: ) = orib(JX,JY,:)
                 sst%jb ( nr ,: ) = orjb(JX,JY,:)

                 sst%pq( nr,: ) =  orpb(JX,JY,:)

                 sst%inc( nr ) = 0._wp
                 sst%bac( nr )  = osum(sst%pq(nr,:),&
                 & dsst(sst%ib(nr,:)-nimppt(narea)+1,sst%jb(nr,:)-njmppt(narea)+1))
                 sst%bia( nr ) = 0._wp
                 sst%err( nr ) = OERR(JX, JY)
                 sst%res( nr ) = sst%val( nr ) - sst%bac( nr )
                 sst%b_a( nr ) = 0._wp

!                 IF( ln_prepbc ) THEN
!                    DO jpred=1,nprd_sst
!                    sst%bcp(nr,jpred) = osum( &
!                    & sst%pq(nr,:),preds(sst%ib(nr,:),sst%jb(nr,:),jpred) )
!                    sst%bcp(nr,jpred) = sst%bcp(nr,jpred) / &
!                    & REAL(n_ssttscnt, KIND = wp)
!                    ENDDO
!                 ELSE
                    sst%bcp(nr,:) = 0._wp
!                 ENDIF

                 sst%track( nr ) = 0
                 sst%tim( nr ) = 0._wp

             ELSE

                 nn_ice = nn_ice + 1

             ENDIF

           ENDIF

        ENDDO
      ENDDO

      sst%no = nr
      sst%nc = nr

      WRITE(CFILE,'(A,I4.4,A)') 'SST_TAB_',narea-1,'.DAT'
      OPEN(913,FILE=CFILE,STATUS='NEW')
      WRITE(913,*) nr
      CLOSE(913)

      WRITE(numuni,*) ' NUMBER OF SST data valid     = ', nr
      WRITE(numuni,*) ' NUMBER OF SST data under ice = ', nn_ice

   END SUBROUTINE Arrange_ReynSST

   SUBROUTINE obsmis_term
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obsmis_term  ***
      !!
      !! ** Purpose :   terminate misfits computation
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      CHARACTER(LEN=99) :: cfd
      INTEGER, PARAMETER :: ufd = 7696

      IF( ln_noaa_sst) THEN
         Call Arrange_ReynSST
      ENDIF

      WRITE(numuni,* ) ' Entering obmis_term'
      CALL flush(numuni)

      ! Write observations
      CALL writeobs

      IF( ALLOCATED(itstep ) ) DEALLOCATE(itstep )
      IF( ALLOCATED(ntobs  ) ) DEALLOCATE(ntobs  )
      IF( ALLOCATED(ntstart) ) DEALLOCATE(ntstart)
      IF( ALLOCATED(ntend  ) ) DEALLOCATE(ntend  )

      WRITE(numuni,* ) '*** OBSMISFITS END ***'
      CALL flush(numuni)
      CLOSE( numuni )

      DEALLOCATE ( south, north )
      DEALLOCATE ( west , east  )
      DEALLOCATE ( tsn_n, tsn_e )
      DEALLOCATE ( sla_n, sla_e )

      IF( ln_smooth_noaa_sst ) DEALLOCATE( sst_wgt )

      DEALLOCATE( dsst, e1e2_i )
      DEALLOCATE( amts_m, amuv_m, amh_m, amtq_m )

   END SUBROUTINE obsmis_term

   SUBROUTINE Read_ReynSST
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE   Read ReynSST    ***
      !!
      !! ** Purpose :   read NOAA Reynolds Sea Surface Temperature (REYNSST)
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE netcdf
      USE lib_mpp
      IMPLICIT NONE

      INTEGER :: stat, ncid, idvar, ierr, jj, myid
      CHARACTER( LEN=30 ) :: csstfile
      REAL( wp ) :: rerr( nreyn_x, nreyn_y)

      myid = narea-1

      WRITE(csstfile,'(A,I4.4,A)') 'reynsst_',myid,'.nc'
      stat = NF90_OPEN(csstfile,NF90_NOWRITE,ncid)
      IF(stat.NE.0) THEN
         WRITE(*,*) ' FILE ', TRIM(csstfile), ' not found'
         CALL ABORT
      ENDIF
      stat = NF90_INQ_VARID (ncid, 'sst', IDVAR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_GET_VAR (ncid, IDVAR, OSST)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      OSST = OSST / REAL(100, wp)
      stat = NF90_INQ_VARID (ncid, 'err', IDVAR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_GET_VAR (ncid, IDVAR, OERR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      OERR = OERR / REAL(100, wp)
      stat = NF90_INQ_VARID (ncid, 'ice', IDVAR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_GET_VAR (ncid, IDVAR, OICE)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_CLOSE(ncid)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)


      WRITE(csstfile,'(A,I4.4,A)') 'reyngrid_',myid,'.nc'
      stat = NF90_OPEN(csstfile,NF90_NOWRITE,ncid)
      IF(stat.NE.0) THEN
         WRITE(*,*) ' FILE ', TRIM(csstfile), ' not found'
         CALL ABORT
      ENDIF
      stat = NF90_INQ_VARID (ncid, 'ib', IDVAR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_GET_VAR (ncid, IDVAR, ORIB)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)

      stat = NF90_INQ_VARID (ncid, 'jb', IDVAR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_GET_VAR (ncid, IDVAR, ORJB)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)

      stat = NF90_INQ_VARID (ncid, 'pb', IDVAR)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
      stat = NF90_GET_VAR (ncid, IDVAR, ORPB)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)

      stat = NF90_CLOSE(ncid)
      IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)

      if( ln_addsstre ) THEN
        WRITE(csstfile,'(A,I4.4,A)') 'sst_reperr_',myid,'.nc'
        stat = NF90_OPEN(csstfile,NF90_NOWRITE,ncid)
        IF(stat.NE.0) THEN
           WRITE(*,*) ' FILE ', TRIM(csstfile), ' not found'
           CALL ABORT
        ENDIF
        stat = NF90_INQ_VARID (ncid, 'sst_reperr', IDVAR)
        IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
        stat = NF90_GET_VAR (ncid, IDVAR, RERR)
        IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)
        stat = NF90_CLOSE(ncid)
        IF (stat /= NF90_NOERR) CALL HANDLE_ERR(stat)

        ! Sum Up SST ERROR VARIANCES
        OERR = SQRT( OERR*OERR + RERR*RERR )

      ENDIF

   END SUBROUTINE Read_ReynSST

   SUBROUTINE readmdt(no_sla)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE   readmdt    ***
      !!
      !! ** Purpose :   read mean dynamic topography from netcdf
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER :: no_sla
      INTEGER :: inum, id1
      CHARACTER(LEN=60) :: cmdtfile

      WRITE(cmdtfile,'(A)') 'MDT'
      ALLOCATE(MDT(jpi,jpj))
      MDT(:,:) = 0._wp
      IF( no_sla > 0 ) THEN
        CALL iom_open(cmdtfile, inum )
        CALL iom_get(inum,jpdom_auto,'mdt',mdt)
        CALL iom_close( inum )
      ENDIF

   END SUBROUTINE readmdt

   SUBROUTINE HANDLE_ERR(ERRCODE)
      USE lib_mpp
      USE NETCDF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ERRCODE
      INTEGER :: IERR
      IF(ERRCODE /= NF90_NOERR) THEN
        PRINT *, 'ERROR: ', TRIM(NF90_STRERROR(ERRCODE))
        CALL ABORT
      ENDIF
   END SUBROUTINE HANDLE_ERR

   SUBROUTINE prepobs
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE    prepobs   ***
      !!
      !! ** Purpose :   prepare time-slots for misfits computation
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ,ALLOCATABLE :: itstep(:)
      REAL(wp),ALLOCATABLE :: tstdiff(:)
      INTEGER ,ALLOCATABLE :: ktstep(:)
      INTEGER :: kny, knm, knd , nmaxobs
      INTEGER :: nobs, ntst, jt, kk, kt, ndtmp
      REAL(wp):: zjulday,zjul1950

      kny  =   ndastp / 10000
      knm  = ( ndastp - (nyear * 10000) ) / 100
      knd  =   ndastp - (nyear * 10000) - ( nmonth * 100 )
      CALL ymds2ju( kny, knm, knd, 0._wp, zjulday )
      CALL ymds2ju( 1950, 1, 1, 0._wp, zjul1950 )

      zjulday = zjulday - zjul1950

      ! Prepare timestep data
      ntst=INT((nitend-nit000+1)/nobsm_freq)+1
      ALLOCATE(itstep(ntst))
      ALLOCATE(ntobs(ntst,NoFams))
      ALLOCATE(ntstart(ntst,NoFams))
      ALLOCATE(ntend(ntst,NoFams))

      ntobs = 0

      WRITE(numuni,*) ' Starting date is ',ndastp
      WRITE(numuni,*) ' Julian day is ',zjulday
      WRITE(numuni,*) ' Timestep in s is ',rdt
      WRITE(numuni,*) ' timestep start is ',nit000
      WRITE(numuni,*) ' timestep end is ',nitend
      WRITE(numuni,*) ' obsmisfits tsteps amount is ',ntst
      CALL flush(numuni)

      nmaxobs = MAXVAL( (/ins%no,sla%no,sst%no,sss%no/) )

      WRITE(numuni,*) ' Max obs per type = ', nmaxobs
      CALL flush(numuni)

      IF( nmaxobs .GT. 0 ) THEN

        WRITE(numuni,*) ' Preparing '
        CALL flush(numuni)

        ALLOCATE(tstdiff(nmaxobs),ktstep(nmaxobs))

        kk = 0
        DO jt=nit000,nitend,nobsm_freq
           kk = kk + 1
           itstep(kk) = jt
        END DO

        ! Insitu part
        nobs=ins%no
        WRITE(numuni,*) ' Treating INS: ',nobs
        CALL flush(numuni)
        IF (nobs.GT.0) THEN
          WRITE(numuni,*) ' MIN/MAX INS obs time:', ins%tim(1) ,ins%tim(nobs)
          CALL flush(numuni)
          tstdiff(1:nobs) = ins%tim(1:nobs) - zjulday
          ktstep(1:nobs)  = NINT(tstdiff(1:nobs) * 86400. / rdt)
          WHERE( ktstep(1:nobs) .EQ. 0 ) ktstep(1:nobs) = 1
          IF( ANY(ktstep(1:nobs).LT. 1 ) .OR. &
              ANY(ktstep(1:nobs).GT. nitend-nit000+1 ) ) THEN
             WRITE(numuni,*) ' INS Obs are out of time range '
             CALL mppstop
          ENDIF
          ktstep(1:nobs)  =  NINT(REAL(ktstep(1:nobs))/REAL(nobsm_freq))
          where(ktstep.eq.0) ktstep = 1
          where(ktstep.gt.ntst) ktstep = ntst

          DO jt=1,nobs
             kt = ktstep(jt)
             ntobs( kt, NnINS ) = ntobs( kt, NnINS ) + 1
          ENDDO
          DO jt=2,ntst
             ntstart(jt, NnINS) = SUM( ntobs(1:jt-1, NnINS) ) + 1
             ntend  (jt, NnINS) = SUM( ntobs(1:jt, NnINS) )
          END DO
          ntstart(1, NnINS) = 1
          ntend  (1, NnINS) = ntobs(1, NnINS)
        ENDIF

        ! SLA part
        nobs=sla%no
        WRITE(numuni,*) ' Treating SLA: ',nobs
        CALL flush(numuni)
        IF (nobs.GT.0) THEN
          WRITE(numuni,*) ' MIN/MAX SLA obs time:', sla%tim(1) ,sla%tim(nobs)
          CALL flush(numuni)
          tstdiff(1:nobs) = sla%tim(1:nobs) - zjulday
          ktstep(1:nobs)  = NINT(tstdiff(1:nobs) * 86400. / rdt)
          WHERE( ktstep(1:nobs) .EQ. 0 ) ktstep(1:nobs) = 1
          IF( ANY(ktstep(1:nobs).LT. 1 ) .OR. &
              ANY(ktstep(1:nobs).GT. nitend-nit000+1 ) ) THEN
             WRITE(numuni,*) ' SLA Obs are out of time range '
             CALL mppstop
          ENDIF
          ktstep(1:nobs)  =  NINT(REAL(ktstep(1:nobs))/REAL(nobsm_freq))
          where(ktstep.eq.0) ktstep = 1
          where(ktstep.gt.ntst) ktstep = ntst

          DO jt=1,nobs
             kt = ktstep(jt)
             ntobs( kt, NnSLA ) = ntobs( kt, NnSLA ) + 1
          ENDDO
          DO jt=2,ntst
             ntstart(jt, NnSLA) = SUM( ntobs(1:jt-1, NnSLA) ) + 1
             ntend  (jt, NnSLA) = SUM( ntobs(1:jt, NnSLA) )
          END DO
          ntstart(1, NnSLA) = 1
          ntend  (1, NnSLA) = ntobs(1, NnSLA)
        ENDIF

        ! SST part
        nobs=sst%no
        WRITE(numuni,*) ' Treating SST: ',nobs
        CALL flush(numuni)
        IF (nobs.GT.0 .AND. .NOT. ln_noaa_sst) THEN
          WRITE(numuni,*) ' MIN/MAX SST obs time:', sst%tim(1) ,sst%tim(nobs)
          CALL flush(numuni)
          tstdiff(1:nobs) = sst%tim(1:nobs) - zjulday
          ktstep(1:nobs)  = NINT(tstdiff(1:nobs) * 86400. / rdt)
          WHERE( ktstep(1:nobs) .EQ. 0 ) ktstep(1:nobs) = 1
          IF( ANY(ktstep(1:nobs).LT. 1 ) .OR. &
              ANY(ktstep(1:nobs).GT. nitend-nit000+1 ) ) THEN
             WRITE(numuni,*) ' SST Obs are out of time range '
             CALL mppstop
          ENDIF
          ktstep(1:nobs)  =  NINT(REAL(ktstep(1:nobs))/REAL(nobsm_freq))
          where(ktstep.eq.0) ktstep = 1
          where(ktstep.gt.ntst) ktstep = ntst

          DO jt=1,nobs
             kt = ktstep(jt)
             ntobs( kt, NnSST ) = ntobs( kt, NnSST ) + 1
          ENDDO
          DO jt=2,ntst
             ntstart(jt, NnSST) = SUM( ntobs(1:jt-1, NnSST) ) + 1
             ntend  (jt, NnSST) = SUM( ntobs(1:jt, NnSST) )
          END DO
          ntstart(1, NnSST) = 1
          ntend  (1, NnSST) = ntobs(1, NnSST)
        ENDIF

        ! SSS part
        nobs=sss%no
        WRITE(numuni,*) ' Treating SSS: ',nobs
        CALL flush(numuni)
        IF (nobs.GT.0) THEN
          WRITE(numuni,*) ' MIN/MAX SSS obs time:', sss%tim(1) ,sss%tim(nobs)
          CALL flush(numuni)
          tstdiff(1:nobs) = sss%tim(1:nobs) - zjulday
          ktstep(1:nobs)  = NINT(tstdiff(1:nobs) * 86400. / rdt)
          WHERE( ktstep(1:nobs) .EQ. 0 ) ktstep(1:nobs) = 1
          IF( ANY(ktstep(1:nobs).LT. 1 ) .OR. &
              ANY(ktstep(1:nobs).GT. nitend-nit000+1 ) ) THEN
             WRITE(numuni,*) ' SSS Obs are out of time range '
             CALL mppstop
          ENDIF
          ktstep(1:nobs)  =  NINT(REAL(ktstep(1:nobs))/REAL(nobsm_freq))
          where(ktstep.eq.0) ktstep = 1
          where(ktstep.gt.ntst) ktstep = ntst

          DO jt=1,nobs
             kt = ktstep(jt)
             ntobs( kt, NnSSS ) = ntobs( kt, NnSSS ) + 1
          ENDDO
          DO jt=2,ntst
             ntstart(jt, NnSSS) = SUM( ntobs(1:jt-1, NnSSS) ) + 1
             ntend  (jt, NnSSS) = SUM( ntobs(1:jt, NnSSS) )
          END DO
          ntstart(1, NnSSS) = 1
          ntend  (1, NnSSS) = ntobs(1, NnSSS)
        ENDIF

        WRITE(numuni,*)
        WRITE(numuni,*) ' -----------------------------'
        WRITE(numuni,*) ' REPORT FROM PREPOBS '
        WRITE(numuni,*) ' -----------------------------'
        WRITE(numuni,*)
        WRITE(numuni,*) ' OTS  MTS     INS     SLA      SST     SSS'

        kk = 0
        DO jt=nit000,nitend,nobsm_freq
           kk = kk + 1
           CALL calc_date( nit000, jt, ndastp, ndtmp )
           WRITE(numuni,*) kk,ndtmp,itstep(kk), ntobs(kk,nnins),&
           & ntobs(kk,nnsla),ntobs(kk,nnsst),ntobs(kk,nnsss)
        END DO
        WRITE(numuni,*) ' -----------------------------'
        WRITE(numuni,*) 'TOTAL     ', sum(ntobs(:,nnins)),&
         sum(ntobs(:,nnsla)),sum(ntobs(:,nnsst)),sum(ntobs(:,nnsss))

        DEALLOCATE(tstdiff,ktstep)

        WRITE(numuni,*)
        IF( SUM( ntobs(:,nnins) ) .GT. 0 ) THEN
          WRITE(numuni,*) ' Minimum, Maximum INS ib',&
          & minval(ins%ib(1:ins%no,:)), maxval(ins%ib(1:ins%no,:))
          WRITE(numuni,*) ' Minimum, Maximum INS jb',&
          & minval(ins%jb(1:ins%no,:)), maxval(ins%jb(1:ins%no,:))
          WRITE(numuni,*) ' Minimum, Maximum INS kb',&
          & minval(ins%kb(1:ins%no)), maxval(ins%kb(1:ins%no))
        ENDIF
        IF( SUM( ntobs(:,nnsla) ) .GT. 0 ) THEN
          WRITE(numuni,*) ' Minimum, Maximum SLA ib',&
          & minval(sla%ib(1:sla%no,:)), maxval(sla%ib(1:sla%no,:))
          WRITE(numuni,*) ' Minimum, Maximum SLA jb',&
          & minval(sla%jb(1:sla%no,:)), maxval(sla%jb(1:sla%no,:))
        ENDIF
        IF( SUM( ntobs(:,nnsst) ) .GT. 0 ) THEN
          WRITE(numuni,*) ' Minimum, Maximum SST ib',&
          & minval(sst%ib(1:sst%no,:)), maxval(sst%ib(1:sst%no,:))
          WRITE(numuni,*) ' Minimum, Maximum SST jb',&
          & minval(sst%jb(1:sst%no,:)), maxval(sst%jb(1:sst%no,:))
        ENDIF
        IF( SUM( ntobs(:,nnsss) ) .GT. 0 ) THEN
          WRITE(numuni,*) ' Minimum, Maximum SSS ib',&
          & minval(sss%ib(1:sss%no,:)), maxval(sss%ib(1:sss%no,:))
          WRITE(numuni,*) ' Minimum, Maximum SSS jb',&
          & minval(sss%jb(1:sss%no,:)), maxval(sss%jb(1:sss%no,:))
        ENDIF

      ENDIF

   END SUBROUTINE prepobs

   SUBROUTINE calc_date( kit000, kt, kdate0, kdate )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE calc_date  ***
      !!
      !! ** Purpose : Compute the calendar date YYYYMMDD at a given time step.
      !!
      !! ** Method  : Compute the calendar date YYYYMMDD at a given time step.
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kit000  ! Initial time step
      INTEGER, INTENT(IN) :: kt      ! Current time step referenced to kit000
      INTEGER, INTENT(IN) :: kdate0  ! Initial date
      INTEGER, INTENT(OUT) :: kdate  ! Current date reference to kdate0
      !
      INTEGER :: iyea0    ! Initial year
      INTEGER :: imon0    ! Initial month
      INTEGER :: iday0    ! Initial day
      INTEGER :: iyea     ! Current year
      INTEGER :: imon     ! Current month
      INTEGER :: iday     ! Current day
      INTEGER :: idaystp  ! Number of days between initial and current date
      INTEGER :: idaycnt  ! Day counter

      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of curr year

      !-----------------------------------------------------------------------
      ! Compute the calendar date YYYYMMDD
      !-----------------------------------------------------------------------

      ! Initial date
      iyea0 =   kdate0 / 10000
      imon0 = ( kdate0 - ( iyea0 * 10000 ) ) / 100
      iday0 =   kdate0 - ( iyea0 * 10000 ) - ( imon0 * 100 )

      ! Check that kt >= kit000 - 1
      IF ( kt < kit000 - 1 ) CALL ctl_stop( ' kt must be >= kit000 - 1')

      ! If kt = kit000 - 1 then set the date to the restart date
      IF ( kt == kit000 - 1 ) THEN

         kdate = ndastp
         RETURN

      ENDIF

      ! Compute the number of days from the initial date
      idaystp = INT( REAL( kt - kit000 ) * rdt / 86400. )

      iday    = iday0
      imon    = imon0
      iyea    = iyea0
      idaycnt = 0

      CALL calc_month_len( iyea, imonth_len )

      DO WHILE ( idaycnt < idaystp )
         iday = iday + 1
         IF ( iday > imonth_len(imon) )  THEN
            iday = 1
            imon = imon + 1
         ENDIF
         IF ( imon > 12 ) THEN
            imon = 1
            iyea = iyea + 1
            CALL calc_month_len( iyea, imonth_len )  ! update month lengths
         ENDIF
         idaycnt = idaycnt + 1
      END DO
      !
      kdate = iyea * 10000 + imon * 100 + iday
      !
   END SUBROUTINE calc_date

   SUBROUTINE calc_month_len( iyear, imonth_len )
      !! ** Purpose : Compute the number of days in a months given a year.
      INTEGER, DIMENSION(12) ::   imonth_len
      INTEGER :: iyear         !: year
      IF ( nleapy < 2 ) THEN
         imonth_len(:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
         IF ( nleapy == 1 ) THEN   ! we are using calendar with leap years
            IF ( MOD(iyear,4) == 0 .AND. ( MOD(iyear,400) == 0 .OR. MOD(iyear,100) /= 0 ) ) THEN
               imonth_len(2) = 29
            ENDIF
         ENDIF
      ELSE
         imonth_len(:) = nleapy   ! all months with nleapy days per year
      ENDIF
      !
   END SUBROUTINE calc_month_len


   INTEGER FUNCTION max_val(in,ia)
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: in, ia(in)
      INTEGER :: i
      max_val = ia(1)
      DO  i=2,in
          max_val = MAX(max_val,ia(i))
      END DO
   END FUNCTION max_val

SUBROUTINE READOBS

!.. Read Obs, in binary or NetCDF format

IMPLICIT NONE

INTEGER :: NPROCS, JP, myproc, kk
INTEGER,ALLOCATABLE :: count(:,:)
INTEGER :: iaux, ist, ien, ier
CHARACTER(LEN=30) :: CFILE

NPROCS = jpnij
MYPROC = NAREA

ALLOCATE(count(nprocs,NoFams))
 WRITE(numuni,*)
 WRITE(numuni,*) ' GLOBAL OBS AMOUNT'
 WRITE(numuni,*) '-------------------------------------'

 CFILE='OBS_TAB.DAT'
 OPEN(913,FILE=CFILE,STATUS='OLD')
 do jp=1,nprocs
    READ(913,*) iaux,count(jp,nnINS),count(jp,nnSLA),&
                   & count(jp,nnSST),count(jp,nnSSS)
    WRITE(numuni,*) iaux, count(jp,:)
 enddo
 CLOSE(913)

 LL_SLAOBS = ( SUM( count(:,nnSLA) ) .GT. 0 )

 WRITE(numuni,*) '-------------------------------------'
 WRITE(numuni,*)
 CALL FLUSH(numuni)

  ins%no = count(myproc,nnINS)
  WRITE(numuni,*) ' INS ',ins%no
  CALL FLUSH(numuni)

  IF(ins%no .GT. 0 ) THEN

     ins%nc = count(myproc,nnINS)
     allocate(ins%ino(ins%no),ins%flg(ins%no),ins%flc(ins%no),ins%par(ins%no))
     allocate(ins%lon(ins%no),ins%lat(ins%no),ins%dpt(ins%no),ins%tim(ins%no))
     allocate(ins%val(ins%no),ins%bac(ins%no),ins%inc(ins%no),ins%tdist(ins%no))
     allocate(ins%bia(ins%no),ins%err(ins%no),ins%eve(ins%no),ins%kty(ins%no))
     allocate(ins%res(ins%no),ins%b_a(ins%no),ins%otype(ins%no),ins%prof(ins%no))
     allocate(ins%ib(ins%no,npq),ins%jb(ins%no,npq),ins%kb(ins%no),ins%plno(ins%no))
     allocate(ins%rb(ins%no),ins%inst(ins%no))
     allocate(ins%pq(ins%no,2*npq) )
     allocate(ins%moi(ins%no,npq) )
     allocate(ins%moj(ins%no,npq) )
     allocate(ins%sd1(ins%no) )
     allocate(ins%sd2(ins%no) )
     allocate(ins%prind(ins%no))
     allocate(ins%nind(ins%no))

  ENDIF

  sla%no=count(myproc,nnSLA)
  WRITE(numuni,*) ' SLA ',sla%no
  CALL FLUSH(numuni)

  IF(SLA%no .GT. 0 ) THEN

        sla%nc=count(myproc,nnSLA)
        allocate(sla%lon(sla%no),&
        &sla%lat(sla%no),&
        &sla%val(sla%no),&
        &sla%tdist(sla%no),&
        &sla%flg(sla%no),&
        &sla%flc(sla%no),&
        &sla%ksat(sla%no),&
        &sla%eve(sla%no),&
        &sla%bot(sla%no),&
        &sla%inc(sla%no),&
        &sla%ino(sla%no),&
        &sla%nind(sla%no),&
        &sla%bac(sla%no),&
        &sla%bia(sla%no),&
        &sla%err(sla%no),&
        &sla%res(sla%no),&
        &sla%b_a(sla%no),&
        &sla%ib (sla%no,npq),&
        &sla%jb (sla%no,npq),&
        &sla%pq(sla%no,npq),&
        &sla%track(sla%no),&
        &sla%dpt(sla%no),&
        &sla%tb(sla%no,jpk),&
        &sla%sb(sla%no,jpk),&
        &sla%bcp(sla%no,nprd_sla),&
        &sla%tim(sla%no))
        allocate(sla%moi(sla%no,npq) )
        allocate(sla%moj(sla%no,npq) )
        allocate(sla%sd1(sla%no) )
        allocate(sla%sd2(sla%no) )

  ENDIF

  sst%no=count(myproc,nnSST)
  WRITE(numuni,*) ' SST ',sst%no
  CALL FLUSH(numuni)

  IF(SST%no .GT. 0 .AND. .NOT. ln_noaa_sst) THEN

        sst%nc=count(myproc,nnSST)

        allocate(sst%lon(sst%no),&
        &sst%lat(sst%no),&
        &sst%val(sst%no),&
        &sst%tdist(sst%no),&
        &sst%flg(sst%no),&
        &sst%nind(sst%no),&
        &sst%flc(sst%no),&
        &sst%ksat(sst%no),&
        &sst%eve(sst%no),&
        &sst%inc(sst%no),&
        &sst%bac(sst%no),&
        &sst%bia(sst%no),&
        &sst%err(sst%no),&
        &sst%res(sst%no),&
        &sst%b_a(sst%no),&
        &sst%ib (sst%no,npq),&
        &sst%jb (sst%no,npq),&
        &sst%pq(sst%no,npq),&
        &sst%bcp(sst%no,nprd_sst),&
        &sst%track(sst%no),&
        &sst%tim(sst%no))
        allocate(sst%moi(sst%no,npq) )
        allocate(sst%moj(sst%no,npq) )
        allocate(sst%sd1(sst%no) )
        allocate(sst%sd2(sst%no) )

  ENDIF

  sss%no=count(myproc,nnSSS)
  WRITE(numuni,*) ' SSS ',sss%no
  WRITE(numuni,*) ' wih preds :',nprd_sss
  CALL FLUSH(numuni)

  IF(SSS%no .GT. 0 ) THEN

        sss%nc=count(myproc,nnSSS)

        allocate(sss%lon(sss%no),&
        &sss%lat(sss%no),&
        &sss%val(sss%no),&
        &sss%tdist(sss%no),&
        &sss%flg(sss%no),&
        &sss%flc(sss%no),&
        &sss%ksat(sss%no),&
        &sss%eve(sss%no),&
        &sss%inc(sss%no),&
        &sss%nind(sss%no),&
        &sss%bac(sss%no),&
        &sss%bia(sss%no),&
        &sss%err(sss%no),&
        &sss%res(sss%no),&
        &sss%b_a(sss%no),&
        &sss%ib (sss%no,npq),&
        &sss%jb (sss%no,npq),&
        &sss%pq(sss%no,npq),&
        &sss%bcp(sss%no,nprd_sss),&
        &sss%track(sss%no),&
        &sss%tim(sss%no))
        allocate(sss%moi(sss%no,npq) )
        allocate(sss%moj(sss%no,npq) )
        allocate(sss%sd1(sss%no) )
        allocate(sss%sd2(sss%no) )

  ENDIF

  jp=myproc
  WRITE(numuni,*) ' About to read file'
  CALL FLUSH(numuni)
  CALL FLUSH(numuni)

  ist=1
  ien=ins%no

  WRITE(numuni,*) ' Insitu part'
  CALL FLUSH(numuni)

  IF(ins%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'INSOBS_',jp-1,'.NC'
     WRITE(numuni,*) ' File ',TRIM(cfile),' ist,ien ',ist,ien
     CALL FLUSH(numuni)

     CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=COUNT(JP,NNINS),&
       NLEVS=jpk,NPREDS=0,NPQ=NPQ,NPQ2=2*NPQ,&
       FLC=INS%FLC(IST:IEN),&
       NIND=INS%NIND(IST:IEN),&
       INO=INS%INO(IST:IEN),&
       OTYPE=INS%OTYPE(IST:IEN),&
       PAR=INS%PAR(IST:IEN),&
       PLNO=INS%PLNO(IST:IEN),&
       INST=INS%INST(IST:IEN),&
       EVE=INS%EVE(IST:IEN),&
       KTY=INS%KTY(IST:IEN),&
       PRIND=INS%PRIND(IST:IEN),&
       PROF=INS%PROF(IST:IEN),&
       TDIST=INS%TDIST(IST:IEN),&
       LON=INS%LON(IST:IEN),&
       LAT=INS%LAT(IST:IEN),&
       DPT=INS%DPT(IST:IEN),&
       TIM=INS%TIM(IST:IEN),&
       VAL=INS%VAL(IST:IEN),&
       BIA=INS%BIA(IST:IEN),&
       IB=INS%IB(IST:IEN,:),&
       JB=INS%JB(IST:IEN,:),&
       KB=INS%KB(IST:IEN),&
       RB=INS%RB(IST:IEN),&
       PQ=INS%PQ(IST:IEN,:),&
       SD1=INS%SD1(IST:IEN),&
       SD2=INS%SD2(IST:IEN),&
       MOI=INS%MOI(IST:IEN,:),&
       MOJ=INS%MOJ(IST:IEN,:) )

   WRITE(numuni,*) ' File read'
   CALL FLUSH(numuni)

  ENDIF

  ist=1
  ien=sla%no

  IF(sla%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'SLAOBS_',jp-1,'.NC'
     WRITE(numuni,*) ' File ',TRIM(cfile),' ist,ien ',ist,ien
     CALL FLUSH(numuni)

     CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=COUNT(JP,NNSLA),&
     NLEVS=jpk,NPREDS=0,NPQ=NPQ,NPQ2=NPQ,&
     FLC=sla%flc(ist:ien),&
     NIND=sla%nind(ist:ien),&
     INO=sla%ino(ist:ien),&
     TRACK=sla%track(ist:ien),&
     EVE=sla%eve(ist:ien),&
     KSAT=sla%ksat(ist:ien),&
     BOT=sla%bot(ist:ien),&
     TDIST=sla%tdist(ist:ien),&
     LON=sla%lon(ist:ien),&
     LAT=sla%lat(ist:ien),&
     TIM=sla%tim(ist:ien),&
     VAL=sla%val(ist:ien),&
     BIA=sla%bia(ist:ien),&
     IB=sla%ib(ist:ien,:),&
     JB=sla%jb(ist:ien,:),&
     PQ=sla%pq(ist:ien,:),&
     SD1=SLA%SD1(IST:IEN),&
     SD2=SLA%SD2(IST:IEN),&
     MOI=SLA%MOI(IST:IEN,:),&
     MOJ=SLA%MOJ(IST:IEN,:),&
     DPT=sla%dpt(ist:ien) )

  ENDIF

  ist=1
  ien=sst%no

  IF(sst%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'SSTOBS_',jp-1,'.NC'
     WRITE(numuni,*) ' File ',TRIM(cfile),' ist,ien ',ist,ien
     CALL FLUSH(numuni)

     CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=COUNT(JP,NNSST),&
     NLEVS=jpk,NPREDS=0,NPQ=NPQ,NPQ2=NPQ,&
     FLC=sst%flc(ist:ien),&
     TRACK=sst%track(ist:ien),&
     EVE=sst%eve(ist:ien),&
     NIND=sst%nind(ist:ien),&
     KSAT=sst%ksat(ist:ien),&
     TDIST=sst%tdist(ist:ien),&
     LON=sst%lon(ist:ien),&
     LAT=sst%lat(ist:ien),&
     TIM=sst%tim(ist:ien),&
     VAL=sst%val(ist:ien),&
     BIA=sst%bia(ist:ien),&
     IB=sst%ib(ist:ien,:),&
     JB=sst%jb(ist:ien,:),&
     SD1=SST%SD1(IST:IEN),&
     SD2=SST%SD2(IST:IEN),&
     MOI=SST%MOI(IST:IEN,:),&
     MOJ=SST%MOJ(IST:IEN,:),&
     PQ=sst%pq(ist:ien,:) )

  ENDIF

  ist=1
  ien=sss%no

  IF(sss%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'SSSOBS_',jp-1,'.NC'
     WRITE(numuni,*) ' File ',TRIM(cfile),' ist,ien ',ist,ien
     CALL FLUSH(numuni)

     CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=COUNT(JP,NNSSS),&
     NLEVS=jpk,NPREDS=NPRD_SSS,NPQ=NPQ,NPQ2=NPQ,&
     NIND=sss%nind(ist:ien),&
     FLC=sss%flc(ist:ien),&
     TRACK=sss%track(ist:ien),&
     EVE=sss%eve(ist:ien),&
     KSAT=sss%ksat(ist:ien),&
     TDIST=sss%tdist(ist:ien),&
     LON=sss%lon(ist:ien),&
     LAT=sss%lat(ist:ien),&
     TIM=sss%tim(ist:ien),&
     VAL=sss%val(ist:ien),&
     BIA=sss%bia(ist:ien),&
     IB=sss%ib(ist:ien,:),&
     JB=sss%jb(ist:ien,:),&
     BCP=sss%bcp(ist:ien,:),&
     SD1=SSS%SD1(IST:IEN),&
     SD2=SSS%SD2(IST:IEN),&
     MOI=SSS%MOI(IST:IEN,:),&
     MOJ=SSS%MOJ(IST:IEN,:),&
     PQ=sss%pq(ist:ien,:))

  ENDIF

END SUBROUTINE READOBS

SUBROUTINE WRITEOBS

!.. Write Obs, in binary or NetCDF format

IMPLICIT NONE

INTEGER :: NPROCS, JP, NFAULT, i, j, i2, myproc, kk
CHARACTER(LEN=30) :: CFILE

  MYPROC = NAREA

  jp=myproc

  IF(ins%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'INSMIS_',jp-1,'.NC'

     CALL IO_OBS(CFILE=CFILE,IOPT=1,NOBS=ins%no,&
     NLEVS=jpk,NPREDS=0,NPQ=NPQ,NPQ2=2*NPQ,&
     flc=ins%flc(1:ins%no),&
     nind=ins%nind(1:ins%no),&
     ino=ins%ino(1:ins%no),&
     otype=ins%otype(1:ins%no),&
     par=ins%par(1:ins%no),&
     plno=ins%plno(1:ins%no),&
     inst=ins%inst(1:ins%no),&
     eve=ins%eve(1:ins%no),&
     kty=ins%kty(1:ins%no),&
     prind=ins%prind(1:ins%no),&
     prof=ins%prof(1:ins%no),&
     tdist=ins%tdist(1:ins%no),&
     lon=ins%lon(1:ins%no),&
     lat=ins%lat(1:ins%no),&
     dpt=ins%dpt(1:ins%no),&
     tim=ins%tim(1:ins%no),&
     val=ins%val(1:ins%no),&
     bia=ins%bia(1:ins%no),&
     ib=ins%ib(1:ins%no,:),&
     jb=ins%jb(1:ins%no,:),&
     kb=ins%kb(1:ins%no),&
     rb=ins%rb(1:ins%no),&
     pq=ins%pq(1:ins%no,:),&
     bac=ins%bac(1:ins%no),&
     res=ins%res(1:ins%no) )

  ENDIF

  IF(sla%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'SLAMIS_',jp-1,'.NC'

     CALL IO_OBS(CFILE=CFILE,IOPT=1,NOBS=sla%no,&
     NLEVS=jpk,NPREDS=NPRD_SLA,NPQ=NPQ,NPQ2=NPQ,&
     flc=sla%flc(1:sla%no),&
     nind=sla%nind(1:sla%no),&
     ino=sla%ino(1:sla%no),&
     track=sla%track(1:sla%no),&
     eve=sla%eve(1:sla%no),&
     ksat=sla%ksat(1:sla%no),&
     bot=sla%bot(1:sla%no),&
     tdist=sla%tdist(1:sla%no),&
     lon=sla%lon(1:sla%no),&
     lat=sla%lat(1:sla%no),&
     tim=sla%tim(1:sla%no),&
     val=sla%val(1:sla%no),&
     bia=sla%bia(1:sla%no),&
     ib=sla%ib(1:sla%no,:),&
     jb=sla%jb(1:sla%no,:),&
     pq=sla%pq(1:sla%no,:),&
     dpt=sla%dpt(1:sla%no),&
     bac=sla%bac(1:sla%no),&
     tb=sla%tb(1:sla%no,1:jpk),&
     sb=sla%sb(1:sla%no,1:jpk),&
     bcp=sla%bcp(1:sla%no,1:nprd_sla),&
     res=sla%res(1:sla%no) )

  ENDIF

  IF(sst%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'SSTMIS_',jp-1,'.NC'

     CALL IO_OBS(CFILE=CFILE,IOPT=1,NOBS=sst%no,&
     NLEVS=jpk,NPREDS=NPRD_SST,NPQ=NPQ,NPQ2=NPQ,&
     flc=sst%flc(1:sst%no),&
     nind=sst%nind(1:sla%no),&
     track=sst%track(1:sst%no),&
     eve=sst%eve(1:sst%no),&
     ksat=sst%ksat(1:sst%no),&
     tdist=sst%tdist(1:sst%no),&
     lon=sst%lon(1:sst%no),&
     lat=sst%lat(1:sst%no),&
     tim=sst%tim(1:sst%no),&
     val=sst%val(1:sst%no),&
     bia=sst%bia(1:sst%no),&
     ib=sst%ib(1:sst%no,:),&
     jb=sst%jb(1:sst%no,:),&
     pq=sst%pq(1:sst%no,:),&
     bcp=sst%bcp(1:sst%no,1:nprd_sst),&
     err=sst%err(1:sst%no),&
     bac=sst%bac(1:sst%no),&
     res=sst%res(1:sst%no) )

  ENDIF

  IF(sss%no .GT. 0) THEN

     Write(cfile,'(A,I4.4,A)') 'SSSMIS_',jp-1,'.NC'

     CALL IO_OBS(CFILE=CFILE,IOPT=1,NOBS=sss%no,&
     NLEVS=jpk,NPREDS=NPRD_SSS,NPQ=NPQ,NPQ2=NPQ,&
     flc=sss%flc(1:sss%no),&
     nind=sss%nind(1:sla%no),&
     track=sss%track(1:sss%no),&
     eve=sss%eve(1:sss%no),&
     ksat=sss%ksat(1:sss%no),&
     tdist=sss%tdist(1:sss%no),&
     lon=sss%lon(1:sss%no),&
     lat=sss%lat(1:sss%no),&
     tim=sss%tim(1:sss%no),&
     val=sss%val(1:sss%no),&
     bia=sss%bia(1:sss%no),&
     bcp=sss%bcp(1:sss%no,:),&
     ib=sss%ib(1:sss%no,:),&
     jb=sss%jb(1:sss%no,:),&
     pq=sss%pq(1:sss%no,:),&
     bac=sss%bac(1:sss%no),&
     res=sss%res(1:sss%no) )

  ENDIF

END SUBROUTINE WRITEOBS

#include "io_obs.h90"

  SUBROUTINE CHECK(STATUS)
    USE mpi
    IMPLICIT NONE
    INTEGER, INTENT ( IN) :: STATUS
    INTEGER :: IERR

    IF(STATUS /= NF90_NOERR) THEN
      WRITE(numuni,*) 'NETCDF ERROR : ',TRIM(NF90_STRERROR(STATUS))
      WRITE(numuni,*) 'STOPPING PROGRAM'
      CALL FLUSH(numuni)
      CALL mppstop
    END IF
  END SUBROUTINE CHECK

! OSUM : AUXILIARY FUNCTION TO INTERPOLATE
!        OVER MODEL POINTS
!        ZR : WEIGTHS VECTOR
!        ZR2: DIAGONAL MATRIX OF MODEL VALUES
!              ONLY DIAGONAL ELEMENTS ARE RELEVANT

    REAL(R8) FUNCTION OSUM(ZR,ZR2)
       IMPLICIT NONE
       REAL(R8), INTENT(IN) :: ZR(NPQ)
       REAL(R8), INTENT(IN) :: ZR2(NPQ,NPQ)
       INTEGER :: J
       OSUM=0._R8
       DO J=1,NPQ
          OSUM = OSUM + ZR(J)*ZR2(J,J)
       ENDDO
   END FUNCTION OSUM

! Debugging version

    REAL(wp) FUNCTION OSUM2(ZR,ZR2)
       IMPLICIT NONE
       REAL(wp), INTENT(IN) :: ZR(NPQ)
       REAL(wp), INTENT(IN) :: ZR2(NPQ,NPQ)
       INTEGER :: J
       OSUM2=0._R8
       DO J=1,NPQ
          OSUM2 = OSUM2 + ZR(J)*ZR2(J,J)
          WRITE(1919,*) J, ZR(J),ZR2(J,J),OSUM2
       ENDDO
   END FUNCTION OSUM2

   REAL(wp) FUNCTION mpp_misfits_ts2(sd1,sd2,ib,jb,kp,kl,pq,mi,mj,ind)

   USE MPI

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: sd1, sd2, ib(npq),jb(npq),kp,kl,mi(npq),mj(npq), ind
   REAL(wp), INTENT(IN):: pq(:)

   REAL(wp) :: zt, zff(2)
   INTEGER  :: jp, jp2, kl2, stat, ierr, nundb, i0,j0, iopt, ncc, nrr
   LOGICAL  :: lmpp, lins, lldbg
   LOGICAL, PARAMETER :: lldebug = .true.

   lins = ( SIZE(pq(:)) .EQ. 2*npq )
   kl2 = kl + 1
   zt = 0._wp
   lldbg = .false.
   nundb = 7000+narea
   IF(narea.eq.sd1 .and. lldebug) lldbg = .true.

   IF(lldbg) WRITE(nundb,*)  '>>> ',sd1,sd2

   DO jp=1,npq

      lmpp=.false.
      DO jp2=1,npq
         if( ib(jp) .eq. mi(jp2) .AND. jb(jp).eq. mj(jp2) ) lmpp=.true.
      ENDDO

      IF(lmpp) Then
         iopt = 0
         if (sd2 - sd1 .eq. noea - nproc ) iopt=1    ! Fetch from east
         if (sd2 - sd1 .eq. nono - nproc ) iopt=2 ! Fetch from north
         if (iopt .eq. 0) THEN
            WRITE(*,*) 'This not possible ', sd1, sd2
            call mppstop
         Endif
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      *** ',jp,ib(jp),jb(jp),kl,kl2,kp,lmpp,iopt

      IF(lmpp) then
          if (iopt .eq. 1 ) Then
               ncc = nimppt(narea)+nlcit(narea)-nimppt(noea+1)
               if ( ncc .ge. jpiglo ) ncc = 3
               i0 = ib(jp)-nimppt(noea+1)-ncc+2
               j0 = jb(jp)-njmpp+1
               IF(lldbg) WRITE(nundb,*)  '         --- ',jb(jp),j0,kl,kp,njmpp
               zff(1) = tsn_e(i0,j0,kl,kp)
               if(lins) zff(2) = tsn_e(i0,j0,kl2,kp)
          Else
               nrr = njmppt(narea)+nlcjt(narea)-njmppt(nono+1)
               i0 = ib(jp)-nimpp+1
               j0 = jb(jp)-njmppt(nono+1)-nrr+2
               IF(lldbg) WRITE(nundb,*)  '         --- ',ib(jp),i0,kl,kp,nimpp
               zff(1) = tsn_n(i0,j0,kl,kp)
               if(lins) zff(2) = tsn_n(i0,j0,kl2,kp)
          Endif
      ELSE
          i0 = ib(jp)-nimpp+1
          j0 = jb(jp)-njmpp+1
          IF(lldbg) WRITE(nundb,*)  '         --- ',i0,j0,kl,kp
          zff(1) = amts_m(i0,j0,kl,kp)
          if (lins) zff(2) = amts_m(i0,j0,kl2,kp)
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      +++ ',zff(1:2)
      zt = zt + pq(jp)*zff(1)
      if (lins) zt = zt + pq(npq+jp)*zff(2)

   ENDDO

   mpp_misfits_ts2 = zt

   IF(lldbg) WRITE(nundb,*)  '      === ',zt
   IF(lldbg) CALL FLUSH(nundb)

   END FUNCTION mpp_misfits_ts2

   REAL(wp) FUNCTION mpp_misfits_ts(sd1,sd2,ib,jb,kp,kl,pq,mi,mj,ind)

   USE MPI

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: sd1, sd2, ib(npq),jb(npq),kp,kl,mi(npq),mj(npq), ind
   REAL(wp), INTENT(IN):: pq(:)

   REAL(wp) :: zt, zff(2)
   INTEGER  :: jp, jp2, kl2, stat, ierr, nundb, i0,j0, iopt, ncc, nrr
   LOGICAL  :: lmpp, lins, lldbg
   LOGICAL, PARAMETER :: lldebug = .false.

   lins = ( SIZE(pq(:)) .EQ. 2*npq )
   kl2 = kl + 1
   zt = 0._wp
   lldbg = .false.
   nundb = 2000+narea
   IF(narea.eq.sd1 .and. lldebug) lldbg = .true.

   IF(lldbg) WRITE(nundb,*)  '>>> ',sd1,sd2

   DO jp=1,npq

      lmpp=.false.
      DO jp2=1,npq
         if( ib(jp) .eq. mi(jp2) .AND. jb(jp).eq. mj(jp2) ) lmpp=.true.
      ENDDO

      IF(lmpp) Then
         iopt = 0
         if (sd2 - sd1 .eq. noea - nproc ) iopt=1    ! Fetch from east
         if (sd2 - sd1 .eq. nono - nproc ) iopt=2 ! Fetch from north
         if (iopt .eq. 0) THEN
            WRITE(*,*) 'This not possible  :: TS ', sd1, sd2, ib, jb, narea, ind
            call mppstop
         Endif
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      *** ',jp,ib(jp),jb(jp),kl,kl2,kp,lmpp,iopt

      IF(lmpp) then
          if (iopt .eq. 1 ) Then
               ncc = nimppt(narea)+nlcit(narea)-nimppt(noea+1)
               if ( ncc .ge. jpiglo ) ncc = 3
               i0 = ib(jp)-nimppt(noea+1)-ncc+2
               j0 = jb(jp)-njmpp+1
               IF(lldbg) WRITE(nundb,*)  '         --- ',jb(jp),j0,kl,kp,njmpp
               zff(1) = tsn_e(i0,j0,kl,kp)
               if(lins) zff(2) = tsn_e(i0,j0,kl2,kp)
          Else
               nrr = njmppt(narea)+nlcjt(narea)-njmppt(nono+1)
               i0 = ib(jp)-nimpp+1
               j0 = jb(jp)-njmppt(nono+1)-nrr+2
               IF(lldbg) WRITE(nundb,*)  '         --- ',ib(jp),i0,kl,kp,nimpp
               zff(1) = tsn_n(i0,j0,kl,kp)
               if(lins) zff(2) = tsn_n(i0,j0,kl2,kp)
          Endif
      ELSE
          i0 = ib(jp)-nimpp+1
          j0 = jb(jp)-njmpp+1
          IF(lldbg) WRITE(nundb,*)  '         --- ',i0,j0,kl,kp
          zff(1) = amts_m(i0,j0,kl,kp)
          if (lins) zff(2) = amts_m(i0,j0,kl2,kp)
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      +++ ',zff(1:2)
      zt = zt + pq(jp)*zff(1)
      if (lins) zt = zt + pq(npq+jp)*zff(2)

   ENDDO

   mpp_misfits_ts = zt

   IF(lldbg) WRITE(nundb,*)  '      === ',zt
   IF(lldbg) CALL FLUSH(nundb)

   END FUNCTION mpp_misfits_ts

   REAL(wp) FUNCTION mpp_misfits_tq(sd1,sd2,ib,jb,kp,pq,mi,mj,ind)

   USE MPI

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: sd1, sd2, ib(npq),jb(npq),kp,mi(npq),mj(npq), ind
   REAL(wp), INTENT(IN):: pq(:)

   REAL(wp) :: zt, zff(2)
   INTEGER  :: jp, jp2, stat, ierr, nundb, i0,j0, iopt, ncc, nrr
   LOGICAL  :: lmpp, lins, lldbg
   LOGICAL, PARAMETER :: lldebug = .false.

   zt = 0._wp
   lldbg = .false.
   nundb = 2000+narea
   IF(narea.eq.sd1 .and. lldebug) lldbg = .true.

   IF(lldbg) WRITE(nundb,*)  '>>> ',sd1,sd2

   DO jp=1,npq

      lmpp=.false.
      DO jp2=1,npq
         if( ib(jp) .eq. mi(jp2) .AND. jb(jp).eq. mj(jp2) ) lmpp=.true.
      ENDDO

      IF(lmpp) Then
         iopt = 0
         if (sd2 - sd1 .eq. noea - nproc ) iopt=1    ! Fetch from east
         if (sd2 - sd1 .eq. nono - nproc ) iopt=2 ! Fetch from north
         if (iopt .eq. 0) THEN
            WRITE(*,*) 'This not possible  :: TS ', sd1, sd2, ib, jb, narea, ind
            call mppstop
         Endif
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      *** ',jp,ib(jp),jb(jp),kp,lmpp,iopt

      IF(lmpp) then
          if (iopt .eq. 1 ) Then
               ncc = nimppt(narea)+nlcit(narea)-nimppt(noea+1)
               if ( ncc .ge. jpiglo ) ncc = 3
               i0 = ib(jp)-nimppt(noea+1)-ncc+2
               j0 = jb(jp)-njmpp+1
               IF(lldbg) WRITE(nundb,*)  '         --- ',jb(jp),j0,kp,njmpp
               zff(1) = tq2_e(i0,j0,kp)
          Else
               nrr = njmppt(narea)+nlcjt(narea)-njmppt(nono+1)
               i0 = ib(jp)-nimpp+1
               j0 = jb(jp)-njmppt(nono+1)-nrr+2
               IF(lldbg) WRITE(nundb,*)  '         --- ',ib(jp),i0,kp,nimpp
               zff(1) = tq2_n(i0,j0,kp)
          Endif
      ELSE
          i0 = ib(jp)-nimpp+1
          j0 = jb(jp)-njmpp+1
          IF(lldbg) WRITE(nundb,*)  '         --- ',i0,j0,kp
          zff(1) = amtq_m(i0,j0,kp)
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      +++ ',zff(1:2)
      zt = zt + pq(jp)*zff(1)

   ENDDO

   mpp_misfits_tq = zt

   IF(lldbg) WRITE(nundb,*)  '      === ',zt
   IF(lldbg) CALL FLUSH(nundb)

   END FUNCTION mpp_misfits_tq


   REAL(wp) FUNCTION mpp_misfits_sla(sd1,sd2,ib,jb,pq,mi,mj,ind,KKnn)

   USE MPI

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: sd1, sd2, ib(npq),jb(npq),mi(npq),mj(npq), ind,KKnn
   REAL(wp), INTENT(IN):: pq(:)

   REAL(wp) :: zt, zff
   INTEGER  :: jp, jp2, stat, ierr, nundb, i0,j0, iopt, ncc, nrr
   LOGICAL  :: lmpp, lldbg
   LOGICAL, PARAMETER :: lldebug = .false.

   zt  = 0._wp
   zff = 0._wp
   lldbg = .false.
   nundb = 2000+narea
   IF(narea.eq.sd1 .and. lldebug) lldbg = .true.

   IF(lldbg) WRITE(nundb,*)  '>>> ',sd1,sd2

   DO jp=1,npq

      lmpp=.false.
      DO jp2=1,npq
         if( ib(jp) .eq. mi(jp2) .AND. jb(jp).eq. mj(jp2) ) lmpp=.true.
      ENDDO

      IF(lmpp) Then
         iopt = 0
         if (sd2 - sd1 .eq. noea - nproc ) iopt=1    ! Fetch from east
         if (sd2 - sd1 .eq. nono - nproc ) iopt=2 ! Fetch from north
         if (iopt .eq. 0) THEN
            WRITE(*,*) 'This not possible - SLA', sd1, sd2, ib, jb, narea, ind
            call mppstop
         Endif
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      *** ',jp,ib(jp),jb(jp),lmpp,iopt

      IF(lmpp) then
          if (iopt .eq. 1 ) Then
               IF(lldbg) WRITE(nundb,*)  ' === ib(jp), jb(jp), nimppt(narea), nlcit(narea), nimppt(noea+1)   ', &
                              &     ib(jp), jb(jp), nimppt(narea), nlcit(narea), nimppt(noea+1)
               ncc = nimppt(narea)+nlcit(narea)-nimppt(noea+1)
               if ( ncc .ge. jpiglo ) ncc = 3
               i0 = ib(jp)-nimppt(noea+1)-ncc+2
               j0 = jb(jp)-njmpp+1
               IF(lldbg) WRITE(nundb,*)  '         --- ',jb(jp),j0,njmpp
               zff = sla_e(i0,j0)
               IF(lldbg) WRITE(nundb,*)  '      === zff east is   ', zff
          Else
               IF(lldbg) WRITE(nundb,*)  ' === ib(jp), jb(jp), njmppt(narea), nlcjt(narea), njmppt(nono+1):   ', &
                              &     ib(jp), jb(jp), njmppt(narea), nlcjt(narea), njmppt(nono+1)
               nrr = njmppt(narea)+nlcjt(narea)-njmppt(nono+1)
               i0 = ib(jp)-nimpp+1
               j0 = jb(jp)-njmppt(nono+1)-nrr+2
               IF(lldbg) WRITE(nundb,*)  '         --- ',ib(jp),i0,nimpp
               zff = sla_n(i0,j0)
               IF(lldbg) WRITE(nundb,*)  '      === zff north is   ', zff
          Endif
      ELSE
          i0 = ib(jp)-nimpp+1
          j0 = jb(jp)-njmpp+1
          IF(lldbg) WRITE(nundb,*)  '         --- ',i0,j0
          zff = ssh(i0,j0,KKnn)-mdt(i0,j0)
      ENDIF

      IF(lldbg) WRITE(nundb,*)  '      +++ ',zff
      zt = zt + pq(jp)*zff

   ENDDO

   mpp_misfits_sla = zt

   IF(lldbg) WRITE(nundb,*)  '      === ',zt
   IF(lldbg) CALL FLUSH(nundb)

   END FUNCTION mpp_misfits_sla

   SUBROUTINE prep_mppmisfits(KKnn)

      USE MPI

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: KKnn
      INTEGER  :: k, ji, jj, jk, jp, narea2, nc, nr
      INTEGER  :: stat, ierr, tag1, tag2
      INTEGER  :: isend_w, irecv_e, isend_s, irecv_n
      INTEGER  :: istatus(mpi_status_size)
      LOGICAL  :: llmpdebug

      tag1=1000+narea
      tag2=2000+narea
      llmpdebug = l_mfc
      IF ( l_mfc ) THEN
           WRITE(numuni,*) ' narea, nproc, noso, nowe, nono, noea, ln_sf, ln_nf, ln_wf, ln_ef'
           WRITE(numuni,*)   narea, nproc, noso, nowe, nono, noea, ln_sf, ln_nf, ln_wf, ln_ef
      ENDIF
      CALL FLUSH(numuni)

      ! Decide which row/column must be sent

      ! Column To West
      if( .not. ln_wf ) Then
        narea2=nowe+1
        nc = nimppt(narea2)+nlcit(narea2)-nimpp
        ! Handle the periodicity
        if ( nc .ge. jpiglo ) nc = nc - (jpiglo-2) +1
      else
        nc = 2
      endif

      ! Row To South
      if( .not. ln_sf ) Then
        narea2=noso+1
        nr = njmppt(narea2)+nlcjt(narea2)-njmpp
      else
        nr = 2
      endif

      IF( llmpdebug ) THEN
        WRITE(numuni,*) ' nc, nr ', nc, nr, &
        & size(west),size(amts_m,1),size(amts_m,2),size(amts_m,3)
        CALL FLUSH(numuni)
      ENDIF

      if( .not. ln_wf ) THEN
        k=0
        do jp=1,2
          do jk=1,jpk
            do jj=1,jpj
              do ji=1,3
                k=k+1
                west (k) = amts_m(nc+ji-1,jj,jk,jp)
              enddo
            enddo
          enddo
        enddo
        do jj=1,jpj
          do ji=1,3
            k=k+1
            west (k) = ssh(nc+ji-1,jj,KKnn)-mdt(nc+ji-1,jj)
          enddo
        enddo
        do jp=1,2
            do jj=1,jpj
              do ji=1,3
                k=k+1
                west (k) = amtq_m(nc+ji-1,jj,jp)
              enddo
            enddo
        enddo
      endif

      if( .not. ln_wf ) THEN
      IF( llmpdebug ) THEN
       WRITE(numuni,*) ' Sending west to ',nowe ,' global index is ',nc+nimpp-1 ; CALL FLUSH(numuni)
      ENDIF
      CALL mpi_isend( west, imwe, mpi_double_precision, nowe, tag1, mpi_comm_oce, isend_w, ierr)
      IF( llmpdebug ) THEN
       WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
      ENDIF
      endif

      if( .not. ln_ef ) THEN
      IF( llmpdebug ) THEN
       WRITE(numuni,*) ' Receiving east from ',noea ; CALL FLUSH(numuni)
      ENDIF
      CALL mpi_irecv(east,  imwe, mpi_double_precision, noea, mpi_any_tag, mpi_comm_oce, irecv_e, ierr)
      IF( llmpdebug ) THEN
       WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
      ENDIF
      endif

      If( .not. ln_sf ) Then
        k=0
        do jp=1,2
          do jk=1,jpk
            do jj=1,3
              do ji=1,jpi
                k=k+1
                south(k) = amts_m(ji,nr+jj-1,jk,jp)
              enddo
            enddo
          enddo
        enddo
        do jj=1,3
          do ji=1,jpi
            k=k+1
            south (k) = ssh(ji,nr+jj-1,KKnn)-mdt(ji,nr+jj-1)
          enddo
        enddo

        do jp=1,2
            do jj=1,3
              do ji=1,jpi
                k=k+1
                south(k) = amtq_m(ji,nr+jj-1,jp)
              enddo
            enddo
        enddo

      Endif

      if( .not. ln_wf ) Then
        IF( llmpdebug ) THEN
         WRITE(numuni,*) ' Waiting send west ' ; CALL FLUSH(numuni)
        ENDIF
        CALL mpi_wait(isend_w, istatus, ierr)
        IF( llmpdebug ) THEN
         WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
        ENDIF
      endif

      if( .not. ln_ef ) Then
        IF( llmpdebug ) THEN
         WRITE(numuni,*) ' Waiting recv east ' ; CALL FLUSH(numuni)
        ENDIF
        CALL mpi_wait(irecv_e, istatus, ierr)
        IF( llmpdebug ) THEN
         WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
        ENDIF
      endif

      if( .not. ln_sf ) THEN
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Sending south to ',noso, ' global index is ',nr+njmpp-1 ; CALL FLUSH(numuni)
         ENDIF
         CALL mpi_isend(south, imsn, mpi_double_precision, noso, tag2, mpi_comm_oce, isend_s,ierr)
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
         ENDIF
      endif

      if( .not. ln_nf ) THEN
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Receiving north from ',nono ; CALL FLUSH(numuni)
         ENDIF
         CALL mpi_irecv(north,  imsn, mpi_double_precision, nono, mpi_any_tag, mpi_comm_oce, irecv_n, ierr)
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
         ENDIF
      endif

      if( .not. ln_ef ) Then
        k=0
        do jp=1,2
          do jk=1,jpk
            do jj=1,jpj
              do ji=1,3
                k=k+1
                tsn_e(ji,jj,jk,jp) = east (k)
              enddo
            enddo
          enddo
        enddo
        do jj=1,jpj
          do ji=1,3
            k=k+1
            sla_e(ji,jj) = east (k)
          enddo
        enddo

        do jp=1,2
            do jj=1,jpj
              do ji=1,3
                k=k+1
                tq2_e(ji,jj,jp) = east (k)
              enddo
            enddo
        enddo

      endif

      if( .not. ln_sf ) THEN
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Waiting send south' ; CALL FLUSH(numuni)
         ENDIF
         CALL mpi_wait(isend_s, istatus, ierr)
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
         ENDIF
      endif

      if( .not. ln_nf ) THEN
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Waiting recv north' ; CALL FLUSH(numuni)
         ENDIF
         CALL mpi_wait(irecv_n, istatus, ierr)
         IF( llmpdebug ) THEN
          WRITE(numuni,*) ' Done'  ; CALL FLUSH(numuni)
         ENDIF
      endif

      if( .not. ln_nf ) Then
        k=0
        do jp=1,2
          do jk=1,jpk
            do jj=1,3
              do ji=1,jpi
                k=k+1
                tsn_n(ji,jj,jk,jp) = north(k)
              enddo
            enddo
          enddo
        enddo
        do jj=1,3
          do ji=1,jpi
            k=k+1
            sla_n(ji,jj) = north (k)
          enddo
        enddo

        do jp=1,2
            do jj=1,3
              do ji=1,jpi
                k=k+1
                tq2_n(ji,jj,jp) = north(k)
              enddo
            enddo
        enddo

      Endif

      IF ( l_mfc ) THEN
           WRITE(numuni,*) ' Exiting from first prep_mppmisfits'
           CALL FLUSH(numuni)
           l_mfc = .FALSE.
      ENDIF

   END SUBROUTINE prep_mppmisfits

   REAL(KIND=wp) FUNCTION dep_to_p( p_dep, p_lat )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION dep_to_p  ***
      !!
      !! ** Purpose : Compute depth from pressure and latitudes
      !!
      !! ** Method  : The expression used in p_to_dep is inverted.
      !!
      !! ** Action  :
      !!
      !! References : P.M Saunders
      !!              Pratical conversion of pressure to depth
      !!              Journal of physical oceanography Vol 11, 1981, pp 573-574
      !!
      !! History :
      !!        !  07-05  (K. Mogensen) Original code
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=wp), INTENT(IN) :: p_dep    ! Depth in meters
      REAL(KIND=wp), INTENT(IN) :: p_lat    ! Latitude in degrees

      !! * Local declarations
      REAL(KIND=wp) :: z_x
      REAL(KIND=wp) :: z_c1
      REAL(KIND=wp) :: z_c2
      REAL(KIND=wp) :: z_d

      z_x = SIN( p_lat / 57.29578 )
      z_x = z_x * z_x
      z_c1 = ( 5.92  + 5.25 * z_x ) * 1e-3
      z_c2 = 2.21e-6
      z_d = ( z_c1 - 1 ) * ( z_c1 - 1  ) - 4 * z_c2 * p_dep
      dep_to_p = (( 1 - z_c1 ) - SQRT( z_d )) / ( 2 * z_c2 )

   END FUNCTION dep_to_p

  REAL(KIND=wp) FUNCTION potemp( ps, pt, pp, ppr )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION potemp  ***
      !!
      !! ** Purpose : Compute potential temperature
      !!
      !! ** Method  : A regression formula is used.
      !!
      !! ** Action  : The code is kept as close to the F77 code as possible
      !!              Check value: potemp(35,20,2000,0) = 19.621967
      !!
      !! References : T. J. Mcdougall, D. R. Jackett, D. G. Wright
      !!              and R. Feistel
      !!              Accurate and computationally efficient algoritms for
      !!              potential temperatures and density of seawater
      !!              Journal of atmospheric and oceanic technology
      !!              Vol 20, 2003, pp 730-741
      !!
      !!
      !! History :
      !!        !  07-05 (K. Mogensen) Original code
      !!----------------------------------------------------------------------

      !! * Arguments

      REAL(KIND=wp), INTENT(IN) :: ps
      REAL(KIND=wp), INTENT(IN) :: pt
      REAL(KIND=wp), INTENT(IN) :: pp
      REAL(KIND=wp), INTENT(IN) :: ppr

      !! * Local declarations
      REAL(KIND=wp) :: zpol
      REAL(KIND=wp), PARAMETER :: a1 =  1.067610e-05
      REAL(KIND=wp), PARAMETER :: a2 = -1.434297e-06
      REAL(KIND=wp), PARAMETER :: a3 = -7.566349e-09
      REAL(KIND=wp), PARAMETER :: a4 = -8.535585e-06
      REAL(KIND=wp), PARAMETER :: a5 =  3.074672e-08
      REAL(KIND=wp), PARAMETER :: a6 =  1.918639e-08
      REAL(KIND=wp), PARAMETER :: a7 =  1.788718e-10

      zpol = a1 + a2 * ps + a3 * ( pp + ppr ) + a4 * pt &
         & + a5 * ps * pt + a6 * pt * pt + a7 * pt * ( pp + ppr )

      potemp = pt + ( pp - ppr ) * zpol

   END FUNCTION potemp

   SUBROUTINE save_background(KKnn)

   INTEGER, INTENT(IN) :: KKnn
   REAL(wp), ALLOCATABLE :: buff(:,:,:,:)
   INTEGER, PARAMETER :: ufd = 9439, nff=2
   CHARACTER(LEN=99) :: cfd

   ALLOCATE( buff(jpi,jpj,jpk,nff) )

   buff(:,:,:,1) = ts(:,:,:,jp_tem,KKnn)
   buff(:,:,:,2) = ts(:,:,:,jp_sal,KKnn)

   WRITE(cfd,'(A,I4.4,A)') 'background_',narea-1,'.bin'
   OPEN(ufd,FILE=cfd,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',&
   & STATUS='NEW')
   WRITE(ufd) jpi, jpj, jpk, nff
   WRITE(ufd) REAL(buff,KIND=sp)
   CLOSE(ufd)

   DEALLOCATE( buff )

   END SUBROUTINE save_background

END MODULE obsmisfits
