!INPUT:
! 1-file name 2-Date

      include './netcdf.inc'
      parameter(nlns=10000000)
      integer :: start(1),count(1),startb(1),countb(1),dimen(1)
      integer :: lend,lendb,dimidpb,dimidp,dimidot,idslab,dimidit
      integer ::  l, ld, np, rejcount, nobs, noqc, rr
      integer :: tmeid, spdid,spdid2, dateb, i, jjj
      real rval
      real , dimension (:), allocatable :: dx, dy, dist, rdac,rocti
      real, dimension (:), allocatable :: sla, dac, flagb, slaf,octi
      real, dimension (:), allocatable :: rinti, inti
      real, dimension (:), allocatable :: slapdac, slafb,slapdacpot
      double precision, dimension (:), allocatable :: data, datao
      integer, dimension(:), allocatable :: lat, lon, lonb, latb
      real, dimension(:), allocatable :: alonb, alatb
      integer, parameter:: zzr=345, ttr=129
      real, dimension(zzr,ttr) :: navlonr, navlatr, ref
      integer refs(zzr,ttr)
      real, dimension(zzr) :: lonr
      real, dimension(ttr) :: latr
      integer  indexsat, tt, zz
      INTEGER, DIMENSION(1) :: ar
      logical nosub
      character*256 :: indir, outdir,infile, listfile, infileb
      character*1 cindexsat
      character*2 cinm, cind, chrs,cwin
      character*4 ciny
      real rdy(nlns), aln(nlns), alt(nlns), numv(nlns), rvl(nlns), err(nlns)
      integer   ino(nlns)
      integer   nam, idrefr, ncir, scd
      integer yyyy, mm, dd, hh, mn
      character*3 sat
      character*8 date
      character*3 satin
      logical :: file_exists
!-----------------------------------------------------
      if(iargc().ne.5)stop 'Stop wrong number of arguments'
      call getarg(1,indir)
      call getarg(2,outdir)
      call getarg(3,infile)
      call getarg(4,date)
      call getarg(5,satin)

      nobs = 0
      noqc = 0
      print*,'read_data'
      !!!!!!!!!Cambiare riga dopo dipende da path del file
      sat=satin(1:3)
      print*,sat
      if(sat=="s3a") then
         indexsat=1
      else if(sat=="c2 ") then
         indexsat=2
      else if(sat=="c2n") then
         indexsat=2
      else if(sat=="j3 ") then
         indexsat=3
      else if(sat=="j3n") then
         indexsat=3
      else if(sat=="j1 ") then
         indexsat=4
      else if(sat=="j2 ") then
         indexsat=5
      else if(sat=="j2n") then
         indexsat=5
      else if(sat=="j2g") then
         indexsat=5
      else if(sat=="al ") then
         indexsat=6
      else if(sat=="alg") then
         indexsat=6
      else if(sat=="s3b") then
         indexsat=7
      else if(sat=="s6a") then
         indexsat=8
      else if(sat=="h2a") then
         indexsat=9
      else
         indexsat=10
      endif

      open(33,file=trim(outdir)//'/SLA_PPREJC.dat.'//date, &
           form='formatted',position='append')
      open(34,file=trim(outdir)//'/SLA_LISTRJ.dat.'//date, &
           form='formatted',position='append')

      print*,trim(indir)//trim(infile),len(trim(indir)//trim(infile))

         ist = nf_open(trim(indir)//"/"//trim(infile), nf_nowrite, ncit)
         call handle_err(ist)
         print*,"TAPAS"
         ist = nf_inq_varid (ncit,'dac', idslab)
         call handle_err(ist)
!         ist = nf_inq_varid (ncit,'time', idbdb)
!         call handle_err(ist)
         ist = nf_inq_varid (ncit,'latitude', idlatb)
         call handle_err(ist)
         ist = nf_inq_varid (ncit,'longitude', idlonb)
         call handle_err(ist)
         ist = nf_inq_varid (ncit,'flag', idflagb)
         call handle_err(ist)
         ist = nf_inq_dimid(ncit,'time',dimidpb)
         call handle_err(ist)
         ist = nf_inq_dimlen(ncit,dimidpb,lendb)
         call handle_err(ist)
         ist = nf_inq_varid (ncit,'sla_filtered',idslaf)
         call handle_err(ist)
         ist = nf_inq_varid (ncit,'ocean_tide',dimidot)
         call handle_err(ist)
         ist = nf_inq_varid (ncit,'internal_tide',dimidit)
         call handle_err(ist)
         ist = nf_inq_varid (ncit,'time', idbd)
         call handle_err(ist)

         !!! DAC data
         allocate ( dac(lendb) ) !! For each cycle the number of point
         startb(1)=1
         countb(1)=lendb
         ist = nf_get_vara_real (ncit,idslab,startb,countb,dac)
         call handle_err(ist)
         !!!OCEAN TIDE
         allocate (octi(lendb) )
         ist = nf_get_vara_real (ncit,dimidot,startb,countb,octi)
         call handle_err(ist)
         !!!INTERNAL TIDE
         allocate (inti(lendb) )
         ist = nf_get_vara_real (ncit,dimidit,startb,countb,inti)
         call handle_err(ist)
         !!! SLA FILTERED data
         allocate ( slaf(lendb) ) !! For each cycle the number of point
         ist = nf_get_vara_real (ncit,idslaf,startb,countb,slaf)
         call handle_err(ist)
         !!! TIME
         allocate ( data(lendb) ) !! For each cycle the number of point
         ist = nf_get_vara_double (ncit,idbd,startb,countb,data)
         call handle_err(ist)
         !!! SLA data flag
         allocate ( flagb(lendb) ) !! For each cycle the number of point
         ist = nf_get_vara_real (ncit,idflagb,startb,countb,flagb)
         call handle_err(ist)
         do jjj=1,lendb
            flagb(jjj)=4
         enddo
         !!! Latitude and Longitude data
         allocate ( latb(lendb) )
         allocate ( lonb(lendb) )
         allocate ( alatb(lendb) )
         allocate ( alonb(lendb) )
         allocate ( dx(lendb) )
         allocate ( dy(lendb) )
         allocate ( rdac(lendb) )
         allocate ( dist(lendb) )
         allocate ( slafb(lendb) )
         allocate ( slapdac(lendb) )
         allocate ( rocti(lendb) )
         allocate ( rinti(lendb) )
         allocate ( slapdacpot(lendb) )
         ist = nf_get_vara_int (ncit,idlatb,startb,countb,latb)
         call handle_err(ist)
         ist = nf_get_vara_int (ncit,idlonb,startb,countb,lonb)
         call handle_err(ist)
         !!! MISSING VALUE
         stat = nf_get_att_real (ncit,idslab,'_FillValue',dacmis)
         stat = nf_get_att_real (ncit,idslaf,'_FillValue',slafmis)
         stat = nf_get_att_real (ncit,dimidot,'_FillValue',otmis)
         stat = nf_get_att_real (ncit,dimidit,'_FillValue',itmis)

         stat = nf_get_att_real (ncit,idbd,'_FillValue',timmis)
         stat = nf_get_att_int (ncid,idlatb,'_FillValue',latmis)
         stat = nf_get_att_int (ncid,idlonb,'_FillValue',lonmis)
         !!! Scale factor for lat,lon,sla
         stat = nf_get_att_real (ncit,idlatb,'scale_factor',sclatb)
         stat = nf_get_att_real (ncit,idlonb,'scale_factor',sclonb)
         stat = nf_get_att_real (ncit,idslab,'scale_factor',scslab)
         stat = nf_get_att_real (ncit,idslaf,'scale_factor',scslafb)
         stat = nf_get_att_real (ncit,dimidot,'scale_factor',scsocb)
         stat = nf_get_att_real (ncit,dimidit,'scale_factor',scsitb) 
         slapdac = slafmis
         slapdacpot = slafmis
         print*,"File read"
!         scd=1
         do ld=1,lendb
             if (.not.((lonb(ld).eq.lonmis).or.(latb(ld).eq.latmis))) then
!                 if (lonb(ld)*sclonb.gt.50) then
!                    lonb(ld)=lonb(ld)-(360*1000000)
!                 endif
                 alonb(ld)=lonb(ld)*sclonb
                 alatb(ld)=latb(ld)*sclatb
             endif
             if ((alonb(ld).ge.-6).and.(alonb(ld).le.36.25).and. &
                  (alatb(ld).ge.30.1875).and.(alatb(ld).le.45.9375)) then
                 if (.not.((alonb(ld).ge.-6).and.(alonb(ld).le.0.).and. &
                     (alatb(ld).ge.43).and.(alatb(ld).le.45.9375))) then
                     if (.not.((alonb(ld).ge.26.5).and.(alonb(ld).le.42.).and. &
                         (alatb(ld).ge.41.).and.(alatb(ld).le.50.))) then
!                     print*,scd
                         if ((dac(ld).eq.dacmis).or.(octi(ld).eq.otmis).or.(lonb(ld).eq.lonmis)&
                             .or.(inti(ld).eq.itmis).or.(latb(ld).eq.latmis).or.(slaf(ld).eq.slafmis)&
                             .or.(data(ld).eq.timmis)) then
!                             if ( scd.eq.1 ) then
                                 noqc=noqc+1
                                 write(34,11)date,',',lonb(ld),',',latb(ld),',', &
                                 slaf(ld),',',indexsat
                                 flagb(ld)=4
!                                 scd=0
!                             else
!                                 noqc=noqc+1
!                                 write(34,11)date,',',lonb(ld),',',latb(ld),',', &
!                                 slaf(ld),',',indexsat
!                                 flagb(ld)=4
!                                 scd=1
!                             endif
                         else
!                         if (lonb(ld)*sclonb.gt.50) then
!                             lonb(ld)=lonb(ld)-(360*1000000)
!                         endif
!                         alonb(ld)=lonb(ld)*sclonb
!                         alatb(ld)=latb(ld)*sclatb
                             rdac(ld)=dac(ld)*scslab
                             slafb(ld)=slaf(ld)*scslafb
                             rocti(ld)=octi(ld)*scsocb
                             rinti(ld)=inti(ld)*scsitb
                             nobs=nobs+1
!                             if ( scd.eq.1 ) then
                                 slapdac(ld)=(slafb(ld)+rdac(ld))/scslafb
                                 slapdacpot(ld)=(slafb(ld)+rdac(ld)+rinti(ld))/scslafb
!ADANI                                 slapdacpot(ld)=(slafb(ld)+rdac(ld)+rocti(ld)+rinti(ld))/scslafb
!ADANI Do not add ocean tides
!ADANI
                                 flagb(ld)=1
!                                 scd=0
!                             else
!                                  scd=1
!                             endif
                         endif
                     endif
                 endif
             endif
         enddo
      write(33,20)noqc,',',nobs
      print*,trim(outdir)//"/"//trim(infile)
      ist = nf_open(trim(outdir)//"/"//trim(infile),nf_write, ncio)
      call handle_err(ist)
      ist = nf_redef(ncio)
      call handle_err(ist)
      ist = nf_inq_varid (ncio,'time', tmeid)
      call handle_err(ist)

      dimen(1)=tmeid

      ist=nf_def_var(ncio, 'SLApDAC', nf_short, 1,dimen, spdid)
      call handle_err(ist)
      ist=nf_def_var(ncio, 'SLApDACpOT', nf_short, 1,dimen, spdid2)
      call handle_err(ist)
      ist=nf_copy_att(ncit, idslaf, 'scale_factor', ncio, spdid)
      call handle_err(ist)
      ist=nf_copy_att(ncit, idslaf, '_FillValue', ncio, spdid)
      call handle_err(ist)
      ist=nf_copy_att(ncit, idslaf, 'scale_factor', ncio, spdid2)
      call handle_err(ist)
      ist=nf_copy_att(ncit, idslaf, '_FillValue', ncio, spdid2)
      call handle_err(ist)
      ist = nf_enddef(ncio)
      call handle_err(ist)
      ist=nf_put_var_real(ncio, spdid, slapdac)
      call handle_err(ist)
      ist=nf_put_var_real(ncio, spdid2, slapdacpot)
      call handle_err(ist)
      ist=nf_put_var_real(ncio, idflagb, flagb)
      call handle_err(ist)
      print*,"Variable written"
      ist = nf_close(ncio)
      ist = nf_close(ncit)

      close(33)
      close(34)
!      deallocate ( sla, data, lon, lat )

1 format(i8)
11 format(a8,a1,i10,a1,i10,a1,f10.4,a1,i2)
20 format(i6,a1,i6)
      stop
      end
!----------------------------------------------------------------------
      SUBROUTINE HANDLE_ERR(stat)
      include './netcdf.inc'
      INTEGER stat
      IF (stat .NE. NF_NOERR) THEN
      PRINT *, NF_STRERROR(stat)
      STOP 'Stopped'
      ENDIF
      END SUBROUTINE HANDLE_ERR
!----------------------------------------------------------------------
        SUBROUTINE ref_inter(refld,navlon,navlat,longitude,latitude,tosub)
        real*4 :: latitude, longitude
        integer, parameter::imt=345, jmt=129, nc=2
        integer :: idnla, idnlo, iddt, idtempv,idsalv
        real*4 ::refldso, refldno,refldse, refldne
        real*4, dimension (imt,jmt):: navlon, navlat, lonm,latm,diff
        real*4, dimension (nc,nc):: refldvt, c
        real*4, dimension (imt,jmt):: refld
        real*4 :: dx, dy, xx, yy, lonv, lonl, latv, latl, p, q
        integer :: ii, jj, nx, ny, nni, i, j, z
        integer ist, nlonv, nlonl, nlatv, nlatl, cc, t
        lonm(:,:)=navlon(:,:)-longitude
        latm(:,:)=navlat(:,:)-latitude
        diff(:,:)=abs(lonm(:,:)) + abs(latm(:,:))
        do ii=2,imt-1
           do jj=2,jmt-1
           if (diff(ii,jj).lt.diff(ii,jj+1).and.diff(ii,jj).lt.diff(ii,jj-1)) then
           if (diff(ii,jj).lt.diff(ii+1,jj).and.diff(ii,jj).lt.diff(ii-1,jj)) then
               nx=ii
               ny=jj
           endif
           endif
           enddo
        enddo
        xx=345-(37-longitude)/0.125
        dx=xx-nx
        if (dx.lt.0) then
           nlonv=nx-1
           nlonl=nx
        else if (dx.eq.0) then
           nlonv=nx
           nlonl=nx+1
        else
           nlonv=nx
           nlonl=nx+1
        endif
        yy=129-(46-latitude)/0.125
        dy=yy-ny
        if (dy.lt.0) then
           nlatv=ny-1
           nlatl=ny
        else if (dy.eq.0) then
           nlatv=ny
           nlatl=ny+1
        else
           nlatv=ny
           nlatl=ny+1
        endif
        lonv=navlon(nlonv,nlatv)
        lonl=navlon(nlonl,nlatl)
        latv=navlat(nlonv,nlatv)
        latl=navlat(nlonl,nlatl)
        refldvt=refld(nlonv:nlonl,nlatv:nlatl)
        refldso=refld(nlonv,nlatv)
        refldno=refld(nlonv,nlatl)
        refldse=refld(nlonl,nlatv)
        refldne=refld(nlonl,nlatl)
        p=(longitude-lonv)/(lonl-lonv)
        q=(latitude-latv)/(latl-latv)
        c=refldvt
        do i=1,2
           do j=1,2
              if (c(i,j).eq.0.) then
                 c(i,j)=0
              else
                 c(i,j)=1
              endif
           enddo
        enddo
        if ((c(1,1)+c(1,2)+c(2,1)+c(2,2)).gt.0.) then
           if (c(1,1).eq.0.) then
              cc=(1-c(2,1))*(1-c(1,2))
              refldso=(refldse*c(2,1)+refldno*c(1,2)+refldne*cc)/(c(2,1)+c(1,2)+cc)
           endif
           if (c(1,2).eq.0.) then
              cc=(1-c(1,1))*(1-c(2,2))
              refldno=(refldso*c(1,1)+refldne*c(2,2)+refldse*cc)/(c(1,1)+c(2,2)+cc)
           endif
           if (c(2,1).eq.0.) then
              cc=(1-c(1,1))*(1-c(2,2))
              refldse=(refldso*c(1,1)+refldne*c(2,2)+refldno*cc)/(c(1,1)+c(2,2)+cc)
           endif
           if (c(2,2).eq.0.) then
              cc=(1-c(2,1))*(1-c(1,2))
              refldne=(refldse*c(2,1)+refldno*c(1,2)+refldso*cc)/(c(2,1)+c(1,2)+cc)
           endif
           tosub=(1-q)*((1-p)*refldso+p*refldse)+q*((1-p)*refldno+p*refldne)
        else
           tosub=-99999
        endif
        return
        end
!--------------------------------------------------------------------
       subroutine conv_date_jul(iiday,iimon,iiyear,iijul)

       dimension idmn(12)
       data idmn/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

         iijul = 0

       if(iiyear.lt.1950) stop 'wrong input year'

         do k=1950,iiyear-1
           iijul = iijul + 365
          if(mod(k,4).eq.0)  iijul = iijul + 1
         enddo

       if(iimon.gt.1)then
         do k=1,iimon-1
          iijul = iijul + idmn(k)
          if(k.eq.2 .and. mod(iiyear,4).eq.0)  iijul = iijul + 1
         enddo
       endif

          iijul = iijul + iiday -1

       return
       end
!----------------------------------------------------------------------
      subroutine conv_jul_date(iiday,iimon,iiyear,hour,minutes,iijul)

      dimension idmn(12)
      data idmn/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      double precision  iijul,iijuld
      real  restday, resthour
      integer hour, minutes

      iiyear = 1950
      if(iijul.gt.365)then
         do while(iijul.gt.365.or.iijul.eq.365)
            iijul=iijul-365
            if(mod(iiyear,4).eq.0) iijul=iijul-1
            iiyear=iiyear+1
        enddo
     endif

     if(iijul.lt.0) then
        iiyear=iiyear-1
        iijul=iijul+366
     endif
        iijuld = iijul
        iijul=ceiling(iijul)
        iimon = 1
     if((iijul).gt.idmn(iimon))then
        mond = idmn(iimon)
        do while((iijul).gt.mond)
           iijul = iijul - mond
           iijuld = iijuld - mond
           iimon=iimon+1
           mond = idmn(iimon)
           if(iimon.eq.2 .and. mod(iiyear,4).eq.0) mond = 29
        enddo
     endif

     iiday = int(iijul)
     restday = iijuld-int(iijuld)
     hour = int(restday * 24)
     resthour = (restday * 24) - hour
     minutes = int(resthour * 60)

!!      return

      end
!----------------------------------------------------------------------

