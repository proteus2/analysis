c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        READSOUND       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine ReadSound (pres_S, zz_S, theta_S, qv_S, uu_S, vv_S, 
     >  snd_pres, snd_zz, pmbCB, tcCB, dz, stnelev,psfc,
     >  sndfile, nmax, n, ifCB, plot_wind, gempak_format, snd_wind,
     >  stid,stnm,slat,slon, iyr,imon,iday,ihr,imin)
c
!	1994/11/14  Added code to keep input qv >= 0.
!	1995/02/15  New GEMPAK 5.0 has extra line. Both fmts supported here
!	1997/03/06  Change GEMPAK format istnm (INTEGER) to stnm (CHAR).
!	1997/04/22  Heights are now MSL, not AGL.
!	1999/10/11  Added slat,slon. [RLC]
c
c  Input: 
c     sndfile
c
c  Output:
c     pres_S, zz_S, theta_S, qv_S = Pres, height, theta, mixing ratio
c     snd_pres	= True if sndfile has pressure
c     snd_zz	= True if sndfile has heights
c     ifCB	= 1 if sndfile contains cloud base information
c     pmbCB, tcCB = Pres (mb) and T (C) of cloud base (if ifCB = 1)
c     dz	= Spacing between levels, if known and const
c     stnelev	= Station elevation
c     psfc	= Pres of surface
c     nmax	= Max size of arrays
c     n		= number of points in sndfile
c     plot_wind	= True if sndfile has winds and winds are non-zero.
c     snd_wind	= True if sndfile has winds
c     gempak_format = True if sndfile is in Gempak format
c     stid,stnm,iyr...imin = Read in for Gempak format
c     slat,slno = stn lat and lon (GEMPAK only)
c
c
c
c  Sounding File Formats
c
c  A. Non-Spiffy formats
c     1. 'SAM' format (theta, qv)
c     2a. 'Raw' format (p, T, Td)
c     2b. 'Gempak' format (p, T, Td)
c     2c. 'FSL' format (p, T, Td)
c
c  B. Spiffy Formats
c     01. p, T, Td
c     02. z, theta, qv
c
c
      Implicit None
      Include	'thermo.consts'
      Integer	nmax, ifCB, n, k, iyr,imon,iday,ihr,imin,
     >  nwords,nwordmax, j,j1,j2,jformat, nwrk, k1,k2,kinc, iter
      Parameter (nwordmax=4, nwrk=1000)
      Real	pmbCB, tcCB, dz, stnelev,psfc, xx, temp, rh,
     > tv,tvkm1,tvavg,tdew, slat,slon
      Real	pres_S(nmax), zz_S(nmax), theta_S(nmax), qv_S(nmax),
     >  uu_S(nmax), vv_S(nmax), tempc_S(nwrk), tdewc_S(nwrk)
      Real      wspd,wdir,spcvt,dtr
      Real kts2ms
      Parameter (kts2ms=0.514444)
      Logical	snd_pres, snd_zz, plot_wind, gempak_format, snd_wind,
     >  minus_z
      Logical   fsl_format
      Integer idata(6)
      Integer iunit,ios
      Integer lintyp,iline,sonde,wmo,wban,ielev,ipsfc,itime
      Integer hydro,mxwd,tropl,lines,nlines,tindex,source
      Integer ktop,kbot,kk,kn,nlevs
      Real    dtdz,dtddz,dudz,dvdz
      Logical foundt,foundw

      Character*4 monlist(12)
      Data monlist /' JAN',' FEB',' MAR',' APR',' MAY',' JUN',
     :              ' JUL',' AUG',' SEP',' OCT',' NOV',' DEC'/

      Character*(*) sndfile,stid,stnm
      Character*49  line2
      Character*64  line,word(nwordmax)
      Character*24  cheight,ctemp,cmoist
      Character*4   chmon
      Character*3   staid
      Character*2   wsunits
      Character*1   ns,ew
      Include	'thermo.stfunc'
c
c----------------------------------------------------------------------=
c
      dtr=acos(-1.)/180.
      ifCB	= 0
      jformat	= 83       ! original : -1
      snd_pres	= .False.
      snd_zz	= .False.
      gempak_format = .False.
      fsl_format = .False.
      stnelev	= 0.
c
      Call ChkExist (sndfile,'*** ')
      Print *
      Print *, 'Reading from: ', TRIM(sndfile)
      iunit=31
      Open (iunit, File=sndfile, Status='Old')
c
c
c...Read the first line, Check for my format
c
c
      Print *, ' Reading first line'
      Read (iunit,'(a64)') line
c
      If (line(1:1).EQ.' ') Then
      jformat     = 83
      END IF

      If (line(1:1).EQ.'%') Then
	Print '(a)', line
	j1	= 1
	j2	= Index(line(2:),'%') + j1
	Read (line(2:),'(I2)') jformat
	Print '(a,i2.2)', '% Sounding format: ', jformat
	nwords	= 0
	If (j2.GT.4) Call Parse (line(4:j2-1),word,nwords,nwordmax)
	Print *, nwords, line(4:j2-1)
	Do j=1,nwords
	  Call Uprc (word(j))
	  Print *, 'Sounding flag: ', word(j)
	  If (word(j).EQ.'W') Then
	    snd_wind	= .True.
	  Else If (word(j).EQ.'-Z') Then
	    minus_z	= .True.
	  Else If (word(j).EQ.'CB') Then
	    ifCB	= 1
	  Else 
	    Print *, 'Unrecognized flag: ', word(j)
	    Stop
	  End If
	End Do
c
      Else If (line(1:1).EQ.'&') Then
	Print *, 'ARPS input format'
	jformat	= -99
	snd_wind	= .True.
	minus_z		= .True.
	Read (iunit,'(a)') line
	print *, line
	Read (iunit,'(a)') line
	print *, line
	Read (iunit,'(a)') line
	print *, line
	Read (iunit,'(a)') line
	print *, line
	Read (iunit,*) cheight,ctemp,cmoist
	Call UPRC (cheight)
	Call UPRC (ctemp)
	Call UPRC (cmoist)
c
	If (cheight(1:1).EQ.'H') Then
	  snd_zz = .True.
	Else If (cheight(1:1).EQ.'P') Then
	  snd_pres = .True.
	Else
	  Print *, 'Unknown vertical coord: ', cheight
	  Stop
	End If
	Print *, 'Vertical coord: ', cheight
c
c  1=Pot.temp., 2=Temp
c
	If (ctemp(1:1).NE.'P' .AND. ctemp(1:1).NE.'T') Then
	  Print *, 'Unknown temperature variable: ', ctemp
	  Stop
	End If
	Print *, 'Temperature variable: ', ctemp
c
c  1=Mixing ratio, 2=RH, 3=Td.
c
	If (cmoist(1:1).NE.'M' .AND. cmoist(1:1).NE.'R'
     >   		       .AND. cmoist(1:1).NE.'D') Then
	  Print *, 'Unknown moisture variable: ', cmoist
	  Stop
	End If
	Print *, 'Moisture variable: ', cmoist
c
	Read (iunit,*) stnelev, psfc
	Read (iunit,*) n
	Read (iunit,'(a)') line
      End If
c
c
c----------------------------------------------------------------------=
c
c...So far, the first two lines have been read
c
c  Older or Odd-ball formats
c
c----------------------------------------------------------------------=
c     ARPS Input Format
      If (jformat.EQ.-99) Then
c----------------------------------------------------------------------=
c
	Print *, 'Reading ARPS format sounding'
	k1	= 1
	k2	= n
	kinc	= 1
	If (minus_z) Then
	  k1	= n
	  k2	= 1
	  kinc	= -1
	End If
	Do k=k1,k2,kinc
	  Read(iunit,*) zz_S(k),theta_S(k),qv_S(k),uu_S(k),vv_S(k)
	  qv_S(k) = MAX( 0.0, qv_S(k) )
	  print *, k,zz_s(k),theta_S(k),qv_S(k)
	End Do
c
c...Compute 'first guess' pressure if given Z (ignoring moisture)
c
	If (snd_zz) Then
	  Print *, 'Computing first guess pressure'
	  k	= 1
          pres_S(k)	= psfc
	  temp	= theta_S(k)
	  If (ctemp(1:1).EQ.'P') temp=theta_S(k)*(pres_S(k)*p00inv)**rcp
	  tv		= temp
c
	  Do k=2,n
	    tvkm1	= tv
            dz		= zz_S(k) - zz_S(k-1)
            pres_S(k) = pres_S(k-1) * Exp (-grav*dz*rdi/tvkm1)
c
c  estimate Tv at k+1/2 to get more accurate integration of pres.
c
	    Do iter=1,2
	      temp	= theta_S(k)
	      If (ctemp(1:1).EQ.'P')
     >             temp=theta_S(k)*(pres_S(k)*p00inv)**rcp
              tv	= temp
              tvavg	= .5 * (tv + tvkm1)
              pres_S(k) = pres_S(k-1) * Exp (-grav*dz*rdi/tvavg)
	    End Do
	  End Do
	Else
          pres_S(:)	= zz_S(:)
	End If
c
c...Compute theta if given temp
c
	If (ctemp(1:1).EQ.'T') Then
	  Print *, 'Computing theta from ', ctemp
	  Do k=1,n
	    temp	= theta_S(k)
	    theta_S(k)	= temp * (p00/pres_S(k)) ** rcp
	  End Do
	End If
c
c...Compute Qv if given RH or Td
c
	If (cmoist(1:1).EQ.'R') Then
	  Print *, 'Computing mixing ratio from ', cmoist
	  Do k=1,n
	    rh		= qv_S(k)
	    temp	= theta_S(k) * (pres_S(k)/p00) ** rcp
	    qv_S(k)	= rh * Fmixrat(pres_S(k),Fsvpres(temp-tfrz))
	    qv_S(k) = MAX( 0.0, qv_S(k) )
	  End Do
	Else If (cmoist(1:1).EQ.'D') Then
	  Print *, 'Computing mixing ratio from ', cmoist
	  Do k=1,n
	    tdew	= qv_S(k)
	    temp	= theta_S(k) * (pres_S(k)/p00) ** rcp
	    qv_S(k)	= Fsvpres(tdew-tfrz)
	    qv_S(k) = MAX( 0.0, qv_S(k) )
	  End Do
	End If
c
c       do k=1,n
cprint '(4f10.2)', pres_s(k)/100., zz_s(k), theta_S(k), qv_S(k)*1000.
cend do
c
	If (.NOT.snd_pres) Then
	Call HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S)
	End If
c
c----------------------------------------------------------------------=
      Else If (jformat.LE.0) Then
c----------------------------------------------------------------------=
c
c...Read the second line, check for GEMPAK format
c
      Print *, 'Non-spiffy sounding format'
c
      If (sndfile.EQ.'skew.sound') Then
	Print *, 'Special format for file ', sndfile
	Read(iunit,*)
	Do j=1,142
	Read(iunit,*,End=9800) k, zz_S(k), pres_S(k), theta_S(k),
     >  xx,  qv_S(k), uu_S(k), vv_S(k)
	qv_S(k) = MAX( 0.0, qv_S(k) )
	End Do
 9800   Continue
        n         = j - 1
c
	stnelev   = zz_S(1)
        psfc      = pres_S(1)
        plot_wind = .True.
	snd_pres= .True.
	snd_zz	= .True.
	!Do k=1,n
	!  zz_S(k)	= zz_S(k) - stnelev
	!End Do

	Return
c
      End If

      print *,' Checking contents of second line'
      Read(iunit,'(a)') line2
      print *, ' line2: ',line2
c
      If (line2(:9).EQ.' SNPARM =') Then
        Print *, '% Sounding file is GEMPAK format.'
        ifCB	= -99
        pmbCB	= 0.
        tcCB	= 0.
      ELSE
        print *,' Checking for FSL format'
        read(line2,'(16x,i5,f7.2,1x,f6.2)',iostat=ios) wmo,slat,slon
        IF(ios .eq. 0) THEN
          print *, ' wmo,slat,slon:',wmo,slat,slon
          IF( slat .gt. -90.1 .and. slat .lt. 90.1 .and.
     >        slon .gt. -180.1 .and. slon .lt. 180.1 ) THEN
            Print *, '% Sounding file is FSL format.'
            ifCB = -199 
            pmbCB = 0.
            tcCB = 0.
          END IF
        END IF 
      End If
c
c
c   Read the third line if not GEMPAK or FSL
c
      If (ifCB.GT.-99) Read(iunit,*) ifCB, pmbCB, tcCB
c
c  Model format
c
      If (ifCB.GT.2) Then
	snd_zz	= .True.
        ifCB	= 0
        pmbCB	= 0.
        tcCB	= 0.
        Backspace (1)
        Read(iunit,*) n,dz,stnelev,psfc, ifCB, pmbCB, tcCB
        n		= n * 2
        dz	= dz * .5
        Print '(1x,a,i3,3f9.2)', 
     >    'Model format sounding: n,dz,stnelev,psfc:', n,dz,stnelev,psfc
        Read(iunit,*) (theta_S(k),k=1,n)
        Read(iunit,*) (qv_S(k),k=1,n)
        Read(iunit,*) (uu_S(k),k=1,n)
        Read(iunit,*) (vv_S(k),k=1,n)
        Close (1)
        Do k=1,n
	  zz_S(k) = (k-1) * dz
	  qv_S(k) = MAX( 0.0, qv_S(k) )
          If (uu_S(k).NE.0. .OR. vv_S(k).NE.0.) plot_wind = .True.
        End Do
c
	Call HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S)
c
c
c  Raw format
c
      Else
        snd_pres= .True.
        snd_zz	= .False.
	plot_wind = .False.
	IF (snd_wind) plot_wind = .TRUE.
	gempak_format = .False.
c
c...Gempak format
!	New GEMPAK 5.0 has extra line. Both fmts supported here (1995/02/15)
c
        If (ifCB.EQ.-99) Then
          print *, ' Reading GEMPAK sounding'
	  snd_zz	= .True.
	  plot_wind	= .True.
	  gempak_format = .True.
          Read(iunit,*)
cSTID = OUN        STNM =    72357   TIME = 921119/ 0 0
c23456789012345678901234567890123456789012345678901234567890123456789012
          Read(iunit,9987) stid,stnm,iyr,imon,iday,ihr,imin
 9987     Format (8x,a3,17x,A6,10x,3i2,1x,2i2)
	    Read(iunit,1972) slat,slon 
 1972     Format (11x,f5.2,11x,f7.2)
          Read(iunit,*)
          Read(iunit,'(A)') line
	  IF (line(7:10).NE.'PRES') READ (1,*) line
	  iyr = iyr + 1900
	  IF (iyr .LT. 1950) iyr = iyr + 100
        Else if (ifCB .EQ. -199) THEN
          print *, ' Reading FSL-formatted sounding'
	  snd_zz	= .True.
	  plot_wind	= .True.
          fsl_format = .True.
          Read(line,'(3i7,6x,a4,i7)') lintyp,ihr,iday,chmon,iyr
          Read(line2,'(3i7,f7.2,a1,f6.2,a1,i6,i7)')
     >     lintyp,wban,wmo,slat,ns,slon,ew,ielev,itime
          Read(iunit,'(7i7)') 
     >         lintyp,hydro,mxwd,tropl,lines,tindex,source
          Read(iunit,'(i7,11x,a3,14x,i7,5x,a2)') 
     >         lintyp,staid,sonde,wsunits
          DO imon=1,11
             IF(chmon .eq. monlist(imon)) GO TO 91
          End Do
  91      CONTINUE
          print *, ' chmon, imon =',chmon,imon

          IF(wsunits.eq.'kt') THEN
            spcvt=kts2ms
          ELSE
            spcvt=1.0
          END IF
          IF(ns .eq. 'S') slat=-slat
          IF(ew .eq. 'W') slon=-slon
        End If
c
c
c   Read T(C),Td(C), convert to theta, qv.
c   For Gempak, also convert DDSS to U,V
c
        If (gempak_format) Then
          Do k=1,nmax
            Read(iunit,*,end=1000) pres_S(k),zz_S(k),
     >        tempc_S(k),tdewc_S(k),uu_S(k),vv_S(k)
	    If (uu_S(k).LT.0. .OR. vv_S(k).LE.0) Then
	      uu_S(k)	= 0.
	      vv_S(k)	= 0.
	    End If
	    pres_S(k)	= pres_S(k) * 100.
	    If (tdewc_S(k).LT.-199.) tdewc_S(k) = -199.
          End Do
 1000     CONTINUE
        Else If (fsl_format) Then
          k=0
          ipsfc=999999
          nlines=lines-4
          DO iline=1,nlines
            Read(iunit,'(7i7)',end=1001) lintyp,(idata(j),j=1,6)
            print *, lintyp,(idata(j),j=1,6)
            IF(idata(2).ge.ielev .AND. idata(2).lt.99990 .AND.
     :        idata(1).lt.ipsfc ) THEN
              k=k+1
              IF(idata(1).lt.99990) THEN
                pres_S(k)=10.0*idata(1)
              ELSE
                pres_S(k)=-199.
              END IF
              IF(k .eq. 1) ipsfc=idata(1)
              zz_S(k)=float(idata(2))
              IF(idata(3).lt.99990) THEN
                tempc_S(k)=idata(3)*0.1
              ELSE
                tempc_S(k)=-199.
              END IF
              IF(idata(4).lt.9990) THEN
                tdewc_S(k)=idata(4)*0.1
              ELSE
                tdewc_S(k)=-199.
              END IF
              IF(idata(5).lt.99990 .and. idata(6).lt.200) THEN
                wdir=dtr*float(idata(5))
                wspd=spcvt*idata(6)
                uu_S(k)=-wspd*sin(wdir)
                vv_S(k)=-wspd*cos(wdir)
              ELSE
                uu_S(k)=-199.
                vv_S(k)=-199.
              END IF
            END IF
          End Do
 1001     CONTINUE
          nlevs=k
!
! Interpolate missing temperatures
!
          DO k=2,nlevs
            IF(tempc_S(k) .lt. -190.) THEN
              kbot=k-1
              foundt=.false.
              DO kn=k+1,nlevs
                IF(tempc_S(kn) .gt. -190.) THEN
                  foundt=.true.
                  EXIT
                END IF
              END DO
              IF(foundt) THEN
                ktop=kn
                dtdz=(tempc_S(ktop)-tempc_S(kbot))/     
     >               (zz_S(ktop)-zz_S(kbot))
                DO kk=k,ktop-1
                  tempc_S(kk)=tempc_S(kbot)+dtdz*(zz_S(kk)-zz_S(kbot))
                END DO
              END IF
            END IF
          END DO
!
! Interpolate missing dew points
!
          DO k=2,nlevs
            IF(tdewc_S(k) .lt. -190.) THEN
              kbot=k-1
              foundt=.false.
              DO kn=k+1,nlevs
                IF(tdewc_S(kn) .gt. -190.) THEN
                  foundt=.true.
                  EXIT
                END IF
              END DO
              IF(foundt) THEN
                ktop=kn
                dtddz=(tdewc_S(ktop)-tdewc_S(kbot))/
     >                (zz_S(ktop)-zz_S(kbot))
                DO kk=k,ktop-1
                  tdewc_S(kk)=tdewc_S(kbot)+dtddz*(zz_S(kk)-zz_S(kbot))
                END DO
              END IF
            END IF
          END DO
!
! Interpolate missing winds
!
          DO k=2,nlevs
            IF(uu_S(k) .lt. -190.) THEN
              kbot=k-1
              foundw=.false.
              DO kn=k+1,nlevs
                IF(uu_S(kn) .gt. -190.) THEN
                  foundw=.true.
                  EXIT
                END IF
              END DO
              IF(foundw) THEN
                ktop=kn
                dudz=(uu_S(ktop)-uu_S(kbot))/
     >               (zz_S(ktop)-zz_S(kbot))
                dvdz=(vv_S(ktop)-vv_S(kbot))/
     >               (zz_S(ktop)-zz_S(kbot))
                DO kk=k,ktop-1
                  uu_S(kk)=uu_S(kbot)+dudz*(zz_S(kk)-zz_S(kbot))
                  vv_S(kk)=vv_S(kbot)+dvdz*(zz_S(kk)-zz_S(kbot))
                END DO
              END IF
            END IF
          END DO
      
        Else If (snd_wind) Then
          Do k=1,nmax
            Read(iunit,*,End=1002) pres_S(k),tempc_S(k),tdewc_S(k),
     >        uu_S(k),vv_S(k)
	  pres_S(k)	= pres_S(k) * 100.
	  If (tdewc_S(k).LT.-199.) tdewc_S(k) = -199.
          End Do
 1002     CONTINUE
        ELSE
          Do k=1,nmax
            Read(iunit,*,End=1003) pres_S(k), tempc_S(k), tdewc_S(k)
	    pres_S(k)	= pres_S(k) * 100.
	    If (tdewc_S(k).LT.-199.) tdewc_S(k) = -199.
          End Do
 1003     Continue
        END IF
c
        Close (1)
        n		= k - 1
	If (gempak_format) Then
	  stnelev	= zz_S(1)
	  psfc		= pres_S(1)
	  !Do k=1,n
	  !  zz_S(k) = zz_S(k) - stnelev
	  !End Do
	  Call DDSS2UV (n, uu_S, vv_S)
	End If
c
	Call ConvertSound (n, pres_S, zz_S, theta_S, qv_S, 
     >    tempc_S, tdewc_S, psfc, snd_pres, snd_zz)
c
c	Call HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S) !*** test ***
c
c9977 Format (i3,f8.2,f7.0,3f7.2)
c       Do k=1,n
c  	print 9977, k,pres_S(k)/100., zz_S(k), theta_S(k), qv_S(k)*1000.
c	end do
      End If
c
c
c----------------------------------------------------------------------=
c     Spiffy New Format
      Else If (jformat.NE.-99) Then
c----------------------------------------------------------------------=
c
c...Format 01: P, T, Td
c
cyhstart::::::::::::::::::::::::::::::::::::::::::::
      If (jformat.EQ.83) Then

        READ (iunit,*) k,zz_S(k),pres_S(k),theta_S(k), xx, qv_S(k),
     >      uu_S(k),vv_S(k)
 
        n = k - 1
 
        DO k = n,1,-1
          READ(iunit,*) xx,zz_S(k),pres_S(k),theta_S(k), xx, qv_S(k),
     >      uu_S(k),vv_S(k)
        END DO

        DO k = n,1,-1
          print*, k, zz_s(k),pres_s(k),theta_s(k)
        ENDDO

cyhend:::::::::::::::::::::::::::::::::::::::::::;:

      Else If (jformat.EQ.01) Then
	Read(iunit,*)
	Read(iunit,*)
	Read(iunit,*)
	Read(iunit,*) n, dz, stnelev, psfc
	If (ifCB.GT.0) Read(iunit,*) ifCB, pmbCB, tcCB
	Do k=1,nmax
	  If (snd_wind) Then
	    Read(iunit,*,End=1111) zz_S(k), tempc_S(k), tdewc_S(k),
     >        uu_S(k), vv_S(k)
	  Else 
	    Read(iunit,*,End=1111) zz_S(k), tempc_S(k), tdewc_S(k)
	  End If
	End Do
 1111   Continue
	n	= k - 1
	!Do k=1,n
	!  zz_S(k) = zz_S(k) - stnelev
	!End Do
c
c
c
c...Format 02: Z, Theta, Qv, U, V
c
      Else If (jformat.LE.02) Then
	Read(iunit,*)
	Read(iunit,*)
	Read(iunit,*)
	Read(iunit,*) n, dz, stnelev, psfc
	If (ifCB.GT.0) Read(iunit,*) ifCB, pmbCB, tcCB
	k1	= 1
	k2	= n
	kinc	= 1
	If (minus_z) Then
	  k1	= n
	  k2	= 1
	  kinc	= -1
	End If
	Do k=k1,k2,kinc
	  If (snd_wind) Then
	    Read(iunit,*) zz_S(k),theta_S(k),qv_S(k),uu_S(k),vv_S(k)
	  Else 
	    Read(iunit,*) zz_S(k),theta_S(k),qv_S(k)
	  End If
	  qv_S(k) = MAX( 0.0, qv_S(k) )
	  print *, k, zz_s(k), theta_S(k), qv_S(k)
	End Do
	!Do k=1,n
	!  zz_S(k)	= zz_S(k) - stnelev
	!End Do

c
	Call HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S)
c
c
c
c
c...Unrecognized Format 
c
      Else
	Print *, 'Unrecognized sounding format: ', jformat
	Stop
      End If
c
c
c----------------------------------------------------------------------=
      End If
c----------------------------------------------------------------------=
c
c
c...Done reading the sounding file
c
       If (snd_wind) plot_wind = .True.
c
c
      Print'(i5,2a)',n,' points read in from ', TRIM(sndfile)
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        DDSS2UV         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine DDSS2UV (n, uu, vv)
      Implicit None
      Integer	k,n
      Real	uu(n), vv(n), dir, spd, deg2rad
      Parameter (deg2rad=.0174532925)
c
c  On input, uu is dir, vv is spd
c  On output, uu, vv are Cartesian speeds
c
      Do k=1,n
	dir	= uu(k)
        spd	= vv(k)
        uu(k)	= - spd * Sin(deg2rad*dir)
        vv(k)	= - spd * Cos(deg2rad*dir)
      End Do
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########      CONVERTSOUND      ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine ConvertSound (n, pres, zz, theta, qv, tempc, tdewc, 
     >    psfc, snd_pres, snd_zz)
      Implicit None
      Include 'thermo.consts'
      Integer	k, n
      Real	pres(n), zz(n), theta(n), qv(n), tempc(n), tdewc(n)
      Real	psfc, dz, tvavg, tvkp1, tvkm1, tv, temp
      Logical	snd_pres, snd_zz
      Include 'thermo.stfunc'
c
c
c  Convert sounding containing (p|z,T,Td) to a standard set of variables
c  (p,z,theta,qv).
c
c
c...Has Z but not P -- compute hydrostatically
c
      If (.NOT.snd_pres) Then
      Print *, '% ConvertSound: computing (p,theta,qv) from (z,T,Td)'
c
      k		= 1
      pres(k)	= psfc
      temp	= tempc(k) + tfrz
      theta(k)	= Ftheta(pres(k),temp)
      qv(k)	= Fmixrat(pres(k),Fsvpres(tdewc(k)))
      tv	= Ftvirtnowl(temp,qv(k))
c
      Do k=1,n-1
        dz	= zz(k+1) - zz(k)
        pres(k+1) = pres(k) * Exp (-grav*dz*rdi/tv)
c
c  estimate Tv at k+1 to get more accurate integration of pres.
c
        temp	= tempc(k+1) + tfrz
        theta(k+1) = Ftheta(pres(k+1),temp)
        qv(k+1)	= Fmixrat(pres(k+1),Fsvpres(tdewc(k+1)))
        tvkp1	= Ftvirtnowl(temp,qv(k+1))
        tvavg	= .5 * (tv + tvkp1)
        pres(k+1) = pres(k) * Exp (-grav*dz*rdi/tvavg)
c
c  recompute Th,Qv
c
        theta(k+1) = Ftheta(pres(k+1),temp)
        qv(k+1)	= Fmixrat(pres(k+1),Fsvpres(tdewc(k+1)))
      End Do
c
c
c...Has pressure and heights
c
      Else If (snd_zz) Then
      Print *, '% ConvertSound: computing (Th,Qv) from (p,z,T,Td)'
c
      Do k=1,n
        theta(k)= Ftheta(pres(k),tempc(k)+tfrz)
        qv(k)	= Fmixrat(pres(k),Fsvpres(tdewc(k)))
      End Do
c
c
c...Has pressure but not heights
c
      Else If (.NOT.snd_zz) Then
      Print *, '% ConvertSound: computing (z,Th,Qv) from (p,T,Td)'
c
      Do k=1,n
        theta(k)= Ftheta(pres(k),tempc(k)+tfrz)
        qv(k)	= Fmixrat(pres(k),Fsvpres(tdewc(k)))
      End Do
c
      zz(1)	= 0.
      tvkm1	= Ftvirtnowl(tempc(1)+tfrz,qv(1))
c
      Do k=2,n
        tv	= Ftvirtnowl(tempc(k)+tfrz,qv(k))
        tvavg	= .5 * (tvkm1 + tv)
	dz	= rd/grav*tvavg * Log(pres(k-1)/pres(k))
	zz(k)	= zz(k-1) + dz
	tvkm1	= tv
      End Do
c
c
      End If
c
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        HYDROPRES       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine HydroPres (n, psfc, pres, zz, theta, qv)
      Implicit None
      Include 'thermo.consts'
      Integer	k,n, iter
      Real	pres(n), zz(n), theta(n), qv(n)
      Real	dz, tv, tvkm1, tvavg, temp, psfc
      Include 'thermo.stfunc'
c
c
c
      Print *, '% HydroPres: compute p from (z,Th,Qv)'
c
      k		= 1
      pres(1)	= psfc
      temp	= theta(k) * (pres(k) * p00inv) ** rcp
      tv	= Ftvirtnowl(temp,qv(k))
c
      Do k=2,n
	tvkm1	= tv
        dz	= zz(k) - zz(k-1)
        pres(k) = pres(k-1) * Exp (-grav*dz*rdi/tvkm1)
c	print '(2i3,f7.2)', 0, k, pres(k)/100.
c
c  estimate Tv at k+1/2 to get more accurate integration of pres.
c
	Do iter=1,2
          temp	= theta(k) * (pres(k) * p00inv) ** rcp
          tv	= Ftvirtnowl(temp,qv(k))
          tvavg	= .5 * (tv + tvkm1)
          pres(k) = pres(k-1) * Exp (-grav*dz*rdi/tvavg)
c	  print '(2i3,f7.2)', iter, k, pres(k)/100.
	End Do
c
      End Do

c     do k=1,n
c     tempc	= theta(k) * (pres(k) * p00inv) ** rcp - tfrz
c     tdewc	= Ftdewc(Fvpres(pres(k),qv(k)))
c     print 9999, k, zz(k), pres(k)/100., tempc, tdewc,
c    >  theta(k), qv(k)*1000.
c     end do
c9999 Format (i4,f7.0,f8.2, 4f7.2)
c
      End
!
!
!
!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########      WRITEARPSSND      ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################
!
!
      SUBROUTINE WriteARPSSnd (file, zrefsfc, psfc, n, 
     &  pres, z, theta, temp, qv, rh, tdew, u, v)
      IMPLICIT NONE

!     RLC 1997/04/22
!     All units MKS

      INTEGER n
      REAL zrefsfc, psfc
      REAL pres(n), z(n), theta(n), temp(n), qv(n), rh(n), tdew(n),
     &  u(n), v(n)
      CHARACTER*(*) file

      INTEGER k
      CHARACTER*16 hgtstr, tempstr, mststr, windstr
      DATA hgtstr/"'height'"/, tempstr/"'temp.'"/,
     &   mststr /"'rel. humidity'"/, windstr/"'uv'"/

!----------------------------------------------------------------------=

      PRINT *, '% WriteARPSSnd: ', file
      OPEN (UNIT=1, FILE=file)

      WRITE (1,9000) 'Line 1'
      WRITE (1,9000) 'Line 2'
      WRITE (1,9000) 'Line 3'
      WRITE (1,9000) 'Line 4'
      WRITE (1,9000) 'Line 5'
      WRITE(1,'(12(A," "))') hgtstr, tempstr, mststr, windstr
      WRITE (1,'(2F12.1)') zrefsfc, psfc
      WRITE (1,*) n
      WRITE (1,*) '-----------------------------------------'

      DO k=n,1,-1
	WRITE (1,9010) z(k), temp(k), rh(k), u(k), v(k)
      END DO

      CLOSE (1)

 9000 FORMAT (16A)
 9010 FORMAT (2F12.2, F10.4, 2F12.3)
      END
