      PROGRAM rout
c     Routing algorithm developed by D. Lohmann.
c     Code first developed by the University of Washington. See WA Hydrology Homepage for more details.
c     Modified by the Resilient Water Systems Group/Singapore University of Technology and Design to account for reservoir operations.
c     Reservoir presentation and operations are incorporated by the following steps
c     1. We first calculate the time lag from each cell to a reservoir / the basin outlet
c     2. We then calculate the flow to each reservoirs and the outlet
c     3. Track the reservoir network
c     4. Reservoir operations are modelled based on rule curves (1)/ operating rules (2)/ predefined time-series data (3)
      IMPLICIT NONE

C************************************************************************************************************************************************************************************
C     Declare variables
c     RCS ID STRING
      CHARACTER*50 RCSID
c     DATA RCSID/"$Id: read_routines.f,v1.0 2019/04/17$"/
      INTEGER   IARGC
      INTEGER   isaleap
      EXTERNAL  isaleap
C     Change dimensions here:
C     1. nrow and ncol should be larger than the grid
C     2. nyr should equal run length yrs+1
C     3. nreach should be equal to the number of stations considered in stations.txt
      INTEGER NROW, NCOL, DAYS, NYR, NREACH, NUMRES
C      INTEGER NMAXREACH
C      PARAMETER (NMAXREACH=NUMRES)
      PARAMETER (NROW = 100, NCOL = 100, NREACH=463, NUMRES=109)
      PARAMETER (NYR = 23)
C     -------------- No changes after here -------------------------------------------------------
C     -------------- Note: NUMRES is the maximum number of reservoirs. Change if more reservoirs added
C     Unit-hydrograph parameters
      INTEGER   KE, LE, TMAX, UH_DAY, PMAX
      REAL      DT
      PARAMETER (DAYS=NYR*366)
      PARAMETER (KE   = 12)
      PARAMETER (LE   = 48)
      PARAMETER (DT   = 3600.0)
      PARAMETER (UH_DAY = 96)
      PARAMETER (TMAX = UH_DAY*24)
      PARAMETER (PMAX = 2500)
      REAL      UH_BOX(PMAX,KE),UHM(NCOL,NROW,LE)
      REAL      UH_S(PMAX,KE+UH_DAY-1,NUMRES)
      REAL      UH_DAILY(PMAX,UH_DAY)
      REAL      FR(TMAX,2)
      CHARACTER*80 UH_STRING(NREACH)

C     Routing
      REAL      AREA(NREACH)
      REAL      FACTOR_SUM
      REAL      XC,YC,SIZE
      REAL      FDUM
      REAL      VELO(NCOL,NROW), DIFF(NCOL,NROW)
      REAL      XMASK(NCOL,NROW), FRACTION(NCOL,NROW)
      REAL      BASE(DAYS,NREACH), RUNO(DAYS,NREACH), FLOW(DAYS,NREACH)
      INTEGER   DIREC(NCOL,NROW,2)
      INTEGER   IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER   MO(12*NYR),YR(12*NYR)
      INTEGER   NO_OF_BOX(NUMRES,NREACH)
      INTEGER   CATCHIJ(PMAX,2,NUMRES,NREACH)
C      INTEGER   REDUCEI(NREACH,15000),REDUCEJ(NREACH,15000)
      INTEGER   H(NCOL,NROW)
      INTEGER   PI(NREACH),PJ(NREACH),NR(NREACH)
      INTEGER   PII,PJJ
      INTEGER   IROW,ICOL
      INTEGER   LP,M,Y,J,I,K
      INTEGER   DPREC,FINAL
      INTEGER   NDAY, CLEN
      INTEGER   NMONTHS
      LOGICAL   TORF
C     Reservoir parameters
      INTEGER   RESER(NCOL,NROW)
      INTEGER   RESORDER(NREACH,NUMRES)
      REAL      VRESER(NUMRES,NREACH),VOL(NUMRES,2,NREACH),QFLUSHED(NUMRES,NREACH),RESFLOWS(NUMRES,DAYS,NREACH),FLOWOUT(NUMRES,NREACH)
      REAL      HHO(NUMRES,2),LUONGDIEN(NUMRES,DAYS),HTK(NUMRES)
      INTEGER   NORESERVOIRS(NREACH)
      INTEGER   NSTATIONS
      INTEGER   N,D
      INTEGER   RES_DIRECT(NUMRES,3,NREACH)
      REAL      RES_EVAPORATION(NUMRES,DAYS)
      REAL      TRVLTIME(NUMRES)
C     Filename
      CHARACTER*21 NAME(NREACH)
      CHARACTER*7  NAMERS51
      CHARACTER*7  NAMERS5(NREACH,NUMRES)
      CHARACTER*20 NAMERS(NREACH,NUMRES)
      CHARACTER*7  NAME5(NREACH)
      CHARACTER*72 FILE_INPUT,FILENAME,RFILENAME
      CHARACTER*72 INPATH,OUTPATH,RPATH
C     Variables for monthly means
      INTEGER   DAYS_IN_MONTH(12)
      DATA      DAYS_IN_MONTH /31,28,31,30,31,30,31,31,30,31,30,31/
      INTEGER   START_YEAR,STOP_YEAR,FIRST_YEAR,LAST_YEAR
      INTEGER   START_MO,STOP_MO,FIRST_MO,LAST_MO
      INTEGER   PORT
      REAL      MONTHLY(12*NYR)
      REAL      MONTHLY_mm(12*NYR)
      REAL      YEARLY(12)
      REAL      YEARLY_mm(12)
      NORESERVOIRS = 0
      NSTATIONS = 0

      

C************************************************************************************************************************************************************************************
C     OPEN NECESSARY FILES
C************************************************************************************************************************************************************************************
c     Process commandline args
      IF(IARGC().NE.1)THEN
           PRINT*, 'USAGE:  rout <infile>'
           STOP
      ENDIF
      CALL GETARG(1,FILE_INPUT)
      OPEN(1,FILE=FILE_INPUT,STATUS='OLD',ERR=9001)
      READ(1,'(//A)') FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)
c     Process velocity file
      READ(1,*)
      READ(1,*) TORF
      IF(TORF)THEN
           READ(1,'(A)') FILENAME
           CALL READ_VELO(VELO,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
           READ(1,*) FDUM
           CALL INIT_ARRAY(VELO,NCOL,NROW,FDUM)
      ENDIF
c     Process diffusion file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
          READ(1,'(A)') FILENAME
          CALL READ_DIFF(DIFF,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
          READ(1,*) FDUM
          CALL INIT_ARRAY(DIFF,NCOL,NROW,FDUM)
      ENDIF
c     Process xmask file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
          READ(1,'(A)') FILENAME
          CALL READ_XMASK(XMASK,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
          READ(1,*) FDUM
          CALL INIT_ARRAY(XMASK,NCOL,NROW,FDUM)
      ENDIF
c     Read fraction file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
          READ(1,'(A)') FILENAME
          CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
          READ(1,*) FDUM
          CALL INIT_ARRAY(FRACTION,NCOL,NROW,FDUM)
      ENDIF
c     Read station file
      READ(1,'(/A)')FILENAME
      OPEN(10,FILE=FILENAME)
c     Read input path and precision of VIC filenames
      READ(1,'(/A)')INPATH
      READ(1,*)DPREC
c     Read output pathname
      READ(1,'(/A)')OUTPATH
c     Read input path of reservoir information
      READ(1,'(/A)') RFILENAME
      READ(1,'(/A)') RPATH
c     Read input file name of reservoir locations
      CALL READ_RESE(RESER,ICOL,IROW,NCOL,NROW,RFILENAME)
c     Number of days to process
      READ(1,*)
c     Start and end year/month from VIC simulation
      READ(1,*) START_YEAR, START_MO, STOP_YEAR, STOP_MO
c     Calculate number of days & months in simulation
      M=START_MO
      Y=START_YEAR
      NMONTHS = 0
      NDAY=0
      DO J=START_MO,12*(STOP_YEAR-START_YEAR)+STOP_MO
        IF(M.EQ.2) THEN
           LP=isaleap(Y)
        ELSE
           LP=0
        ENDIF
        NDAY = NDAY+DAYS_IN_MONTH(M)+LP
        NMONTHS = NMONTHS + 1
        MO(NMONTHS) = M
        YR(NMONTHS) = Y
        M = M + 1
        IF (M .GT. 12) THEN
            M = 1
            Y  = Y + 1
        ENDIF
      END DO
      IF(NDAY.GT.DAYS) THEN
         PRINT*, 'IN ROUT.F RESET DAYS TO ', NDAY
         STOP
      ENDIF
      PRINT*,'NDAY = ',NDAY, ' NMONTHS = ',NMONTHS
c     Read start and end year/month for writing output
      READ(1,*) FIRST_YEAR, FIRST_MO, LAST_YEAR, LAST_MO
c     Read uh file
      READ(1,'(/A)')FILENAME
      READ(1,*)PORT
c      PORT = 31416


C************************************************************************************************************************************************************************************
C     START MODELLING
C************************************************************************************************************************************************************************************

c     Calculate impulse response function for grid cells
      CALL MAKE_UHM(UHM,VELO,DIFF,XMASK,NCOL,NROW,LE,DT,IROW,ICOL)
c     loop over required stations
      I=1
 100  CONTINUE
      READ(10,*,END=110) NR(I),NAME(I),PI(I),PJ(I),AREA(I)
      READ(10,'(A80)',END=110) UH_STRING(I)
      IF (NR(I) .EQ. 1) THEN
c            WRITE(*,'(I2,2X,A,I4,I4,G12.6)')
c     &             NR(I), NAME(I), PI(I), PJ(I)
c            PRINT*, 'Routing station: ', NAME(I)
            PI(I)=ICOL+1-PI(I)						!note, the arrays are flipped left to right
            NAME5(I) = NAME(I)
            NSTATIONS=NSTATIONS+1
c     Look for cells, contributing to the outlet                
            PII=PI(I)
            PJJ=PJ(I)
c            print*, 'Searching catchment...'
            CALL SEARCH_WHOLECATCHMENT
     &           (PI(I),PJ(I),DIREC,NCOL,NROW,
     &            NO_OF_BOX(:,I),CATCHIJ(:,:,:,I),PMAX,IROW,ICOL,
     &            NORESERVOIRS(I), RES_DIRECT(:,:,I), RESER,NUMRES)
c            print*, 'Reading grid_UH...'

c     Read a pre-defined UH grid
            CALL READ_GRID_UH
     &          (UH_BOX,KE,PMAX,NO_OF_BOX(:,I), CATCHIJ(:,:,:,I),FILENAME,NORESERVOIRS(I),NUMRES)
c     Make UH grid for reservoir catchments
c            print*, 'Making grid UH for reservoirs...'
            D=1
            DO N = 1,NO_OF_BOX(NORESERVOIRS(I),I)
               
               IF ((RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I))
     &              .GT.0) .AND. (RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I))
     &              .NE.9999)) THEN					! the second condition can be removed
                 WRITE(NAMERS(I,D),*) RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I))
                 NAMERS51 = trim(adjustl(NAMERS(I,D)))//"_"//trim(adjustl(NAME5(I)))
                 CALL SEARCH_CATCHMENTRS(CATCHIJ(N,1,NORESERVOIRS(I),I),
     &                 CATCHIJ(N,2,NORESERVOIRS(I),I),DIREC,NCOL,NROW,NO_OF_BOX(:,I),CATCHIJ(:,:,:,I),PMAX,
     &                 IROW,ICOL,NORESERVOIRS(I),RES_DIRECT(:,:,I),RESER,0,SIZE,VELO,PI(I),PJ(I),NUMRES)
                 CALL MAKE_GRID_UH
     &                 (DIREC, NO_OF_BOX(:,I), UH_DAY, TMAX, PI(I), PJ(I), LE, UH_DAILY, KE,
     &                 CATCHIJ(:,:,:,I),UHM, FR, PMAX, NCOL, NROW, UH_BOX, UH_S,
     &                 UH_STRING(I),NAMERS51,NORESERVOIRS(I),RESER,
     &                 RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I)),
     &                 RES_DIRECT(:,:,I),RESORDER(I,D),0,NUMRES)
                  NAMERS5(I,RESORDER(I,D)) = trim(adjustl(NAMERS(I,D)))//"_"//trim(adjustl(NAME5(I)))
                  D=D+1
                  
                  
               END IF               
            END DO
c     Make UH grid for the rest of the basin to the basin outlet
c            print*, 'making grid UH...'
            
            CALL SEARCH_CATCHMENTRS(PI(I),PJ(I),
     &            DIREC,NCOL,NROW,NO_OF_BOX(:,I),CATCHIJ(:,:,:,I),PMAX,
     &            IROW,ICOL,NORESERVOIRS(I),RES_DIRECT(:,:,I),RESER,1,SIZE,VELO,PI(I),PJ(I),NUMRES)
            CALL MAKE_GRID_UH
     &           (DIREC, NO_OF_BOX(:,I), UH_DAY, TMAX, PI(I), PJ(I), LE, UH_DAILY, KE,
     &            CATCHIJ(:,:,:,I),UHM, FR, PMAX, NCOL, NROW, UH_BOX, UH_S,
     &            UH_STRING(I),NAME5(I),NORESERVOIRS(I),RESER,NORESERVOIRS(I),RES_DIRECT(:,:,I),RESORDER(I,D),1,NUMRES)

C            DO N = 1,NORESERVOIRS(I)
C                        DO J= 1,NO_OF_BOX(RES_DIRECT(N,1,I),I)
C                              REDUCEI(I,J)=CATCHIJ(J,1,RES_DIRECT(N,1,I))
C                              REDUCEJ(I,J)=CATCHIJ(J,2,RES_DIRECT(N,1,I))
C                        ENDDO
C            ENDDO
            I=I+1
      ENDIF    
      GOTO 100
      
 110  CONTINUE
      
c Internal check to verify that the number of stations is equal to NREACH
      IF (NSTATIONS .EQ. NREACH) THEN
c            print*, 'Number of stations equal number of reaches'
      ELSE
            print*, 'ERROR: Number of stations considered in stations.txt is different from number of reaches in rout.f'
            STOP
      ENDIF

c     Flow generation for the required station step by step
c               print*, 'making convolution...'
               CALL MAKE_CONVOLUTIONRS
     &               (RPATH,RESER,NCOL, NROW, NO_OF_BOX, PMAX, DAYS,
     &               CATCHIJ, BASE, RUNO, FLOW, KE, UH_DAY, FRACTION,
     &               FACTOR_SUM,XC,YC,SIZE,DPREC,INPATH,ICOL,NDAY,
     &               IDAY,IMONTH,IYEAR, START_YEAR, START_MO, MO, YR, NYR, VOL, QFLUSHED,
     &               RESFLOWS, FLOWOUT, HHO, LUONGDIEN,HTK,DIREC,IROW,
     &               PI,PJ,NORESERVOIRS,RES_DIRECT,RES_EVAPORATION,TRVLTIME,RESORDER,NAMERS5,NAME5,NREACH,UH_S,VRESER,PORT, NUMRES)   


c     Writing data into files
               print*, 'writing data...'
c     Writing data (only daily discharges) for all the stations into the file 'output'
               CLEN = INDEX(OUTPATH,' ')-1
               OPEN(80, FILE = OUTPATH(1:CLEN)//'OUTPUT.day')
               OPEN(81, FILE = OUTPATH(1:CLEN)//'OUTPUT.day_mm')
               WRITE(80,*) '       YEAR','        MONTH','        DAY ',NAME(:)
               WRITE(81,*) '       YEAR','        MONTH','        DAY ',NAME(:)
               DO I = 1,NDAY
                  WRITE(80,*) IYEAR(I),IMONTH(I),IDAY(I),FLOW(I,:)
                  WRITE(81,*) IYEAR(I),IMONTH(I),IDAY(I),FLOW(I,:) / FACTOR_SUM
               END DO
               CLOSE(80)
               CLOSE(81)
c     Writing data for each station into a specific file
c               DO I=1,NREACH
c                   CALL WRITE_DATA
c     &                 (FLOW(:,I), NDAY, NAME(I), FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,
c     &                  VOL(:,:,I),FLOWIN, FLOWOUT(:,I), HHO, LUONGDIEN,HTK,RESER,NCOL,NROW,
c     &                  ICOL,IROW, RPATH, NORESERVOIRS(I),RES_DIRECT(:,:,I),RES_EVAPORATION,NO_OF_BOX(:,I),NUMRES)
c                   CALL WRITE_MONTH
c     &                 (DAYS_IN_MONTH,START_YEAR, STOP_YEAR, FIRST_YEAR,
c     &                  LAST_YEAR, START_MO, STOP_MO, FIRST_MO,LAST_MO,
c     &                  NAME(I), DAYS, FLOW(:,I), FACTOR_SUM, MONTHLY, MONTHLY_mm,
c     &                  YEARLY,YEARLY_mm,OUTPATH,NDAY,IMONTH,IYEAR,MO,YR,NMONTHS,NYR)      
c               ENDDO 
               CALL WRITE_RES
     &                 (FLOW(:,NREACH), NDAY, NAME(NREACH), FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,
     &                  VOL(:,:,NREACH),QFLUSHED(:,NREACH), RESFLOWS(:,:,NREACH), FLOWOUT(:,NREACH), HHO, LUONGDIEN,HTK,RESER,NCOL,NROW,
     &                  ICOL,IROW, RPATH, NORESERVOIRS(NREACH),RES_DIRECT(:,:,NREACH),RES_EVAPORATION,NO_OF_BOX(:,NREACH),VRESER(:,NREACH),NUMRES)
                  


c MODIFICA A WRITING
C               DO I=1,NREACH
C                   CALL WRITE_DATA1
C     &                 (FLOW(:,I), NDAY, NAME(I), FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR)
C     
C                   CALL WRITE_MONTH
C     &                 (DAYS_IN_MONTH,START_YEAR, STOP_YEAR, FIRST_YEAR,
C     &                  LAST_YEAR, START_MO, STOP_MO, FIRST_MO,LAST_MO,
C     &                  NAME(I), DAYS, FLOW(:,I), FACTOR_SUM, MONTHLY, MONTHLY_mm,
C     &                  YEARLY,YEARLY_mm,OUTPATH,NDAY,IMONTH,IYEAR,MO,YR,NMONTHS,NYR)      
C               ENDDO 
C               CALL WRITE_RESERVOIR  
C     &                 (FLOW(:,NREACH), NDAY, NAME(NREACH), FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,  
C     &                  VOL(:,:,NREACH),FLOWIN, FLOWOUT, HHO, LUONGDIEN,HTK,RESER,NCOL,NROW,
C     &                  ICOL,IROW, RPATH, NORESERVOIRS(NREACH),RES_DIRECT(:,:,NREACH),RES_EVAPORATION,NO_OF_BOX(:,NREACH),NUMRES)


C      DO I=1,NREACH
C            DO J=1,NORESERVOIRS(I)
C                  DO N = 1,NO_OF_BOX(RESORDER(I,J),I)
C                        WRITE (*,*)(UH_S(N,K,J,I),  K = 1,5)
C                  ENDDO
C                  PRINT*,'//'
C            ENDDO
C            PRINT *,'//'
C      ENDDO
C      PRINT *, 'LEGGE DA UH_S'
C      DO I = 1,NREACH
C            DO J = 1,NORESERVOIRS(I)-1
C                  OPEN(98, file = trim(adjustl(NAMERS5(I,J)))//'.uh_s', status='old')
C                  DO N = 1,NO_OF_BOX(RESORDER(I,J),I)
C                        READ(98, *) (UH_S(N,K,J), K = 1,5)  
C                  ENDDO
C                  CLOSE(98) 
C            ENDDO
C            J=NORESERVOIRS(I)
C            print*,J
C            print*,trim(adjustl(NAME5(I)))
C            OPEN(98, file = trim(adjustl(NAME5(I)))//'.uh_s', status='old')
C            DO N= 1,NO_OF_BOX(RESORDER(I,J),I)
C                  READ(98, *) (UH_S(N,K,J), K = 1,5)
C            ENDDO
C            CLOSE(98)
C            print*, UH_S(1,1:5,2)
C      ENDDO  
      
C      DO I=1,NREACH
C            DO J=1,NORESERVOIRS(I)
C                  DO N = 1,NO_OF_BOX(RESORDER(I,J),I)
C                        WRITE (*,*)(UH_S(N,K,J,I),  K = 1,5)
C                  ENDDO
C                  PRINT*,'//'
C            ENDDO
C            PRINT *,'//'
C      ENDDO
      


      STOP     
 9001 WRITE(*,*) 'CANNOT OPEN: ', FILE_INPUT
      END
C************************************************************************************************************************************************************************************
c     FUNCTION  ISALEAP
      integer function isaleap( iyr )
c     return 1 if a leap yr else 0
      if( (mod(iyr,4) .eq. 0 .and. mod(iyr,100) .ne.0)
     $                       .or. mod(iyr,400) .eq. 0) then
         isaleap = 1
      else
         isaleap = 0
      endif
      end
C     END OF FILE
C************************************************************************************************************************************************************************************