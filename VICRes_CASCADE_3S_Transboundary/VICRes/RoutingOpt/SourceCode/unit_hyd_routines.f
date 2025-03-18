c     SUBROUTINES RELATED TO UNIT HYDROGRAPH
c     Calculate impulse response function for grid cells
c     using equation (15) from Lohmann, et al. (1996)  Tellus article
      SUBROUTINE MAKE_UHM(UH,VELO,DIFF,XMASK,NCOL,NROW,LE,DT,
     $        IROW,ICOL)
      CHARACTER*50 RCSID
      REAL    UH(NCOL,NROW,LE), INTE
      INTEGER I, J, K, LE
      REAL    T, DT, POT, VELO(NCOL,NROW), DIFF(NCOL,NROW)
      REAL    XMASK(NCOL,NROW), H, PI
      PI = 4.0*atan(1.0)
      DO I = 1, ICOL
         DO J =1, IROW
            T = 0.0
            DO K = 1,LE
               T = T + DT
               IF (VELO(I,J) .GT. 0.0) THEN
                  POT = ((VELO(I,J)*T-XMASK(I,J))**2.0)/
     &                  (4.0*DIFF(I,J)*T)
                  IF (POT .GT. 69.0) THEN
                     H = 0.0
                     ELSE
                     H = 1.0/(2.0*SQRT(PI*DIFF(I,J))) *
     &                   XMASK(I,J)/(T**1.5) * EXP(-POT)
                  END IF
               ELSE
                  H = 0.0
               END IF
               UH(I,J,K) = H
            END DO
         END DO
      END DO
      DO I = 1, ICOL
         DO J = 1, IROW
            INTE = 0.0
            DO K = 1,LE
               INTE = INTE + UH(I,J,K)
            END DO
            IF (INTE .GT. 0.0) THEN
               DO K = 1,LE
                  UH(I,J,K) = UH(I,J,K)/INTE
               END DO
            END IF
         END DO
      END DO
      RETURN
      END

c     Solve for
      SUBROUTINE MAKE_GRID_UH
     & (DIREC, NOB, UH_DAY, TMAX, PI, PJ, LE, UH_DAILY, KE,
     &  CATCHIJ, UHM, FR, PMAX, NCOL, NROW, UH_BOX,
     &  UH_S, UH_STRING, NAME5,NORESERVOIRS, RESER,
     &	RESERNAME,RES_DIRECT,RESORDER,FINAL,NUMRES)
      IMPLICIT NONE
      INTEGER UH_DAY, TMAX, PMAX, KE, NUMRES
      INTEGER NOB(NUMRES)
      INTEGER PI, PJ, LE, NCOL, NROW,FINAL
      INTEGER DIREC(NCOL,NROW,2)
      REAL    UH_DAILY(PMAX,UH_DAY)
      INTEGER CATCHIJ(PMAX,2,NUMRES)
      REAL    UHM(NCOL,NROW,LE), FR(TMAX,2)
      REAL    UH_BOX(PMAX,KE)
      REAL    UH_S(PMAX,KE+UH_DAY-1,NUMRES)
      INTEGER N, I, J, K, L, T, II, JJ, TT, U
      INTEGER 	RESER(NCOL,NROW)
      REAL    SUM
      INTEGER stat
      INTEGER 	RES_DIRECT(NUMRES,3)
      CHARACTER*80 UH_STRING       !new, AW
      CHARACTER*7  NAME5
      INTEGER NORESERVOIRS,RESORDER, RESERNAME
      RESORDER = NORESERVOIRS
      DO I = 1, (NORESERVOIRS-1)
        IF (RESERNAME .EQ. RES_DIRECT(I,1)) THEN
            RESORDER = I
        END IF
   	END DO
      IF (FINAL .EQ. 1) THEN
            RESORDER = NORESERVOIRS
      END IF
C         IF (UH_STRING(1:4) .ne. 'NONE') THEN       ! read UH_S grid, not make it
C           print*, 'reading UH_S grid from file'
C           open(98, file=UH_STRING, status='old')
C           DO N = 1,NOB(RESORDER)
C             READ(98, *) (UH_S(N,K,RESORDER), K = 1,KE+UH_DAY-1)
C           END DO
C         ELSE				         ! make UH_S grid, and save it
C           print*, 'making UH_S grid...it takes a while...'
C           print*, 'NOTE:  your new UH_S grid file will be written in the'
C           print*, '       directory you run from, and will be called'
C           print*, '            '//NAME5//'.uh_s'
C           print*, '       save this file and specify it in your station'
C           print*, '       location file to avoid this step in the future'
C   		open(98, iostat=stat, file = NAME5//'.uh_s', status='replace')
C   		if (stat .NE. 0) open(98, file = NAME5//'.uh_s', status='new')
C   		DO N = 1, NOB(RESORDER)
C             print*, 'grid cell', N,' out of', NOB(RESORDER)
C             DO K = 1,UH_DAY
C               UH_DAILY(N,K) = 0.0
C             END DO
C             I = CATCHIJ(N,1,RESORDER)
C             J = CATCHIJ(N,2,RESORDER)
C             DO K = 1,24
C               FR(K,1) = 1.0 / 24.0
C               FR(K,2) = 0.0
C             END DO
C             DO K = 25,TMAX
C               FR(K,1) = 0.0
C               FR(K,2) = 0.0
C             END DO
C   100      CONTINUE
C             IF ((I .NE. PI) .OR. (J .NE. PJ)) THEN
C               DO T = 1, TMAX
C                 DO L = 1, LE
C                    IF ((T-L) .GT. 0) THEN
C                      FR(T,2) = FR(T,2) + FR(T-L,1)*UHM(I,J,L)
C                    END IF
C                  END DO
C               END DO
C               II = DIREC(I,J,1)
C               JJ = DIREC(I,J,2)
C               I = II
C               J = JJ
C               DO T = 1, TMAX
C                  FR(T,1) = FR(T,2)
C                  FR(T,2) = 0.0
C               END DO
C             END IF
C             IF ((I .NE. PI) .OR. (J .NE. PJ)) THEN
C               GOTO 100
C             END IF
C             DO T = 1,TMAX
C               TT = (T+23)/24
C               UH_DAILY(N,TT) = UH_DAILY(N,TT) + FR(T,1)
C             END DO
C           END DO
C           DO N = 1,NOB(RESORDER)
C             DO K = 1, KE+UH_DAY-1
C               UH_S(N,K,RESORDER) = 0.0
C             END DO
C           END DO
C           DO N = 1,NOB(RESORDER)
C             DO K = 1,KE
C               DO U =1,UH_DAY
C                 UH_S(N,K+U-1,RESORDER) = UH_S(N,K+U-1,RESORDER) +
C      &                         UH_BOX(N,K) * UH_DAILY(N,U)
C               END DO
C             END DO
C             SUM = 0
C             DO K = 1,KE+UH_DAY-1
C               SUM = SUM + UH_S(N,K,RESORDER)
C             END DO
C             DO K = 1,KE+UH_DAY-1
C               UH_S(N,K,RESORDER) = UH_S(N,K,RESORDER) / SUM
C             END DO
C           END DO
C c   write out the grid for future reference...
C           PRINT *, 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO', NOB(RESORDER)
C           DO N = 1,NOB(RESORDER)
C             WRITE(98, *) (UH_S(N,K,RESORDER), K = 1,KE+UH_DAY-1)
C           END DO
C         END IF
C        close(98)
      RETURN
      END
