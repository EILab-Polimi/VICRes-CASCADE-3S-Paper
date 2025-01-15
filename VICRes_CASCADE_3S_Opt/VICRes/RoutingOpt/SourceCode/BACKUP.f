

        DO I = 1, NDAY
            DO J = 1, NO_STAS-1
                CURRENTYEAR = START_YEAR + INT(I/365)		! approximate, does not consider leap years
            IF (CURRENTYEAR>=OPEYEAR(J)) THEN
                VRESER(J) = VRESERTHAT(J)
                REALHEAD(J) = HYDRAUHEAD(J)
                QRESER(J) = QRESERTHAT(J)
            ELSE
                VRESER(J) = 0.001
                REALHEAD(J) = 0.001
                QRESER(J) = 999999
            END IF
            FLOWIN(J,I) = RESFLOWS(J,I)
c           Calculate the designed water level
c           Note RULE = 1: simplified rule curve - 2: rule curve - 3: operating rules: - 4 pre-defined time-series data - 5 12 month operating rule
            CRTDATE = INT(1.0* mod(I,365)+(START_MO-1)*30)						! approximate
            IF ((RULE(J) .EQ. 1) .or. (RULE(J) .EQ. 2)) THEN   ! (Options 1 and 2: rule curves)
                IF (CURRENTYEAR<OPEYEAR(J)) THEN
                    FLOWOUT(J,I) = FLOWIN(J,I)
                    VOL(J,I+1) = VOL(J,I)
                    GOTO 123
                END IF
                IF (RULE(J) .EQ. 1) THEN
                    IF (OP1(J,1)>OP1(J,2)) THEN
                        TEMPO = OP1(J,1)
                        OP1(J,1) = OP1(J,2)
                        OP1(J,2) = TEMPO
                    END IF
                    ! Caculate target water level
                    IF ((CRTDATE .GT. OP1(J,1)) .and. (CRTDATE .LT. OP1(J,2))) THEN
                        IF (HMAX(J) .GE. HMIN(J)) THEN
                            DESIGNWL=(CRTDATE-OP1(J,1))/(OP1(J,2)-OP1(J,1))
     &                      *(HMAX(J)-HMIN(J))
                        ELSE
                            DESIGNWL=(HMIN(J)-HMAX(J))-(CRTDATE-OP1(J,1))/(OP1(J,2)-OP1(J,1))
     &                      *(HMIN(J)-HMAX(J))
                        END IF
                    ELSE IF (CRTDATE .GE. OP1(J,2)) THEN
                        IF (HMAX(J) .GE. HMIN(J)) THEN 
                           DESIGNWL=(HMAX(J)-HMIN(J))
     &                     -(CRTDATE-OP1(J,2))/(365-OP1(J,2)+OP1(J,1))*
     &                     (HMAX(J)-HMIN(J))
                        ELSE
                           DESIGNWL=(CRTDATE-OP1(J,2))/(365-OP1(J,2)+OP1(J,1))*
     &                     (HMIN(J)-HMAX(J))
                        END IF
                    ELSE
                        IF (HMAX(J) .GE. HMIN(J)) THEN 
                            DESIGNWL=(HMAX(J)-HMIN(J))
     &                      -(CRTDATE+365-OP1(J,2))/(365-OP1(J,2)+OP1(J,1))*
     &                     (HMAX(J)-HMIN(J))
                        ELSE
                        DESIGNWL=(CRTDATE+365-OP1(J,2))/(365-OP1(J,2)+OP1(J,1))*
     &                  (HMIN(J)-HMAX(J))
                    END IF
                END IF
                DESIGNWL = DESIGNWL + (HMIN(J) - HRESERMIN(J))
            ELSE
                DESIGNWL = RESLV(J,CAL_MONTH(CRTDATE))
                DESIGNWL = DESIGNWL - H0(J)
            END IF
            HTK(J,I) = DESIGNWL + H0(J)													! water head
            CURRENTWL = VOL(J,I) * (HRESERMAX(J)-H0(J))/VRESER(J)
            IF (CURRENTWL>=DESIGNWL) THEN												! Zone 3
                IF ((VOL(J,I)+FLOWIN(J,I)*24*3.6 -QRESER(J)*24*3.6)						! Case 2
     &          >(DESIGNWL* VRESER(J)) /(HRESERMAX(J)-H0(J))) THEN
                    VOL(J,I+1) = VOL(J,I) + FLOWIN(J,I)*24*3.6-QRESER(J)*24*3.6
                    FLOWOUT(J,I) = QRESER(J)
                    VOL(J,I+1)=(DESIGNWL*VRESER(J))/(HRESERMAX(J)-H0(J))
                ELSE																	! Case 1
                    VOL(J,I+1)=(DESIGNWL*VRESER(J))/(HRESERMAX(J)-H0(J))
                    FLOWOUT(J,I)=(VOL(J,I)-VOL(J,I+1))/24/3.6 + FLOWIN(J,I)
                END IF
                GOTO 123
            ELSE																		! Zone 2
                IF ((VOL(J,I)+FLOWIN(J,I)*24*3.6)>((DESIGNWL* VRESER(J))				! Case 2
     &           /(HRESERMAX(J)-H0(J)))) THEN
                    VOL(J,I+1)=(DESIGNWL * VRESER(J))/(HRESERMAX(J)-H0(J))
                    FLOWOUT(J,I)=FLOWIN(J,I)-(VOL(J,I+1)-VOL(J,I))/24/3.6
                    IF (FLOWOUT(J,I)>QRESER(J)) THEN
                        VOL(J,I+1)=VOL(J,I+1)+(FLOWOUT(J,I)-QRESER(J))*24*3.6
                        FLOWOUT(J,I) = QRESER(J)
                    END IF
                    GOTO 123
                ELSE																	! Case 1 (this case covers Zone 1 + Zone 2 case 1)
                    VOL(J,I+1) = VOL(J,I) + FLOWIN(J,I)*24*3.6
                    GOTO 123
                END IF
            END IF
        ELSE IF (RULE(J) .EQ. 3) THEN ! opearting rule (Option 3)
        ! Note x1 and x4 in radian (0 to pi/2), not degree
            IF (VOL(J,I) < VDEAD(J)) THEN ! below dead water level
                FLOWOUT(J,I) = 0																					! case 1
            ELSE IF (VOL(J,I) .LE. X2(J)) THEN ! hedging
                IF ((VOL(J,I)-VDEAD(J)+FLOWIN(J,I)*24*3.6) .LE. (Demand(J)+(VOL(J,I)-X2(J))*tan(X1(J))*24*3.6)) THEN		! discharge more than the water amount in reservoir (case 2)
                    FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6+FLOWIN(J,I)											! discharge all of water in the reservoir
                ELSE	! (case 3)
                    FLOWOUT(J,I) = Demand(J)+(VOL(J,I)-X2(J))*tan(X1(J))
                END IF
            ELSE IF (VOL(J,I).GE. X3(J)) THEN 	! spilling
                IF ((Demand(J) + (VOL(J,I)-X3(J))*tan(X4(J)))*24*3.6>(VOL(J,I))-VDEAD(J)) THEN		! discharge more than the water in reservoir (just to make sure)
                    FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6										! discharge all of water in the reservoir
                ELSE IF ((Demand(J) + (VOL(J,I)-X3(J))*tan(X4(J)))>QRESER(J)) THEN
                    FLOWOUT(J,I) = QRESER(J)
                ELSE	 ! case 5
                    FLOWOUT(J,I) = Demand(J) + (VOL(J,I)-X3(J))*tan(X4(J))
                END IF
            ELSE ! releasing
                IF (Demand(J)*24*3.6>(VOL(J,I)-VDEAD(J))) THEN
                    FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6
                ELSE
                    FLOWOUT(J,I) = Demand(J)			! case 4
                END IF
            END IF
            IF (FLOWOUT(J,I)>QRESER(J)) THEN			! double check just in case users chose unrealistic x1 and x4
                FLOWOUT(J,I) = QRESER(J)
            END IF
            VOL(J,I+1) = VOL(J,I) + (FLOWIN(J,I)-FLOWOUT(J,I))*24*3.6
            GOTO 123
        ELSE IF (RULE(J) .EQ. 4) THEN! pre-defined time series (Option 4)
            IF ((OP2(J,I+STARTDAY(J)) .GT. QRESER(J)) .AND. (VOL(J,I) .LT. VRESER(J))) THEN
                OP2(J,I+STARTDAY(J)) = QRESER(J)
            END IF
            IF (OP2(J,I+STARTDAY(J))*24*3.6 .GT. ((VOL(J,I)-VDEAD(J)))+FLOWIN(J,I)*24*3.6) THEN
                FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6 + FLOWIN(J,I)
            ELSE
                FLOWOUT(J,I) = OP2(J,I+STARTDAY(J))
            END IF
            VOL(J,I+1) = VOL(J,I) + (FLOWIN(J,I)-FLOWOUT(J,I))*24*3.6
            GOTO 123
        ELSE IF (RULE(J) .EQ. 5) THEN ! note: this option is similar to OP3 but for a periodic demand
            ! Note x1 and x4 in radian (0 to pi/2), not degree, this part can be shorthen
            X1(J) = OP5X1(J,CAL_MONTH(CRTDATE))
            X2(J) = OP5X2(J,CAL_MONTH(CRTDATE))
            X3(J) = OP5X3(J,CAL_MONTH(CRTDATE))
            X4(J) = OP5X4(J,CAL_MONTH(CRTDATE))
            Demand(J) = DEMAND5(J,CAL_MONTH(CRTDATE))
            IF (VOL(J,I) < VDEAD(J)) THEN ! below dead water level
                FLOWOUT(J,I) = 0																					! case 1
            ELSE IF (VOL(J,I) .LE. X2(J)) THEN ! hedging
                IF ((VOL(J,I)-VDEAD(J)+FLOWIN(J,I)*24*3.6) .LE. (Demand(J)+(VOL(J,I)-X2(J))*tan(X1(J))*24*3.6)) THEN		! discharge more than the water amount in reservoir (case 2)                    
                    FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6+FLOWIN(J,I)											! discharge all of water in the reservoir
                ELSE	! (case 3)
                    FLOWOUT(J,I) = Demand(J)+(VOL(J,I)-X2(J))*tan(X1(J))
                END IF
            ELSE IF (VOL(J,I).GE. X3(J)) THEN 	! spilling
                IF ((Demand(J) + (VOL(J,I)-X3(J))*tan(X4(J)))*24*3.6>(VOL(J,I))-VDEAD(J)) THEN		! discharge more than the water in reservoir (just to make sure)
                    FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6										! discharge all of water in the reservoir
                ELSE IF ((Demand(J) + (VOL(J,I)-X3(J))*tan(X4(J)))>QRESER(J)) THEN
                    FLOWOUT(J,I) = QRESER(J)
                ELSE	 ! case 5
                    FLOWOUT(J,I) = Demand(J) + (VOL(J,I)-X3(J))*tan(X4(J))
                END IF
            ELSE ! releasing
                IF (Demand(J)*24*3.6>(VOL(J,I)-VDEAD(J))) THEN
                    FLOWOUT(J,I) = (VOL(J,I)-VDEAD(J))/24/3.6
                ELSE
                    FLOWOUT(J,I) = Demand(J)			! case 4
                END IF
            END IF
            IF (FLOWOUT(J,I)>QRESER(J)) THEN			! double check just in case users chose unrealistic x1 and x4
                FLOWOUT(J,I) = QRESER(J)
            END IF
            VOL(J,I+1) = VOL(J,I) + (FLOWIN(J,I)-FLOWOUT(J,I))*24*3.6
            GOTO 123
 123    END IF
        ! Check if there are any negative values
        IF (ENERGYPRO(J,I)<0) THEN
            ENERGYPRO(J,I)=0
        END IF
        IF (FLOWOUT(J,I)<0) THEN
            FLOWOUT(J,I)=0
        END IF
        IF (VOL(I,J)<0)THEN ! Not allow dropping below the minimum water level (mostly due to evaporation)
            VOL(I,J)=0
        END IF
c       Remote water for irrigation
        IF (FLOWOUT(J,I)>=IRRIGATION(J,I)) THEN
            FLOWOUT(J,I) = FLOWOUT(J,I) - IRRIGATION(J,I)
        ELSE
            FLOWOUT(J,I) = 0
        END IF
        ! Calculate energy production
        IF (VOL(J,I+1)<=VRESER(J)) THEN
            ENERGYPRO(J,I) = 0.9 * 9.81 * FLOWOUT(J,I)
        ELSE
            ENERGYPRO(J,I) = 0.9 * 9.81 * QRESER(J)
        END IF
        !IF (J .EQ. 3) THEN
        !    WRITE(*,*) RES_DIRECT(J,1),' - ',ENERGYPRO(J,I),'  ',HHO(J,I),'  ',HHO(J,I+1),'  ',FLOWOUT(J,I),'  ',HRESERMAX(J),'  ',REALHEAD(J)
        !END IF
        ! Update water losses due to seepage and infiltration
        ! Infiltration is permanent losses + water seepage is added to outflow
        ! Note seepage occurs until there is no water left (considering dead volume also)
        IF (VOL(J,I+1) - (SEEPAGE(J)+INFILTRATION(J))*24*3.6 .GT. 0) THEN
            VOL(J,I+1) = VOL(J,I+1) - (SEEPAGE(J)+INFILTRATION(J))*24*3.6
            !Note that seepage does not contribute to energy production, so we add here
            FLOWOUT(J,I) = FLOWOUT(J,I) + SEEPAGE(J)
        ELSE
            VOL(J,I+1) = 0
        END IF
c       Check the neccesity to spill water
        IF (VOL(J,I+1)>VRESER(J)) THEN
            FLOWOUT(J,I) = FLOWOUT(J,I)+(VOL(J,I+1)-VRESER(J))/24/3.6
            VOL(J,I+1) = VRESER(J)
        END IF
        HHO(J,I)=VOL(J,I)/VRESER(J)*(HRESERMAX(J)-H0(J))+H0(J)
        HHO(J,I+1)=VOL(J,I+1)/VRESER(J)* (HRESERMAX(J)-H0(J))+H0(J)
        ! Note: hydraulic head calculated from the maximum water level
        ENERGYPRO(J,I) = ENERGYPRO(J,I) *
     &  (( HHO(J,I)+HHO(J,I+1))/2-(HRESERMAX(J)-REALHEAD(J)))/1000			! this part is for hydropower production estimation, ignore if work with irrigation reservoirs
        ! Propagate water to the downstream reservoir, considering the time lag
        RESORDER = NO_STAS
        DO K = 1, NO_STAS
            IF ((RES_DIRECT(K,1) .EQ. RES_DIRECT(J,2))) THEN
                RESORDER = K
            END IF
        END DO
        IF (FLOWOUT(J,I) .LT. 0) THEN
            FLOWOUT(J,I) = 0
        END IF
        IF (RES_DIRECT(J,2) .EQ. 0) THEN
            RESFLOWS(NO_STAS,I+1+INT(TRVLTIME(J))) = RESFLOWS(NO_STAS,I+1+INT(TRVLTIME(J))) + FLOWOUT(J,I)
        ELSE
            RESFLOWS(RESORDER,I+1+INT(TRVLTIME(J))) = RESFLOWS(RESORDER,I+1+INT(TRVLTIME(J))) + FLOWOUT(J,I)
        END IF
      !   WRITE(*,*) '-------------------------------------------------'
      !   WRITE(*,*) 'Hdesign ',HTK(J,I)
      !   WRITE(*,*) 'CRTDATE',CRTDATE,'Reservoir No.',J,' FLOWIN ',FLOWIN(J,I),
     &!  'FLOWOUT',FLOWOUT(J,I),' VOL ',VOL(J,I),' ELEC', ENERGYPRO(J,I)
        END DO
        IF (RESFLOWS(NO_STAS,I) .LT. 0) THEN
             RESFLOWS(NO_STAS,I) = 0
        END IF
        FLOW(I) = RESFLOWS(NO_STAS,I)
c       We prepare the array to pass to cascade (for the moment only one reach)
        a(1,1) = I*1.0
        a(1,2) = FLOW(I)
        a(1,3) = 3.23
C       IN REALTà INUTILE RIDEFINIRE QUESTI DUE OGNI LOOP 
        sh = shape(a)
        do j=1,4   
c           reads from the socket one message header: STATUS, INIT, SENDDATA, GETDATA 
            call readbuffer(socket, header, msglen)
c           All the write are commented. We can decomment them if we need to control what is happening between client and server
c            write(*,*) "@ message from server: ", trim(header)

c           if header is STATUS the wrapper is asking the client what we are doing: we need inizialization (just for the first step), we have data to pass to cascade or we are waiting to receive the cascade results    
            if (ADJUSTL(trim(header)) .eq. "STATUS") then   ! the wrapper is inquiring on what we are doing
               if (.not. isinit) then
                  call writebuffer(socket, "NEEDINIT    ", msglen)   ! signals that we need initialization
c                  write(*,*) "@ message to server: NEEDINIT"
               elseif (hasdata) then
                  call writebuffer(socket, "HAVEDATA    ", msglen)   ! signals that we are done computing and can data
c                  write(*,*) "@ message to server: HAVEDATA"
               else
                  call writebuffer(socket, "READY       ", msglen)   ! we are idling and eager to compute something
c                  write(*,*) "@ message to server: READY"
               endif

            elseif (ADJUSTL(trim(header)) .eq. "INIT") then     ! the driver is kindly sending a string for initialization
c               write(*,*) " Initializing system from server"
               isinit=.true.   ! we actually do nothing with this string, thanks anyway. could be used to pass some information 

            elseif (ADJUSTL(trim(header)) .eq. "SENDDATA") then  ! Server wants to send data to the driver
               if (.not. isinit) then
                     write(*,*) "Driver not iniliasied."
               elseif (hasdata) then
                    write(*,*) "Driver has data to send back to server"
               else   ! Driver is ready to receive data
c                 for the moment we do not know what b will be (let's assume same shape of a but unlikely)               
                  call readbufferreal(socket, msgbuffer, 4*n_reach)
                  b = reshape(msgbuffer, shape(a))
c                 print the array received, decomented now to see that everything works but to coment when we will use the acutal model
                  write(*,*) "Received array from server:"
                  do li=1,sh(1)
                     write(*,*) b(1,:)
                  end do
                  hasdata = .true.
               end if

            elseif (ADJUSTL(trim(header)) .eq. "GETDATA") then   ! Server signaling driver to send data
               if (.not. isinit) then
                   write(*,*) "Driver not iniliasied."
               elseif (.not. hasdata) then
                  write(*,*) "Driver does not have data to send"
               else
                  call writebuffer(socket, "DATAREADY   ", msglen)
c                  write(*,*) "@ message to server: DATAREADY"

                  msgbuffer = reshape(a, [n_reach])
                  call writebufferreal(socket, msgbuffer, 4*n_reach)   ! writing data
c                 print the array sent, decomented now to see that everything works but needs to be commented when we will use the acutal model
                  write(*,*) "Sent array to server:"
                  do li=1,sh(1)
                     write(*,*) a(1,:)
                  end do                  
                  hasdata = .false.
               end if

            else
               write(*,*) " unexpected header ", header
               stop "ended"
            endif
        enddo        
      END DO
      RETURN
 9001 WRITE(*,*) 'Error reading UH data'
 9002 WRITE(*,*) 'Error in reading reservoir data'
 9003 WRITE(*,*) 'Error in reading irrigation data'
      END
C     END OF FILE
C************************************************************************************************************************************************************************************

C************************************************************************************************************************************************************************************
C       Convert from day to month
C************************************************************************************************************************************************************************************
        INTEGER FUNCTION CAL_MONTH(CRTDATE)
        INTEGER CRTDATE
        IF (CRTDATE .LE. 31) THEN !jan
            CAL_MONTH = 1
        ELSE IF (CRTDATE .LE. 59) THEN ! feb
            CAL_MONTH = 2
        ELSE IF (CRTDATE .LE. 90) THEN ! mar
            CAL_MONTH = 3
        ELSE IF (CRTDATE .LE. 120) THEN ! apr
            CAL_MONTH = 4
        ELSE IF (CRTDATE .LE. 151) THEN ! may
            CAL_MONTH = 5
        ELSE IF (CRTDATE .LE. 181) THEN ! jun
            CAL_MONTH = 6
        ELSE IF (CRTDATE .LE. 212) THEN ! jul
            CAL_MONTH = 7
        ELSE IF (CRTDATE .LE. 243) THEN ! aug
            CAL_MONTH = 8
        ELSE IF (CRTDATE .LE. 273) THEN ! sep
            CAL_MONTH = 9
        ELSE IF (CRTDATE .LE. 304) THEN ! oct
            CAL_MONTH = 10
        ELSE IF (CRTDATE .LE. 334) THEN ! nov
            CAL_MONTH = 11
        ELSE ! dec
            CAL_MONTH = 12
        END IF
        RETURN
        END