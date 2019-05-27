     SUBROUTINE COLLELAS(E,DMFP)
     use PARAMETER_2
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)

     COMMON /NO4/ENGY(0:201),DIELA(0:201,0:90)
     COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI

     RA=40D0*DLOG10(E)+1d0
     K1=RA; K2=K1+1
     RATIO=(E-ENGY(K1))/(ENGY(K2)-ENGY(K1))
     
     DTHT=PAI/90d0
     CALL DRND(RAN)

   !---エネルギー ENGY(K1) の時---!
     R1=RAN*DIELA(K1,90)
     CALL BINSEARCH(R1,K1)
     ELAS1=DBLE(N1)+(R1-DIELA(K1,N1))/(DIELA(K1,N2)-DIELA(K1,N1))
   !---エネルギー ENGY(K2) の時---!
     R2=RAN*DIELA(K2,90)
     CALL BINSEARCH(R2,K2)
     ELAS2=DBLE(N1)+(R2-DIELA(K2,N1))/(DIELA(K2,N2)-DIELA(K2,N1))

     CALL SCATANGLE
     
     DIMFP=DIELA(K1,90)+RATIO*(DIELA(K2,90)-DIELA(K1,90))
     DMFP=1d0/DIMFP

!----------------内部サブルーチン----------------------------------------!
     CONTAINS
      !---二分探査を行う---!
       SUBROUTINE BINSEARCH(R0,K0)
          N1=0; N2=90
          DO WHILE(N2-N1 > 1)
             N0=(N2+N1)/2
             IF(R0 > DIELA(K0,N0)) THEN
                N1=N0
             ELSE
                N2=N0
             ENDIF
          END DO
       END SUBROUTINE BINSEARCH

      !---散乱角の決定---!
       SUBROUTINE SCATANGLE
         OMG=(ELAS1+RATIO*(ELAS2-ELAS1))*DTHT
         IF(OMG > PAI) OMG=PAI
         SOMG=DSIN(OMG)
         COMG=DCOS(OMG)
         CALL DRND(RAN)
         SPHI=DSIN(RAN*2D0*PAI)
         CPHI=DCOS(RAN*2D0*PAI)
       END SUBROUTINE SCATANGLE
!------------------------------------------------------------------------!

     END
