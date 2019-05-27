!**********************************************************************!
     SUBROUTINE INTRND
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)
     COMMON /RANDOM/MT_R(0:623),MTI_R,N_R,M_R,ISEED,IMTRXA,  &
&              IUPMSK,ILWMSK,ITMSKB,ITMSKC

!     ISEED=4338
!     ISEED=4385
     ISEED=4356
     N_R=624; M_R=397
     MTI_R=N_R+1
     IMTRXA='9908B0DF'X
     IUPMSK='80000000'X
     ILWMSK='7FFFFFFF'X
     ITMSKB='9D2C5680'X
     ITMSKC='EFC60000'X

     END
!**********************************************************************!
     SUBROUTINE SGENRND
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)
     COMMON /RANDOM/MT_R(0:623),MTI_R,N_R,M_R,ISEED,  &
&            IMTRXA,IUPMSK,ILWMSK,ITMSKB,ITMSKC

     MT_R(0)=IAND(ISEED,'FFFFFFFF'X)
     DO MTI_R=1,N_R-1
       MT_R(MTI_R)=IAND(69069*MT_R(MTI_R-1),'FFFFFFFF'X)
     END DO

     END
!**********************************************************************!
     SUBROUTINE DRND(RND)
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)
     COMMON /RANDOM/MT_R(0:623),MTI_R,N_R,M_R,ISEED,  &
&            IMTRXA,IUPMSK,ILWMSK,ITMSKB,ITMSKC
     DIMENSION MAG01(0:1)

 10  MAG01(0)='0'X
     MAG01(1)=IMTRXA
     IF(MTI_R >= N_R) THEN
       IF(MTI_R == N_R+1) CALL SGENRND
         DO KK=0,N_R-M_R-1
           IY=IOR(IAND(MT_R(KK),IUPMSK),IAND(MT_R(KK+1),ILWMSK))
           MT_R(KK)=IEOR(IEOR(MT_R(KK+M_R),ISHFT(IY,-1))  &
&                                        ,MAG01(IAND(IY,'1'X)))
         END DO
         DO KK=N_R-M_R,N_R-2
           IY=IOR(IAND(MT_R(KK),IUPMSK),IAND(MT_R(KK+1),ILWMSK))
           MT_R(KK)=IEOR(IEOR(MT_R(KK+(M_R-N_R)),ISHFT(IY,-1))  &
&                                        ,MAG01(IAND(IY,'1'X)))
         END DO
         IY=IOR(IAND(MT_R(N_R-1),IUPMSK),IAND(MT_R(0),ILWMSK))
         MT_R(N_R-1)=IEOR(IEOR(MT_R(M_R-1),ISHFT(IY,-1))  &
&                                        ,MAG01(IAND(IY,'1'X)))
         MTI_R=0
       ENDIF

     IY=MT_R(MTI_R)
     MTI_R=MTI_R+1
     IY=IEOR(IY,ISHFT(IY,-11))
     IY=IEOR(IY,IAND(ISHFT(IY,7),ITMSKB))
     IY=IEOR(IY,IAND(ISHFT(IY,15),ITMSKC))
     IY=IEOR(IY,ISHFT(IY,-18))
     IF(IY >= 0) THEN
       RND=DBLE(IY)/4294967295D0
!       IF(RND <= 0D0) GOTO 10
     ELSE
       IX=ISHFT(IY,1)
       IX=ISHFT(IX,-1)
       RND=( 2147483648D0+DBLE(IX))/4294967295D0
!       IF(1D0-RND <= 0D0) GOTO 10
     ENDIF

     END
