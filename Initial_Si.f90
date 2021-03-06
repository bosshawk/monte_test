     SUBROUTINE INITIAL_SI
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)

     COMMON /NO1/PAI,E,C,DM0,ENG0,A0,AVNUM
     COMMON /NO2/EZERO,NEL

     COMMON /SI_1/ANUM,AWEIT,DENST
     COMMON /SI_3/EF,EP,FAI
     COMMON /SI_5/OMEGA(10),AUGEDE(10)
     COMMON /SI_7/V_SI,ESTOP_SI
     COMMON /KYOUKAI/ZKYOUKAI,V_KYOUKAI

     COMMON /PMMA_3/EF_PM,EP_PM,FAI_PM

     RS=(47.132D0/EP)**(2D0/3D0)
     EF=50.1123D0/(RS*RS)
     ESTOP_SI=EF+FAI

!	 V_SI=EF+FAI          !-- Si -> SINKUU --!
     V_SI=EF+FAI-FAI_PM   !-- Si -> PMMA --!
!PRINT *,'POT.   SI->PM,   SI->VAC=',V_SI,EF+FAI ; PAUSE
     V_KYOUKAI=V_SI

!   オージェ電子 か 特性Ｘ線 の比
!   FLUORESCENCE YIELD --- OMEGA --- FROM EPMA BOOK P.142
!     OMA=-0.0217D0; OMB=0.0332D0; OMC=-1.14D-6
!     OMEG=(OMA+(OMB+OMC*ANUM**2)*ANUM)**4
!     OMEGA(1)=OMEG/(1D0+OMEG)
!     OMEGA(3)=ANUM**4D0/(1.02D8+ANUM**4)!
!     OMEGA(2)=OMEGA(3) !??
!     OMEGA(4)=OMEGA(3) !??

!   断面積の計算
     CALL MAKEELAS;  PRINT *, '  MAKE_ELAS  OK'
!     CALL MAKEPLAS;  PRINT *, '  MAKE_PLAS  OK'
!     CALL MAKETUNG;  PRINT *, '  MAKE_TUNG  OK'
!     CALL MAKEGRYZ;  PRINT *, '  MAKE_GRYZ  OK'

     END
