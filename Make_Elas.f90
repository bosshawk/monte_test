     SUBROUTINE MAKEELAS
    
     use PARAMETER_2
    
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)
     DIMENSION DCSEC(201,90)
     COMMON /NO4/ENGY(0:201),DIELA(0:201,0:90)

     OPEN(1,FILE='mott_data_Cu.txt',STATUS='OLD')

     DTHT=2.*PAI/180d0
     DO I=1,201
        READ(1,10) (DCSEC(I,J),J=1,89)
  10    FORMAT(7D12.4)
        DCSEC(I,90)=DCSEC(I,89)
        TOT=0D0
        DO J=1,90
           TOT=TOT+DCSEC(I,J)*DSIN(DBLE(J)*DTHT)*2.0D0*PAI*DTHT
           DIELA(I,J)=TOT*ROU*Avo/Ai
        END DO
     END DO

     CLOSE(1)

     END
