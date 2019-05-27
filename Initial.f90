! 乱数の初期化、Mott断面積(エネルギー分割、断面積読み込み)
     SUBROUTINE INITIAL
     
     use PARAMETER_2
     
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	 
     COMMON /NO4/ENGY(0:201),DIELA(0:201,0:90)

!   乱数の初期化
     CALL INTRND

!   エネルギーの初期化
     DO I=1,201
        ENGY(I)=10D0**((DBLE(I)-1D0)/40D0)
     END DO

!   Cuの材料定数定義・断面積読み込み
     CALL MAKEELAS;  PRINT *, '  MAKE_ELAS  OK'

     END
