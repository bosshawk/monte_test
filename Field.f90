      SUBROUTINE Field(E,x,y,z,x1,y1,z1,x2,y2,z2)
    
      use PARAMETER_1
    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI 
      
      
      
      EEE=E
      
      E=E-V0
      
      DT=1d-11				!---真空中 固定ステップ時間 0.1ns　
      YE=1.7588196D+11
      Q=1.602d-19
      Efield=-Vs/WD    
      Eacc=Ye*Efield

      DRDT=DSQRT(E*YE*2.D0)
      DXDT=DRDT*ST1*CP1
      DYDT=DRDT*ST1*SP1
      DZDT=DRDT*CT1     

  10  x1=x;y1=y;z1=z
      x=x+dxdt*dt
      y=y+dydt*dt
      z=z+dzdt*dt+0.5*Eacc*dt*dt

      dzdt=dzdt+Eacc*dt
       
      x2=x;y2=y;z2=z
      
if(view==1)WRITE(4,*) X,Y,Z	

	  if(z > 0d0 .or. z < WD ) goto 20	! 試料または対物に突入した(Fieldの無駄毛はkyoukaiで処理)
	  
      goto 10

  20  RRR=DSQRT(DXDT**2+DYDT**2+DZDT**2)
      RR=DSQRT(DXDT**2+DYDT**2)
      IF(RR.NE.0) THEN
         CT1=DZDT/RRR
         ST1=RR/RRR
         CP1=DXDT/RR
         SP1=DYDT/RR
      ELSE
         CT1=DZDT/RRR
         ST1=RR/RRR
         CP1=1d0
         SP1=0d0
      ENDIF

      DRDT2=DXDT**2+DYDT**2+DZDT**2
   
      E=EEE!+V0

       RETURN
      END
