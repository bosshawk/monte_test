! モンテカルロシミュレーション ！ 距離：ｍ単位　エネルギー：eV単位
    
    !Monte_Mott_Fog_Field_Nishino4
    
	Program Monte_Mott_Fog_Field
    
    use PARAMETER_1
    use PARAMETER_2
    
	IMPLICIT DOUBLEPRECISION (A-H,O-Z)
	
    COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI 
    COMMON /NO4/ENGY(0:201),DIELA(0:201,0:90)
    common/Iend/Iend
    common/countPE1/count_BSEyld,BSEyld
    common/countSE3/count_SEyld,SEyld,SE50yld,count_SE50yld
    
    common/rad/fogPE(0:100000)
    common/rad2/fogSE(0:100000),fogSE2(0:100000)
    common/count1/PEsample,PElens
    common/count2/SEsample,SElens
    
    
	OPEN(2,FILE='MONTE.txt')
    OPEN(3,FILE='monte2.txt')
    OPEN(4,FILE='Field.txt')
    open(5,file='yield.csv')
    
    open(6,file='count.csv')
    open(7,file='count2.csv')
    
    CALL initial  !乱数・エネルギー分割・材料定数・Mott断面積読み込み

    !Do JJ=80,200,5
    !E0=ENGY(JJ)
    
    !do aa=-200,200,5
    !   Vs=aa
    
    BSEyld=0;SEyld=0;SE50yld=0
    PElens=0;PEsample=0
    SElens=0;SEsample=0
    
    
    do reset=0,100000
    fogPE(reset)=0
    enddo
    
    do reset=0,100000
    fogSE(reset)=0
    enddo
  
!-------------------------------------------------------------------------------------------    
    DO I=1,imax
    print*, I,E0
    
	E=E0+Vs                             !　単位：eV
    Ebend=PEstop                        !　１次電子計算終了エネルギー
9   X=0;Y=0;Z=0+1d-18			        !　単位：m
    count_goto10=0                      !  初期入射でCT1にマイナスが選択されるともう一度計算
    Iend=0                              !  Rmin
    count_BSEyld=0                      !  試料から出た一発目のみをyieldとしてカウント
    
	CT1=1;   ST1=0;   CP1=1;   SP1=0	!垂直入射
    
10 CALL COLLELAS(E,HJK)
    count_goto10=count_goto10+1
	CALL DRND(R)
	Ds=-HJK*LOG(R)
	xb=x; yb=y; zb=z
    Eb=E
	Call xyz(ds,x,y,z)
    
    if(count_goto10==1)then
        if(z<0)goto 9
    endif
    
    call kyoukaiP(xb,yb,zb,x,y,z,E,Ds,Ebend)
        if(Iend==1)cycle
        
        if(E > Ebend)then
            goto 10
        else
            call countPE(x,y,z)
            cycle
        endif
        
    END DO                         
!-------------------------------------------------------------------------------------------    
	BSEyield=real(BSEyld+SE50yld)/real(imax)
	SEyield=real(SEyld)/real(imax)

    write(5,100) E0,BSEyield,E0,SEyield
100 format (f0.0,',',f0.3,',',f0.0,',',f0.3)
    
    do rad=0,300,1
    WRITE(6,200)rad/10,fogPE(rad),fogSE(rad),fogPE(rad)+fogSE(rad),E0,VS
    enddo
    
200 format (f0.5,',',f0.5,',',f0.5,',',f0.5,',',f0.5,',',f0.5)
    
    !do rad=0,300,1    
    !if(rad==0)then
    !WRITE(7,200)rad/10,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) / (pai*(1d-4)**2)  )*1d-3,   E0,VS
    !!print*,(1d-4)**2,(((rad*1d-4)+1d-4)*((rad*1d-4)+1d-4))
    !!print*,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) / (pai*(1d-4)**2)  )*1d-3,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) /(pai*((rad*1d-4)+1d-4)*((rad*1d-4)+1d-4)))*1d-3
    !else    
    !WRITE(7,200)rad/10,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) / (pai*((rad*1d-4)+1d-4)**2)-(pai*(rad*1d-4)**2))*1d-3,   E0,VS
    !!print*,(((rad*1d-4)+1d-4)**2)-((rad*1d-4)**2),((((rad*1d-4)+1d-4)*((rad*1d-4)+1d-4))-((rad*1d-4)*(rad*1d-4)))
    !!print*,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) / (pai*((rad*1d-4)+1d-4)**2)-(pai*(rad*1d-4)**2))*1d-3,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) /((pai*((rad*1d-4)+1d-4)*((rad*1d-4)+1d-4))-(pai*(rad*1d-4)*(rad*1d-4))))*1d-3
    !endif
    !enddo 
    do rad=0,300,1
    if(rad==0)then
    WRITE(7,200)rad/10,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) /(pai*((rad*1d-4)+1d-4)*((rad*1d-4)+1d-4)))*1d-3   
    else
    WRITE(7,200)rad/10,((  real(fogPE(rad)+fogSE(rad))/real(Imax)  ) /((pai*((rad*1d-4)+1d-4)*((rad*1d-4)+1d-4))-(pai*(rad*1d-4)*(rad*1d-4))))*1d-3
    endif
    enddo
    
    !WRITE(6,200) E0,real(PEsample)/real(imax),real(PElens)/real(imax),real(SEsample)/real(imax),real(SElens)/real(imax)
    !enddo
    END

!**********************************************************************!
! エネルギー損失 ！
	SUBROUTINE ELOSS(DEDS,E)
    use PARAMETER_2
	IMPLICIT DOUBLEPRECISION (A-H,O-Z)
	DJi=(9.76D0*Zi+58.5D0*Zi**(-0.19D0))  ! eV単位
	IF(E<=6.338d0*DJi)  THEN   !*1d5は[keV/cm]→[eV/m]の為、1d-6はAi,ROUの単位の為
        DEDS=7.85d4/(1.26D0*SQRT(E*DJi*1d-6))*Zi/Ai*ROU *1d5*1d-6
	ELSE
		DEDS=7.85d4/(E*1d-3)*Zi/Ai*LOG(1.166D0*E/DJI)*ROU *1d5*1d-6
	END IF
	END
	
	
!**********************************************************************!
! 座標 X,Y,Z ！
    Subroutine xyz(s,x,y,z)
	implicit doubleprecision(a-h,o-z)
    
    COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI 

	CT2=CT1;  ST2=ST1; CP2=CP1; SP2=SP1
	IF(abs(ST2) < 1d-15) THEN						!垂直入射
		CT1=COMG;		SP1=SPHI
		ST1=SOMG;		CP1=CPHI
        
	ELSE
		CT1=CT2*COMG-ST2*SOMG*CPHI
		ST1=SQRT(1-CT1**2D0)
		AA=(COMG-CT2*CT1)/(ST2*ST1)
		BB=(SPHI*SOMG)/ST1
		SP1=AA*SP2+BB*CP2
		CP1=AA*CP2-BB*SP2
	END IF
	X=X+s*ST1*CP1
	Y=Y+s*ST1*SP1
	Z=Z+s*CT1
    end
    
    Subroutine countPE(x,y,z)
    use PARAMETER_1
	IMPLICIT DOUBLEPRECISION (A-H,O-Z)
    common/rad/fogPE(0:100000)
    common/count1/PEsample,PElens
    if(z>0d0) then      
!-------------------------------0.1mmに10プロット-----------------------------------    
    do rad=0,300,1
    rad0=dsqrt(X**2+Y**2)
    if(rad0>rad*1d-4 .and. rad0<(rad*1d-4)+1d-4) fogPE(rad)=fogPE(rad)+1
    !print*,rad*1d-4,(rad*1d-4)+1d-4
    enddo    
!-----------------------------------------------------------------------------------        
          
    rad=dsqrt(X**2+Y**2)
    
    if(rad>Rmin1)PEsample=PEsample+1       
    
    elseif(z<WD)then 
        
    PElens=PElens+1          
        
!-----------------------------------------------------------------------------------
    
    endif

    
    end
    
	

