	Subroutine SE(Ep,DE,xpb,ypb,zpb,xp,yp,zp)	!--- 二次電子軌道 ---
    use PARAMETER_1
    use PARAMETER_2
	implicit doubleprecision(a-h,o-z)

    
    COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI 
    common/countSE3/count_SEyld,SEyld,SE50yld,count_SE50yld
	parameter (NSE=10000)
	dimension Es(NSE),xs(NSE),ys(NSE),zs(NSE),st1s(NSE),ct1s(NSE),sp1s(NSE),cp1s(NSE)
	
	
	Phi=4.45d0; Ef=5.0d0  !銅のフェルミエネルギーを5eV、仕事関数を4.45eVとする。
    Esend=V0			  ! 打ち切りエネルギー    
!		Ese=10.0	! Ese：SE1つを作るための平均エネルギー　ISはそのstepでの二次電子発生数
!		IS=int(DE/Ese)
!		CALL RANDOM_NUMBER (R)
!		if(R  <  (DE/Ese-dble(IS)) ) IS=IS+1
!		do ii=1,IS
    Er=DE                           !エネルギー保存の法則を満たす
10  CALL DRND(R)
   
    if(Er < Esend) return	
	A=(EP-Ef)/(EP-Ef-Phi)                                            
    ESE=(R*Ef-A*(Ef+Phi))/(R-A)     !二次電子のエネルギー分布　DEKIより抜粋                                   
    if(ESE > Er) goto 10
    Er=Er-ESE
    IF(ESE < Esend) goto 10  
    
    CALL DRND(R)
	X=Xpb+(Xp-Xpb)*R;  Y=Ypb+(Yp-Ypb)*R;  Z=Zpb+(Zp-Zpb)*R	!二次電子の「種」発生場所
    x1=x; y1=y; z1=z
    
    if(z>0 .and. dsqrt(X**2+Y**2)>Rmax1)then
        goto 10
    elseif(z<WD .and. (dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2))then
        goto 10
    endif
    
    E=ESE
!******************************!   
    call countSE1(x,y,z,E)
!******************************!   
	st1=0;ct1=1;sp1=0;cp1=1
    CALL DRND(R)
      COMG=1.0-2.0*R                 ! コサイン分布で電子放出                                  
      SOMG=DSQRT(1.0-COMG*COMG)                                         
    CALL DRND(R)
      SPHI=DSIN(R*2.0*PAI)           ! コサイン分布で電子放出             
      CPHI=DCOS(R*2.0*PAI)

	call SEstep
	call xyz(s,x,y,z)
   
    count_SEyld=0;count_SE50yld=0
    call kyoukaiS(x1,y1,z1,x,y,z,E,s,Esend)
    
    if(z>0 .and. dsqrt(X**2+Y**2)>Rmax1)then
        goto 10
    elseif(z<WD .and. (dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2))then
        goto 10
    endif

!******************************!    
    call countSE2(x,y,z,E)
!******************************!      
	n1=0; n2=0					! n1:生成したSE数, n2:計算済みSE数
 50	nn=0			! ---------------増殖スタート--------------
 60 E2=E; X2=X; Y2=Y; Z2=Z  ! 電離場所
    
		st12=st1; ct12=ct1; sp12=sp1; cp12=cp1 ! 種電子の進行方向
    CALL DRND(R)
		E1=E*dsqrt(R)				! 散乱電子
 		comg=dsqrt(dsqrt(R));   somg=dsqrt(1d0-comg*comg)
    CALL DRND(R)
		sphi=dsin(R*2d0*pai);  cphi=dcos(R*2d0*pai)	! Eに依存した発生方向
		E=E1 

		x=x2; y=y2; z=z2; st1=st12; ct1=ct12; sp1=sp12; cp1=cp12
	call SEstep
	call xyz(s,x,y,z)
    x1=x; y1=y; z1=z; st11=st1; ct11=ct1; sp11=sp1; cp11=cp1  ! 散乱電子到着場所
    
    call kyoukaiS(x2,y2,z2,x,y,z,E,s,Esend)
    
        x1=x; y1=y; z1=z; st11=st1; ct11=ct1; sp11=sp1; cp11=cp1
    
    if(z>0 .and. dsqrt(X**2+Y**2)>Rmax1)then
        E1=1
        goto 65
    elseif(z<WD .and. (dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2))then
        E1=1
        goto 65
    endif

!******************************!    
    call countSE2(x,y,z,E)
!******************************!  
    
 65 continue
	nn=nn+1						    ! 電離電子
	E=E2*somg*somg
	comg=somg; somg=dsqrt(1d0-comg*comg); sphi=-sphi; cphi=-cphi
	x=x2; y=y2; z=z2; 
    st1=st12; ct1=ct12; sp1=sp12; cp1=cp12
	call SEstep
	call xyz(s,x,y,z)
!	    x3=x; y3=y; z3=z; st13=st1; ct13=ct1; sp13=sp1; cp13=cp1  ! 電離電子到着場所
    
    
    call kyoukaiS(x2,y2,z2,x,y,z,E,s,Esend)
    
    if(z>0 .and. dsqrt(X**2+Y**2)>Rmax1)then
        E=1
        goto 66
    elseif(z<WD .and. (dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2))then
        E=1
        goto 66
    endif

!******************************!    
    call countSE2(x,y,z,E)
!******************************!      
    
66    Es(n1+nn)=E
	  xs(n1+nn)=x; ys(n1+nn)=y; zs(n1+nn)=z 
	  st1s(n1+nn)=st1; ct1s(n1+nn)=ct1
	  sp1s(n1+nn)=sp1; cp1s(n1+nn)=cp1		! 電離電子のデータ収納完了	
!E3=E
    
	E=E1; x=x1; y=y1; z=z1; st1=st11; ct1=ct11
	sp1=sp11; cp1=cp11						! 散乱電子のデータに戻す
    
	if(E < Esend) goto 70		

	goto 60				

 70		n1=n1+nn
80      if(n2 == n1) goto 10
		n2=n2+1					! 電離電子のデータをロード
		E=Es(n2)
		x=xs(n2); y=ys(n2); z=zs(n2)
		st1=st1s(n2); ct1=ct1s(n2); sp1=sp1s(n2); cp1=cp1s(n2)
	if(E < Esend) goto 80		
	goto 50
		
		
	CONTAINS	!------- 内部サブルーチン -------!
	SUBROUTINE SEstep		!--- SEステップ長 ---
    
		!RAMD=5d-10											! 20A
		!IF (E < 25.0) RAMD=10.0**(-2.6*DLOG10(E)+4.3)*1d-10		! A unit
  !      CALL DRND(R)
		!S=-1.0*RAMD*DLOG(R)
		!if(S > 5d-10) S=4d-10
!----------------------------------------------------------!        
        RAMD=15d-10
        
        CALL DRND(R)
		S=-1.0*RAMD*DLOG(R)
        
!----------------------------------------------------------! 
  !      RAMD=5d-10											! 20A
		!IF (E < 25.0) RAMD=10.0**(-2.6*DLOG10(E)+4.3)*0.8d-10		! A unit
  !      !S=ramd
  !     ! print*,E,S;pause
  !      CALL DRND(R)
		!S=-1.0*RAMD*DLOG(R)
  !      !print*,E,S;pause

        
	END SUBROUTINE SEstep
	

	SUBROUTINE EJECT		!--- SE　真空へ放出 ---
	common /SEYeld/ K
		A=-1.0*CT1
		B=dSQRT(V0/E)
		IF (A < B) return  ! 表面エネルギー障壁を超えられない！
		E=E-V0				! 放出エネルギー
		if(E > 50.) then
            ibse=ibse+1		
		else
		    ise=ise+1
        endif
    END subroutine eject

    END
    
    SUBROUTINE countSE1(x,y,z,E)
    use PARAMETER_1
    IMPLICIT DOUBLEPRECISION (A-H,O-Z)
    common/rad2/fogSE(0:100000),fogSE2(0:100000)
    common/count2/SEsample,SElens
    
    rad0=dsqrt(X**2+Y**2)
    
    if(z > 0d0)then
        
    if(E>V0.and.rad0>Rmin1)then
       SEsample=SEsample-1 
       
    do rad=0,300,1
    rad0=dsqrt(X**2+Y**2)
    if(rad0>rad*1d-4 .and. rad0<(rad*1d-4)+1d-4) fogSE(rad)=fogSE(rad)-1
    enddo 
    
    
    endif
    
    elseif(z < WD)then
        
    if(E>V0)then
        SElens=SElens-1
    endif
    
    endif
    

    end
    
    SUBROUTINE countSE2(x,y,z,E)
    use PARAMETER_1
    IMPLICIT DOUBLEPRECISION (A-H,O-Z)
    common/rad2/fogSE(0:100000),fogSE2(0:100000)
    common/count2/SEsample,SElens
    
    rad0=dsqrt(X**2+Y**2)
    
    if(z > 0d0)then
     
    if(E>V0.and.rad0>Rmin1)then
        
        SEsample=SEsample-1
        
    do rad=0,300,1
    rad0=dsqrt(X**2+Y**2)
    if(rad0>rad*1d-4 .and. rad0<(rad*1d-4)+1d-4) fogSE(rad)=fogSE(rad)-1
    enddo   
    
    elseif(E<V0.and.rad0>Rmin1)then
       
        SEsample=SEsample+1 
        
    do rad=0,300,1
    rad0=dsqrt(X**2+Y**2)
    if(rad0>rad*1d-4 .and. rad0<(rad*1d-4)+1d-4) fogSE(rad)=fogSE(rad)+1
    enddo    
        
    endif
    elseif(z < WD)then
        
    if(E>V0)then
        SElens=SElens-1
    elseif(E<V0)then
        SElens=SElens+1
    endif         
    endif
   
    end