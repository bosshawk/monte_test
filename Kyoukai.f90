    Subroutine KyoukaiP(xb,yb,zb,x,y,z,E,Ds,Ebend)
    
    use PARAMETER_1
    
    implicit doubleprecision(A-H,O-Z)
    
    COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI  
    common/Iend/Iend
    common/countPE1/count_BSEyld,BSEyld
    if(E < Ebend) return
!------------------------------------------------------------------------------------------
	if( (zb>0d0 .and. z>0d0) .OR. (zb<WD .and. z<WD) ) then 
        
        if(zb>0d0 .and. z>0d0)then
            if(dsqrt(X**2+Y**2)>Rmax1)then
            Iend=Iend+1
            x=xb; y=yb; z=zb
            return   
            endif
        elseif(zb<WD .and. z<WD)then
            if(dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2)then
            Iend=Iend+1
            x=xb; y=yb; z=zb
            return   
            endif
        endif

        if(view==1)CALL kidou(xb,yb,zb,x,y,z)         !write
        CALL ELOSS(DEDS,E); DE=DEDS*Ds                !エネルギー損失
        ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                  !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                        !次のエネルギー
	    return   
!------------------------------------------------------------------------------------------
    elseif(zb>0d0 .and. z<0d0) then                        

	    xb0=xb; yb0=yb; zb0=zb; Ds0=Ds
	    Dsratio=(0.-zb)/(z-zb); Ds=Ds0*Dsratio
	    x0=xb+(x-xb)*Dsratio; y0=yb+(y-yb)*Dsratio; z0=1d-15            
	    x=x0; y=y0; z=z0
       
        if(view==1)CALL kidou_mudage(xb,yb,zb,x,y,z)   !write
	    CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
        
	    ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                   !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                         !次のエネルギー
	    if(E < Ebend) return    
         
!        if(DABS(CT1) < DSQRT(V0/E) )then            
!            
!         CT1=-CT1
!         Dsres=1d0-Dsratio;Ds=Ds0*Dsres
!         X=X0+Ds*ST1*CP1;Y=Y0+Ds*ST1*SP1;Z=Z0+Ds*CT1
!         if(view==1)CALL kidou_mudage(x0,y0,z0,x,y,z)   !write
!         
!         CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
!	    
!         ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
!	     CALL SE(E,DE,x0,y0,z0,x,y,z)                   !２次電子
!         ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
!         E=E-DE
!        
!        else

        if(count_BSEyld==0)BSEyld=BSEyld+1
        count_BSEyld=count_BSEyld+1

        Call Field(E,x,y,z,x1,y1,z1,x2,y2,z2)
           xb=x1;yb=y1;zb=z1;x=x2;y=y2;z=z2
        
        if(zb>0d0 .and. z>0d0)then
            if(dsqrt(X**2+Y**2)>Rmax1)then
            Iend=Iend+1
            x=x0; y=y0; z=z0
            return   
            endif
        elseif(zb<WD .and. z<WD)then
            if(dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2)then
            Iend=Iend+1
            x=x0; y=y0; z=z0
            return   
            endif
        endif
        
        
        if(z<WD) then
            
        xWD=xb+(x-xb)*(WD-zb)/(z-zb); yWD=yb+(y-yb)*(WD-zb)/(z-zb); zWD=WD-1d-15
        xb=xWD;yb=yWD;zb=zWD
        Dsres=1d0-Dsratio
        Ds=Ds0*Dsres
		x=xb+(x0-xb0)*Dsres/Dsratio; y=yb+(y0-yb0)*Dsres/Dsratio; z=zb+(z0-zb0)*Dsres/Dsratio

        if(view==1)CALL kidou_mudage(xb,yb,zb,x,y,z)   !write
        CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
	    
        ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                   !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                         !次のエネルギー
        
        return
        
        elseif(z>0)then
            
        x00=x+(xb-x)*(0-z)/(zb-z); y00=y+(yb-y)*(0-z)/(zb-z); z00=1d-15
        xb=x00;yb=y00;zb=z00
        Dsres=1d0-Dsratio; Ds=Ds0*Dsres
        x=xb+(x0-xb0)*Dsres/Dsratio; y=yb+(y0-yb0)*Dsres/Dsratio; z=zb+(zb0-z0)*Dsres/Dsratio!zb0-z0
        if(view==1)CALL kidou_mudage(xb,yb,zb,x,y,z)   !write
        CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
	    
        ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                   !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                         !次のエネルギー
        return
        
        elseif( z<0d0  .and.  z>WD ) then
            return
        endif
        
        !endif


	    return

!------------------------------------------------------------------------------------------
    elseif( zb<WD .and. z>WD ) then                       
        
	    xb0=xb; yb0=yb; zb0=zb; Ds0=Ds
	    Dsratio=(WD-zb)/(z-zb); Ds=Ds0*Dsratio
	    xWD=xb+(x-xb)*Dsratio; yWD=yb+(y-yb)*Dsratio; zWD=WD-1d-15      
	    x=xWD; y=yWD; z=zWD
        if(view==1)CALL kidou_mudage(xb,yb,zb,x,y,z)   !write
        CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
	    
        ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                   !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                         !次のエネルギー
	    if(E < Ebend) return 

        
!        if(DABS(CT1)<DSQRT(V0/E))then            
!            
!         CT1=-CT1
!         Dsres=1d0-Dsratio;Ds=Ds0*Dsres
!         X=xWD+Ds*ST1*CP1;Y=yWD+Ds*ST1*SP1;Z=zWD+Ds*CT1
!         if(view==1)CALL kidou_mudage(xWD,yWD,zWD,x,y,z)   !write
!         
!         CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
!	    
!         ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
!	     CALL SE(E,DE,xWD,yWD,zWD,x,y,z)                   !２次電子
!         ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
!         E=E-DE
!        
!        else
 
       Call Field(E,x,y,z,x1,y1,z1,x2,y2,z2)
        xb=x1;yb=y1;zb=z1;x=x2;y=y2;z=z2
        
        if(zb>0d0 .and. z>0d0)then
            if(dsqrt(X**2+Y**2)>Rmax1)then
                
            Iend=Iend+1
            x=xWD; y=yWD; z=zWD
            return   
            endif
        elseif(zb<WD .and. z<WD)then
            if(dsqrt(X**2+Y**2)<Rmin2.or.dsqrt(X**2+Y**2)>Rmax2)then
            
            Iend=Iend+1
            x=xWD; y=yWD; z=zWD
            return   
            endif
        endif
        
        if(z<WD) then!負バイアス
            
        xWD2=xb+(x-xb)*(WD-zb)/(z-zb); yWD2=yb+(y-yb)*(WD-zb)/(z-zb); zWD2=WD-1d-15   
        xb=xWD2;yb=yWD2;zb=zWD2
        Dsres=1d0-Dsratio; Ds=Ds0*Dsres
		x=xb+(xWD-xb0)*Dsres/Dsratio; y=yb+(yWD-yb0)*Dsres/Dsratio; z=zb+(zb0-zWD)*Dsres/Dsratio!zb0-zWD
        if(view==1)CALL kidou_mudage(xb,yb,zb,x,y,z)   !write
        CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
	    
        ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                   !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                         !次のエネルギー
        return
        
        elseif(z>0)then
            
        x0=x+(xb-x)*(0-z)/(zb-z); y0=y+(yb-y)*(0-z)/(zb-z); z0=1d-15
        
        xb=x0;yb=y0;zb=z0
        Dsres=1d0-Dsratio; Ds=Ds0*Dsres
		x=xb+(xWD-xb0)*Dsres/Dsratio; y=yb+(yWD-yb0)*Dsres/Dsratio; z=zb+(zWD-zb0)*Dsres/Dsratio
        if(view==1)CALL kidou_mudage(xb,yb,zb,x,y,z)   !write
        CALL ELOSS(DEDS,E); DE=DEDS*Ds                 !エネルギー損失
	    
        ST11=ST1;CT11=CT1;SP11=SP1;CP11=CP1
	    CALL SE(E,DE,xb,yb,zb,x,y,z)                   !２次電子
        ST1=ST11;CT1=CT11;SP1=SP11;CP1=CP11
        E=E-DE                                         !次のエネルギー
        return
        
        elseif( z<0d0  .and.  z>WD ) then
            return
        endif 

        !endif
        
	    return

    endif
!------------------------------------------------------------------------------------------
    
    end


    Subroutine KyoukaiS(xb,yb,zb,x,y,z,E,Ds,Esend)
    use PARAMETER_1
    implicit doubleprecision(A-H,O-Z)
    
    COMMON /NO7/ST1,CT1,SP1,CP1,SOMG,COMG,SPHI,CPHI 
    
    common/countSE3/count_SEyld,SEyld,SE50yld,count_SE50yld

!------------------------------------------------------------------------------------------
	if( (zb>0d0 .and. z>0d0) .OR. (zb<WD .and. z<WD) ) then
        
        if(view==1)CALL SE_kidou(xb,yb,zb,x,y,z)            !write
        return   
!------------------------------------------------------------------------------------------
    elseif(zb>0d0 .and. z<0d0) then                         

	    xb0=xb; yb0=yb; zb0=zb; Ds0=Ds
	    Dsratio=(0.-zb)/(z-zb); Ds=Ds0*Dsratio
	    x0=xb+(x-xb)*Dsratio; y0=yb+(y-yb)*Dsratio; z0=1d-15            
	    x=x0; y=y0; z=z0
        if(view==1)CALL SE_kidou_mudage(xb,yb,zb,x,y,z)   !write
        if((E-Esend) < 0)then
            return
        endif
        if(DABS(CT1) < DSQRT(V0/E) )then            
            
         CT1=-CT1
         Dsres=1d0-Dsratio;Ds=Ds0*Dsres
         X=X0+Ds*ST1*CP1;Y=Y0+Ds*ST1*SP1;Z=Z0+Ds*CT1
         if(view==1)CALL SE_kidou_mudage(x0,y0,z0,x,y,z)   !write
         
        else
            
        if(count_SEyld==0)then
        if(E>50)then
            SE50yld=SE50yld+1
        else
            SEyld=SEyld+1
        endif
        endif
        
        if(CT1<0) then
        Call Field(E,x,y,z,x1,y1,z1,x2,y2,z2)
        endif

        if(z<WD) then
    
        xWD=x1+(x2-x1)*(WD-z1)/(z2-z1); yWD=y1+(y2-y1)*(WD-z1)/(z2-z1); zWD=WD-1d-15
        !xb=xWD;yb=yWD;zb=zWD
        Dsres=1d0-Dsratio;
        Ds=Ds0*Dsres
		x=xWD+(x0-xb0)*Dsres/Dsratio; y=yWD+(y0-yb0)*Dsres/Dsratio; z=zWD+(z0-zb0)*Dsres/Dsratio
        if(view==1)CALL SE_kidou_mudage(xWD,yWD,zWD,x,y,z)   !write
        return
        
        elseif(z>0)then
   
        x00=x1+(x2-x1)*(0-z1)/(z2-z1); y00=y1+(y2-y1)*(0-z1)/(z2-z1); z00=1d-15
        Dsres=1d0-Dsratio; Ds=Ds0*Dsres
		x=x00+(x0-xb0)*Dsres/Dsratio; y=y00+(y0-yb0)*Dsres/Dsratio; z=z00+(zb0-z0)*Dsres/Dsratio!zb0-z0
        if(view==1)CALL SE_kidou_mudage(x00,y00,z00,x,y,z)   !write
        return
        
        elseif( z<0d0  .and.  z>WD ) then
            return
        endif
        
        endif

        return

!------------------------------------------------------------------------------------------
    elseif( zb<WD .and. z>WD ) then                        
        
	    xb0=xb; yb0=yb; zb0=zb; Ds0=Ds
	    Dsratio=(WD-zb)/(z-zb); Ds=Ds0*Dsratio
	    xWD=xb+(x-xb)*Dsratio; yWD=yb+(y-yb)*Dsratio; zWD=WD-1d-15      
	    x=xWD; y=yWD; z=zWD
        if(view==1)CALL SE_kidou_mudage(xb,yb,zb,x,y,z)   !write
        if((E-Esend) < 0)then
            return
        endif

        if(DABS(CT1)<DSQRT(V0/E))then            
            
         CT1=-CT1
         Dsres=1d0-Dsratio;Ds=Ds0*Dsres
         X=xWD+Ds*ST1*CP1;Y=yWD+Ds*ST1*SP1;Z=zWD+Ds*CT1
         if(view==1)CALL SE_kidou_mudage(xWD,yWD,zWD,x,y,z)   !write
         
         else

        count_SEyld=count_SEyld+1
        count_SE50yld=count_SE50yld+1

        if(CT1>0) then
        Call Field(E,x,y,z,x1,y1,z1,x2,y2,z2)
        endif

        if(z<WD) then!負バイアス
            
        xWD2=x1+(x2-x1)*(WD-z1)/(z2-z1); yWD2=y1+(y2-y1)*(WD-z1)/(z2-z1); zWD2=WD-1d-15   
        Dsres=1d0-Dsratio; Ds=Ds0*Dsres
		x=xWD2+(xWD-xb0)*Dsres/Dsratio; y=yWD2+(yWD-yb0)*Dsres/Dsratio; z=zWD2+(zb0-zWD)*Dsres/Dsratio!zb0-zWD
        if(view==1)CALL SE_kidou_mudage(xWD2,yWD2,zWD2,x,y,z)   !write
        return
        
        elseif(z>0)then
            
        x0=x2+(x1-x2)*(0-z2)/(z1-z2); y0=y2+(y1-y2)*(0-z2)/(z1-z2); z0=1d-15
        Dsres=1d0-Dsratio; Ds=Ds0*Dsres
		x=x0+(xWD-xb0)*Dsres/Dsratio; y=y0+(yWD-yb0)*Dsres/Dsratio; z=z0+(zWD-zb0)*Dsres/Dsratio
        if(view==1)CALL SE_kidou_mudage(x0,y0,z0,x,y,z)   !write
        return
        
        elseif( z<0d0  .and.  z>WD ) then
            return
        endif
        
        endif
        
        return

    endif
!------------------------------------------------------------------------------------------
    
    end
    
Subroutine kidou(xb,yb,zb,x,y,z)
implicit doubleprecision(A-H,O-Z)

do kkk=1,100; x10=xb+kkk*(x-xb)/100;y10=yb+kkk*(y-yb)/100;z10=zb+kkk*(z-zb)/100;write(2,*) x10,y10,z10;enddo

endsubroutine

Subroutine kidou_mudage(xb,yb,zb,x,y,z)
implicit doubleprecision(A-H,O-Z)

do kkk=1,10; x10=xb+kkk*(x-xb)/10;y10=yb+kkk*(y-yb)/10;z10=zb+kkk*(z-zb)/10;write(2,*) x10,y10,z10;enddo

endsubroutine


   
Subroutine SE_kidou(xb,yb,zb,x,y,z)
implicit doubleprecision(A-H,O-Z)

do kkk=1,50; x10=xb+kkk*(x-xb)/50;y10=yb+kkk*(y-yb)/50;z10=zb+kkk*(z-zb)/50;write(3,*) x10,y10,z10;enddo

endsubroutine

Subroutine SE_kidou_mudage(xb,yb,zb,x,y,z)
implicit doubleprecision(A-H,O-Z)

do kkk=1,10; x10=xb+kkk*(x-xb)/10;y10=yb+kkk*(y-yb)/10;z10=zb+kkk*(z-zb)/10;write(3,*) x10,y10,z10;enddo

endsubroutine
