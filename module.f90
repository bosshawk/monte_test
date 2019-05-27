    MODULE PARAMETER_1                    
    IMPLICIT DOUBLEPRECISION (A-H,O-Z)
        
        
        integer,        PARAMETER :: view=0                        !�O���`��
        integer,        PARAMETER :: Imax=10000                   !�v�Z��
        DOUBLEPRECISION,PARAMETER :: E0=3000                       !���˃G�l���M�[
        DOUBLEPRECISION,PARAMETER :: WD=-5d-3                     !���[�L���O�f�B�X�^���X
        DOUBLEPRECISION,PARAMETER :: Vs=50                          

        
        
        DOUBLEPRECISION,PARAMETER :: PEstop=50.                   
        DOUBLEPRECISION,PARAMETER :: dR=1.0d-3                     !�~�d�ɂ̕�
        
        
        DOUBLEPRECISION,PARAMETER :: Rmax1=100d-3
        DOUBLEPRECISION,PARAMETER :: Rmin1=0 
        DOUBLEPRECISION,PARAMETER :: Rmin2=0                      
        DOUBLEPRECISION,PARAMETER :: Rmax2=11d-3                             
        DOUBLEPRECISION,PARAMETER :: V0=9.45d0

        
                
    END MODULE PARAMETER_1
    
    MODULE PARAMETER_2                    
    IMPLICIT DOUBLEPRECISION (A-H,O-Z)
    
    !Ai=63.5463d-3;   Zi=29.d0;   ROU=8.94d3
    
    DOUBLEPRECISION,PARAMETER :: Ai=63.5463d-3
    DOUBLEPRECISION,PARAMETER :: Zi=29.d0
    DOUBLEPRECISION,PARAMETER :: ROU=8.94d3
    
    
    DOUBLEPRECISION,PARAMETER :: PAI=3.1415927
    DOUBLEPRECISION,PARAMETER :: Avo=6.022D23
    
    END MODULE PARAMETER_2