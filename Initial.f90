! �����̏������AMott�f�ʐ�(�G�l���M�[�����A�f�ʐϓǂݍ���)
     SUBROUTINE INITIAL
     
     use PARAMETER_2
     
     IMPLICIT DOUBLEPRECISION(A-H,O-Z)
	 
     COMMON /NO4/ENGY(0:201),DIELA(0:201,0:90)

!   �����̏�����
     CALL INTRND

!   �G�l���M�[�̏�����
     DO I=1,201
        ENGY(I)=10D0**((DBLE(I)-1D0)/40D0)
     END DO

!   Cu�̍ޗ��萔��`�E�f�ʐϓǂݍ���
     CALL MAKEELAS;  PRINT *, '  MAKE_ELAS  OK'

     END
