!#######################################################################################
!Esta subrotina faz a transformação de coordenadas para os pontos de gauss num dominio
!circular
!#######################################################################################
SUBROUTINE T_COORD_CIRCULO(X,Y,QSI,ETA,TETA2,TETA1,XT,YT,RAIOi,PQSI,PETA,JAC)
IMPLICIT NONE

REAL*8::QSI,ETA,TETA2,TETA1,XT,YT,X,Y,RAIOi,JAC,PQSI,PETA
    
!calculo das coordenadas dos pontos de Gauss no sistema cartesiano
	XT=X+0.5d0*RAIOi*(1.d0+QSI)*DCOS(TETA2*0.5d0*(1.d0-ETA)+TETA1*0.5d0*(1.d0+ETA))                !QSI=XG(i)   ETA=XG(j)
	YT=Y+0.5d0*RAIOi*(1.d0+QSI)*DSIN(TETA2*0.5d0*(1.d0-ETA)+TETA1*0.5d0*(1.d0+ETA)) 

!calculo do jacobiano
    JAC=DABS(0.5D0*RAIOi*(1.D0+QSI))*DABS(0.25D0*RAIOi*(TETA1-TETA2))*PQSI*PETA
    
    
RETURN
END SUBROUTINE T_COORD_CIRCULO