!##################################################################################################
!esta subrotina faz a transfomação de coordenadas para integração no contorno
!e calcula o jacobiano já multiplicado pelo pelo peso do ponto de gauss
!##################################################################################################
    
SUBROUTINE T_COORD_CONTORNO(QSI,XV1,YV1,XV2,YV2,PQSI,XT,YT,JAC)
IMPLICIT NONE

REAL*8::QSI,XV1,YV1,XV2,YV2,PQSI,XT,YT,JAC

XT=0.5d0*(XV1+XV2+QSI*(XV2-XV1))	
YT=0.5d0*(YV1+YV2+QSI*(YV2-YV1))	
JAC=PQSI*0.5d0*DSQRT(((XV2-XV1)**2.d0)+(YV2-YV1)**2.d0)

RETURN
ENDSUBROUTINE T_COORD_CONTORNO