!###########################################################################
!esta subrotina faz a transfomação de coordenadas para o triangulo 
!e calcula o jacobiano já multiplicado pelo pelo peso do ponto de gauss
!###########################################################################
SUBROUTINE T_COORD_TRI(X,Y,QSI,ETA,TETA2,TETA1,PQSI,PETA,DPLINT,TETADPL,XT,YT,JAC)
IMPLICIT NONE

REAL*8::X,Y,QSI,ETA,TETA2,TETA1,XT,YT,PQSI,PETA,JAC,DPLINT,TETADPL,AUX1,AUX2

real*8:: a1,a2,a3,a4


XT=X+(DPLINT/(2*DCOS(TETADPL-0.5D0*TETA1*(1-ETA)-0.5D0*TETA2*(1+ETA))))*(1+QSI)*DCOS(0.5D0*TETA1*(1-ETA)+0.5D0*TETA2*(1+ETA))
YT=Y+(DPLINT/(2*DCOS(TETADPL-0.5D0*TETA1*(1-ETA)-0.5D0*TETA2*(1+ETA))))*(1+QSI)*DSIN(0.5D0*TETA1*(1-ETA)+0.5D0*TETA2*(1+ETA))
a1=(2*DCOS(TETADPL-0.5D0*TETA1*(1-ETA)-0.5D0*TETA2*(1+ETA)))
a2=(TETADPL-0.5D0*TETA1*(1-ETA)-0.5D0*TETA2*(1+ETA))


!calculo do jacobiano
AUX1=DABS((DPLINT/(2.D0*DCOS(TETADPL-0.5D0*TETA1*(1.D0-ETA)-0.5D0*TETA2*(1.D0+ETA))))*(1.D0+QSI))
AUX2=DABS((DPLINT/(4.D0*DCOS(TETADPL-0.5D0*TETA1*(1.D0-ETA)-0.5D0*TETA2*(1.D0+ETA))))*(TETA2-TETA1))
JAC=AUX1*AUX2*PQSI*PETA

RETURN
END SUBROUTINE T_COORD_TRI