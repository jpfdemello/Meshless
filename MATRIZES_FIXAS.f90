!#################################################################################
!esta subrotina monta as matrizes que não variam na aplicação do método
!#################################################################################
SUBROUTINE MATRIZES_FIXAS(D,Ni,MATRIZD1,T,CARGA,MATRIZR1,MATRIZR2,MATRIZD2,MATRIZD3,MATRIZD4,MATRIZD5)
IMPLICIT NONE

REAL*8::D,Ni,T,CARGA,LAMB
REAL*8,DIMENSION(2,2)::MATRIZD2
REAL*8,DIMENSION(3,3)::MATRIZD1
REAL*8,DIMENSION(3)::MATRIZR1,MATRIZR2,MATRIZD5
REAL*8,DIMENSION(3,2)::MATRIZD3,MATRIZD4

!montagem da matriz D1
MATRIZD1=0.d0
MATRIZD1(1,1)=D
MATRIZD1(2,2)=D
MATRIZD1(1,2)=D*NI
MATRIZD1(2,1)=MATRIZD1(1,2)
MATRIZD1(3,3)=D*(1.d0-NI)/2.d0

!calculo da constante caracteristica da teoria de Reissner - Lambda
LAMB=(DSQRT(10.D0))/T

!montagem da matriz R1
MATRIZR1=0.D0
MATRIZR1(1)=Ni*CARGA/((1.D0-Ni)*LAMB**2)
MATRIZR1(2)=Ni*CARGA/((1.D0-Ni)*LAMB**2)

!montagem da matriz R1
MATRIZR2=0.D0
MATRIZR2(1)=Ni*CARGA/((1.D0-Ni)*LAMB**2)

!montagem da matriz D2
MATRIZD2=0.D0
MATRIZD2(1,1)=(D*(1.d0-NI)*LAMB**2)/2.d0
MATRIZD2(2,2)=(D*(1.d0-NI)*LAMB**2)/2.d0

!montagem da matriz D3
MATRIZD3=0.D0
MATRIZD3(1,1)=D
MATRIZD3(1,2)=D*Ni

!montagem da matriz D4
MATRIZD4=0.D0
MATRIZD4(2,1)=D*(1.d0-NI)/2.d0
MATRIZD4(2,2)=D*(1.d0-NI)/2.d0

!montagem da matriz D5
MATRIZD5=0.D0
MATRIZD5(3)=(D*(1.d0-NI)*LAMB**2)/2.d0

RETURN
END SUBROUTINE MATRIZES_FIXAS