!##################################################################################################
!esta subrotina monta as matrizes que contem componente normais ao contorno
!##################################################################################################
    
SUBROUTINE MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
IMPLICIT NONE


REAL*8::xv1,yv1,xv2,yv2,S1,S2,AUX1,AUX2,N1,N2
REAL*8,DIMENSION(3,3)::MATRIZN3
REAL*8,DIMENSION(2,2)::MATRIZN4,MATRIZN5

!calculo das direções tangenciais
AUX1=(XV2-XV1)**2.d0
AUX2=(YV2-YV1)**2.d0
S1=(XV2-XV1)/(DSQRT((AUX1)+(AUX2)))
S2=(YV2-YV1)/(DSQRT((AUX1)+(AUX2)))

!obtendo as normais
N1=S2
N2=-S1

!montagem das matrizes N3,N4 e N5
MATRIZN3=0.d0
MATRIZN4=0.d0
MATRIZN5=0.d0
MATRIZN3(1,1)=S2
MATRIZN3(1,2)=S1
MATRIZN3(2,1)=-S1
MATRIZN3(2,2)=S2
MATRIZN3(3,3)=1.d0
MATRIZN4(1,1)=S2
MATRIZN4(1,2)=-S1
MATRIZN5(2,1)=S1
MATRIZN5(2,2)=S2

RETURN
ENDSUBROUTINE MATRIZES_NORMAIS
