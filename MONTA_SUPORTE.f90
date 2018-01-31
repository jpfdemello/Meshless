!##################################################################################
!Esta subrotina monta o suporte para o ponto de gauss
!##################################################################################    
    
SUBROUTINE MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
USE TIPO
IMPLICIT NONE

INTEGER::i,NPT,PSUPG
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::XT,YT,AUX1,AUX2
REAL*8,DIMENSION(NPT)::D,X,Y
REAL*8,DIMENSION(PSUPG)::DG
TYPE(GROUP),DIMENSION(NPT)::DORD
DO i=1,NPT
    AUX1=(XT-X(i))**2
    AUX2=(YT-Y(i))**2
    !calcula a distância do ponto de Gauss até todos os pontos do problema.
    D(i)=DSQRT(AUX1+AUX2)
    DORD(i)%VALOR=D(i)
    DORD(i)%ORDEM=NO(i)
    DORD(i)%ORDEM2=i
ENDDO

CALL QuickSort(DORD,NPT)

DO i=1,PSUPG
    FG(i)=DORD(i)%ORDEM2
    DG(i)=DORD(i)%VALOR
ENDDO
   

END SUBROUTINE MONTA_SUPORTE
    
