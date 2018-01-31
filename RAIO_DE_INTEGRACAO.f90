!################################################################################
!SUBROTINA QUE CALCULA O TAMANHO DO RAIO DE INTEGRAÇÃO DO PONTO i
!################################################################################
SUBROUTINE RAIO_DE_INTEGRACAO(i,NPT,X,Y,RAIOi)
IMPLICIT NONE

INTEGER::i,J,NPT
REAL*8::AUX1,AUX2,RAIOi
REAL*8,DIMENSION(NPT)::X,Y,DiST

RAIOi=1.d308
DO J=1,NPT
    AUX1=(X(i)-X(J))**2
    AUX2=(Y(i)-Y(J))**2
    DiST(J)=DSQRT(AUX1+AUX2)
    IF(DiST(J).GT.0)THEN
        IF(RAIOi.GT.DiST(J))THEN
            RAIOi=DiST(J)
        ENDIF
    ENDIF
ENDDO

ENDSUBROUTINE RAIO_DE_INTEGRACAO
