    
 
!##################################################################################
!Subrotina para gravar os resultados no arquivo de saida
!##################################################################################
SUBROUTINE OUTPUT(ARQS,NCANT,iDC,XC,YC,Cadj,NNOC,NDOM,NPG,ADP,NPT,DESL,NO)
IMPLICIT NONE

CHARACTER*20:: ARQS
INTEGER::MON,NCANT,i,NNOC,NDOM,NPG,ADP,NPT
INTEGER,DIMENSION(NCANT)::iDC
INTEGER,DIMENSION(NCANT,2)::Cadj
INTEGER  IANO,IMES,IDIA,IHORA,IMIN,ISEC,I100TH
REAL*8,DIMENSION(NCANT)::XC,YC
REAL*8,DIMENSION(3*NPT)::DESL
INTEGER,DIMENSION(NPT)::NO

 WRITE (2,41) NNOC,NDOM,NCANT,NPG,ADP
 41 FORMAT (/,10X,'N DE PONTOS NO CONTORNO =',I3,//,10X,'N DE PONTOS &
 INTERNOS =',I4,//,10X,'N DE CANTOS =',I3,//,10X,'N DE PONTOS DE GAUSS =',I3,//,10X,'INICIAR COM PONTOS ADICIONAIS NO SUPORTE =',&
 I3,///,25X,'GEOMETRIA DA PLACA',/,5X,'CANTO',5X,'COORD X',5X,'COORD Y',5X,'CANTO &
 ANTERIOR',5X,'CANTO POSTERIOR')

DO i=1,NCANT						!impressao da conectividade dos cantos da placa e calculo dos versores dos lados da placa.
 		WRITE(2,42) iDC(i),XC(i),YC(i),Cadj(i,1),Cadj(i,2)
ENDDO
42 FORMAT (2X,i5,5X,F8.3,4X,F8.3,9X,I5,14X,I5)

WRITE(2,150)
150 FORMAT(//,74('*'),/,25X,'DESLOCAMENTOS NODAIS',//,4X,'NO',3X,'ROTACAO 1',3X,'ROTACAO 2',3X,'DESL. TRANV.',&
    3X,'M11',12X,'M22',11X,'M12',10X,'Q1',10X,'Q2')

 DO i=1,NPT
     WRITE(2,201) NO(i),DESL(3*i-2),DESL(3*i-1),DESL(3*i)
 ENDDO

201 FORMAT (2X,I4.2,2X,E12.5,4X,E12.5,4X,E12.5) !,4X,ES9.2E2,4X,ES9.2E2,4X,ES9.2E2,4X,ES9.2E2,4X,ES9.2E2)

CALL GETDAT(IANO,IMES,IDIA)
CALL GETTIM(IHORA,IMIN,ISEC,I100TH)


write(2,500)IDIA,IMES,IANO,IHORA,IMIN,ISEC

500 format(//,5X,74('*'),/,5X,1('*'),52X,1('*'),2X,'DATA: ',I2.2,1H-,I2.2,1H-,I4.4,1X,1H*,/,5X,1H*,10X, 'Danilo &
Hiroshi Konda',21X,1H*,19X,1H*,/,5X,1H*10X,'Prof. Jose Antonio Fontes Santiago',8X, 1H*,2X,'HORA: &
',I2.2,1H:,I2.2,1H:,I2.2,2X,1H*,/,5X,1H*,52X, 1H*,19X,1H*,/,5X,74('*'))
 
    
RETURN
ENDSUBROUTINE OUTPUT
    
    