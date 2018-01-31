!###################################################################
!ESTA SUBROTINA CALCULA OS CANTOS ADJACENTES A UM CANTO E TAMBÉM
!OS VERSORES DOS LADOS DA PLACA (SENTIDO ANTI-HORÁRIO POSITIVO)
!###################################################################

SUBROUTINE CANTOS(NCANT,NNOC,C,NO,NPT,XC,YC,X,Y,VCANT,iDC,CA,iDCadj,Cadj)
IMPLICIT NONE

INTEGER::i,J,NCANT,NNOC,NPT
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(NNOC)::C,CA
INTEGER,DIMENSION(NCANT)::iDC
INTEGER,DIMENSION(NCANT,2)::iDCadj
REAL*8::AUX1,AUX2,AUX3
REAL*8,DIMENSION(NCANT)::XC,YC
REAL*8,DIMENSION(NPT)::X,Y
INTEGER,DIMENSION(NCANT,2)::Cadj
REAL*8, DIMENSION(NCANT,2)::Xadj,Yadj
REAL*8,DIMENSION(NCANT,2)::VCANT

J=0
DO i=1,NNOC
	IF (C(I).EQ.1) THEN
		!contador para identificar a quantidade de cantos
        J=J+1
		IF(J.GT.NCANT)THEN
			WRITE(*,*) 'ERRO NO LANCAMETO DOS CANTOS DA PLACA'
			READ(*,*)
        ENDIF
        !identificando se o Nó(i) é nó de canto
		iDC(J)=NO(I)
   	    !salvando as coordenadas dos nós de canto
        XC(J)=X(i)
		YC(J)=Y(i)
        !Salvando o canto anterior como canto adjacente 1.
		Cadj(J,1)=CA(i)
	ENDIF
ENDDO

DO i=1,NCANT		 !encontrando o canto posterior a um canto
	DO J=1,NCANT
		IF(Cadj(i,1).EQ.iDC(J)) THEN
            !salvando o canto posterior como canto adjacente 2
			Cadj(J,2)=iDC(i)
			!salvando a identidade dos cantos adjacentes
            iDCadj(i,1)=J
			iDCadj(J,2)=i
		ENDIF
	ENDDO																  
ENDDO
DO i=1,NCANT		 
		DO J=1,NCANT
			IF(iDC(i).EQ.Cadj(J,1)) THEN
				!identificando e salvando as coordenadas dos cantos adjacentes 1
                Xadj(J,1)=XC(i)
				Yadj(J,1)=YC(i)
            ELSE IF (iDC(i).EQ.Cadj(J,2)) THEN
                !identificando e salvando as coordenadas dos cantos adjacentes 2
				Xadj(J,2)=XC(i)
				Yadj(J,2)=YC(i)
			ENDIF			
        ENDDO
ENDDO

DO i=1,NCANT						!impressao da conectividade dos cantos da placa e calculo dos versores dos lados da placa.
    AUX1=(Xadj(i,2)-XC(i))**2
    AUX2=((Yadj(i,2))-YC(i))**2
    AUX3=DSQRT(AUX1+AUX2)
    !calculando os versores dos lados da placa. Orientação positiva anti-horária
	VCANT(i,1)=(Xadj(i,2)-XC(i))/AUX3
	VCANT(i,2)=	((Yadj(i,2))-YC(i))/AUX3
ENDDO

RETURN
ENDSUBROUTINE CANTOS    
