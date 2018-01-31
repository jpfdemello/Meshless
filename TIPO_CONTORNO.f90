!##################################################################################################
!esta subrotina monta as matrizes que contem componente normais ao contorno
!##################################################################################################
    
SUBROUTINE TIPO_CONTORNO(TIPOSUP,NNOC,NPT,NO,iDLMENOR,NCANT,NPRES,NOPRES,VPRES,TIPOPRES,iDC,iDCadj,iDCMENOR,&
    CPRES1,CPRES2,VSUPPRES1,VSUPPRES2,VALORP)
IMPLICIT NONE


INTEGER::TIPOSUP,i,J,K,NPT,NCANT,NPRES,NNOC,iDCMENOR
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(2)::iDLMENOR
INTEGER,DIMENSION(3)::CPRES1,CPRES2
INTEGER,DIMENSION(NCANT)::iDC
INTEGER,DIMENSION(NCANT,2)::iDCadj
INTEGER,DIMENSION(NNOC,3)::TIPOPRES
INTEGER,DIMENSION(NPRES)::NOPRES
REAL*8,DIMENSION(3)::VSUPPRES1,VSUPPRES2
REAL*8,DIMENSION(NNOC,3)::VALORP
REAL*8,DIMENSION(NPRES,3)::VPRES

!correlacionando o valor prescrito com o respectivo nó

!somente valores prescritos diferentes de zero precisam ser declarados
VALORP=0.D0

!por enquanto não vou utilizar valorp para todos os pontos, pois estou utilizando a prescrição do no do canto para toda a lateral da placa
DO i=1,NNOC
    DO J=1,NPRES
        IF (NO(i).EQ.NOPRES(J)) THEN
			!correlacionando a prescrição com os nós do contorno
            VALORP(i,1)=VPRES(J,1)
			VALORP(i,2)=VPRES(J,2)
			VALORP(i,3)=VPRES(J,3)
		ENDIF
	ENDDO
ENDDO


!se não tiver contorno associado, o valor de cpres é 2
CPRES1=2
CPRES2=2



IF(TIPOSUP.EQ.1)THEN
    RETURN
ELSEIF(TIPOSUP.EQ.2)THEN
    DO J=1,NNOC
		!encontrando o nó de canto que irá determinar a condição de contorno
        IF (iDC(iDLMENOR(1)).EQ.NO(J)) THEN			 						
			DO K=1,3
				CPRES1(K)=TIPOPRES(J,K)
				VSUPPRES1(K)=VALORP(J,K)
			ENDDO
		ENDIF
       ENDDO
ELSE        
			!definindo a condição de contorno para o trecho do contorno envolvido pelo subdominio local
    DO J=1,NNOC
	    IF (iDC(iDCadj(iDCMENOR,1)).EQ.NO(J)) THEN			 						!AQUI ACHEI QUAL NO DO CONTORNO VAI DETERMINAR A CONDICAO DE CONTORNO	 .
		    DO K=1,3
		    	CPRES1(K)=TIPOPRES(J,K)
		    	VSUPPRES1(K)=VALORP(J,K)
		    ENDDO
    	ENDIF
	    IF(iDC(iDCMENOR).EQ.NO(J)) THEN									 !AQUI ACHEI QUAL NO DO CONTORNO VAI DETERMINAR A CONDICAO DE CONTORNO
		    DO K=1,3
			    CPRES2(K)=TIPOPRES(J,K)
			    VSUPPRES2(K)=VALORP(J,K)
	    	ENDDO
	    ENDIF
    ENDDO                   
ENDIF
    
RETURN
ENDSUBROUTINE TIPO_CONTORNO