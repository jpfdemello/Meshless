SUBROUTINE CONDICAO_DE_CONTORNO(ALFA,ALFA2,i,NNOC,C,NCANT,iDLMENOR,NO,NPT,iDC,X,Y,NPSUP,MON,TMON,KGLOBAL,FGLOBAL,VALORP,TIPOPRES,VCANT)
    IMPLICIT NONE

INTEGER:: i,J,K,NNOC,NCANT,NPT,NPSUP,MON,TMON
INTEGER,DIMENSION(NNOC)::C,CA
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(NCANT)::iDC
INTEGER,DIMENSION(2)::iDLMENOR
INTEGER,DIMENSION(NPSUP)::FG
INTEGER,DIMENSION(NNOC,3)::TIPOPRES

REAL*8::UM,ALFA,ALFA2
REAL*8,DIMENSION(1)::XB,YB
REAL*8,DIMENSION(2)::NORMC
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(NCANT,2)::VCANT
REAL*8,DIMENSION(NPSUP)::DG,Fi
REAL*8,DIMENSION(TMON,NPSUP)::AiB,dAiXB,dAiYB,AidBX,AidBY
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY
REAL*8,DIMENSION(3*NPT,3*NPT)::KGLOBAL
REAL*8,DIMENSION(3*NPT)::FGLOBAL
REAL*8,DIMENSION(NNOC,3)::VALORP


    IF(ABS(ALFA).LT.ALFA2)THEN
        
        IF(i.LE.NNOC)THEN
        
            IF(C(i).NE.1)THEN
                NORMC(1)=VCANT(iDLMENOR(1),2)
                NORMC(2)=-VCANT(iDLMENOR(1),1)
            ELSE !c é igual a 1 e o nó é canto
                DO J=1,NCANT
                    IF(iDC(J).EQ.NO(i)) THEN            !AQUI EU POSSO ESCOLHER AS NORMAIS DOS CANTOS
                        NORMC(1)=VCANT(J,2)
                        NORMC(2)=-VCANT(J,1)
                    ENDIF
                ENDDO            
            ENDIF
!if(i.eq.1)then
!    normc(1)=-1.
!    normc(2)=0.
!endif
!if(i.eq.17)then
!    normc(1)=1.
!    normc(2)=0.
!endif

            CALL MONTA_SUPORTE(X(i),Y(i),X,Y,NPT,NPSUP,NO,FG,DG)
        
            CALL MINIMOS_QUADRADOS(i,X(i),Y(i),X,Y,NPSUP,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
        
            XB(1)=X(i)
            YB(1)=Y(i)
        
            CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XB,YB,M,DMX,DMY)

            !obtendo a função fi para o ponto base
            Fi=0.D0
            UM=0.d0
            DO J=1,NPSUP
                DO K=1,TMON
                    !calculo de Fi
                    Fi(J)=M(1,K)*AiB(K,J)+Fi(J)
                ENDDO
            UM=fi(j)+UM
            ENDDO
            UM=UM-1.
            IF(ABS(UM).GT.0.1)THEN
                WRITE(*,*) 'PROBLEMA NA FI DO PONTO',I
                WRITE(*,*) 'AUMENTE OS PONTOS ADICIONAIS NO SUPORTE - cod z'
                PAUSE
            ENDIF

            IF(C(i).EQ.1)THEN   !AQUI MODIFICO A MATRIZ K e F NOS NÓS DE CANTOS
                DO J=1,NNOC
                    IF(NO(J).EQ.CA(i))THEN
                        IF(TIPOPRES(i,1).EQ.0.AND.TIPOPRES(J,1).EQ.1)THEN
                            DO K=1,3*NPT
                                KGLOBAL(3*i-2,K)=0.D0
                            ENDDO
                            DO K=1,NPSUP
                                KGLOBAL(3*i-2,FG(K)*3-2)=Fi(K)*NORMC(1)
                                KGLOBAL(3*i-2,FG(K)*3-1)=Fi(K)*NORMC(2)
                            ENDDO
                            FGLOBAL(3*i-2)=VALORP(J,1)
                        ENDIF
                        IF(TIPOPRES(i,2).EQ.0.AND.TIPOPRES(J,2).EQ.1)THEN
                            DO K=1,3*NPT
                                KGLOBAL(3*i-1,K)=0.D0
                            ENDDO
                            DO K=1,NPSUP
                                KGLOBAL(3*i-1,FG(K)*3-1)=Fi(K)*NORMC(1)
                                KGLOBAL(3*i-1,FG(K)*3-2)=-Fi(K)*NORMC(2)
                            ENDDO
                            FGLOBAL(3*i-1)=VALORP(J,2)
                        ENDIF

                        IF(TIPOPRES(i,3).EQ.0.AND.TIPOPRES(J,3).EQ.1)THEN
                            DO K=1,3*NPT
                                KGLOBAL(3*i,K)=0.D0
                            ENDDO
                            DO K=1,NPSUP
                                KGLOBAL(3*i,FG(K)*3)=Fi(K)
                            ENDDO
                            FGLOBAL(3*i)=VALORP(J,3)
                        ENDIF                    
                    
                    ENDIF
                ENDDO
            ENDIF        
 
            IF(TIPOPRES(i,1).EQ.1)THEN
                !zerar a respectiva linha da matriz KGLOBAL
                DO J=1,3*NPT
			        KGLOBAL(3*i-2,J)=0.d0
                ENDDO
                DO J=1,NPSUP
			        KGLOBAL(3*i-2,FG(J)*3-2)=Fi(J)*NORMC(1)						
			        KGLOBAL(3*i-2,FG(J)*3-1)=Fi(J)*NORMC(2)                 
                ENDDO
                FGLOBAL(3*i-2)=VALORP(i,1)
            ENDIF
            IF(TIPOPRES(i,2).EQ.1)THEN
                !zerar a respectiva linha da matriz KGLOBAL
                DO J=1,3*NPT
			        KGLOBAL(3*i-1,J)=0.d0
                ENDDO
                DO J=1,NPSUP
			        KGLOBAL(3*i-1,FG(J)*3-1)=Fi(J)*NORMC(1)						
			        KGLOBAL(3*i-1,FG(J)*3-2)=-Fi(J)*NORMC(2)                 
                ENDDO
                FGLOBAL(3*i-1)=VALORP(i,2)
            ENDIF
            IF(TIPOPRES(i,3).EQ.1)THEN
                !zerar a respectiva linha da matriz KGLOBAL
                DO J=1,3*NPT
			        KGLOBAL(3*i,J)=0.d0
                ENDDO
                DO J=1,NPSUP
			        KGLOBAL(3*i,FG(J)*3)=Fi(J)						
                ENDDO
                FGLOBAL(3*i)=VALORP(i,3)
            ENDIF
        
        ENDIF
        
    ENDIF

END SUBROUTINE CONDICAO_DE_CONTORNO