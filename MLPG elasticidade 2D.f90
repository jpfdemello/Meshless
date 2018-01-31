!  MLPGversao20.f90 
!
!  FUNCTIONS:
!  MLPGversao20 - Entry point of console application.
!
!****************************************************************************
!
!  PROGRAM: MLPGversao20
!
!  PURPOSE:  Entry point for the console application.
!   
!****************************************************************************

PROGRAM MLPGversao20
USE ALOCAVEIS
IMPLICIT NONE
!não esquece de corrigir o problema dos angulos

CHARACTER*20::ARQS
CHARACTER::TRANS
INTEGER:: I,J,K,L,NNOC,NDOM,MAXSUP,NCANT,MON,NPG,NPRES,ADP,NPT,TIPOSUP,PSUPG,TMON,iDCMENOR,NPSUP,INFO
REAL*8::  E,NI,T,CARGA,RAIOi,TETA2,TETA1,XT,YT,D,JAC,TETAi,TETAF,DPLINTAUX,TETADPLAUX,TETA3,TETA4
REAL*8::  X1,Y1,X2,Y2,X3,Y3,X4,Y4,XV1,YV1,XV2,YV2,N1,N2,UM,ALFA,teste
INTEGER,DIMENSION(3)::CPRES1,CPRES2,CPRES
INTEGER,DIMENSION(2)::IDLMENOR
INTEGER,DIMENSION(:),ALLOCATABLE::iDC,FG,IPIV !C,NO,CA,NOPRES,
INTEGER,DIMENSION(:,:),ALLOCATABLE::iDCadj,Cadj !TIPOPRES,
REAL*8,DIMENSION(:),ALLOCATABLE::XC,YC,XG,WG,DG,FGLOBAL,XGLOBAL,Fi,DESL,dFidx,dFidy !X,Y,
REAL*8,DIMENSION(:,:),ALLOCATABLE::VCANT,KGLOBAL,M,DMX,DMY,VALORP,Stress !VPRES,
!variáveis calculadas na rotina minimos quadrados
REAL*8,DIMENSION(:,:),ALLOCATABLE::AiB,dAiXB,dAiYB,AidBX,AidBY
!matrizes constantes utilizadas na integração
REAL*8,DIMENSION(2,2)::MATRIZD2,MATRIZN4,MATRIZN5
REAL*8,DIMENSION(3,3)::MATRIZD1,MATRIZN3
REAL*8,DIMENSION(3)::MATRIZR1,MATRIZR2,VSUPPRES1,VSUPPRES2,MATRIZD5,VSUPPRES
REAL*8,DIMENSION(2)::DPLINT,TETADPL,NORMC
REAL*8,DIMENSION(3,2)::MATRIZD3,MATRIZD4
REAL*8,DIMENSION(1)::XB,YB
REAL*8,PARAMETER::Pi=4.d0*datan(1.d0) !3.14159265358979323846d0
Real*8,Parameter::alfa2=1.e-5


CALL INPUT(NNOC,NDOM,MAXSUP,NCANT,MON,NPG,E,NI,T,CARGA,NPRES,ADP,NPT,PSUPG,TMON,D,ARQS,NPSUP,ALFA) !,X,Y,C,NO,CA,NOPRES,TIPOPRES,VPRES,

!subrotina input verificada - ok!

CALL OPEN_OUTPUT(ARQS,MON)

ALLOCATE(iDC(NCANT),XC(NCANT),YC(NCANT),VCANT(NCANT,2),iDCadj(NCANT,2),VALORP(NNOC,3),Cadj(NCANT,2))

CALL CANTOS(NCANT,NNOC,C,NO,NPT,XC,YC,X,Y,VCANT,iDC,CA,iDCadj,Cadj)
!vcant,iDCadj,Cadj ok

ALLOCATE(XG(NPG),WG(NPG),KGLOBAL(3*NPT,3*NPT),FGLOBAL(3*NPT))

    
CALL Calcula_Pontos_Gauss(NPG,XG,WG)

CALL MATRIZES_FIXAS(D,Ni,MATRIZD1,T,CARGA,MATRIZR1,MATRIZR2,MATRIZD2,MATRIZD3,MATRIZD4,MATRIZD5)

!iniciando a matriz KGLOBAL e FGLOBAL
KGLOBAL=0.D0
FGLOBAL=0.D0


DO i=1,NPT

    ALLOCATE(FG(PSUPG),DG(PSUPG),AiB(TMON,PSUPG),dAiXB(TMON,PSUPG),dAiYB(TMON,PSUPG),AidBX(TMON,PSUPG),AidBY(TMON,PSUPG))
    
    CALL RAIO_DE_INTEGRACAO(i,NPT,X,Y,RAIOi)

raioi=1.d0*raioi

    CALL CLASSIFICA_PONTO(i,NCANT,NPT,X,Y,XC,YC,VCANT,RAIOi,iDCadj,NNOC,C,TIPOSUP,TETA1,TETA2,TETA3,TETA4,DPLINT,TETADPL,&
        X1,Y1,X2,Y2,X3,Y3,X4,Y4,iDCMENOR,iDLMENOR)
     
    CALL TIPO_CONTORNO(TIPOSUP,NNOC,NPT,NO,iDLMENOR,NCANT,NPRES,NOPRES,VPRES,TIPOPRES,iDC,iDCadj,iDCMENOR,&
    CPRES1,CPRES2,VSUPPRES1,VSUPPRES2,VALORP)
 
    DO J=1,NPG
        DO K=1,NPG
            
            !a integral de domínio de teta2 à teta1 é feita para todos os TIPOSUP
            INFO=1
            CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)

            CALL T_COORD_CIRCULO(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,XT,YT,RAIOi,WG(J),WG(K),JAC)
              
            CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)

            CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                
            CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)
      
            IF(TIPOSUP.EQ.2)THEN

                !integrar o triangulo
                !para tiposup2 a integral é feita de teta1 à teta2
                
                INFO=2
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
            
                DPLINTAUX=DPLINT(1)
                TETADPLAUX=TETADPL(1)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
                
                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
           
                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)
               
                !integrais de contorno apenas com o laço J 
                IF(K.EQ.1)THEN
                
                    XV1=X1
                    YV1=Y1
                    XV2=X2
                    YV2=Y2
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
                
                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                    
                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES1(L)
			            VSUPPRES(L)=VSUPPRES1(L)
		            ENDDO
                    
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAiXB,dAiYB,AidBX,AidBY,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,&
                        MATRIZD5,ALFA,MATRIZD1)
                    
                ENDIF
                
            ELSEIF(TIPOSUP.EQ.3) THEN

                !integrar o primeiro triangulo
                !Para tiposup3 o primeiro triangulo é integrado de teta1 à teta3

                INFO=3
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
                
                DPLINTAUX=DPLINT(1)
                TETADPLAUX=TETADPL(1)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)

                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
           
                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)                

                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X1
                    YV1=Y1
                    XV2=X3
                    YV2=Y3
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
                
                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                    
                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES1(L)
			            VSUPPRES(L)=VSUPPRES1(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAiXB,dAiYB,AidBX,AidBY,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,&
                        MATRIZD5,ALFA,MATRIZD1)
                
                ENDIF

                !integrar o segundo triangulo
                !Para tiposup3 o segundo triangulo é integrado de teta3 à teta2
                INFO=31
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
                
                DPLINTAUX=DPLINT(2)
                TETADPLAUX=TETADPL(2)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
                
                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
           
                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)                
          
                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X3
                    YV1=Y3
                    XV2=X2
                    YV2=Y2
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
                
                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                    
                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES2(L)
			            VSUPPRES(L)=VSUPPRES2(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAiXB,dAiYB,AidBX,AidBY,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,&
                        MATRIZD5,ALFA,MATRIZD1)
                                
                ENDIF
        
            ELSEIF(TIPOSUP.EQ.4)THEN
       
                !integrar o setor circular
                !a integral é feita de teta3 à teta4
                INFO=4
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
                
                CALL T_COORD_CIRCULO(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,XT,YT,RAIOi,WG(J),WG(K),JAC)
                
                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)

                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                
                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)
                             
                !integrar o primeiro triangulo
                !Para tiposup4 o primeiro triangulo é integrado de teta1 à teta3                
                INFO=41
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
                
                DPLINTAUX=DPLINT(1)
                TETADPLAUX=TETADPL(1)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
                
                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
           
                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)                

                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X1
                    YV1=Y1
                    XV2=X3
                    YV2=Y3
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
                
                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                    
                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES1(L)
			            VSUPPRES(L)=VSUPPRES1(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAiXB,dAiYB,AidBX,AidBY,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,&
                        MATRIZD5,ALFA,MATRIZD1)
                                          
                ENDIF

                !integrar o segundo triangulo
                !Para tiposup4 o segundo triangulo é integrado de teta4 à teta2                
                INFO=42
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
    
                DPLINTAUX=DPLINT(2)
                TETADPLAUX=TETADPL(2)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
                
                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
           
                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
                    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)    

                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X4
                    YV1=Y4
                    XV2=X2
                    YV2=Y2
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
                
                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)
                    
                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES2(L)
			            VSUPPRES(L)=VSUPPRES2(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAiXB,dAiYB,AidBX,AidBY,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,&
                        MATRIZD5,ALFA,MATRIZD1)
                                          
                ENDIF
            ENDIF
        ENDDO
    ENDDO    
    
    CALL CONDICAO_DE_CONTORNO(ALFA,ALFA2,i,NNOC,C,NCANT,iDLMENOR,NO,NPT,iDC,X,Y,NPSUP,MON,TMON,KGLOBAL,FGLOBAL,VALORP,TIPOPRES,VCANT)
    !aplicar o MMQM no nós do contorno para  obter as funções de forma para impor as condicoes de contorno
    
    DEALLOCATE(FG,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)

    !aplicar condições de contorno (somente para nós do contorno)

    CALL EVO(i,NPT)
ENDDO
 
!CALL IMPOSICAO_CC(TIPOPRES,NNOC,NPT,KGLOBAL,FGLOBAL)

ALLOCATE(XGLOBAL(3*NPT),IPIV(3*NPT))
XGLOBAL=0.D0

J=2 !J=1 usa o solver do lapack... j=2 utiliza o meu solver ou o svd
IF(J.EQ.1)THEN
    CALL DGETRF(3*NPT,3*NPT,KGLOBAL,3*NPT,IPIV,INFO)
    IF(INFO.EQ.0)THEN
        CONTINUE
    ELSE
        WRITE(*,*) 'PROBLEMA EM DGETRF'
        PAUSE
    ENDIF
    TRANS='N'
    XGLOBAL=FGLOBAL
    CALL DGETRS(TRANS,3*NPT,1,KGLOBAL,3*NPT,IPIV,XGLOBAL,3*NPT,INFO)
    IF(INFO.EQ.0)THEN
        CONTINUE
    ELSE
        WRITE(*,*) 'PROBLEMA EM DGETRS'
        PAUSE
    ENDIF
ELSEIF(J.EQ.2)THEN
    i=999 
    CALL r8mat_fs (3*NPT,KGLOBAL,FGLOBAL,XGLOBAL,i)
!    i=998
    IF(i.EQ.998)THEN
       
        CALL SVDSOLVE(kglobal,fglobal,3*npt,1,XGLOBAL)
      

    ENDIF
ENDIF

!interpolar os valores ficticios para obter os valores reais
ALLOCATE(DESL(3*NPT),STRESS(3,3*NPT))
ALLOCATE(FG(NPSUP),DG(NPSUP),AiB(TMON,NPSUP),dAiXB(TMON,NPSUP),dAiYB(TMON,NPSUP),AidBX(TMON,NPSUP),AidBY(TMON,NPSUP))
ALLOCATE(Fi(NPSUP),M(1,TMON),DMX(1,TMON),DMY(1,TMON),dFidx(NPSUP),dFidy(NPSUP))
DESL=0.D0
Stress=0.d0 
DO i=1,NPT
     
    CALL MONTA_SUPORTE(X(i),Y(i),X,Y,NPT,NPSUP,NO,FG,DG)                
    
    !aplicar o MMQM no nós do contorno para  obter as funções de forma para impor as condicoes de contorno
    CALL MINIMOS_QUADRADOS(i,X(i),Y(i),X,Y,NPSUP,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)

    XB(1)=X(i)
    YB(1)=Y(i)
    CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XB,YB,M,DMX,DMY)

    !obtendo a função fi para o ponto base
    Fi=0.D0
    dfidx=0.d0
    dfidy=0.d0
    UM=0.D0
    DO J=1,NPSUP
        DO K=1,TMON
            !calculo de Fi
            Fi(J)=M(1,K)*AiB(K,J)+Fi(J)
            !Cálculo de dPhidx e dPhidy
            !1º Termo
            dFidx(j)=dMx(1,k)*AiB(k,j)+dFidx(j)
            dFidy(j)=dMy(1,k)*Aib(k,j)+dFidy(j)
            !2º Termo
            dFidx(j)=M(1,k)*dAiXB(k,j)+dFidx(j)
            dFidy(j)=M(1,k)*dAiYB(k,j)+dFidy(j)
            !3º Termo
            dFidx(j)=M(1,k)*AidBx(k,j)+dFidx(j)
            dFidy(j)=M(1,k)*AidBy(k,j)+dFidy(j)
        ENDDO
        UM=UM+FI(J)
    ENDDO    
    UM=UM-1.
    IF(ABS(UM).GT.0.1)THEN
        WRITE(*,*) 'PROBLEMA NA FI DO PONTO',I
        WRITE(*,*) 'AUMENTE OS PONTOS ADICIONAIS NO SUPORTE - cod s'
        PAUSE
    ENDIF

    DO J=1,NPSUP
        DESL(i*3-2)=Fi(J)*XGLOBAL(FG(J)*3-2)+DESL(i*3-2)
        DESL(i*3-1)=Fi(J)*XGLOBAL(FG(J)*3-1)+DESL(i*3-1)
        DESL(i*3)=Fi(J)*XGLOBAL(FG(J)*3)+DESL(i*3)
!        write(*,*) fi(j),xglobal(fg(j)*3-2)
    ENDDO

    do j=1,NPSUP
        Stress(1,3*i-2)=Stress(1,3*i-2)+dFidx(j)*xGlobal(Fg(j)*3-2)
        Stress(2,3*i-1)=Stress(2,3*i-1)+dFidy(j)*xGlobal(Fg(j)*3-1)
        Stress(3,3*i)=Stress(3,3*i)+dFidx(j)*xGlobal(Fg(j)*3-2)+dFidy(j)*xGlobal(Fg(j)*3-1)
    end do
    
ENDDO

CALL OUTPUT(ARQS,NCANT,iDC,XC,YC,Cadj,NNOC,NDOM,NPG,ADP,NPT,DESL,NO)

do k=1,3*NPT
    write(3,"(3x,3E15.5)") (Stress(j,k), j=1,3)
end do

write(*,*) 'chegou ao fim'
read(*,*)



END PROGRAM MLPGversao20