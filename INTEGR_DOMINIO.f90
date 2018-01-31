    
!#############################################################################################
!Esta subrotina faz a integração de domínio para um ponto de gauss
!#############################################################################################
SUBROUTINE INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAiXB,dAiYB,AidBX,AidBY,TMON,&
    PSUPG,MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL)
IMPLICIT NONE

INTEGER::i,J,K,L,NPT,TMON,PSUPG,MON,N
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::RAIOi,XT,YT,AUX1,AUX2,Dist,W,DWX,DWY,CARGA,JAC,UM
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(3,3)::MATRIZW1,MATRIZD1,AUX33
REAL*8,DIMENSION(TMON,PSUPG)::AiB,dAiXB,dAiYB,AidBX,AidBY
REAL*8,DIMENSION(PSUPG)::dFix,dFiy,Fi
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY
REAL*8,DIMENSION(1)::XBM,YBM
REAL*8,DIMENSION(2,3*PSUPG)::MATRIZFI2
REAL*8,DIMENSION(3,3*PSUPG)::MATRIZFi1,INT133N,INT633N
REAL*8,DIMENSION(3)::MATRIZR1
REAL*8,DIMENSION(3,2)::MATRIZW2,AUX32
REAL*8,DIMENSION(2,2)::MATRIZD2
REAL*8,DIMENSION(3*NPT,3*NPT)::KGLOBAL
REAL*8,DIMENSION(3*NPT)::FGLOBAL

!calcular a distancia do ponto de gauss ao ponto base
AUX1=(X(i)-XT)**2.
AUX2=(Y(i)-YT)**2.
Dist=DSQRT(AUX1+AUX2)
!avaliar a função peso(spline de 4 ordem) no ponto de gauss (função peso centrada no ponto base)
W=1.-6.*((Dist/RAIOi)**2)+8.*((Dist/RAIOi)**3)-3.*((Dist/RAIOi)**4.)
DWX=-(12.D0*(X(i)-XT)/RAIOi**2.)*(-1.D0+2.D0*((Dist)/RAIOi)-((Dist)/RAIOi)**2.)
DWY=-(12.D0*(Y(i)-YT)/RAIOi**2.)*(-1.D0+2.D0*((Dist)/RAIOi)-((Dist)/RAIOi)**2.)

!montar a matriz w1(3x3) para o ponto de gauss
MATRIZW1=0.
MATRIZW1(1,1)=DWX
MATRIZW1(1,3)=DWY
MATRIZW1(2,2)=DWY
MATRIZW1(2,3)=DWX

!construindo a base monomial para o ponto de gauss
XBM(1)=XT
YBM(1)=YT
N=1

CALL MONTA_BASE_MONOMIAL(MON,TMON,N,XBM,YBM,M,DMX,DMY)

!obtendo as derivadas da função de forma fi e do próprio fi
dFix=0.D0
dFiy=0.D0
Fi=0.D0
UM=0.D0
DO J=1,PSUPG
    DO K=1,TMON
        !dFi é composto por uma soma de tres termos
        !calculo do primeiro termo
        dFiX(J)=DMX(1,K)*AiB(K,J)+dFiX(J)
        dFiY(J)=DMY(1,K)*AiB(K,J)+dFiY(J)
        !calculo do segundo termo
        dFiX(J)=M(1,K)*dAiXB(K,J)+dFiX(J)
        dFiY(J)=M(1,K)*dAiYB(K,J)+dFiY(J)
        !calculo do terceiro termo
        dFiX(J)=M(1,K)*AidBX(K,J)+dFiX(J)
        dFiY(J)=M(1,K)*AidBY(K,J)+dFiY(J)
        !calculo de Fi
        Fi(J)=M(1,K)*AiB(K,J)+Fi(J)
    ENDDO
    UM=UM+Fi(J)
ENDDO
UM=UM-1.D0

IF(ABS(UM).GT.0.1)THEN
    write(*,*) um
    WRITE(*,*) 'PROBLEMA NA FI DO PONTO',I
    WRITE(*,*) 'AUMENTE OS PONTOS ADICIONAIS NO SUPORTE - cod x'
    PAUSE
ENDIF

!montar a matriz Fi1 (3x3*psupg) 
MATRIZFI1=0.D0
DO J=1,3*PSUPG,3                   
    MATRIZFI1(1,J)=dFiX((J+2)/3)
    MATRIZFI1(2,J+1)=dFiY((J+2)/3)
    MATRIZFI1(3,J)=dFiY((J+2)/3)
    MATRIZFI1(3,J+1)=dFiX((J+2)/3)
ENDDO
!!a primeira e a segunda integral da eq.20 já podem ser resolvidas aqui
!
!!montagem da matriz W2
!MATRIZW2=0.D0
!MATRIZW2(1,1)=W
!MATRIZW2(2,2)=W
!MATRIZW2(3,1)=DWX
!MATRIZW2(3,2)=DWY
!
!!montar a matriz Fi2 (2x3*psupg)
!MATRIZFI2=0.D0
!DO J=1,3*PSUPG,3                
!    MATRIZFI2(1,J)=Fi((J+2)/3)
!    MATRIZFI2(2,J+1)=Fi((J+2)/3)
!    MATRIZFI2(1,J+2)=dFiX((J+2)/3)
!    MATRIZFI2(2,J+2)=dFiY((J+2)/3)
!ENDDO
!
!!a sexta integral da eq.20 já pode ser resolvida aqui. A nona integral também!
!
!!RESOLVENDO AS INTEGRAIS
!! primeira integral
!AUX33=0.D0
!produto da matriz W1 pela matriz D1
AUX33=MATMUL(MATRIZW1,MATRIZD1)

!produto de AUX33(MATRIZW1*MATRIZD1) pela matriz Fi1
INT133N=0.D0
INT133N=MATMUL(AUX33,MATRIZFi1)
!multiplicando o resultado da integral no ponto de gauss pelo jacobiano da transformada e pelo peso do ponto de gauss
INT133N=JAC*INT133N


!salvar o resultado da integral 1 (int133n) na matriz K global
DO J=1,PSUPG
    DO K=1,3
        !primeira linha relativa ao nó i da matriz K global
        KGLOBAL(i*3-2,FG(J)*3-3+K)=INT133N(1,J*3-3+K)+KGLOBAL(i*3-2,FG(J)*3-3+K)
        !segunda linha relativa ao nó i da matriz K global
        KGLOBAL(i*3-1,FG(J)*3-3+K)=INT133N(2,J*3-3+K)+KGLOBAL(i*3-1,FG(J)*3-3+K)        
        !terceira linha relativa ao nó i da matriz K global
        KGLOBAL(i*3,FG(J)*3-3+K)=INT133N(3,J*3-3+K)+KGLOBAL(i*3,FG(J)*3-3+K)       
    ENDDO
ENDDO

!!segunda integral
!FGLOBAL(3*i-2)=-MATRIZW1(1,1)*MATRIZR1(1)*JAC+FGLOBAL(3*i-2)
!FGLOBAL(3*i-1)=-MATRIZW1(2,2)*MATRIZR1(2)*JAC+FGLOBAL(3*i-1)
!
!!sexta integral
!AUX32=0.D0
!!produto de matriz W2 por matriz D2
!AUX32=MATMUL(MATRIZW2,MATRIZD2)
!!produto de AUX32(MATRIZW2*MATRIZD2) pela matriz Fi2
!INT633N=0.D0
!INT633N=MATMUL(AUX32,MATRIZFi2)
!!multiplicando o resultado da integral no ponto de gauss pelo jacobiano da transformada e pelo peso do ponto de gauss
!INT633N=JAC*INT633N
!!salvar o resultado da integral 6 (int633n) na matriz K global 
!DO J=1,PSUPG  !tem erro aqui
!    DO K=1,3
!        !primeira linha relativa ao nó i da matriz K global
!        KGLOBAL(i*3-2,FG(J)*3-3+K)=INT633N(1,J*3-3+K)+KGLOBAL(i*3-2,FG(J)*3-3+K)
!        !segunda linha relativa ao nó i da matriz K global
!        KGLOBAL(i*3-1,FG(J)*3-3+K)=INT633N(2,J*3-3+K)+KGLOBAL(i*3-1,FG(J)*3-3+K)        
!        !terceira linha relativa ao nó i da matriz K global
!        KGLOBAL(i*3,FG(J)*3-3+K)=INT633N(3,J*3-3+K)+KGLOBAL(i*3,FG(J)*3-3+K)       
!    ENDDO
!ENDDO
!
!!nona integral
!FGLOBAL(3*i)=W*CARGA*JAC+FGLOBAL(3*i)
!

END SUBROUTINE INTEGR_DOMINIO

    