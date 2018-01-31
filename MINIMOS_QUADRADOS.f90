!##########################################################################################
!SUBROTINA PARA MONTAR A FUNÇÃO DE FORMA PARA O PONTO BASE XT,YT
!##########################################################################################    
SUBROUTINE  MINIMOS_QUADRADOS(JJ,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,AiB,dAiXB,dAiYB,AidBX,AidBY)  
IMPLICIT NONE

INTEGER::i,J,K,NPT,PSUPG,MON,TMON,JJ
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::XT,YT,RSG
REAL*8,DIMENSION(PSUPG)::XB,YB,DG,W,DWX,DWY
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(PSUPG,TMON)::M,DMX,DMY
REAL*8,DIMENSION(TMON,PSUPG)::B,dBx,dBy,AiB,dAiXB,dAiYB,AidBX,AidBY
REAL*8,DIMENSION(TMON,TMON)::A,Ai,IDENT,dAx,dAy,dAiX,dAiY,AUX1,AUX2

DO i=1,PSUPG
    XB(i)=X(FG(i))
    YB(i)=Y(FG(i))
ENDDO

!montando a matriz P(armazenada em M) para o suporte do ponto XT,YT
CALL MONTA_BASE_MONOMIAL(MON,TMON,PSUPG,XB,YB,M,DMX,DMY)

!calculo do raio da função peso para aplicação do MMQM no suporte do ponto de Gauss
RSG=1.1*DG(PSUPG)
!montando a matriz W (spline de 4 ordem) e suas derivadas
!WRITE(*,*) RSG
DO i=1,PSUPG
    W(i)=1.D0-6.D0*((DG(i)/RSG)**2.)+8.D0*((DG(i)/RSG)**3.)-3.D0*((DG(i)/RSG)**4.)
    DWX(i)=(12.D0*(XT-X(FG(i)))/RSG**2.)*(-1.D0+2.D0*((DG(i))/RSG)-((DG(i))/RSG)**2.)
    DWY(i)=(12.D0*(YT-Y(FG(i)))/RSG**2.)*(-1.D0+2.D0*((DG(i))/RSG)-((DG(i))/RSG)**2.)
ENDDO

!construção da matriz B e suas derivadas
DO i=1,TMON
    DO J=1,PSUPG
        B(i,J)=M(J,i)*W(J)
        dBx(i,J)=M(J,i)*DWX(J)
        dBy(i,J)=M(J,i)*DWY(J)
    ENDDO
ENDDO
!construção da matriz A e suas derivadas
A=0.
dAx=0.
dAy=0.
A=MATMUL(B,M)
dAx=MATMUL(dBX,M)!+MATMUL(DMX,B)
dAy=MATMUL(dBy,M)!+MATMUL(DMY,B)

!determinando a matriz A inversa
DO i=1,TMON
    IDENT=0.
    IDENT(i,i)=1.
    CAll r8mat_fs(TMON,A,IDENT(:,i),Ai(:,i),JJ)
ENDDO
!calculo das derivadas da matriz A inversa
AUX1=0.
AUX2=0.
AUX1=MATMUL(dAx,Ai)
AUX2=MATMUL(dAy,Ai)
dAiX=0.
dAiY=0.
dAiX=MATMUL(-Ai,AUX1)
dAiY=MATMUL(-Ai,AUX2)
!calculo dos parametros que serão utilizados na integração
AiB=0.
dAiXB=0.
dAiYB=0.
AidBX=0.
AidBY=0.

AiB=MATMUL(Ai,B)
dAiXB=MATMUL(dAiX,B)
dAiYB=MATMUL(dAiY,B)
AidBX=MATMUL(Ai,dBX)
AidBY=MATMUL(Ai,dbY)



RETURN
END SUBROUTINE MINIMOS_QUADRADOS