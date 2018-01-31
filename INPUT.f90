SUBROUTINE INPUT(NNOC,NDOM,MAXSUP,NCANT,MON,NPG,E,NI,T,CARGA,NPRES,ADP,NPT,PSUPG,TMON,D,ARQS,NPSUP,ALFA)  
USE ALOCAVEIS !X,Y,C,NO,CA,NOPRES,TIPOPRES,VPRES
IMPLICIT NONE

CHARACTER*20:: ARQUENT,ARQS,Kglobal
LOGICAL:: EXISTE
INTEGER:: I,J,NNOC,NDOM,MAXSUP,NCANT,MON,NPG,NPRES,ADP,NPT,PSUPG,TMON,NPSUP
REAL*8::E,NI,T,CARGA,D,ALFA

!------------------------------------------------------------------------------------
! 			   IMPRESSAO DO TITULO DO PROGRAMA
!------------------------------------------------------------------------------------

1 FORMAT (10X,'PROGRAMA : ANALISE DE PLACA ESPESSA COM O MLPG',//,20X,'ORIENTADOR:JOSE &
ANTONIO FONTES SANTIAGO',/,20X,'ORIENTADOR:JOSE CLAUDIO DE FARIA TELLES',//,20X,'ALUNO: DANILO &
HIROSHI KONDA',//,' DESENVOLVIDO EM 04/2017 - PEC/COPPE/UFRJ',/,' AUTOR: DANILO HIROSHI KONDA',//)
WRITE (*,1)
!---------------------------------------------------------------------------------------------------------------------------------------------
!	ENTRADA E LEITURA DOS NOMES DOS ARQUIVOS DE ENTRADA E SAIDA DE DADOS
!----------------------------------------------------------------------------------------------------------------------------------------------

!WRITE (*,'(//A\)')' DIGITE O NOME DO ARQUIVO DE ENTRADA DE DADOS : '
!READ (*,'(A)') ARQUENT
!WRITE (*,'(//A\)')' DIGITE O NOME DO ARQUIVO DE SAIDA DE DADOS : '
!READ (*,'(A)') ARQS
ARQUENT='anal2121.txt'
ARQS='s2121.txt'
Kglobal='Kglobal5.txt'
INQUIRE (FILE=ARQUENT,EXIST=EXISTE)
IF (.NOT.EXISTE) THEN
    WRITE (*,'(/3A)') ' ARQUIVO DE DADOS',ARQUENT,'NAO ENCONTRADO !!'
	WRITE (*,'(/3A)') ' PROGRAMA TERMINADO!!'
	READ(*,*)
	STOP 
ENDIF

OPEN (1,FILE=ARQUENT)
OPEN (2,FILE=ARQS,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
Open(3,File=Kglobal)
DO i =1, 9
    READ(1,*)
END DO
20 FORMAT (A100)                                        

READ(1,*)  NNOC,NDOM,MAXSUP,NCANT,MON  

DO i =1, 4
    READ(1,*)
END DO

READ(1,*)  NPG,E,NI,T,CARGA,NPRES,ADP  

DO i =1, 4
    READ(1,*)
END DO

READ(1,*)  PSUPG,ALFA


DO i =1, 5
    READ(1,*)
END DO

NPT=NNOC+NDOM
TMON=(MON+1)*(MON+2)/2
ALLOCATE(NO(NPT),X(NPT),Y(NPT),C(NNOC),CA(NNOC),TIPOPRES(NNOC,3))

J=0
C=0
DO i=1,NNOC
	READ(1,*) NO(I),X(I),Y(I),C(I),CA(I),TIPOPRES(i,1),TIPOPRES(i,2),TIPOPRES(i,3)
ENDDO

DO i =1, 4
    READ(1,*)
END DO

ALLOCATE(NOPRES(NPRES),VPRES(NPRES,3))

DO i=1,NPRES
	READ(1,*) NOPRES(i),VPRES(i,1),VPRES(i,2),VPRES(i,3)
ENDDO

DO i =1, 3
    READ(1,*)
END DO

DO I=NNOC+1,NPT
	READ(1,*) NO(I),X(I),Y(I)
    write(*,*) i
ENDDO
!calculo do modulo de rigidez para elasticidade 2D 
D=E/(1.-(NI**2.d0))

NPSUP=TMON+ADP

    ENDSUBROUTINE INPUT