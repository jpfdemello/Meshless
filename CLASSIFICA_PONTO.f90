!##########################################################################################################
!Esta subrotina classifica os pontos do problema para definir a área de integração e calcula os parametros geometricos da
!area de integração. Aqui também são determinadas as condições de contorno para integração no contorno, quando o subdominio
!avança sobre o contorno do problema (Ainda que, as condições de contorno sejam definidas em cada ponto do contorno,
!para evitar algumas dificuldades, cada lado da placa receberá a mesma prescrição do canto em que o lado se inicia)
!##########################################################################################################
SUBROUTINE CLASSIFICA_PONTO(JJ,NCANT,NPT,X,Y,XC,YC,VCANT,RAIOi,iDCadj,NNOC,C,TIPOSUP,TETA1,TETA2,TETA3,TETA4,DPLINT,TETADPL,&
    X1,Y1,X2,Y2,X3,Y3,X4,Y4,iDCMENOR,iDLMENOR)
IMPLICIT NONE

INTEGER::JJ,i,NCANT,NPT,IDCMENOR,TIPOSUP,NNOC
INTEGER,DIMENSION(2)::iDLMENOR
INTEGER,DIMENSION(NCANT,2)::iDCadj
INTEGER,DIMENSION(NNOC)::C
REAL*8::AUX1,AUX2,MENORDPC,RAIOi,TETA1,TETA2,X1,Y1,X2,Y2,X3,Y3,TETA3,X4,Y4,TETA4
REAL*8,DIMENSION(2)::MENORDPL,D1,D2,N1,N2,N3,N4,DPLINT,TETADPL
REAL*8,DIMENSION(NCANT):: DPC,XC,YC,DPL
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(NCANT,2)::VCANT
REAL*8,PARAMETER::Pi=3.14159265358979323846d0

DPLINT=0.D0
TETADPL=0.D0

DO i=1,NCANT
    AUX1=(X(JJ)-XC(i))**2
    AUX2=(Y(JJ)-YC(i))**2
!distancia do ponto base aos cantos da placa
    DPC(i)=DSQRT(AUX1+AUX2)	 
ENDDO

MENORDPC=DPC(1)
iDCMENOR=1

DO i=2,NCANT
    IF(DPC(i).LT.MENORDPC) THEN
        !calculo da distancia do ponto base ao canto mais próximo
		MENORDPC=DPC(i)
		!identificador do canto mais próximo do ponto base
        iDCMENOR=i
	ENDIF
ENDDO

DO i=1,NCANT
    !calculo da distancia do ponto base aos lados da placa
    DPL(i)=DABS((X(JJ)-XC(i))*VCANT(i,2)-(Y(JJ)-YC(i))*VCANT(i,1))
ENDDO

MENORDPL(1)=1.d308
iDLMENOR(1)=1

DO i=1,NCANT
	!salvar a distancia e a identificação dos dois lados mais próximos ao ponto base 
    IF(DPL(i).LT.MENORDPL(1)) THEN
		MENORDPL(2)=MENORDPL(1)			
		MENORDPL(1)=DPL(i)
		iDLMENOR(2)= iDLMENOR(1)
		iDLMENOR(1)=i
	ELSE IF(DPL(i).LT.MENORDPL(2)) THEN
		MENORDPL(2)=DPL(i)
		iDLMENOR(2)=i							
	ENDIF
ENDDO

IF(RAIOi.LT.MENORDPC.AND.RAIOi.LE.MENORDPL(1)) THEN
    !o suporte de integração não toca o contorno e a integração será feita em um circulo completo
    TIPOSUP=1 
	TETA2=0.d0
	TETA1=2.d0*PI
ELSE IF(RAIOi.LE.MENORDPC.AND.RAIOi.GT.MENORDPL(1).AND.RAIOi.LE.MENORDPL(2)) THEN
    TIPOSUP=2
	D1(1)=(DSQRT((RAIOi**2)-(MENORDPL(1)**2)))*VCANT(iDLMENOR(1),1)
	D1(2)=(DSQRT((RAIOi**2)-(MENORDPL(1)**2)))*VCANT(iDLMENOR(1),2)
	D2(1)=(DSQRT((X(JJ)-XC(iDLMENOR(1)))**2+(Y(JJ)-YC(iDLMENOR(1)))**2-(MENORDPL(1))**2))*VCANT(iDLMENOR(1),1)
	D2(2)=(DSQRT((X(JJ)-XC(iDLMENOR(1)))**2+(Y(JJ)-YC(iDLMENOR(1)))**2-(MENORDPL(1))**2))*VCANT(iDLMENOR(1),2)
	X1=XC(iDLMENOR(1))+D2(1)-D1(1)					 					  !estas são as coordenadas para integração no triangulo
	Y1=YC(iDLMENOR(1))+D2(2)-D1(2)
	X2=XC(iDLMENOR(1))+D2(1)+D1(1)
	Y2=YC(iDLMENOR(1))+D2(2)+D1(2)
	N1(1)=(X1-X(JJ))/(DSQRT((X1-X(JJ))**2+(Y1-Y(JJ))**2))
	N1(2)=(Y1-Y(JJ))/(DSQRT((X1-X(JJ))**2+(Y1-Y(JJ))**2))
	N2(1)=(X2-X(JJ))/(DSQRT((X2-X(JJ))**2+(Y2-Y(JJ))**2))
	N2(2)=(Y2-Y(JJ))/(DSQRT((X2-X(JJ))**2+(Y2-Y(JJ))**2))

    IF(N1(2).GE.0) THEN
		TETA1=DACOS(N1(1))
	ELSE 
		TETA1=2.d0*(PI)-DACOS(N1(1))
    ENDIF
	
    IF(N2(2).GE.0) THEN										
		TETA2=DACOS(N2(1))									
	ELSE 
		TETA2=2.d0*(PI)-DACOS(N2(1))					   
    ENDIF

    !para tiposup=2, dplint(1) é a distancia do ponto ao lado mais próximo
    DPLINT(1)=MENORDPL(1)
    !calculando o angulo do dplint
    IF(-VCANT(iDLMENOR(1),2).GE.0.)THEN
        TETADPL(1)=DACOS(VCANT(iDLMENOR(1),2))
    ELSE
        TETADPL(1)=2.D0*Pi-DACOS(VCANT(iDLMENOR(1),2))
    ENDIF
    
ELSE IF(RAIOi.GE.MENORDPC.AND.RAIOi.GT.MENORDPL(2)) THEN
	TIPOSUP=3	
    IF(iDLMENOR(1).EQ.iDCadj(iDCMENOR,1)) THEN
		AUX1=MENORDPL(1)
		AUX2=MENORDPL(2)
	ELSE
		AUX1=MENORDPL(2)
		AUX2=MENORDPL(1)
    ENDIF    

    D1(1)=(DSQRT((RAIOi**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),1)
	D1(2)=(DSQRT((RAIOi**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),2)
	D2(1)=(DSQRT(((X(JJ)-XC(iDCadj(iDCMENOR,1)))**2)+((Y(JJ)-YC(iDCadj(iDCMENOR,1)))**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),1)
   	D2(2)=(DSQRT(((X(JJ)-XC(iDCadj(iDCMENOR,1)))**2)+((Y(JJ)-YC(iDCadj(iDCMENOR,1)))**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),2)
	X1=XC(iDCadj(iDCMENOR,1))+D2(1)-D1(1)
	Y1=YC(iDCadj(iDCMENOR,1))+D2(2)-D1(2)
	N1(1)=(X1-X(JJ))/(DSQRT((X1-X(JJ))**2+(Y1-Y(JJ))**2))
	N1(2)=(Y1-Y(JJ))/(DSQRT((X1-X(JJ))**2+(Y1-Y(JJ))**2))    
 
    AUX1=((X(JJ)-XC(iDCMENOR))**2)+((Y(JJ)-YC(iDCMENOR))**2)
    
 	D1(1)=(DSQRT(AUX1-AUX2**2))*VCANT(iDCMENOR,1)
	D1(2)=(DSQRT(((X(JJ)-XC(iDCMENOR))**2)+((Y(JJ)-YC(iDCMENOR))**2)-AUX2**2))*VCANT(iDCMENOR,2)
	D2(1)=(DSQRT((RAIOi**2)-AUX2**2))*VCANT(iDCMENOR,1)
	D2(2)=(DSQRT((RAIOi**2)-AUX2**2))*VCANT(iDCMENOR,2)
	X2=XC(iDCMENOR)+D1(1)+D2(1)
	Y2=YC(iDCMENOR)+D1(2)+D2(2)
	N2(1)=(X2-X(JJ))/(DSQRT((X2-X(JJ))**2+(Y2-Y(JJ))**2))
	N2(2)=(Y2-Y(JJ))/(DSQRT((X2-X(JJ))**2+(Y2-Y(JJ))**2))
	X3=XC(iDCMENOR)
	Y3=YC(iDCMENOR)   

    IF(JJ.GT.NNOC)THEN
		N3(1)=(X3-X(JJ))/(DSQRT((X3-X(JJ))**2+(Y3-Y(JJ))**2))		
		N3(2)=(Y3-Y(JJ))/(DSQRT((X3-X(JJ))**2+(Y3-Y(JJ))**2))            
    ELSE IF(C(JJ).EQ.0)THEN    
		N3(1)=(X3-X(JJ))/(DSQRT((X3-X(JJ))**2+(Y3-Y(JJ))**2))		
		N3(2)=(Y3-Y(JJ))/(DSQRT((X3-X(JJ))**2+(Y3-Y(JJ))**2))            
    ELSE IF(C(JJ).EQ.1)THEN    
		N3(1)=0. !ESTABELECI UM VALOR QUALQUQER POIS SE O NÓ FOI DE CANTO, NÃO FAZ SENTIDO A DIREÇÃO N3
		N3(2)=0.
    ENDIF

    IF(N1(2).GE.0) THEN
		TETA1=DACOS(N1(1))
	ELSE 
		TETA1=2.d0*(PI)-DACOS(N1(1))
    ENDIF
	
    IF(N2(2).GE.0) THEN
		TETA2=DACOS(N2(1))
	ELSE 
		TETA2=2.d0*(PI)-DACOS(N2(1))
    ENDIF
	
    IF(N3(2).GE.0) THEN
	    TETA3=DACOS(N3(1))
	ELSE 
		TETA3=2.d0*(PI)-DACOS(N3(1))
    ENDIF
 
   !para tiposup=3, dplint(1) é a distancia do ponto ao lado anterior e dplint(2) referente ao lado posterior
    IF(IDLMENOR(1).EQ.IDCMENOR)THEN
        DPLINT(2)=MENORDPL(1)
        DPLINT(1)=MENORDPL(2)
		IF(-VCANT(iDLMENOR(2),1).GE.0.)THEN						  
			TETADPL(1)=DACOS(VCANT(iDLMENOR(2),2))				  
		ELSE                                                    
			TETADPL(1)=2.D0*PI-DACOS(VCANT(iDLMENOR(2),2))		  
		ENDIF													
		IF(-VCANT(iDLMENOR(1),1).GE.0.)THEN						 
			TETADPL(2)=DACOS(VCANT(iDLMENOR(1),2))				 
		ELSE														
			TETADPL(2)=2.D0*PI-DACOS(VCANT(iDLMENOR(1),2))		  
		ENDIF	        
     ELSEIF(IDLMENOR(2).EQ.IDCMENOR)THEN
        DPLINT(2)=MENORDPL(2)
        DPLINT(1)=MENORDPL(1)
        IF(-VCANT(iDLMENOR(1),1).GE.0.)THEN						  
			TETADPL(1)=DACOS(VCANT(iDLMENOR(1),2))				  
		ELSE													
			TETADPL(1)=2.D0*PI-DACOS(VCANT(iDLMENOR(1),2))		
		ENDIF													
		IF(-VCANT(iDLMENOR(2),1).GE.0.)THEN						
			TETADPL(2)=DACOS(VCANT(iDLMENOR(2),2))				
		ELSE													
			TETADPL(2)=2.D0*PI-DACOS(VCANT(iDLMENOR(2),2))		
		ENDIF		
    ENDIF
    
    
ELSE IF(RAIOi.LT.MENORDPC.AND.RAIOi.GT.MENORDPL(2)) THEN
    TIPOSUP=4

    IF(iDLMENOR(1).EQ.iDCadj(iDCMENOR,1)) THEN
		AUX1=MENORDPL(1)
		AUX2=MENORDPL(2)
	ELSE
		AUX1=MENORDPL(2)
		AUX2=MENORDPL(1)
    ENDIF
    
	D1(1)=(DSQRT((RAIOi**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),1)
	D1(2)=(DSQRT((RAIOi**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),2)
	D2(1)=(DSQRT(((X(JJ)-XC(iDCadj(iDCMENOR,1)))**2)+((Y(JJ)-YC(iDCadj(iDCMENOR,1)))**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),1)
   	D2(2)=(DSQRT(((X(JJ)-XC(iDCadj(iDCMENOR,1)))**2)+((Y(JJ)-YC(iDCadj(iDCMENOR,1)))**2)-AUX1**2))*VCANT(iDCadj(iDCMENOR,1),2)
	X1=XC(iDCadj(iDCMENOR,1))+D2(1)-D1(1)
	Y1=YC(iDCadj(iDCMENOR,1))+D2(2)-D1(2)
	N1(1)=(X1-X(JJ))/(DSQRT((X1-X(JJ))**2+(Y1-Y(JJ))**2))
	N1(2)=(Y1-Y(JJ))/(DSQRT((X1-X(JJ))**2+(Y1-Y(JJ))**2))
	X3=XC(iDCadj(iDCMENOR,1))+D2(1)+D1(1)
	Y3=YC(iDCadj(iDCMENOR,1))+D2(2)+D1(2)
	N3(1)=(X3-X(JJ))/(DSQRT((X3-X(JJ))**2+(Y3-Y(JJ))**2))
	N3(2)=(Y3-Y(JJ))/(DSQRT((X3-X(JJ))**2+(Y3-Y(JJ))**2))    

    D1(1)=(DSQRT(((X(JJ)-XC(iDCMENOR))**2)+((Y(JJ)-YC(iDCMENOR))**2)-AUX2**2))*VCANT(iDCMENOR,1)
	D1(2)=(DSQRT(((X(JJ)-XC(iDCMENOR))**2)+((Y(JJ)-YC(iDCMENOR))**2)-AUX2**2))*VCANT(iDCMENOR,2)
	D2(1)=(DSQRT((RAIOi**2)-AUX2**2))*VCANT(iDCMENOR,1)
	D2(2)=(DSQRT((RAIOi**2)-AUX2**2))*VCANT(iDCMENOR,2)
	X2=XC(iDCMENOR)+D1(1)+D2(1)
	Y2=YC(iDCMENOR)+D1(2)+D2(2)
	N2(1)=(X2-X(JJ))/(DSQRT((X2-X(JJ))**2+(Y2-Y(JJ))**2))
	N2(2)=(Y2-Y(JJ))/(DSQRT((X2-X(JJ))**2+(Y2-Y(JJ))**2))
	X4=XC(iDCMENOR)+D1(1)-D2(1)
	Y4=YC(iDCMENOR)+D1(2)-D2(2)
	N4(1)=(X4-X(JJ))/(DSQRT((X4-X(JJ))**2+(Y4-Y(JJ))**2))
	N4(2)=(Y4-Y(JJ))/(DSQRT((X4-X(JJ))**2+(Y4-Y(JJ))**2))

	IF(N1(2).GE.0) THEN
		TETA1=DACOS(N1(1))
	ELSE 
		TETA1=2.d0*(PI)-DACOS(N1(1))
	ENDIF
	IF(N2(2).GE.0) THEN
		TETA2=DACOS(N2(1))
	ELSE 
		TETA2=2.d0*(PI)-DACOS(N2(1))
	ENDIF
	IF(N3(2).GE.0) THEN
		TETA3=DACOS(N3(1))
	ELSE 
		TETA3=2.d0*(PI)-DACOS(N3(1))
	ENDIF
	IF(N4(2).GE.0) THEN
		TETA4=DACOS(N4(1))
	ELSE 
		TETA4=2.d0*(PI)-DACOS(N4(1))
    ENDIF 
    
   !para tiposup=3, dplint(1) é a distancia do ponto ao lado anterior e dplint(2) referente ao lado posterior
    IF(IDLMENOR(1).EQ.IDCMENOR)THEN
        DPLINT(2)=MENORDPL(1)
        DPLINT(1)=MENORDPL(2)
		IF(-VCANT(iDLMENOR(2),1).GE.0.)THEN						 
			TETADPL(1)=DACOS(VCANT(iDLMENOR(2),2))				
		ELSE													
			TETADPL(1)=2.D0*PI-DACOS(VCANT(iDLMENOR(2),2))		  
		ENDIF														
		IF(-VCANT(iDLMENOR(1),1).GE.0.)THEN						 
			TETADPL(2)=DACOS(VCANT(iDLMENOR(1),2))				
		ELSE														  
			TETADPL(2)=2.D0*PI-DACOS(VCANT(iDLMENOR(1),2))		 
		ENDIF	        
    ELSEIF(IDLMENOR(2).EQ.IDCMENOR)THEN
        DPLINT(2)=MENORDPL(2)
        DPLINT(1)=MENORDPL(1)
		IF(-VCANT(iDLMENOR(1),1).GE.0.)THEN						  
			TETADPL(1)=DACOS(VCANT(iDLMENOR(1),2))			
		ELSE														
			TETADPL(1)=2.D0*PI-DACOS(VCANT(iDLMENOR(1),2))		
		ENDIF														
		IF(-VCANT(iDLMENOR(2),1).GE.0.)THEN						 
			TETADPL(2)=DACOS(VCANT(iDLMENOR(2),2))				
		ELSE														
			TETADPL(2)=2.D0*PI-DACOS(VCANT(iDLMENOR(2),2))		
		ENDIF	
    ENDIF    
    
    
ENDIF
    ENDSUBROUTINE CLASSIFICA_PONTO
    
