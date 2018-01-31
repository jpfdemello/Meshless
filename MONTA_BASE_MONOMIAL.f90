!##########################################################################################################
!Esta subrotina monta a base monomial até 3 ordem e suas derivadas (ordem,qtd de pontos,x dos pontos, y dos pontos 
!##########################################################################################################
SUBROUTINE  MONTA_BASE_MONOMIAL(MON,TMON,NPT,X,Y,M,DMX,DMY)
IMPLICIT NONE

INTEGER::i
INTEGER,INTENT(IN)::MON,NPT,TMON
REAL*8,DIMENSION(NPT),INTENT(IN)::X,Y
REAL*8,DIMENSION(NPT,TMON),INTENT(OUT)::M,DMX,DMY

IF (MON.EQ.1) THEN
	DO I=1,NPT
		m(i,1)=1.
		m(i,2)=X(i)
		m(i,3)=y(i)
		dmx(i,1)=0.
		dmx(i,2)=1.
		dmx(i,3)=0.
		dmy(i,1)=0.
		dmy(i,2)=0.
		dmy(i,3)=1.
	ENDDO
ELSE IF (MON.EQ.2) THEN
	DO I=1,NPT
		m(i,1)=1.
		m(i,2)=X(i)
		m(i,3)=(Y(i))
		m(i,4)=(X(i))**2.
		m(i,5)=(X(i))*(Y(i))
		m(i,6)=(Y(i))**2.
		dmx(i,1)=0.
		dmx(i,2)=1.
		dmx(i,3)=0.
		dmx(i,4)=2.*(X(i))
		dmx(i,5)=(Y(i))
		dmx(i,6)=0.
		dmy(i,1)=0.
		dmy(i,2)=0.
		dmy(i,3)=1.
		dmy(i,4)=0.
		dmy(i,5)=(X(i))
		dmy(i,6)=2.*(Y(i))
	ENDDO
ELSE IF (MON.EQ.3) THEN
	DO I=1,NPT
		m(i,1)=1.
		m(i,2)=X(i)
		m(i,3)=(Y(i))
		m(i,4)=(X(i))**2.
		m(i,5)=(X(i))*(Y(i))
		m(i,6)=(Y(i))**2.
		m(i,7)=(X(i))**3.
		m(i,8)=(Y(i))*((X(i))**2.)
		m(i,9)=(X(i))*((Y(i))**2.)
		m(i,10)=(Y(i))**3.
		dmx(i,1)=0.
		dmx(i,2)=1.
		dmx(i,3)=0.
		dmx(i,4)=2.*(X(i))
		dmx(i,5)=(Y(i))
		dmx(i,6)=0.
		dmx(i,7)=3.*((X(i))**2.)
		dmx(i,8)=2.*(X(i))*(Y(i))
		dmx(i,9)=(Y(i))**2.
		dmx(i,10)=0.
		dmy(i,1)=0.
		dmy(i,2)=0.
		dmy(i,3)=1.
		dmy(i,4)=0.
		dmy(i,5)=(X(i))
		dmy(i,6)=2.*(Y(i))
		dmy(i,7)=0.
		dmy(i,8)=(X(i))**2.
		dmy(i,9)=2.*(X(i))*(Y(i))
		dmy(i,10)=3.*((Y(i))**2.)
	ENDDO
ENDIF
RETURN
END SUBROUTINE  MONTA_BASE_MONOMIAL


