!##########################################################################################################
!Esta subrotina calcula as abscissas (XG) e os pesos (WG) para (N) pontos de Gauss. Basta fornecer a quantidade de (N) de 
!pontos que a subrotina devolve as abscissas e os pesos.
!##########################################################################################################
	SUBROUTINE Calcula_Pontos_Gauss(N,XG,WG)
	 IMPLICIT NONE
	INTEGER N,I,J,M
	REAL*8,DIMENSION(N):: XG,WG
	REAL*8 X1,X2,EPS,P1,P2,P3,PP,XL,XM,Z,Z1,XI,XJ,XN,pi
	PARAMETER (EPS=3.E-14)

	pi=atan(1.0d0)*4.0d0
	XN=float(N)
	M=(N+1)/2
	X2=1.d0
	X1=-1.d0
	XM=0.5d0*(X2+X1)
	XL=0.5d0*(X2-X1)
	DO I=1,M
      XI=float(I)
      Z=COS(pi*(XI-0.25d0)/(XN+0.5d0))
      Z1=Z+10.0d0*EPS
      Do While(DABS(Z-Z1).GT.EPS)
        P1=1.0d0
        P2=0.0d0
         DO J=1,N
           XJ=float(J)
           P3=P2
           P2=P1
           P1=((2.0d0*XJ-1.0d0)*Z*P2-(XJ-1.0d0)*P3)/XJ
         ENDDO
        PP=XN*(Z*P1-P2)/(Z*Z-1.0d0)
        Z1=Z
        Z=Z1-P1/PP
      End Do
      XG(I)=XM-XL*Z
      XG(N+1-I)=XM+XL*Z
      WG(I)=2.0d0*XL/((1.0d0-Z*Z)*PP*PP)
      WG(N+1-I)=WG(I)
	ENDDO
	RETURN
    END	SUBROUTINE Calcula_Pontos_Gauss
    