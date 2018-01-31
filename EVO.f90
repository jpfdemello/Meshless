    
!##################################################################
!subrotina para imprimir o andamento do programa
!##################################################################
SUBROUTINE EVO(i,NPT)
IMPLICIT NONE

INTEGER::I,NPT
REAL::J
J=i*100./NPT

WRITE(*,10) J,'%'
10 FORMAT (F7.3,A1)



RETURN
END SUBROUTINE EVO