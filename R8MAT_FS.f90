!#####################################################################
!############SOLVER AX=B #############################################
!#####################################################################
subroutine r8mat_fs ( n, a, b, x,JJ )

!*****************************************************************************80
!
!! R8MAT_FS factors and solves a linear system with one right hand side.
!
!  Discussion:
!
!    This routine uses partial pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the coefficient matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  INTEGER::JJ

  a2(1:n,1:n) = a(1:n,1:n)
  x(1:n) = b(1:n)

     
  
  do k = 1, n
!
!  Find the maximum element in column I.
!
    p = k

    do i = k + 1, n
      if ( abs ( a2(p,k) ) < abs ( a2(i,k) ) ) then
        p = i
      end if
    end do

    if ( a2(p,k) == 0.0D00 ) then
        write(*,*) a2(p,k),p,k
        IF(JJ.EQ.999)THEN
            WRITE(*,*) 'PROBLEMA NA MATRIZ GLOBAL'
            WRITE(*,*) 'SOLUCAO PROSSEGUE COM SVD'
            READ(*,*)
            JJ=998
            RETURN
        ELSE
            WRITE(*,*) 'PROBLEMA NO PONTO',JJ
        ENDIF
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', k
      pause
      stop 1
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##########!!!!!!!!!!!!!!!!!
    if ( ABS(a2(p,k)).lt.1.0E-20 ) then
        write(*,*) a2(p,k)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FS - Fatal error!'
      write ( *, '(a,i8)' ) ' proximo de Zero pivot on step ', k
        IF(JJ.EQ.999)THEN
            WRITE(*,*) 'PROBLEMA NA MATRIZ GLOBAL'
            WRITE(*,*) 'SOLUCAO PROSSEGUE COM SVD'
            READ(*,*)
            JJ=998
            RETURN
        ENDIF 
      pause
      stop 1
    end if

!!!!!!!!!!!!!!!!!!###############!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Switch rows K and P.
!
    if ( k /= p ) then

      row(1:n) = a2(k,1:n)
      a2(k,1:n) = a2(p,1:n)
      a2(p,1:n) = row(1:n)

      t    = x(k)
      x(k) = x(p)
      x(p) = t

    end if
!
!  Scale the pivot row.
!
    a2(k,k+1:n) = a2(k,k+1:n) / a2(k,k)
    x(k) = x(k) / a2(k,k)
    a2(k,k) = 1.0D+00 
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = k + 1, n
      if ( a2(i,k) /= 0.0D+00 ) then
        t = - a2(i,k)
        a2(i,k) = 0.0D+00
        a2(i,k+1:n) = a2(i,k+1:n) + t * a2(k,k+1:n)
        x(i) = x(i) + t * x(k)
      end if
    end do
  
  end do
!
!  Back solve.
!
  do j = n, 2, -1
    x(1:j-1) = x(1:j-1) - a2(1:j-1,j) * x(j)
  end do

  
  return
end	subroutine r8mat_fs