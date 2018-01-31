
    
!######################################################################
!Decomposicação SVD
!######################################################################
subroutine svdsolve(a,FGLOBAL,row,soln,XGLOBAL)
	implicit none

    INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)

   INTEGER, INTENT(IN) :: row,soln
   REAL(DBP),INTENT(IN)::FGLOBAL(ROW)
	REAL(dbp) , INTENT(IN) :: a(row,row)
	REAL(dbp):: ans(row,soln)
    REAL(DBP),INTENT(OUT)::XGLOBAL(ROW)
    REAL(DBP):: b(row,soln)
    
   REAL(dbp), allocatable :: usvd(:,:), wsvd(:), vsvd(:,:), winvutb(:,:)
   INTEGER :: ii,ij
   REAL(dbp) :: wmax, wmin, factor1
   CHARACTER :: ask2*80

   REAL(dbp), parameter :: EPS = 3.0d-15
   REAL(dbp) :: tmpmat(row,row), wmat(row,row), wmatinv(row, row)
    
   ALLOCATE ( usvd(row,row), wsvd(row), vsvd(row,row), winvutb(row,soln) )

   
    
   B(:,1)=FGLOBAL(:)
   ANS(:,1)=XGLOBAL


   usvd(:,:) = a(:,:)
   vsvd(:,:) = 0.0d0
   wsvd(:) = 0.0d0

! Produces a SINGULAR VALUE DECOMPOSITION of matrix a
  
   call SVDCMP(usvd, row,row, row,row, wsvd, vsvd)

   wmat=0.d0 !identity_matrix(row)
do ij=1,row
    wmat(ij,ij)=1.d0
enddo
   
   do ij = 1,row
   	wmat(ij,ij) = wsvd(ij)
   end do
   tmpmat = a - MATMUL(usvd, MATMUL(wmat,TRANSPOSE(vsvd)) )
   WRITE(60,*) 'Here is it'
   do ii = 1, row
   	WRITE(60,333) (tmpmat(ii,ij), ij=1,row)
   ENDDO
333 FORMAT(300g15.3)



! KEY TO SVD - If Diagonal elements are too small then MAKE IT ZERO!!!!!!!
   PRINT *, 'Solution by Singular Value Decomposition'
   PRINT *, 'Is matrix expected to be close to singular (y/n) ?'
   READ(*, *) ask2

   if (ask2 .eq. 'y') then
	   wmax = MAXVAL(dabs(wsvd(:)))
   	PRINT *, 'enter precision for svd (wmax = ',wmax, ')'
	   READ *, factor1
   	wmin = factor1 * eps * wmax          ! 10 times the set precision
! Actually DON'T set w(i,i) = 0, the fact that it's equivalent to zero lies in the
! next few lines in Solving for Ax=b or x=V.[1/w(i,i)].U`.b
		PRINT *, 'w-min threshold   Precision'
   	PRINT *, wmin, eps
   else
      wmin = 0.0d0
   end if


   wmatinv(:,:) = 0.0d0
   do ij = 1,row                       ! modified 9/9/98
     	if ( dabs(wsvd(ij)) .ge.wmin) 	wmatinv(ij,ij) = 1.0/wsvd(ij)
   end do
   winvutb = MATMUL(wmatinv, MATMUL(TRANSPOSE(usvd), b))

   ans(:,:) = 0.0d0
   ans(:,:) = MATMUL(vsvd,winvutb)

   XGLOBAL(:)=ANS(:,1)
   
   deallocate (usvd, wsvd, vsvd, winvutb)

   return
    END subroutine svdsolve