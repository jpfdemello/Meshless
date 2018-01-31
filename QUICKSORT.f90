!##########################################################################################################
!Esta subrotina ordena um conjunto de dados em ordem crescente 
!##########################################################################################################
		
recursive subroutine QuickSort(A,na)
USE TIPO
IMPLICIT NONE 
! DUMMY ARGUMENTS
integer, intent(in) :: nA

!type group
!    integer :: ordem,ordem2    ! original order of unsorted data
!    real*8:: valor       ! values to be sorted by
!end type group


type (group), dimension(nA), intent(in out) :: A

 
! LOCAL VARIABLES
integer :: left, right
real*8 :: random
real*8 :: pivot
type (group) :: temp
integer :: marker

    if (nA > 1) then
 
        call random_number(random)
        pivot = A(int(random*real(nA-1))+1)%valor   ! random pivor (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right)%valor > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%valor < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call QuickSort(A(:marker-1),marker-1)
        call QuickSort(A(marker:),nA-marker+1)
 
    end if
 
end subroutine QuickSort
    