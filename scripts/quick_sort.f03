!===========================================
! Quick sort arrays of various data types.
!===========================================




!===================================
module quick_sort
!===================================

  use constants_and_parameters
  implicit none

  contains

  SUBROUTINE quick_sort_int(list, order, n)

    !------------------------------------------------------------
    ! Sort array of integers (list), rearrange array of integers
    ! (order) in the same way.
    !------------------------------------------------------------
    
    IMPLICIT NONE
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified to sort the second given array by the same rules.


    INTEGER :: n
    INTEGER(i8), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: order


    CALL quick_sort_1_int_int(1, n)

    CONTAINS

      RECURSIVE SUBROUTINE quick_sort_1_int_int(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        !     Local variables
        INTEGER             :: i, j, itemp
        INTEGER(i8)         :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
           ! Use interchange sort for small lists
           CALL interchange_sort_int_int(left_end, right_end)

        ELSE
           ! Use partition ("quick") sort
           reference = list((left_end + right_end)/2)
           i = left_end - 1; j = right_end + 1

          DO
              ! Scan list from left end until element >= reference is found
              DO
                 i = i + 1
                 IF (list(i) >= reference) EXIT
              END DO
              ! Scan list from right end until element <= reference is found
              DO
                 j = j - 1
                 IF (list(j) <= reference) EXIT
              END DO


              IF (i < j) THEN
                 ! Swap two out-of-order elements
                 temp = list(i); list(i) = list(j); list(j) = temp
                 itemp = order(i); order(i) = order(j); order(j) = itemp
              ELSE IF (i == j) THEN
                 i = i + 1
                 EXIT
              ELSE
                 EXIT
              END IF
           END DO
           IF (left_end < j) CALL quick_sort_1_int_int(left_end, j)
           IF (i < right_end) CALL quick_sort_1_int_int(i, right_end)
        END IF

      END SUBROUTINE quick_sort_1_int_int


      SUBROUTINE interchange_sort_int_int(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        !     Local variables
        INTEGER             :: i, j, itemp
        INTEGER(i8)         :: temp

        DO i = left_end, right_end - 1
           DO j = i+1, right_end
              IF (list(i) > list(j)) THEN
                 temp = list(i); list(i) = list(j); list(j) = temp
                 itemp = order(i); order(i) = order(j); order(j) = itemp
              END IF
           END DO
        END DO

      END SUBROUTINE interchange_sort_int_int

  END SUBROUTINE quick_sort_int



  SUBROUTINE quick_sort_real(list, order, n)

    !------------------------------------------------------------
    ! Sort array of integers (list), rearrange array of integers
    ! (order) in the same way.
    !------------------------------------------------------------
    
    IMPLICIT NONE
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified to sort the second given array by the same rules.


    INTEGER :: n
    REAL(dp), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: order


    CALL quick_sort_1_real(1, n)

    CONTAINS

      RECURSIVE SUBROUTINE quick_sort_1_real(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        !     Local variables
        INTEGER             :: i, j, itemp
        REAL(dp)            :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
           ! Use interchange sort for small lists
           CALL interchange_sort_real(left_end, right_end)

        ELSE
           ! Use partition ("quick") sort
           reference = list((left_end + right_end)/2)
           i = left_end - 1; j = right_end + 1

          DO
              ! Scan list from left end until element >= reference is found
              DO
                 i = i + 1
                 IF (list(i) >= reference) EXIT
              END DO
              ! Scan list from right end until element <= reference is found
              DO
                 j = j - 1
                 IF (list(j) <= reference) EXIT
              END DO


              IF (i < j) THEN
                 ! Swap two out-of-order elements
                 temp = list(i); list(i) = list(j); list(j) = temp
                 itemp = order(i); order(i) = order(j); order(j) = itemp
              ELSE IF (i == j) THEN
                 i = i + 1
                 EXIT
              ELSE
                 EXIT
              END IF
           END DO
           IF (left_end < j) CALL quick_sort_1_real(left_end, j)
           IF (i < right_end) CALL quick_sort_1_real(i, right_end)
        END IF

      END SUBROUTINE quick_sort_1_real


      SUBROUTINE interchange_sort_real(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        !     Local variables
        INTEGER             :: i, j, itemp
        REAL(dp)            :: temp

        DO i = left_end, right_end - 1
           DO j = i+1, right_end
              IF (list(i) > list(j)) THEN
                 temp = list(i); list(i) = list(j); list(j) = temp
                 itemp = order(i); order(i) = order(j); order(j) = itemp
              END IF
           END DO
        END DO

      END SUBROUTINE interchange_sort_real


  END SUBROUTINE quick_sort_real

end module quick_sort




