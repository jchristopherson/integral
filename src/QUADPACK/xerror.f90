subroutine xerror ( xmess, nmess, nerr, level )

!***********************************************
!
! XERROR replaces the SLATEC XERROR routine.
!
!  Modified:
!
!    12 September 2015
!

    integer ( kind = 4 ) level
    integer ( kind = 4 ) nerr
    integer ( kind = 4 ) nmess
    character ( len = * ) xmess

    if ( 1 <= LEVEL ) then
        WRITE ( *,'(1X,A)') XMESS(1:NMESS)
        WRITE ( *,'('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)') &
            NERR,LEVEL
    end if

    return
end
