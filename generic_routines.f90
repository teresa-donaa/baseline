MODULE generic_routines
!
! Various generic routines 
!
IMPLICIT NONE
!
CONTAINS
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION convertNumberBase ( n, b, l )
    !
    ! Converts an integer n from base 10 to base b, 
    ! generating a vector of integers of length l
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    ! 
    INTEGER, INTENT(IN) :: n, b, l
    !
    ! Declaring function's type
    !
    INTEGER, DIMENSION(l) :: convertNumberBase
    !
    ! Declaring local variables
    !
    INTEGER :: i, tmp
    !
    ! Beginning execution
    !
    tmp = n
    DO i = 1, l
        !
        convertNumberBase(l-i+1) = MOD(tmp,b)+1
        tmp = tmp/b
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION convertNumberBase
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION ran2 ( idum, iv, iy, idum2 )
    !
    ! Thread safe function that generates U(0,1) random deviates (to be used with OpenMP)
    ! See function RAN2 on p. 272 of NRF77
    ! Long period (> 2 � 10^18) random number generator of L�Ecuyer with Bays-Durham shuffle
    ! and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
    ! of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
    ! alter idum between successive deviates in a sequence. RNMX should approximate the largest
    ! floating value that is less than 1.
    ! Always set idum2 = 123456789, iv = 0, iy = 0 upon initialization
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Declaring function's type
    !
    REAL(8) :: ran2
    !
    ! Declaring local variables and parameters
    !
    INTEGER, PARAMETER :: IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1-1, IA1 = 40014, &
        IA2 = 40692, IQ1 = 53668, IQ2 = 52774, IR1 = 12211, IR2 = 3791, NDIV = 1+IMM1/32
    REAL(8), PARAMETER :: AM = 1.d0/IM1, EPS = 1.2d-7, RNMX = 1.d0-EPS
    INTEGER :: j, k
    !
    ! Beginning execution
    !
    ! Initializing
    !
    IF (idum .LE. 0) THEN
        !
        idum = MAX(-idum,1) 
        idum2 = idum
        DO j = 32+8,1,-1 
            !
            k = idum/IQ1
            idum = IA1*(idum-k*IQ1)-k*IR1
            IF (idum .LT. 0) idum = idum+IM1
            IF (j .LE. 32) iv(j) = idum
            !
        END DO
        !
        iy = iv(1)
        !
    END IF
    !
    ! Start here when not initializing
    !
    k = idum/IQ1 
    idum = IA1*(idum-k*IQ1)-k*IR1 
    IF (idum .LT. 0) idum = idum+IM1 
    k = idum2/IQ2
    idum2 = IA2*(idum2-k*IQ2)-k*IR2 
    IF (idum2 .LT. 0) idum2 = idum2+IM2
    j = 1+iy/NDIV 
    iy = iv(j)-idum2 
    iv(j) = idum
    !
    IF (iy .LT. 1) iy = iy+IMM1
    ran2 = MIN(AM*iy,RNMX) 
    RETURN
    !
    ! Ending execution and returning control
    !
    END function ran2
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! End of execution
!
END MODULE generic_routines