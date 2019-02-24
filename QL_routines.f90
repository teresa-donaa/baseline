MODULE QL_routines
!
USE globals
USE m_median
!
! Various routines used in exp2
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE initQMatrices ( PI, delta, Q, maxValQ, maxLocQ )
    !
    ! Randomly initializing Q matrices
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), DIMENSION(numActions,numAgents), INTENT(IN) :: PI
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    REAL(8), DIMENSION(numStates,numPrices,numAgents), INTENT(OUT) :: Q
    INTEGER, DIMENSION(numStates,numAgents), INTENT(OUT) :: maxLocQ
    REAL(8), DIMENSION(numStates,numAgents), INTENT(OUT) :: maxValQ
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice, p(DepthState,numAgents), iAction, iNash, iDevNash
    INTEGER :: devPrices(numAgents)
    REAL(8) :: den
    !
    ! Beginning execution
    !
    ! 1. Randomizing over the opponents decisions
    !
    DO iAgent = 1, numAgents
        !
        DO iPrice = 1, numPrices
            !
            den = COUNT(indexActions(:,iAgent) .EQ. iPrice)*(1.d0-delta(iAgent))
            Q(:,iPrice,iAgent) = SUM(PI(:,iAgent),MASK = indexActions(:,iAgent) .EQ. iPrice)/den
            !
        END DO
        !
    END DO
    !
    ! 2. Repeated Nash equilibrium
    !
    !! IMPORTANT: 
    !! This only works when the repeated Nash equilibrium belongs to the prices grid
    !
    !iNash = 0
    !DO iAction = 1, numActions
    !    !
    !    IF (ALL(indexActions(iAction,:) .EQ. indexNashPrices)) THEN
    !        !
    !        iNash = iAction
    !        EXIT
    !        !
    !    END IF
    !    !
    !END DO
    !!
    !DO iAgent = 1, numAgents
    !    !
    !    DO iPrice = 1, numPrices
    !        !
    !        IF (iPrice .EQ. indexNashPrices(iAgent)) &
    !            Q(:,iPrice,iAgent) = PI(iNash,iAgent)/(1.d0-delta(iAgent))
    !        IF (iPrice .NE. indexNashPrices(iAgent)) THEN
    !            !
    !            devPrices = indexNashPrices
    !            devPrices(iAgent) = iPrice
    !            DO iAction = 1, numActions
    !                !
    !                IF (ALL(indexActions(iAction,:) .EQ. devPrices)) THEN
    !                    !
    !                    iDevNash = iAction
    !                    EXIT
    !                    !
    !                END IF
    !                !
    !            END DO
    !            Q(:,iPrice,iAgent) = &
    !                PI(iDevNash,iAgent)+delta(iAgent)*PI(iNash,iAgent)/(1.d0-delta(iAgent))
    !            !
    !        END IF
    !        !
    !    END DO
    !    !
    !END DO
    !
    maxValQ = MAXVAL(Q,DIM = 2)
    maxLocQ = MAXLOC(Q,DIM = 2)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initQMatrices
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE initState ( u, p, stateNumber, actionNumber )
    !
    ! Randomly initializing prices 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: u(DepthState,numAgents)
    INTEGER, DIMENSION(DepthState,numAgents), INTENT(OUT) :: p
    INTEGER, INTENT(OUT) :: stateNumber, actionNumber
    !
    ! Beginning execution
    !
    p = 1+INT(numPrices*u)
    stateNumber = computeStateNumber(p)
    actionNumber = computeActionNumber(p(1,:))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initState
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generate_uIniPrice ( numGames, uIniPrice, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: numGames
    REAL(8), INTENT(OUT) :: uIniPrice(DepthState,numAgents,numGames)
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Declaring local variables
    !
    INTEGER :: iGame, iAgent, iDepth
    !
    ! Beginning execution
    !
    ! Generate U(0,1) draws for price initialization
    !
    DO iGame = 1, numGames
        !
        DO iDepth = 1, DepthState
            !
            DO iAgent = 1, numAgents
                !
                uIniPrice(iDepth,iAgent,iGame) = ran2(idum,iv,iy,idum2)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE generate_uIniPrice
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generateUExploration ( uExploration, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(OUT) :: uExploration(2,numAgents)
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Declaring local variables
    !
    INTEGER :: iDecision, iAgent
    !
    ! Beginning execution
    !
    ! Generate U(0,1) draws for price initialization
    !
    DO iDecision = 1, 2
        !
        DO iAgent = 1, numAgents
            !
            uExploration(iDecision,iAgent) = ran2(idum,iv,iy,idum2)
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE generateUExploration
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStateNumber ( p )
    !
    ! Given the price vectors, computes the state number
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(DepthState,numAgents), INTENT(IN) :: p
    !
    ! Declaring function's type
    !
    INTEGER :: computeStateNumber
    !
    ! Declaring local variables
    !
    INTEGER, DIMENSION(LengthStates) :: stateVector
    !
    ! Beginning execution
    !
    IF (DepthState0 .GT. 0) THEN
        !
        stateVector = RESHAPE(TRANSPOSE(p),(/ LengthStates /))
        computeStateNumber = 1+SUM(cStates*(stateVector-1))
        !
    ELSE IF (DepthState0 .EQ. 0) THEN
        !
        computeStateNumber = 1
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStateNumber
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeActionNumber ( p )
    !
    ! Given the prices, computes the action number
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(numAgents), INTENT(IN) :: p
    !
    ! Declaring function's type
    !
    INTEGER :: computeActionNumber
    !
    ! Declaring local variables
    !
    INTEGER, DIMENSION(numAgents) :: tmp
    !
    ! Beginning execution
    !
    tmp = cActions*(p-1)
    computeActionNumber = 1+SUM(tmp)
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeActionNumber
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStatesCodePrint ( )
    !
    ! Compute the states code in printable format (with '.')
    !
    IMPLICIT NONE
    !
    ! Declaring function's type
    !
    CHARACTER(len = lengthStatesPrint) :: computeStatesCodePrint(numStates)
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, indexState(LengthStates)
    CHARACTER(len = lengthFormatActionPrint) :: tmp
    CHARACTER(len = lengthStatesPrint) :: labelState
    !
    ! Beginning execution
    !
    DO i = 1, numStates
        !
        indexState = convertNumberBase(i-1,numPrices,LengthStates)
        !
        DO j = 1, LengthStates
            !
            WRITE(tmp,'(I0.<lengthFormatActionPrint>)') indexState(j)
            IF (j .EQ. 1) THEN 
                !
                labelState = TRIM(tmp)   
                !
            ELSE IF (MOD(j,numAgents) .NE. 1) THEN
                !
                labelState = TRIM(labelState) // '.' // TRIM(tmp)   
                !
            ELSE
                !
                labelState = TRIM(labelState) // '-' // TRIM(tmp)   
                !
            END IF
            !
        END DO
        !
        computeStatesCodePrint(i) = labelState
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStatesCodePrint
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStrategyNumber ( maxLocQ )
    !
    ! Given the maxLocQ vectors, computes the lengthStrategies-digit strategy number
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(numStates,numAgents), INTENT(IN) :: maxLocQ
    !
    ! Declaring function's type
    !
    INTEGER :: computeStrategyNumber(lengthStrategies)
    !
    ! Declaring local variables
    !
    INTEGER :: i, il, iu
    !
    ! Beginning execution
    !
    iu = 0
    DO i = 1, numAgents
        !
        il = iu+1
        iu = iu+numStates
        computeStrategyNumber(il:iu) = maxLocQ(:,i)
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStrategyNumber
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeIndicators ( iModel, converged, timeToConvergence, &
        numGamesConverged, meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence ) 
    !
    ! Computes the output indicators from the resilts of a simulation experiment
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel
    INTEGER, INTENT(IN) :: converged(numGames)
    REAL(8), DIMENSION(numGames), INTENT(IN) :: timeToConvergence
    INTEGER, INTENT(OUT) :: numGamesConverged
    REAL(8), INTENT(OUT) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
    !
    ! Declaring local variables and parameters
    !
    INTEGER :: i, j, h, l, iAgent
    !
    ! Beginning execution
    !
    ! Preliminaries
    !
    numGamesConverged = SUM(converged)
    maskConverged = (converged .EQ. 1)
    meanNashProfit = SUM(NashProfits)/numAgents
    meanCoopProfit = SUM(CoopProfits)/numAgents
    !
    ! Time to convergence
    !
    meanTimeToConvergence = SUM(timeToConvergence,MASK = maskConverged)/numGamesConverged
    seTimeToConvergence = &
        SQRT(SUM(timeToConvergence**2,MASK = maskConverged)/numGamesConverged-meanTimeToConvergence**2)
    medianTimeToConvergence = median(timeToConvergence)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeIndicators
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! End of execution
!
END MODULE QL_routines