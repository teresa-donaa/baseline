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
    SUBROUTINE initQMatrices ( iGames, uRandomSampling, PI, delta, Q, maxValQ, maxLocQ )
    !
    ! Randomly initializing Q matrices
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iGames
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: uRandomSampling
    REAL(8), DIMENSION(numActions,numAgents), INTENT(IN) :: PI
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    REAL(8), DIMENSION(numStates,numPrices,numAgents), INTENT(OUT) :: Q
    INTEGER, DIMENSION(numStates,numAgents), INTENT(OUT) :: maxLocQ
    REAL(8), DIMENSION(numStates,numAgents), INTENT(OUT) :: maxValQ
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice, iState, i 
    INTEGER :: status
    REAL(8) :: den
    CHARACTER(len = 225) :: QFileName
    CHARACTER(len = 5) :: iChar
    CHARACTER(len = 5) :: codModelChar
    CHARACTER(len = 200) :: QFileFolderName
    !
    ! Beginning execution
    !
    DO iAgent = 1, numAgents
        !
        IF (typeQInitialization(iAgent) .EQ. 0) THEN
            !
            ! Randomizing over the opponents decisions
            !
            DO iPrice = 1, numPrices
                !
                den = COUNT(indexActions(:,iAgent) .EQ. iPrice)*(1.d0-delta(iAgent))
                Q(:,iPrice,iAgent) = SUM(PI(:,iAgent),MASK = indexActions(:,iAgent) .EQ. iPrice)/den
                !
            END DO
            !
        ELSE IF (typeQInitialization(iAgent) .GT. 0) THEN
            !
            ! Start from a randomly drawn Q matrix at convergence 
            ! on model "typeQInitialization(iAgent)"
            !
            WRITE(codModelChar,'(I0.5)') typeQInitialization(iAgent)
            i = 1+INT(DBLE(numGames)*uRandomSampling(iAgent))
            WRITE(iChar,'(I0.5)') i
            QFileName = 'Q_' // codModelChar // '_' // iChar // '.txt'
            IF (iAgent .EQ. 1) THEN
                !
                QFileFolderName = Q1FileFolderName
                !
            ELSE IF (iAgent .EQ. 2) THEN
                !
                QFileFolderName = Q2FileFolderName
                !
            ELSE
                !
                PRINT*, 'Reading external Q matrix only implemented for 2 agents in QL_routines.f90'
                PAUSE
                !
            END IF
            QFileName = TRIM(QFileFolderName) // TRIM(QFileName)
            !
            ! Write on Q matrices to file
            !
            OPEN(UNIT = iGames,FILE = QFileName,READONLY,RECL = 10000,IOSTAT = status)
            IF (iAgent .GT. 1) READ(iGames,100)
100         FORMAT(<(iAgent-1)*numStates-1>(/))
            DO iState = 1, numStates
                !
                READ(iGames,*) Q(iState,:,iAgent)
                !
            END DO
            CLOSE(UNIT = iGames)
            !
        END IF
        !
    END DO
    !
    ! Find initial optimal strategy
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
!@SP Sanity check    
!p(1,1) = 14
!@SP Sanity check    
    stateNumber = computeStateNumber(p)
    actionNumber = computeActionNumber(p(1,:))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initState
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generate_uIniPrice ( uIniPrice, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
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
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generateURandomSampling ( uRandomSampling, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(OUT) :: uRandomSampling(numGames,numAgents)
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Declaring local variables
    !
    INTEGER :: iGame, iAgent
    !
    ! Beginning execution
    !
    ! Generate U(0,1) draws for price initialization
    !
    DO iGame = 1, numGames
        !
        DO iAgent = 1, numAgents
            !
            uRandomSampling(iGame,iAgent) = ran2(idum,iv,iy,idum2)
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE generateURandomSampling
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
    SUBROUTINE computeIndicators ( iModel, converged, timeToConvergence ) 
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
    !
    ! Declaring local variables 
    !
    INTEGER :: i, j, h, l, iAgent
    LOGICAL :: maskConverged(numGames)
    INTEGER :: numGamesConverged
    REAL(8) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
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
    ! Print output
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10002,891) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents)
891         FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I0.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I0.2, ' ')), &
            '   numConv ', &
            '    avgTTC      seTTC     medTTC ')
        !
    END IF
    !
    WRITE(10002,991) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        numGamesConverged, &
        meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
991 FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
        <numPrices*numAgents>(F10.7, 1X), &
        I10, 1X, &
        <3>(F10.2, 1X))
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