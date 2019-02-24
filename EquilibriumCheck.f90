MODULE EquilibriumCheck
!
USE globals
USE QL_routines
!
! Computes check for best response and equilibrium in all states and for all agents
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeEqCheck ( iModel )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel
    !
    ! Declaring local variable
    !
    INTEGER, PARAMETER :: numThresPathCycleLength = 10
    INTEGER, PARAMETER :: ThresPathCycleLength(numThresPathCycleLength) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
    INTEGER :: numPeriods, iGame, iAgent, iState, iPrice, iPeriod, iThres, i, j, PreCycleLength, CycleLength, &
        OptimalStrategy(numStates,numAgents), p(DepthState,numAgents), pPrime(numAgents), VisitedStates(numStates+1), &
        OptimalPrice, PathCycleStates(numStates), lastObservedStateNumber, PathCycleLength(numGames), &
        numStatesBRAll(numAgents,numGames), numStatesBRPathCycle(numAgents,numGames), &
        numStatesEqAll(numGames), numStatesEqPathCycle(numGames), &
        IsBestReply(numStates,numAgents), &
        numImprovedPrices, ImprovedPrices(numStates), numPathCycleLength(numThresPathCycleLength)
    REAL(8) :: DiscountFactors(0:numStates,numAgents), VisitedProfits(numStates+1), PreCycleProfit, CycleProfit, &
        StateValueFunction(numStates,numPrices), MaxStateValueFunction, &
        freqStatesBRAll(numAgents,numGames), freqStatesBRPathCycle(numAgents,numGames), &
        freqStatesEqAll(numGames), freqStateseqPathCycle(numGames), &
        avgFreqStatesBRAll(0:numAgents,numThresPathCycleLength), avgFreqStatesBRPathCycle(0:numAgents,numThresPathCycleLength), &
        avgFreqStatesEqAll(numThresPathCycleLength), avgFreqStatesEqPathCycle(numThresPathCycleLength)
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    !
    ! Beginning execution
    !
    PRINT*, 'Computing equilibrium checks'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    numPeriods = numStates+1        ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
    DiscountFactors = TRANSPOSE(RESHAPE((/ (delta**iPeriod, iPeriod = 0, numPeriods-1) /),(/ numAgents,numPeriods /)))
    !
    numStatesBRAll = 0
    numStatesBRPathCycle = 0
    numStatesEqAll = 0
    numStatesEqPathCycle = 0
    PathCycleLength = 0
    !
    ! Reading strategies and states at convergence from file
    !
    OPEN(UNIT = 998,FILE = FileNameIndexStrategies,STATUS = "OLD")    ! Open indexStrategies file
    DO i = 1, lengthStrategies
        !
        IF (MOD(i,10000) .EQ. 0) PRINT*, 'Read ', i, ' lines of indexStrategies'
        READ(998,21) (indexStrategies(i,iGame), iGame = 1, numGames)
    21  FORMAT(<numGames>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    CLOSE(UNIT = 998)                   ! Close indexStrategies file
    OPEN(UNIT = 999,FILE = FileNameIndexLastState,STATUS = "OLD")     ! Open indexLastState file
    DO iGame = 1, numGames
        !
        READ(999,22) indexLastState(:,iGame)
    22  FORMAT(<LengthStates>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    PRINT*, 'Read indexLastState'
    CLOSE(UNIT = 999)                   ! Close indexLastState file
    !
    ! Beginning loop over games
    !
    !$omp parallel do &
    !$omp private(OptimalStrategy,lastObservedStateNumber,IsBestReply,iAgent,StateValueFunction,PathCycleStates, &
    !$omp   iState,p,pPrime,OptimalPrice,iPrice,VisitedStates,VisitedProfits,iPeriod, &
    !$omp   PreCycleLength,CycleLength,PreCycleProfit,CycleProfit,OptimalStrategyVec,LastStateVec,i, &
    !$omp   MaxStateValueFunction,numImprovedPrices,ImprovedPrices) &
    !$omp firstprivate(numGames,PI)
    DO iGame = 1, numGames                  ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        !$omp critical
        OptimalStrategyVec = indexStrategies(:,iGame)
        LastStateVec = indexLastState(:,iGame)
        !$omp end critical
        !
        optimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        lastObservedStateNumber = computeStateNumber(RESHAPE(LastStateVec, (/ DepthState,numAgents /) ))
        !
        IsBestReply = 0
        DO iAgent = 1, numAgents            ! Start of loop over agents
            !
            StateValueFunction = 0.d0
            PathCycleStates = 0
            !
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute state value function for the optimal strategy in all states and actions
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            DO iState = 1, numStates        ! Start of loop over states
                !
                p = RESHAPE(convertNumberBase(iState-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
                pPrime = OptimalStrategy(iState,:)
                OptimalPrice = pPrime(iAgent)
                !
                DO iPrice = 1, numPrices    ! Start of loop over prices
                    !
                    ! First period deviation from optimal strategy
                    !
                    pPrime(iAgent) = iPrice
                    !
                    ! 1. Compute pre-cycle and cycle data
                    !
                    VisitedStates = 0
                    VisitedProfits = 0.d0
                    DO iPeriod = 1, numPeriods
                        !
                        IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                        p(1,:) = pPrime
                        VisitedStates(iPeriod) = computeStateNumber(p)
                        VisitedProfits(iPeriod) = PI(computeActionNumber(pPrime),iAgent)
                        !
                        ! Check if the state has already been visited
                        !
                        IF ((iPeriod .GE. 2) .AND. (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) THEN
                            !
                            PreCycleLength = MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
                            CycleLength = iPeriod-PreCycleLength
                            EXIT
                            !
                        END IF
                        !
                        ! After period 1, every agent follows the optimal strategy
                        !
                        pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
                        !
                    END DO
                    !
                    ! 2. Compute state value function for the optimal strategy
                    !
                    IF ((iPrice .EQ. OptimalPrice) .AND. (iState .EQ. lastObservedStateNumber)) THEN
                        !
                        PathCycleStates(:CycleLength) = VisitedStates(PreCycleLength+1:iPeriod)
                        PathCycleLength(iGame) = CycleLength
                        !
                    END IF
                    PreCycleProfit = SUM(DiscountFactors(0:PreCycleLength-1,iAgent)*VisitedProfits(1:PreCycleLength))
                    CycleProfit = SUM(DiscountFactors(0:CycleLength-1,iAgent)*VisitedProfits(PreCycleLength+1:iPeriod))
                    StateValueFunction(iState,iPrice) = &
                        PreCycleProfit+delta(iAgent)**PreCycleLength*CycleProfit/(1.d0-delta(iAgent)**CycleLength)
                    !
                END DO                      ! End of loop over prices
                !
            END DO                          ! End of loop over initial states
            !
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Check best reply condition on all states and on the observed cycle only
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            DO iState = 1, numStates        ! Start of loop over states
                !
                pPrime = OptimalStrategy(iState,:)
                OptimalPrice = pPrime(iAgent)
                MaxStateValueFunction = MAXVAL(StateValueFunction(iState,:))
                numImprovedPrices = 0
                ImprovedPrices = 0
                DO iPrice = 1, numPrices    ! Start of loop over prices to find optimal price(s)
                    !
                    IF (ABS(StateValueFunction(iState,iPrice)-MaxStateValueFunction) .LE. EPSILON(MaxStateValueFunction)) THEN
                        !
                        numImprovedPrices = numImprovedPrices+1                        
                        ImprovedPrices(numImprovedPrices) = iPrice
                        !
                    END IF
                    !
                END DO
                !
                IF (ANY(ImprovedPrices(:numImprovedPrices) .EQ. OptimalPrice)) THEN
                    !
                    IsBestReply(iState,iAgent) = 1
                    !
                    ! For all states
                    !
                    numStatesBRAll(iAgent,iGame) = numStatesBRAll(iAgent,iGame)+1
                    IF (ANY(PathCycleStates(:PathCycleLength(iGame)) .EQ. iState)) THEN
                        !
                        ! This state belongs to the observed cycle
                        !
                        numStatesBRPathCycle(iAgent,iGame) = numStatesBRPathCycle(iAgent,iGame)+1
                        !
                    END IF
                    !
                END IF
                !
            END DO                          ! End of loop over states
            freqStatesBRAll(iAgent,iGame) = DBLE(numStatesBRAll(iAgent,iGame))/DBLE(numStates)
            freqStatesBRPathCycle(iAgent,iGame) = DBLE(numStatesBRPathCycle(iAgent,iGame))/DBLE(PathCycleLength(iGame))
            !
        END DO                              ! End of loop over agents
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Check equilibrium condition on all states and on the observed cycle only
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iState = 1, numStates        ! Start of loop over states
            !
            IF (ALL(IsBestReply(iState,:) .EQ. 1)) THEN
                !
                ! For all states
                !
                numStatesEqAll(iGame) = numStatesEqAll(iGame)+1
                IF (ANY(PathCycleStates(:PathCycleLength(iGame)) .EQ. iState)) THEN
                    !
                    ! This state belongs to the observed cycle
                    !
                    numStatesEqPathCycle(iGame) = numStatesEqPathCycle(iGame)+1
                    !
                END IF
                !
            END IF
            !
        END DO                          ! End of loop over states
        freqStatesEqAll(iGame) = DBLE(numStatesEqAll(iGame))/DBLE(numStates)
        freqStatesEqPathCycle(iGame) = DBLE(numStatesEqPathCycle(iGame))/DBLE(PathCycleLength(iGame))
        !
    END DO                                  ! End of loop over games
    !$omp end parallel do
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Computing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Frequencies for optimal cycles' lengths
    !
    avgFreqStatesBRAll = 0.d0
    avgFreqStatesBRPathCycle = 0.d0
    numPathCycleLength = 0
    DO iThres = 1, numThresPathCycleLength
        !
        IF (iThres .LT. numThresPathCycleLength) THEN
            !
            numPathCycleLength(iThres) = COUNT(PathCycleLength .EQ. ThresPathCycleLength(iThres))
            DO iAgent = 1, numAgents
                !
                IF (numPathCycleLength(iThres) .GT. 0) THEN
                    !
                    avgFreqStatesBRAll(iAgent,iThres) = &
                        SUM(freqStatesBRAll(iAgent,:),MASK = PathCycleLength .EQ. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                    avgFreqStatesBRPathCycle(iAgent,iThres) = &
                        SUM(freqStatesBRPathCycle(iAgent,:),MASK = PathCycleLength .EQ. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                    !
                ELSE
                    !
                    avgFreqStatesBRAll(iAgent,iThres) = 0.d0
                    avgFreqStatesBRPathCycle(iAgent,iThres) = 0.d0
                    !
                END IF
                !
            END DO
            !
            IF (numPathCycleLength(iThres) .GT. 0) THEN
                !
                avgFreqStatesBRAll(0,iThres) = &
                    SUM(freqStatesBRAll,MASK = SPREAD(PathCycleLength,DIM = 1,NCOPIES = numAgents) .EQ. ThresPathCycleLength(iThres))/DBLE(numAgents*numPathCycleLength(iThres))
                avgFreqStatesBRPathCycle(0,iThres) = &
                    SUM(freqStatesBRPathCycle,MASK = SPREAD(PathCycleLength,DIM = 1,NCOPIES = numAgents) .EQ. ThresPathCycleLength(iThres))/DBLE(numAgents*numPathCycleLength(iThres))
                avgFreqStatesEqAll(iThres) = &
                    SUM(freqStatesEqAll,MASK = PathCycleLength .EQ. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                avgFreqStatesEqPathCycle(iThres) = &
                    SUM(freqStatesEqPathCycle,MASK = PathCycleLength .EQ. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                !
            ELSE
                !
                avgFreqStatesBRAll(0,iThres) = 0.d0
                avgFreqStatesBRPathCycle(0,iThres) = 0.d0
                avgFreqStatesEqAll(iThres) = 0.d0
                avgFreqStatesEqPathCycle(iThres) = 0.d0
                !
            END IF
            !
        ELSE
            !
            numPathCycleLength(iThres) = COUNT(PathCycleLength .GE. ThresPathCycleLength(iThres))
            DO iAgent = 1, numAgents
                !
                IF (numPathCycleLength(iThres) .GT. 0) THEN
                    !
                    avgFreqStatesBRAll(iAgent,iThres) = &
                        SUM(freqStatesBRAll(iAgent,:),MASK = PathCycleLength .GE. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                    avgFreqStatesBRPathCycle(iAgent,iThres) = &
                        SUM(freqStatesBRPathCycle(iAgent,:),MASK = PathCycleLength .GE. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                    !
                ELSE
                    !
                    avgFreqStatesBRAll(iAgent,iThres) = 0.d0
                    avgFreqStatesBRPathCycle(iAgent,iThres) = 0.d0
                    !
                END IF
                !
            END DO
            !
            IF (numPathCycleLength(iThres) .GT. 0) THEN
                !
                avgFreqStatesBRAll(0,iThres) = &
                    SUM(freqStatesBRAll,MASK = SPREAD(PathCycleLength,DIM = 1,NCOPIES = numAgents) .GE. ThresPathCycleLength(iThres))/DBLE(numAgents*numPathCycleLength(iThres))
                avgFreqStatesBRPathCycle(0,iThres) = &
                    SUM(freqStatesBRPathCycle,MASK = SPREAD(PathCycleLength,DIM = 1,NCOPIES = numAgents) .GE. ThresPathCycleLength(iThres))/DBLE(numAgents*numPathCycleLength(iThres))
                avgFreqStatesEqAll(iThres) = &
                    SUM(freqStatesEqAll,MASK = PathCycleLength .GE. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                avgFreqStatesEqPathCycle(iThres) = &
                    SUM(freqStatesEqPathCycle,MASK = PathCycleLength .GE. ThresPathCycleLength(iThres))/DBLE(numPathCycleLength(iThres))
                !
            ELSE
                !
                avgFreqStatesBRAll(0,iThres) = 0.d0
                avgFreqStatesBRPathCycle(0,iThres) = 0.d0
                avgFreqStatesEqAll(iThres) = 0.d0
                avgFreqStatesEqPathCycle(iThres) = 0.d0
                !
            END IF
            !
        END IF
        !
    END DO
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10004,1) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (i, i, i, i = 1, numAgents), &
            (iThres, iThres, iThres, iThres, iThres, (iThres, iAgent, iThres, iAgent, iThres, iAgent, iAgent = 1, numAgents), iThres = 1, numThresPathCycleLength)
1           FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>(' MExplPar', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I2.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            'TotNum TotAvgFreqBRAll TotAvgFreqBRPath TotAvgFreqEqAll TotAvgFreqEqPath ', &
            <numAgents>('Ag', I1, 'Num Ag', I1, 'AvgFreqBRAll Ag', I1, 'AvgFreqBRPath '), &
            <numThresPathCycleLength>('PathLen', I2.2, 'Num ', &
                                      'PathLen', I2.2, 'AvgFreqBRAll PathLen', I2.2, 'AvgFreqBRPath ', &
                                      'PathLen', I2.2, 'AvgFreqEqAll PathLen', I2.2, 'AvgFreqEqPath ', &
                <numAgents>('PathLen', I2.2, 'Ag', I1, 'Num PathLen', I2.2, 'Ag', I1, 'AvgFreqBRAll PathLen', I2.2, 'Ag', I1, 'AvgFreqBRPath ')) &
            )
        !
    END IF
    !
    WRITE(10004,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        numAgents*numGames, &
            SUM(freqStatesBRAll)/DBLE(numAgents*numGames), SUM(freqStatesBRPathCycle)/DBLE(numAgents*numGames), &
            SUM(freqStatesEqAll)/DBLE(numGames), SUM(freqStatesEqPathCycle)/DBLE(numGames), &
            (numGames, SUM(freqStatesBRAll(iAgent,:))/DBLE(numGames), SUM(freqStatesBRPathCycle(iAgent,:))/DBLE(numGames), iAgent = 1, numAgents), &
        (numAgents*numPathCycleLength(iThres), &
            avgFreqStatesBRAll(0,iThres), avgFreqStatesBRPathCycle(0,iThres), &
            avgFreqStatesEqAll(iThres), avgFreqStatesEqPathCycle(iThres), &
            (numPathCycleLength(iThres), avgFreqStatesBRAll(iAgent,iThres), avgFreqStatesBRPathCycle(iAgent,iThres), iAgent = 1, numAgents), &
            iThres = 1, numThresPathCycleLength)
2   FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.3, 1X), &
        <6*numAgents>(F10.7, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        (I6, 1X, F15.7, 1X, F16.7, 1X, F15.7, 1X, F16.7, 1X), &
            <numAgents>(I6, 1X, F15.7, 1X, F16.7, 1X), &
        <numThresPathCycleLength>(I12, 1X, F21.7, 1X, F22.7, 1X, F21.7, 1X, F22.7, 1X, &
            <numAgents>(I15, 1X, F24.7, 1X, F25.7, 1X)) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeEqCheck
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE EquilibriumCheck
