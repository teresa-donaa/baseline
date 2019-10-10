MODULE PIGapToMaximum
!
USE globals
USE QL_routines
USE EquilibriumCheck
!
! Computes gap in average undiscounted profit values w.r.t. maximum
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeAvgPIGapToMax ( iModel )
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
    INTEGER :: iPeriod, iGame, iState, iAgent, iPrice, iThres, i, j, PathCycleStates(numPeriods), &
        PathCycleLength(numGames), OptimalStrategy(numStates,numAgents), lastObservedStateNumber, &
        p(DepthState,numAgents), pPrime(numAgents), OptimalPrice, &
        VisitedStates(numPeriods), PreCycleLength, GameCycleLength, &
        PossibleActions(numPrices), PriceVector(numPrices)
    INTEGER, DIMENSION(0:numThresPathCycleLength,0:numAgents) :: NumAvgPIGapTot, NumAvgPIGapOnPath, &
        NumAvgPIGapNotOnPath, NumAvgPIGapNotBRAllStates, NumAvgPIGapNotBRonPath, NumAvgPIGapNotEqAllStates, NumAvgPIGapNotEqonPath
    INTEGER, DIMENSION(numPeriods,numGames) :: CycleStates
    INTEGER, DIMENSION(numGames) :: CycleLength    
    INTEGER, DIMENSION(numAgents,numPeriods,numGames) :: CyclePrices
    REAL(8) :: VisitedProfits(numPeriods), MaxProfits(numPeriods), &
        AvgVisitedProfits(numStates,numAgents), AvgMaxProfits(numStates,numAgents), &
        AvgPIGap(numStates,numAgents), AvgProfitGainGap(numStates,numAgents), &
        QTrue(numStates,numPrices,numAgents), QGap(numStates,numAgents), MaxQTrue(numStates,numAgents), tmp
    REAL(8), DIMENSION(0:numThresPathCycleLength,0:numAgents) :: SumAvgPIGapTot, SumAvgPIGapOnPath, &
        SumAvgPIGapNotOnPath, SumAvgPIGapNotBRAllStates, SumAvgPIGapNotBRonPath, SumAvgPIGapNotEqAllStates, SumAvgPIGapNotEqonPath
    REAL(8), DIMENSION(numAgents,numPeriods,numGames) :: CycleProfits
    LOGICAL, DIMENSION(numStates,numAgents) :: IsOnPath, IsBRAllStates, IsBRonPath, IsEqAllStates, IsEqonPath
    LOGICAL, DIMENSION(numAgents) :: IsBR
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    !
    ! Beginning execution
    !
    PRINT*, 'Computing profit gaps'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    PathCycleLength = 0
    !
    SumAvgPIGapTot = 0.d0
    SumAvgPIGapOnPath = 0.d0
    SumAvgPIGapNotOnPath = 0.d0
    SumAvgPIGapNotBRAllStates = 0.d0
    SumAvgPIGapNotBRonPath = 0.d0
    SumAvgPIGapNotEqAllStates = 0.d0
    SumAvgPIGapNotEqonPath = 0.d0
    NumAvgPIGapTot = 0
    NumAvgPIGapOnPath = 0
    NumAvgPIGapNotOnPath = 0
    NumAvgPIGapNotBRAllStates = 0
    NumAvgPIGapNotBRonPath = 0
    NumAvgPIGapNotEqAllStates = 0
    NumAvgPIGapNotEqonPath = 0
    !
    ! Reading strategies and states at convergence from file
    !
    CALL ReadInfoModel(converged,timeToConvergence, & 
        CycleLength,CycleStates,CyclePrices,CycleProfits,indexStrategies)
    !
    ! Beginning loop over games
    !
    !$omp parallel do &
    !$omp private(QTrue,MaxQTrue,QGap,AvgVisitedProfits,AvgMaxProfits,AvgPIGap,AvgProfitGainGap,PathCycleStates, &
    !$omp   IsOnPath,IsBRAllStates,IsBRonPath,IsEqAllStates,IsEqonPath,IsBR,OptimalStrategyVec,LastStateVec, &
    !$omp   OptimalStrategy,lastObservedStateNumber,iState,iAgent, &
    !$omp   p,VisitedProfits,pPrime,OptimalPrice,VisitedStates,MaxProfits,iPeriod, &
    !$omp   PossibleActions,iPrice,PriceVector,PreCycleLength,GameCycleLength,iThres,tmp,i) &
    !$omp firstprivate(numGames,PI) &
    !$omp reduction(+ : SumAvgPIGapTot,SumAvgPIGapOnPath,SumAvgPIGapNotOnPath,SumAvgPIGapNotBRAllStates, &
    !$omp   SumAvgPIGapNotBRonPath,SumAvgPIGapNotEqAllStates,SumAvgPIGapNotEqonPath,NumAvgPIGapTot, &
    !$omp   NumAvgPIGapOnPath,NumAvgPIGapNotOnPath,NumAvgPIGapNotBRAllStates,NumAvgPIGapNotBRonPath, &
    !$omp   NumAvgPIGapNotEqAllStates,NumAvgPIGapNotEqonPath)
    !
    DO iGame = 1, numGames                  ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        ! Read strategy and last observed state at convergence from file
        !
        !$omp critical
        OptimalStrategyVec = indexStrategies(:,iGame)
        !$omp end critical
        !
        OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        lastObservedStateNumber = computeStateNumber(RESHAPE(LastStateVec, (/ DepthState,numAgents /) ))
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute true Q for the optimal strategy for all agents, in all states and actions
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        QTrue = 0.d0
        MaxQTrue = 0.d0
        QGap = 0.d0
        AvgVisitedProfits = 0.d0
        AvgMaxProfits = 0.d0
        AvgPIGap = 0.d0
        AvgProfitGainGap = 0.d0
        PathCycleStates = 0
        !
        DO iState = 1, numStates                ! Start of loop over states
            !
            DO iAgent = 1, numAgents            ! Start of loop over agents
                !
                pPrime = OptimalStrategy(iState,:)
                OptimalPrice = pPrime(iAgent)
                !
                DO iPrice = 1, numPrices        ! Start of loop over prices
                    !
                    ! First period deviation from optimal strategy
                    !
                    ! 1. Compute state value function for the optimal strategy in (iState,iPrice)
                    !
                    CALL computeQCell(OptimalStrategy,iState,iPrice,iAgent,delta, &
                        QTrue(iState,iPrice,iAgent),VisitedStates,PreCycleLength,GameCycleLength)
                    !
                    ! 2. Check if state is on path
                    !
                    IF ((iPrice .EQ. OptimalPrice) .AND. (iState .EQ. lastObservedStateNumber)) THEN
                        !
                        PathCycleStates(:GameCycleLength) = VisitedStates(PreCycleLength+1:PreCycleLength+GameCycleLength)
                        PathCycleLength(iGame) = GameCycleLength
                        !
                    END IF
                    !
                END DO                          ! End of loop over prices
                !
                ! Compute gap in Q function values w.r.t. maximum
                !
                MaxQTrue(iState,iAgent) = MAXVAL(QTrue(iState,:,iAgent))
                QGap(iState,iAgent) = &
                    (MaxQTrue(iState,iAgent)-QTrue(iState,OptimalStrategy(iState,iAgent),iAgent))/ABS(MaxQTrue(iState,iAgent))
                !
            END DO                              ! End of loop over agents
            !
        END DO                                  ! End of loop over initial states
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute true AvgPI for the optimal strategy for all agents, in all states and actions
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iState = 1, numStates                ! Start of loop over states
            !
            DO iAgent = 1, numAgents            ! Start of loop over agents
                !
                p = RESHAPE(convertNumberBase(iState-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
                pPrime = OptimalStrategy(iState,:)
                OptimalPrice = pPrime(iAgent)
                !
                ! 1. Compute pre-cycle and cycle data, for all states and for OptimalStrategy prices
                !
                VisitedStates = 0
                VisitedProfits = 0.d0
                MaxProfits = 0.d0
                !
                DO iPeriod = 1, numPeriods  ! Start of loop over future dates
                    !
                    IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                    p(1,:) = pPrime
                    VisitedStates(iPeriod) = computeStateNumber(p)
                    VisitedProfits(iPeriod) = PI(computeActionNumber(pPrime),iAgent)
                    !
                    ! Compute the highest possible profit given the opponents' prices
                    !
                    PossibleActions = 0
                    DO iPrice = 1, numPrices
                        !
                        PriceVector = pPrime
                        PriceVector(iAgent) = iPrice
                        PossibleActions(iPrice) = computeActionNumber(PriceVector)
                        !
                    END DO
                    MaxProfits(iPeriod) = MAXVAL(PI(PossibleActions,iAgent))
                    !
                    ! Check if the state has already been visited
                    !
                    IF ((iPeriod .GE. 2) .AND. (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) THEN
                        !
                        PreCycleLength = MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
                        GameCycleLength = iPeriod-PreCycleLength
                        EXIT
                        !
                    END IF
                    !
                    ! After period 1, every agent follows the optimal strategy
                    !
                    pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
                    !
                END DO                      ! End of loop over future dates
                !
                ! 2. Compute AvgPI function value for (state,agent) triplet
                !
                AvgVisitedProfits(iState,iAgent) = SUM(VisitedProfits(PreCycleLength+1:iPeriod))/DBLE(GameCycleLength)
                AvgMaxProfits(iState,iAgent) = SUM(MaxProfits(PreCycleLength+1:iPeriod))/DBLE(GameCycleLength)
                !
                ! Compute gap in AvgPI function values w.r.t. maximum
                !
                AvgPIGap(iState,iAgent) = AvgMaxProfits(iState,iAgent)-AvgVisitedProfits(iState,iAgent)
                AvgProfitGainGap(iState,iAgent) = AvgPIGap(iState,iAgent)/(CoopProfits(iAgent)-NashProfits(iAgent))
                !
            END DO                              ! End of loop over agents
            !
        END DO                                  ! End of loop over initial states
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute mask matrices
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        IsOnPath = .FALSE.
        IsBRAllStates = .FALSE.
        IsBRonPath = .FALSE.
        IsEqAllStates = .FALSE.
        IsEqonPath = .FALSE.
        !
        DO iState = 1, numStates                ! Start of loop over states
            !
            IF (ANY(PathCycleStates(:PathCycleLength(iGame)) .EQ. iState)) IsOnPath(iState,:) = .TRUE.
            IsBR = .FALSE.
            !
            DO iAgent = 1, numAgents            ! Start of loop over agents
                !
                OptimalPrice = OptimalStrategy(iState,iAgent)
                IF (ABS(QTrue(iState,OptimalPrice,iAgent)-MaxQTrue(iState,iAgent)) .LE. EPSILON(MaxQTrue(iState,iAgent))) THEN
                    !
                    IsBR(iAgent) = .TRUE.
                    IsBRAllStates(iState,iAgent) = .TRUE.
                    IF (ANY(PathCycleStates(:PathCycleLength(iGame)) .EQ. iState)) IsBRonPath(iState,iAgent) = .TRUE.
                    !
                END IF
                !
            END DO
            !
            DO iAgent = 1, numAgents            ! Start of loop over agents
                !
                OptimalPrice = OptimalStrategy(iState,iAgent)
                IF (ALL(IsBR)) THEN
                    !
                    IsEqAllStates(iState,iAgent) = .TRUE.
                    IF (ANY(PathCycleStates(:PathCycleLength(iGame)) .EQ. iState)) &
                        IsEqonPath(iState,iAgent) = .TRUE.
                    !
                END IF
                !
            END DO                          ! End of loop over agents
            !
        END DO                              ! End of loop over states
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Computing averages and descriptive statistics
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Summing by agent and threshold
        !
        iThres = MIN(PathCycleLength(iGame),ThresPathCycleLength(numThresPathCycleLength))
        DO iAgent = 1, numAgents
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent))
            SumAvgPIGapTot(0,0) = SumAvgPIGapTot(0,0)+tmp
            SumAvgPIGapTot(0,iAgent) = SumAvgPIGapTot(0,iAgent)+tmp
            SumAvgPIGapTot(iThres,0) = SumAvgPIGapTot(iThres,0)+tmp
            SumAvgPIGapTot(iThres,iAgent) = SumAvgPIGapTot(iThres,iAgent)+tmp
            NumAvgPIGapTot(0,0) = NumAvgPIGapTot(0,0)+numStates
            NumAvgPIGapTot(0,iAgent) = NumAvgPIGapTot(0,iAgent)+numStates
            NumAvgPIGapTot(iThres,0) = NumAvgPIGapTot(iThres,0)+numStates
            NumAvgPIGapTot(iThres,iAgent) = NumAvgPIGapTot(iThres,iAgent)+numStates
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent),MASK = IsOnPath(:,iAgent))
            SumAvgPIGapOnPath(0,0) = SumAvgPIGapOnPath(0,0)+tmp
            SumAvgPIGapOnPath(0,iAgent) = SumAvgPIGapOnPath(0,iAgent)+tmp
            SumAvgPIGapOnPath(iThres,0) = SumAvgPIGapOnPath(iThres,0)+tmp
            SumAvgPIGapOnPath(iThres,iAgent) = SumAvgPIGapOnPath(iThres,iAgent)+tmp
            NumAvgPIGapOnPath(0,0) = NumAvgPIGapOnPath(0,0)+COUNT(IsOnPath(:,iAgent))
            NumAvgPIGapOnPath(0,iAgent) = NumAvgPIGapOnPath(0,iAgent)+COUNT(IsOnPath(:,iAgent))
            NumAvgPIGapOnPath(iThres,0) = NumAvgPIGapOnPath(iThres,0)+COUNT(IsOnPath(:,iAgent))
            NumAvgPIGapOnPath(iThres,iAgent) = NumAvgPIGapOnPath(iThres,iAgent)+COUNT(IsOnPath(:,iAgent))
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent),MASK = .NOT.(IsOnPath(:,iAgent)))
            SumAvgPIGapNotOnPath(0,0) = SumAvgPIGapNotOnPath(0,0)+tmp
            SumAvgPIGapNotOnPath(0,iAgent) = SumAvgPIGapNotOnPath(0,iAgent)+tmp
            SumAvgPIGapNotOnPath(iThres,0) = SumAvgPIGapNotOnPath(iThres,0)+tmp
            SumAvgPIGapNotOnPath(iThres,iAgent) = SumAvgPIGapNotOnPath(iThres,iAgent)+tmp
            NumAvgPIGapNotOnPath(0,0) = NumAvgPIGapNotOnPath(0,0)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            NumAvgPIGapNotOnPath(0,iAgent) = NumAvgPIGapNotOnPath(0,iAgent)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            NumAvgPIGapNotOnPath(iThres,0) = NumAvgPIGapNotOnPath(iThres,0)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            NumAvgPIGapNotOnPath(iThres,iAgent) = NumAvgPIGapNotOnPath(iThres,iAgent)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent),MASK = .NOT.(IsBRAllStates(:,iAgent)))
            SumAvgPIGapNotBRAllStates(0,0) = SumAvgPIGapNotBRAllStates(0,0)+tmp
            SumAvgPIGapNotBRAllStates(0,iAgent) = SumAvgPIGapNotBRAllStates(0,iAgent)+tmp
            SumAvgPIGapNotBRAllStates(iThres,0) = SumAvgPIGapNotBRAllStates(iThres,0)+tmp
            SumAvgPIGapNotBRAllStates(iThres,iAgent) = SumAvgPIGapNotBRAllStates(iThres,iAgent)+tmp
            NumAvgPIGapNotBRAllStates(0,0) = NumAvgPIGapNotBRAllStates(0,0)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            NumAvgPIGapNotBRAllStates(0,iAgent) = NumAvgPIGapNotBRAllStates(0,iAgent)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            NumAvgPIGapNotBRAllStates(iThres,0) = NumAvgPIGapNotBRAllStates(iThres,0)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            NumAvgPIGapNotBRAllStates(iThres,iAgent) = NumAvgPIGapNotBRAllStates(iThres,iAgent)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent),MASK = .NOT.(IsBRonPath(:,iAgent)))
            SumAvgPIGapNotBRonPath(0,0) = SumAvgPIGapNotBRonPath(0,0)+tmp
            SumAvgPIGapNotBRonPath(0,iAgent) = SumAvgPIGapNotBRonPath(0,iAgent)+tmp
            SumAvgPIGapNotBRonPath(iThres,0) = SumAvgPIGapNotBRonPath(iThres,0)+tmp
            SumAvgPIGapNotBRonPath(iThres,iAgent) = SumAvgPIGapNotBRonPath(iThres,iAgent)+tmp
            NumAvgPIGapNotBRonPath(0,0) = NumAvgPIGapNotBRonPath(0,0)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            NumAvgPIGapNotBRonPath(0,iAgent) = NumAvgPIGapNotBRonPath(0,iAgent)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            NumAvgPIGapNotBRonPath(iThres,0) = NumAvgPIGapNotBRonPath(iThres,0)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            NumAvgPIGapNotBRonPath(iThres,iAgent) = NumAvgPIGapNotBRonPath(iThres,iAgent)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent),MASK = .NOT.(IsEqAllStates(:,iAgent)))
            SumAvgPIGapNotEqAllStates(0,0) = SumAvgPIGapNotEqAllStates(0,0)+tmp
            SumAvgPIGapNotEqAllStates(0,iAgent) = SumAvgPIGapNotEqAllStates(0,iAgent)+tmp
            SumAvgPIGapNotEqAllStates(iThres,0) = SumAvgPIGapNotEqAllStates(iThres,0)+tmp
            SumAvgPIGapNotEqAllStates(iThres,iAgent) = SumAvgPIGapNotEqAllStates(iThres,iAgent)+tmp
            NumAvgPIGapNotEqAllStates(0,0) = NumAvgPIGapNotEqAllStates(0,0)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            NumAvgPIGapNotEqAllStates(0,iAgent) = NumAvgPIGapNotEqAllStates(0,iAgent)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            NumAvgPIGapNotEqAllStates(iThres,0) = NumAvgPIGapNotEqAllStates(iThres,0)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            NumAvgPIGapNotEqAllStates(iThres,iAgent) = NumAvgPIGapNotEqAllStates(iThres,iAgent)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            !
            tmp = SUM(AvgProfitGainGap(:,iAgent),MASK = .NOT.(IsEqonPath(:,iAgent)))
            SumAvgPIGapNotEqonPath(0,0) = SumAvgPIGapNotEqonPath(0,0)+tmp
            SumAvgPIGapNotEqonPath(0,iAgent) = SumAvgPIGapNotEqonPath(0,iAgent)+tmp
            SumAvgPIGapNotEqonPath(iThres,0) = SumAvgPIGapNotEqonPath(iThres,0)+tmp
            SumAvgPIGapNotEqonPath(iThres,iAgent) = SumAvgPIGapNotEqonPath(iThres,iAgent)+tmp
            NumAvgPIGapNotEqonPath(0,0) = NumAvgPIGapNotEqonPath(0,0)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            NumAvgPIGapNotEqonPath(0,iAgent) = NumAvgPIGapNotEqonPath(0,iAgent)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            NumAvgPIGapNotEqonPath(iThres,0) = NumAvgPIGapNotEqonPath(iThres,0)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            NumAvgPIGapNotEqonPath(iThres,iAgent) = NumAvgPIGapNotEqonPath(iThres,iAgent)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            !
        END DO                              ! End of loop over agents
        !
    END DO                                  ! End of loop over games
    !$omp end parallel do
    !
    ! Averaging
    !
    SumAvgPIGapTot = SumAvgPIGapTot/DBLE(NumAvgPIGapTot)
    WHERE (ISNAN(SumAvgPIGapTot)) SumAvgPIGapTot = -999.999d0
    SumAvgPIGapOnPath = SumAvgPIGapOnPath/DBLE(NumAvgPIGapOnPath)
    WHERE (ISNAN(SumAvgPIGapOnPath)) SumAvgPIGapOnPath = -999.999d0 
    SumAvgPIGapNotOnPath = SumAvgPIGapNotOnPath/DBLE(NumAvgPIGapNotOnPath)
    WHERE (ISNAN(SumAvgPIGapNotOnPath)) SumAvgPIGapNotOnPath = -999.999d0
    SumAvgPIGapNotBRAllStates = SumAvgPIGapNotBRAllStates/DBLE(NumAvgPIGapNotBRAllStates)
    WHERE (ISNAN(SumAvgPIGapNotBRAllStates)) SumAvgPIGapNotBRAllStates = -999.999d0
    SumAvgPIGapNotBRonPath = SumAvgPIGapNotBRonPath/DBLE(NumAvgPIGapNotBRonPath)
    WHERE (ISNAN(SumAvgPIGapNotBRonPath)) SumAvgPIGapNotBRonPath = -999.999d0
    SumAvgPIGapNotEqAllStates = SumAvgPIGapNotEqAllStates/DBLE(NumAvgPIGapNotEqAllStates)
    WHERE (ISNAN(SumAvgPIGapNotEqAllStates)) SumAvgPIGapNotEqAllStates = -999.999d0
    SumAvgPIGapNotEqonPath = SumAvgPIGapNotEqonPath/DBLE(NumAvgPIGapNotEqonPath)
    WHERE (ISNAN(SumAvgPIGapNotEqonPath)) SumAvgPIGapNotEqonPath = -999.999d0
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10007,1) &
            (i, i = 1, numAgents), &
            (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, (j, i, j = 1, 2), i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (i, i, i, i, i, i, i, i = 1, numAgents), &
            (j, j, j, j, j, j, j, &
                (j, i, j, i, j, i, j, i, j, i, j, i, j, i, i = 1, numAgents), j = 1, numThresPathCycleLength)
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), <numAgents>('    delta', I1, ' '), &
            <numAgents>('typeQini', I1, ' ', <2>('par', I1, 'Qini', I1, ' ')), &
            <numDemandParameters>('  DemPar', I0.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            '   AvgPGGapTot AvgPGGaponPath AvgPGGapNotOnPath ', &
                'AvgPGGapNotBRAllSt AvgPGGapNotBROnPath AvgPGGapNotEqAllSt AvgPGGapNotEqOnPath ', &
            <numAgents>('AvgPGGapTotAg', I1, ' AvgPGGaponPathAg', I1, ' AvgPGGapNotOnPathAg', I1, ' ', &
                'AvgPGGapNotBRAllStAg', I1, ' AvgPGGapNotBROnPathAg', I1, ' AvgPGGapNotEqAllStAg', I1, ' AvgPGGapNotEqOnPathAg', I1, ' '), &
            <numThresPathCycleLength>('PL', I2.2, 'AvgPGGapTot PL', I2.2, 'AvgPGGaponPath PL', I2.2, 'AvgPGGapNotOnPath ', &
                'PL', I2.2, 'AvgPGGapNotBRAllSt PL', I2.2, 'AvgPGGapNotBROnPath PL', I2.2, 'AvgPGGapNotEqAllSt PL', I2.2, 'AvgPGGapNotEqOnPath ', &
            <numAgents>('PL', I2.2, 'AvgPGGapTotAg', I1, ' PL', I2.2, 'AvgPGGaponPathAg', I1, ' PL', I2.2, 'AvgPGGapNotOnPathAg', I1, ' ', &
                'PL', I2.2, 'AvgPGGapNotBRAllStAg', I1, ' PL', I2.2, 'AvgPGGapNotBROnPathAg', I1, ' PL', I2.2, 'AvgPGGapNotEqAllStAg', I1, ' PL', I2.2, 'AvgPGGapNotEqOnPathAg', I1, ' ') )&
            )
        !
    END IF
    !
    WRITE(10007,2) codModel, &
        alpha, MExpl, delta, &
        (typeQInitialization(i), parQInitialization(i, :), i = 1, numAgents), &
        DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        SumAvgPIGapTot(0,0), SumAvgPIGapOnPath(0,0), SumAvgPIGapNotOnPath(0,0), &
            SumAvgPIGapNotBRAllStates(0,0), SumAvgPIGapNotBRonPath(0,0), SumAvgPIGapNotEqAllStates(0,0), SumAvgPIGapNotEqonPath(0,0), &
        (SumAvgPIGapTot(0,i), SumAvgPIGapOnPath(0,i), SumAvgPIGapNotOnPath(0,i), &
            SumAvgPIGapNotBRAllStates(0,i), SumAvgPIGapNotBRonPath(0,i), SumAvgPIGapNotEqAllStates(0,i), SumAvgPIGapNotEqonPath(0,i), i = 1, numAgents), &
        (SumAvgPIGapTot(j,0), SumAvgPIGapOnPath(j,0), SumAvgPIGapNotOnPath(j,0), &
            SumAvgPIGapNotBRAllStates(j,0), SumAvgPIGapNotBRonPath(j,0), SumAvgPIGapNotEqAllStates(j,0), SumAvgPIGapNotEqonPath(j,0), &
        (SumAvgPIGapTot(j,i), SumAvgPIGapOnPath(j,i), SumAvgPIGapNotOnPath(j,i), &
            SumAvgPIGapNotBRAllStates(j,i), SumAvgPIGapNotBRonPath(j,i), SumAvgPIGapNotEqAllStates(j,i), SumAvgPIGapNotEqonPath(j,i), i = 1, numAgents), j = 1, numThresPathCycleLength)
2   FORMAT(I5, 1X, &
        <3*numAgents>(F10.5, 1X), &
        <numAgents>(A9, 1X, <2>(F9.2, 1X)), &
        <numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        F14.7, 1X, F14.7, 1X, F17.7, 1X, F18.7, 1X, F19.7, 1X, F18.7, 1X, F19.7, 1X, &
        <numAgents>(F14.7, 1X, F17.7, 1X, F20.7, 1X, F21.7, 1X, F22.7, 1X, F21.7, 1X, F22.7, 1X), &
        <numThresPathCycleLength>(F15.7, 1X, F18.7, 1X, F21.7, 1X, F22.7, 1X, F23.7, 1X, F22.7, 1X, F23.7, 1X, &
            <numAgents>(F18.7, 1X, F21.7, 1X, F24.7, 1X, F25.7, 1X, F26.7, 1X, F25.7, 1X, F26.7, 1X)) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeAvgPIGapToMax
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE PIGapToMaximum
