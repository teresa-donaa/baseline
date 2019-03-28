MODULE QGapToMaximum
!
USE globals
USE QL_routines
!
! Computes gap in Q function values w.r.t. maximum
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeQGapToMax ( iModel )
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
    INTEGER :: numPeriods, iPeriod, iGame, iState, iAgent, iPrice, iThres, i, j, &
        PathCycleStates(numStates+1), PathCycleLength(numGames), &
        OptimalStrategy(numStates,numAgents), lastObservedStateNumber, &
        p(DepthState,numAgents), pPrime(numAgents), OptimalPrice, &
        VisitedStates(numStates+1), PreCycleLength, CycleLength
    INTEGER, DIMENSION(0:numThresPathCycleLength,0:numAgents) :: NumQGapTot, NumQGapOnPath, &
        NumQGapNotOnPath, NumQGapNotBRAllStates, NumQGapNotBRonPath, NumQGapNotEqAllStates, NumQGapNotEqonPath
    REAL(8) :: DiscountFactors(0:numStates,numAgents), &
        QTrue(numStates,numPrices,numAgents), VisitedProfits(numStates+1), &
        PreCycleProfit, CycleProfit, QGap(numStates,numAgents), MaxQTrue(numStates,numAgents), tmp
    REAL(8), DIMENSION(0:numThresPathCycleLength,0:numAgents) :: SumQGapTot, SumQGapOnPath, &
        SumQGapNotOnPath, SumQGapNotBRAllStates, SumQGapNotBRonPath, SumQGapNotEqAllStates, SumQGapNotEqonPath
    LOGICAL, DIMENSION(numStates,numAgents) :: IsOnPath, IsBRAllStates, IsBRonPath, IsEqAllStates, IsEqonPath
    LOGICAL, DIMENSION(numAgents) :: IsBR
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Q gaps'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    numPeriods = numStates+1        ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
    DiscountFactors = TRANSPOSE(RESHAPE((/ (delta**iPeriod, iPeriod = 0, numPeriods-1) /),(/ numAgents,numPeriods /)))
    PathCycleLength = 0
    !
    SumQGapTot = 0.d0
    SumQGapOnPath = 0.d0
    SumQGapNotOnPath = 0.d0
    SumQGapNotBRAllStates = 0.d0
    SumQGapNotBRonPath = 0.d0
    SumQGapNotEqAllStates = 0.d0
    SumQGapNotEqonPath = 0.d0
    NumQGapTot = 0
    NumQGapOnPath = 0
    NumQGapNotOnPath = 0
    NumQGapNotBRAllStates = 0
    NumQGapNotBRonPath = 0
    NumQGapNotEqAllStates = 0
    NumQGapNotEqonPath = 0
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
    !$omp private(QTrue,MaxQTrue,QGap,IsOnPath,IsBRAllStates,IsBRonPath,IsEqAllStates,IsEqonPath, &
    !$omp   IsBR,iThres,tmp,OptimalStrategy,lastObservedStateNumber,iAgent,PathCycleStates, &
    !$omp   iState,p,pPrime,OptimalPrice,iPrice,VisitedStates,VisitedProfits,iPeriod, &
    !$omp   PreCycleLength,CycleLength,PreCycleProfit,CycleProfit,OptimalStrategyVec,LastStateVec,i) &
    !$omp firstprivate(numGames,PI) &
    !$omp reduction(+ : SumQGapTot,SumQGapOnPath,SumQGapNotOnPath,SumQGapNotBRAllStates,SumQGapNotBRonPath, &
    !$omp   SumQGapNotEqAllStates,SumQGapNotEqonPath,NumQGapTot,NumQGapOnPath,NumQGapNotOnPath,NumQGapNotBRAllStates, &
    !$omp   NumQGapNotBRonPath,NumQGapNotEqAllStates,NumQGapNotEqonPath)
    !
    DO iGame = 1, numGames                  ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        ! Read strategy and last observed state at convergence from file
        !
        !$omp critical
        OptimalStrategyVec = indexStrategies(:,iGame)
        LastStateVec = indexLastState(:,iGame)
        !$omp end critical
        !
        optimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        lastObservedStateNumber = computeStateNumber(RESHAPE(LastStateVec, (/ DepthState,numAgents /) ))
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute true Q for the optimal strategy for all agents, in all states and actions
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        QTrue = 0.d0
        MaxQTrue = 0.d0
        QGap = 0.d0
        PathCycleStates = 0
        !
        DO iState = 1, numStates                ! Start of loop over states
            !
            DO iAgent = 1, numAgents            ! Start of loop over agents
                !
                p = RESHAPE(convertNumberBase(iState-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
                pPrime = OptimalStrategy(iState,:)
                OptimalPrice = pPrime(iAgent)
                !
                DO iPrice = 1, numPrices        ! Start of loop over prices
                    !
                    ! First period deviation from optimal strategy
                    !
                    pPrime(iAgent) = iPrice
                    !
                    ! 1. Compute pre-cycle and cycle data
                    !
                    VisitedStates = 0
                    VisitedProfits = 0.d0
                    !
                    DO iPeriod = 1, numPeriods  ! Start of loop over future dates
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
                    END DO                      ! End of loop over future dates
                    !
                    ! 2. Compute Q function value for (state,price,agent) triplet
                    !
                    IF ((iPrice .EQ. OptimalPrice) .AND. (iState .EQ. lastObservedStateNumber)) THEN
                        !
                        PathCycleStates(:CycleLength) = VisitedStates(PreCycleLength+1:iPeriod)
                        PathCycleLength(iGame) = CycleLength
                        !
                    END IF
                    PreCycleProfit = SUM(DiscountFactors(0:PreCycleLength-1,iAgent)*VisitedProfits(:PreCycleLength))
                    CycleProfit = SUM(DiscountFactors(0:CycleLength-1,iAgent)*VisitedProfits(PreCycleLength+1:iPeriod))
                    QTrue(iState,iPrice,iAgent) = &
                        PreCycleProfit+delta(iAgent)**PreCycleLength*CycleProfit/(1.d0-delta(iAgent)**CycleLength)
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
            tmp = SUM(QGap(:,iAgent))
            SumQGapTot(0,0) = SumQGapTot(0,0)+tmp
            SumQGapTot(0,iAgent) = SumQGapTot(0,iAgent)+tmp
            SumQGapTot(iThres,0) = SumQGapTot(iThres,0)+tmp
            SumQGapTot(iThres,iAgent) = SumQGapTot(iThres,iAgent)+tmp
            NumQGapTot(0,0) = NumQGapTot(0,0)+numStates
            NumQGapTot(0,iAgent) = NumQGapTot(0,iAgent)+numStates
            NumQGapTot(iThres,0) = NumQGapTot(iThres,0)+numStates
            NumQGapTot(iThres,iAgent) = NumQGapTot(iThres,iAgent)+numStates
            !
            tmp = SUM(QGap(:,iAgent),MASK = IsOnPath(:,iAgent))
            SumQGapOnPath(0,0) = SumQGapOnPath(0,0)+tmp
            SumQGapOnPath(0,iAgent) = SumQGapOnPath(0,iAgent)+tmp
            SumQGapOnPath(iThres,0) = SumQGapOnPath(iThres,0)+tmp
            SumQGapOnPath(iThres,iAgent) = SumQGapOnPath(iThres,iAgent)+tmp
            NumQGapOnPath(0,0) = NumQGapOnPath(0,0)+COUNT(IsOnPath(:,iAgent))
            NumQGapOnPath(0,iAgent) = NumQGapOnPath(0,iAgent)+COUNT(IsOnPath(:,iAgent))
            NumQGapOnPath(iThres,0) = NumQGapOnPath(iThres,0)+COUNT(IsOnPath(:,iAgent))
            NumQGapOnPath(iThres,iAgent) = NumQGapOnPath(iThres,iAgent)+COUNT(IsOnPath(:,iAgent))
            !
            tmp = SUM(QGap(:,iAgent),MASK = .NOT.(IsOnPath(:,iAgent)))
            SumQGapNotOnPath(0,0) = SumQGapNotOnPath(0,0)+tmp
            SumQGapNotOnPath(0,iAgent) = SumQGapNotOnPath(0,iAgent)+tmp
            SumQGapNotOnPath(iThres,0) = SumQGapNotOnPath(iThres,0)+tmp
            SumQGapNotOnPath(iThres,iAgent) = SumQGapNotOnPath(iThres,iAgent)+tmp
            NumQGapNotOnPath(0,0) = NumQGapNotOnPath(0,0)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            NumQGapNotOnPath(0,iAgent) = NumQGapNotOnPath(0,iAgent)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            NumQGapNotOnPath(iThres,0) = NumQGapNotOnPath(iThres,0)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            NumQGapNotOnPath(iThres,iAgent) = NumQGapNotOnPath(iThres,iAgent)+COUNT(.NOT.(IsOnPath(:,iAgent)))
            !
            tmp = SUM(QGap(:,iAgent),MASK = .NOT.(IsBRAllStates(:,iAgent)))
            SumQGapNotBRAllStates(0,0) = SumQGapNotBRAllStates(0,0)+tmp
            SumQGapNotBRAllStates(0,iAgent) = SumQGapNotBRAllStates(0,iAgent)+tmp
            SumQGapNotBRAllStates(iThres,0) = SumQGapNotBRAllStates(iThres,0)+tmp
            SumQGapNotBRAllStates(iThres,iAgent) = SumQGapNotBRAllStates(iThres,iAgent)+tmp
            NumQGapNotBRAllStates(0,0) = NumQGapNotBRAllStates(0,0)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            NumQGapNotBRAllStates(0,iAgent) = NumQGapNotBRAllStates(0,iAgent)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            NumQGapNotBRAllStates(iThres,0) = NumQGapNotBRAllStates(iThres,0)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            NumQGapNotBRAllStates(iThres,iAgent) = NumQGapNotBRAllStates(iThres,iAgent)+COUNT(.NOT.(IsBRAllStates(:,iAgent)))
            !
            tmp = SUM(QGap(:,iAgent),MASK = .NOT.(IsBRonPath(:,iAgent)))
            SumQGapNotBRonPath(0,0) = SumQGapNotBRonPath(0,0)+tmp
            SumQGapNotBRonPath(0,iAgent) = SumQGapNotBRonPath(0,iAgent)+tmp
            SumQGapNotBRonPath(iThres,0) = SumQGapNotBRonPath(iThres,0)+tmp
            SumQGapNotBRonPath(iThres,iAgent) = SumQGapNotBRonPath(iThres,iAgent)+tmp
            NumQGapNotBRonPath(0,0) = NumQGapNotBRonPath(0,0)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            NumQGapNotBRonPath(0,iAgent) = NumQGapNotBRonPath(0,iAgent)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            NumQGapNotBRonPath(iThres,0) = NumQGapNotBRonPath(iThres,0)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            NumQGapNotBRonPath(iThres,iAgent) = NumQGapNotBRonPath(iThres,iAgent)+COUNT(.NOT.(IsBRonPath(:,iAgent)))
            !
            tmp = SUM(QGap(:,iAgent),MASK = .NOT.(IsEqAllStates(:,iAgent)))
            SumQGapNotEqAllStates(0,0) = SumQGapNotEqAllStates(0,0)+tmp
            SumQGapNotEqAllStates(0,iAgent) = SumQGapNotEqAllStates(0,iAgent)+tmp
            SumQGapNotEqAllStates(iThres,0) = SumQGapNotEqAllStates(iThres,0)+tmp
            SumQGapNotEqAllStates(iThres,iAgent) = SumQGapNotEqAllStates(iThres,iAgent)+tmp
            NumQGapNotEqAllStates(0,0) = NumQGapNotEqAllStates(0,0)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            NumQGapNotEqAllStates(0,iAgent) = NumQGapNotEqAllStates(0,iAgent)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            NumQGapNotEqAllStates(iThres,0) = NumQGapNotEqAllStates(iThres,0)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            NumQGapNotEqAllStates(iThres,iAgent) = NumQGapNotEqAllStates(iThres,iAgent)+COUNT(.NOT.(IsEqAllStates(:,iAgent)))
            !
            tmp = SUM(QGap(:,iAgent),MASK = .NOT.(IsEqonPath(:,iAgent)))
            SumQGapNotEqonPath(0,0) = SumQGapNotEqonPath(0,0)+tmp
            SumQGapNotEqonPath(0,iAgent) = SumQGapNotEqonPath(0,iAgent)+tmp
            SumQGapNotEqonPath(iThres,0) = SumQGapNotEqonPath(iThres,0)+tmp
            SumQGapNotEqonPath(iThres,iAgent) = SumQGapNotEqonPath(iThres,iAgent)+tmp
            NumQGapNotEqonPath(0,0) = NumQGapNotEqonPath(0,0)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            NumQGapNotEqonPath(0,iAgent) = NumQGapNotEqonPath(0,iAgent)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            NumQGapNotEqonPath(iThres,0) = NumQGapNotEqonPath(iThres,0)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            NumQGapNotEqonPath(iThres,iAgent) = NumQGapNotEqonPath(iThres,iAgent)+COUNT(.NOT.(IsEqonPath(:,iAgent)))
            !
        END DO                              ! End of loop over agents
        !
    END DO                                  ! End of loop over games
    !$omp end parallel do
    !
    ! Averaging
    !
    SumQGapTot = SumQGapTot/DBLE(NumQGapTot)
    WHERE (ISNAN(SumQGapTot)) SumQGapTot = -999.999d0
    SumQGapOnPath = SumQGapOnPath/DBLE(NumQGapOnPath)
    WHERE (ISNAN(SumQGapOnPath)) SumQGapOnPath = -999.999d0
    SumQGapNotOnPath = SumQGapNotOnPath/DBLE(NumQGapNotOnPath)
    WHERE (ISNAN(SumQGapNotOnPath)) SumQGapNotOnPath = -999.999d0
    SumQGapNotBRAllStates = SumQGapNotBRAllStates/DBLE(NumQGapNotBRAllStates)
    WHERE (ISNAN(SumQGapNotBRAllStates)) SumQGapNotBRAllStates = -999.999d0
    SumQGapNotBRonPath = SumQGapNotBRonPath/DBLE(NumQGapNotBRonPath)
    WHERE (ISNAN(SumQGapNotBRonPath)) SumQGapNotBRonPath = -999.999d0
    SumQGapNotEqAllStates = SumQGapNotEqAllStates/DBLE(NumQGapNotEqAllStates)
    WHERE (ISNAN(SumQGapNotEqAllStates)) SumQGapNotEqAllStates = -999.999d0
    SumQGapNotEqonPath = SumQGapNotEqonPath/DBLE(NumQGapNotEqonPath)
    WHERE (ISNAN(SumQGapNotEqonPath)) SumQGapNotEqonPath = -999.999d0
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10006,1) (i, i = 1, numAgents) &
            , (i, i = 1, numExplorationParameters), (i, i = 1, numAgents) &
            , (i, i = 1, numDemandParameters) &
            , (i, i = 1, numAgents), (i, i = 1, numAgents) &
            , (i, i = 1, numAgents), (i, i = 1, numAgents) &
            , (i, i = 1, numAgents), (i, i = 1, numAgents) &
            , ((i, j, j = 1, numPrices), i = 1, numAgents) &
            , (i, i, i, i, i, i, i, i = 1, numAgents) &
            , (j, j, j, j, j, j, j, &
                (j, i, j, i, j, i, j, i, j, i, j, i, j, i, i = 1, numAgents), j = 1, numThresPathCycleLength)
1       FORMAT('Model ' &
            , <numAgents>('    alpha', I1, ' ') &
            , <numExplorationParameters>(' MExplPar', I1, ' '), <numAgents>('    delta', I1, ' ') &
            , <numDemandParameters>('  DemPar', I2.2, ' ') &
            , <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' ') &
            , <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' ') &
            , <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' ') &
            , <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')) &
            , '   QGapTot QGaponPath QGapNotOnPath ' &
                , 'QGapNotBRAllSt QGapNotBROnPath QGapNotEqAllSt QGapNotEqOnPath ' &
            , <numAgents>('QGapTotAg', I1, ' QGaponPathAg', I1, ' QGapNotOnPathAg', I1, ' ', &
                'QGapNotBRAllStAg', I1, ' QGapNotBROnPathAg', I1, ' QGapNotEqAllStAg', I1, ' QGapNotEqOnPathAg', I1, ' ') &
            , <numThresPathCycleLength>('   PL', I2.2, 'QGapTot PL', I2.2, 'QGaponPath PL', I2.2, 'QGapNotOnPath ', &
                'PL', I2.2, 'QGapNotBRAllSt PL', I2.2, 'QGapNotBROnPath PL', I2.2, 'QGapNotEqAllSt PL', I2.2, 'QGapNotEqOnPath ', &
            <numAgents>('PL', I2.2, 'QGapTotAg', I1, ' PL', I2.2, 'QGaponPathAg', I1, ' PL', I2.2, 'QGapNotOnPathAg', I1, ' ', &
                'PL', I2.2, 'QGapNotBRAllStAg', I1, ' PL', I2.2, 'QGapNotBROnPathAg', I1, ' PL', I2.2, 'QGapNotEqAllStAg', I1, ' PL', I2.2, 'QGapNotEqOnPathAg', I1, ' ') ) &
            )
        !
    END IF
    !
    WRITE(10006,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        SumQGapTot(0,0), SumQGapOnPath(0,0), SumQGapNotOnPath(0,0), &
            SumQGapNotBRAllStates(0,0), SumQGapNotBRonPath(0,0), SumQGapNotEqAllStates(0,0), SumQGapNotEqonPath(0,0), &
        (SumQGapTot(0,i), SumQGapOnPath(0,i), SumQGapNotOnPath(0,i), &
            SumQGapNotBRAllStates(0,i), SumQGapNotBRonPath(0,i), SumQGapNotEqAllStates(0,i), SumQGapNotEqonPath(0,i), i = 1, numAgents), &
        (SumQGapTot(j,0), SumQGapOnPath(j,0), SumQGapNotOnPath(j,0), &
            SumQGapNotBRAllStates(j,0), SumQGapNotBRonPath(j,0), SumQGapNotEqAllStates(j,0), SumQGapNotEqonPath(j,0), &
        (SumQGapTot(j,i), SumQGapOnPath(j,i), SumQGapNotOnPath(j,i), &
            SumQGapNotBRAllStates(j,i), SumQGapNotBRonPath(j,i), SumQGapNotEqAllStates(j,i), SumQGapNotEqonPath(j,i), i = 1, numAgents), j = 1, numThresPathCycleLength)
2   FORMAT(I5, 1X &
        , <3*numAgents+numDemandParameters>(F10.3, 1X) &
        , <6*numAgents>(F10.3, 1X) &
        <numPrices*numAgents>(F10.5, 1X), &
        , F10.7, 1X, F10.7, 1X, F13.7, 1X, F14.7, 1X, F15.7, 1X, F14.7, 1X, F15.7, 1X &
        , <numAgents>(F10.7, 1X, F13.7, 1X, F16.7, 1X, F17.7, 1X, F18.7, 1X, F17.7, 1X, F18.7, 1X) &
        , <numThresPathCycleLength>(F14.7, 1X, F14.7, 1X, F17.7, 1X, F18.7, 1X, F19.7, 1X, F18.7, 1X, F19.7, 1X &
        ,    <numAgents>(F14.7, 1X, F17.7, 1X, F20.7, 1X, F21.7, 1X, F22.7, 1X, F21.7, 1X, F22.7, 1X)) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeQGapToMax
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE QGapToMaximum
