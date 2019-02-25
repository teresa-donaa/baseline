MODULE ImpulseResponseToNash
!
USE globals
USE QL_routines
!
! Computes analysis of the Impulse Response to a permanent or temporary deviation to the Nash price 
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeIRToNashAnalysis ( iModel )
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
    INTEGER, PARAMETER :: numThresPeriodsLength = 10
    INTEGER, PARAMETER :: numThresPeriodsLength0 = numThresPeriodsLength+1
    INTEGER, PARAMETER :: numShockPeriodsPrint = 25
    INTEGER :: ThresPeriodsLength(numThresPeriodsLength), ThresPeriodsLength0(numThresPeriodsLength0)
    INTEGER :: PeriodsLengthPre, PeriodsLengthShock, PeriodsLengthPost, PunishmentStrategy, &
        visitedStatesPre(numStates+1), visitedStates(MAX(numShockPeriodsPrint,numStates+1)), &
        p(DepthState,numAgents), pPrime(numAgents), &
        iPeriod, iAgent, jAgent, iPrice, tmp1(numAgents), &
        iGame, optimalStrategy(numStates,numAgents), LastObservedPrices(DepthState,numAgents), &
        indexShockState(LengthStates), numPeriods, iThres, i, j
    INTEGER :: FreqPeriodLengthPre(numThresPeriodsLength)
    INTEGER, DIMENSION(numAgents,numThresPeriodsLength) :: FreqPeriodLengthShock, FreqPeriodLengthPost
    INTEGER :: FreqPunishmentStrategy(numAgents,0:numThresPeriodsLength)
    REAL(8) :: visitedPrices(numStates+1,numAgents), visitedProfits(numStates+1,numAgents), &
        pNash
    REAL(8), DIMENSION(numAgents) :: avgPricesPre, avgProfitsPre, avgPricesPreQ, avgProfitsPreQ
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents,numAgents) :: &
        avgPricesShock, avgProfitsShock, avgPricesShockQ, avgProfitsShockQ, &
        avgPricesPercShock, avgProfitsPercShock, avgPricesPercShockQ, avgProfitsPercShockQ
    REAL(8), DIMENSION(numAgents,numAgents) :: &
        avgPricesPost, avgProfitsPost, avgPricesPostQ, avgProfitsPostQ, &
        avgPricesPercPost, avgProfitsPercPost, avgPricesPercPostQ, avgProfitsPercPostQ
    LOGICAL :: FlagReturnedToState
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    REAL(8) :: AggrPricesPre, AggrDevPricesPost, AggrNonDevPricesPost, AggrProfitsPre, AggrDevProfitsPost, AggrNonDevProfitsPost, &
        AggrDevPricesPercPost, AggrNonDevPricesPercPost, AggrDevProfitsPercPost, AggrNonDevProfitsPercPost
    REAL(8), DIMENSION(numShockPeriodsPrint) :: &
        AggrDevPricesShock, AggrNonDevPricesShock, AggrDevProfitsShock, AggrNonDevProfitsShock, &
        AggrDevPricesPercShock, AggrNonDevPricesPercShock, AggrDevProfitsPercShock, AggrNonDevProfitsPercShock
    REAL(8) :: AggrPricesPreQ, AggrDevPricesPostQ, AggrNonDevPricesPostQ, AggrProfitsPreQ, AggrDevProfitsPostQ, AggrNonDevProfitsPostQ, &
        AggrDevPricesPercPostQ, AggrNonDevPricesPercPostQ, AggrDevProfitsPercPostQ, AggrNonDevProfitsPercPostQ
    REAL(8), DIMENSION(numShockPeriodsPrint) :: &
        AggrDevPricesShockQ, AggrNonDevPricesShockQ, AggrDevProfitsShockQ, AggrNonDevProfitsShockQ, &
        AggrDevPricesPercShockQ, AggrNonDevPricesPercShockQ, AggrDevProfitsPercShockQ, AggrNonDevProfitsPercShockQ
    REAL(8), DIMENSION(numAgents) :: PricePre, ProfitPre
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Impulse Response functions to Nash'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !
    ThresPeriodsLength = (/ ( i, i = 1, numThresPeriodsLength ) /)
    ThresPeriodsLength0 = (/ 0, ThresPeriodsLength /)
    !
    numPeriods = numStates+1        ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
    FreqPeriodLengthPre = 0
    FreqPeriodLengthShock = 0
    FreqPeriodLengthPost = 0
    FreqPunishmentStrategy = 0
    avgPricesPre = 0.d0
    avgPricesShock = 0.d0
    avgPricesPost = 0.d0
    avgPricesPercShock = 0.d0
    avgPricesPercPost = 0.d0
    avgProfitsPre = 0.d0
    avgProfitsShock = 0.d0
    avgProfitsPost = 0.d0
    avgProfitsPercShock = 0.d0
    avgProfitsPercPost = 0.d0
    avgPricesPreQ = 0.d0
    avgPricesShockQ = 0.d0
    avgPricesPostQ = 0.d0
    avgPricesPercShockQ = 0.d0
    avgPricesPercPostQ = 0.d0
    avgProfitsPreQ = 0.d0
    avgProfitsShockQ = 0.d0
    avgProfitsPostQ = 0.d0
    avgProfitsPercShockQ = 0.d0
    avgProfitsPercPostQ = 0.d0
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
    !$omp private(OptimalStrategy,LastObservedPrices,visitedStatesPre,visitedPrices, &
    !$omp   visitedProfits,p,pPrime,iPeriod,iAgent,pNash,OptimalStrategyVec,LastStateVec, &
    !$omp   visitedStates,iPrice,tmp1,flagReturnedToState,jAgent,indexShockState,i, &
    !$omp   PricePre,ProfitPre) &
    !$omp firstprivate(PI,PricesGrids) &
    !$omp reduction(+ : FreqPeriodLengthPre, &
    !$omp   avgPricesPre,avgProfitsPre,avgPricesShock,avgProfitsShock,avgPricesPost,avgProfitsPost, &
    !$omp   avgPricesPreQ,avgProfitsPreQ,avgPricesShockQ,avgProfitsShockQ,avgPricesPostQ,avgProfitsPostQ, &
    !$omp   avgPricesPercShock,avgProfitsPercShock,avgPricesPercPost,avgProfitsPercPost, &
    !$omp   avgPricesPercShockQ,avgProfitsPercShockQ,avgPricesPercPostQ,avgProfitsPercPostQ, &
    !$omp   FreqPeriodLengthShock,FreqPunishmentStrategy,FreqPeriodLengthPost)
    DO iGame = 1, numGames        ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        !$omp critical
        OptimalStrategyVec = indexStrategies(:,iGame)
        LastStateVec = indexLastState(:,iGame)
        !$omp end critical
        !
        optimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        IF (DepthState0 .EQ. 0) THEN
            !
            LastObservedPrices = optimalStrategy
            !
        ELSE IF (DepthState0 .GE. 1) THEN
            !
            LastObservedPrices = RESHAPE(LastStateVec, (/ DepthState,numAgents /))
            !
        END IF
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Pre-shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        visitedStatesPre = 0
        visitedPrices = 0.d0
        visitedProfits = 0.d0
        p = LastObservedPrices
        pPrime = optimalStrategy(computeStateNumber(p),:)
        DO iPeriod = 1, numPeriods
            !
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            LastObservedPrices = p
            visitedStatesPre(iPeriod) = computeStateNumber(p)
            DO iAgent = 1, numAgents
                !
                visitedPrices(iPeriod,iAgent) = PricesGrids(pPrime(iAgent),iAgent)
                visitedProfits(iPeriod,iAgent) = PI(computeActionNumber(pPrime),iAgent)
                !
            END DO
            !
            ! Check if the state has already been visited
            !
            IF ((iPeriod .GE. 2) .AND. (ANY(visitedStatesPre(:iPeriod-1) .EQ. visitedStatesPre(iPeriod)))) THEN
                !
                PeriodsLengthPre = &
                    iPeriod-MINVAL(MINLOC((visitedStatesPre(:iPeriod-1)-visitedStatesPre(iPeriod))**2))
                FreqPeriodLengthPre(MIN(numThresPeriodsLength,PeriodsLengthPre)) = &
                    FreqPeriodLengthPre(MIN(numThresPeriodsLength,PeriodsLengthPre))+1
                PricePre = SUM(visitedPrices(iPeriod-PeriodsLengthPre+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPre)
                avgPricesPre = avgPricesPre+PricePre
                avgPricesPreQ = avgPricesPreQ+PricePre**2
                ProfitPre = SUM(visitedProfits(iPeriod-PeriodsLengthPre+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPre)
                avgProfitsPre = avgProfitsPre+ProfitPre
                avgProfitsPreQ = avgProfitsPreQ+ProfitPre**2
                EXIT
                !
            END IF
            pPrime = optimalStrategy(visitedStatesPre(iPeriod),:)
            !
        END DO
        visitedStatesPre(:PeriodsLengthPre) = &
            visitedStatesPre(iPeriod-PeriodsLengthPre+1:iPeriod)
        visitedStatesPre(PeriodsLengthPre+1:) = 0
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iAgent = 1, numAgents        ! Start of loop aver shocking agent 
            !
            visitedStates = 0
            p = LastObservedPrices
            !
            ! Price selection in shock period:
            ! Agent "iAgent" selects the price closest to the Nash one,
            ! The other agents stick to the strategy at convergence
            !
            pPrime = optimalStrategy(computeStateNumber(p),:)
            pNash = MINVAL(MINLOC((pricesGrids(:,iAgent)-NashPrices(iAgent))**2))
            pPrime(iAgent) = pNash
            !
            flagReturnedToState = .FALSE.
            DO iPeriod = 1, MAX(numShockPeriodsPrint,numPeriods)
                !
                IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                p(1,:) = pPrime
                visitedStates(iPeriod) = computeStateNumber(p)
                DO jAgent = 1, numAgents
                    !
                    IF (iPeriod .LE. numShockPeriodsPrint) THEN
                        !
                        avgPricesShock(iPeriod,iAgent,jAgent) = & 
                            avgPricesShock(iPeriod,iAgent,jAgent)+PricesGrids(pPrime(jAgent),jAgent)
                        avgPricesShockQ(iPeriod,iAgent,jAgent) = & 
                            avgPricesShockQ(iPeriod,iAgent,jAgent)+PricesGrids(pPrime(jAgent),jAgent)**2
                        avgPricesPercShock(iPeriod,iAgent,jAgent) = & 
                            avgPricesPercShock(iPeriod,iAgent,jAgent)+(PricesGrids(pPrime(jAgent),jAgent)/PricePre(jAgent))
                        avgPricesPercShockQ(iPeriod,iAgent,jAgent) = & 
                            avgPricesPercShockQ(iPeriod,iAgent,jAgent)+(PricesGrids(pPrime(jAgent),jAgent)/PricePre(jAgent))**2
                        !
                        avgProfitsShock(iPeriod,iAgent,jAgent) = & 
                            avgProfitsShock(iPeriod,iAgent,jAgent)+PI(computeActionNumber(pPrime),jAgent)
                        avgProfitsShockQ(iPeriod,iAgent,jAgent) = & 
                            avgProfitsShockQ(iPeriod,iAgent,jAgent)+PI(computeActionNumber(pPrime),jAgent)**2
                        avgProfitsPercShock(iPeriod,iAgent,jAgent) = & 
                            avgProfitsPercShock(iPeriod,iAgent,jAgent)+(PI(computeActionNumber(pPrime),jAgent)/ProfitPre(jAgent))
                        avgProfitsPercShockQ(iPeriod,iAgent,jAgent) = & 
                            avgProfitsPercShockQ(iPeriod,iAgent,jAgent)+(PI(computeActionNumber(pPrime),jAgent)/ProfitPre(jAgent))**2
                        !
                    END IF
                    !
                END DO
                !
                ! Check if the state has already been visited
                ! Case 1: the state retuns to one of the states in the pre-shock cycle
                !
                IF ((.NOT.(flagReturnedToState)) .AND. &
                    (ANY(visitedStatesPre(:PeriodsLengthPre) .EQ. visitedStates(iPeriod)))) THEN
                    !
                    PeriodsLengthShock = iPeriod
                    PunishmentStrategy = iPeriod
                    FreqPeriodLengthShock(iAgent,MIN(numThresPeriodsLength,PeriodsLengthShock)) = &
                        FreqPeriodLengthShock(iAgent,MIN(numThresPeriodsLength,PeriodsLengthShock))+1
                    FreqPunishmentStrategy(iAgent,MIN(numThresPeriodsLength,PunishmentStrategy)) = &
                        FreqPunishmentStrategy(iAgent,MIN(numThresPeriodsLength,PunishmentStrategy))+1
                    indexShockState = RESHAPE(p,(/ LengthStates /))
                    flagReturnedToState = .TRUE.
                    !
                END IF
                !
                ! Case 2: after some time, the state starts cycling among a new set of states
                !
                IF ((iPeriod .GE. 2) .AND. (.NOT.(flagReturnedToState)) .AND. &
                    (ANY(visitedStates(:iPeriod-1) .EQ. visitedStates(iPeriod)))) THEN
                    !
                    PeriodsLengthShock = &
                        MINVAL(MINLOC((visitedStates(:iPeriod-1)-visitedStates(iPeriod))**2))
                    PunishmentStrategy = 0
                    FreqPeriodLengthShock(iAgent,MIN(numThresPeriodsLength,PeriodsLengthShock)) = &
                        FreqPeriodLengthShock(iAgent,MIN(numThresPeriodsLength,PeriodsLengthShock))+1
                    FreqPunishmentStrategy(iAgent,MIN(numThresPeriodsLength,PunishmentStrategy)) = &
                        FreqPunishmentStrategy(iAgent,MIN(numThresPeriodsLength,PunishmentStrategy))+1
                    indexShockState = RESHAPE(p,(/ LengthStates /))
                    flagReturnedToState = .TRUE.
                    !
                END IF
                pPrime = optimalStrategy(visitedStates(iPeriod),:)
                IF (computeImpulseResponseToNash .EQ. -1) &
                    pPrime(iAgent) = pNash          ! The deviation to the Nash price is permanent
                IF ((computeImpulseResponseToNash .GE. 1) .AND. (iPeriod .LT. computeImpulseResponseToNash)) &
                    pPrime(iAgent) = pNash          ! The deviation to the Nash price is temporary
                !
            END DO
            !
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Post-shock period analysis
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            visitedStates = 0
            visitedPrices = 0.d0
            visitedProfits = 0.d0
            p = RESHAPE(indexShockState, (/ DepthState,numAgents /) )
            pPrime = optimalStrategy(computeStateNumber(p),:)
            IF (computeImpulseResponseToNash .EQ. -1) &            
                pPrime(iAgent) = pNash          ! The deviation to the Nash price is permanent!
            !            
            DO iPeriod = 1, numPeriods
                !
                IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                p(1,:) = pPrime
                visitedStates(iPeriod) = computeStateNumber(p)
                DO jAgent = 1, numAgents
                    !
                    visitedPrices(iPeriod,jAgent) = PricesGrids(pPrime(jAgent),jAgent)
                    visitedProfits(iPeriod,jAgent) = PI(computeActionNumber(pPrime),jAgent)
                    !
                END DO
                !
                ! Check if the state has already been visited
                !
                IF ((iPeriod .GE. 2) .AND. (ANY(visitedStates(:iPeriod-1) .EQ. visitedStates(iPeriod)))) THEN
                    !
                    PeriodsLengthPost = &
                        iPeriod-MINVAL(MINLOC((visitedStates(:iPeriod-1)-visitedStates(iPeriod))**2))
                    FreqPeriodLengthPost(iAgent,MIN(numThresPeriodsLength,PeriodsLengthPost)) = &
                        FreqPeriodLengthPost(iAgent,MIN(numThresPeriodsLength,PeriodsLengthPost))+1
                    !
                    avgPricesPost(iAgent,:) = avgPricesPost(iAgent,:)+ &
                        SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost)
                    avgPricesPostQ(iAgent,:) = avgPricesPostQ(iAgent,:)+ &
                        (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost))**2
                    avgPricesPercPost(iAgent,:) = avgPricesPercPost(iAgent,:)+ &
                        (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost)/PricePre)
                    avgPricesPercPostQ(iAgent,:) = avgPricesPercPostQ(iAgent,:)+ &
                        (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost)/PricePre)**2
                    !
                    avgProfitsPost(iAgent,:) = avgProfitsPost(iAgent,:)+ &
                        SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost)
                    avgProfitsPostQ(iAgent,:) = avgProfitsPostQ(iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost))**2
                    avgProfitsPercPost(iAgent,:) = avgProfitsPercPost(iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost)/ProfitPre)
                    avgProfitsPercPostQ(iAgent,:) = avgProfitsPercPostQ(iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/DBLE(PeriodsLengthPost)/ProfitPre)**2
                    !
                    EXIT
                    !
                END IF
                pPrime = optimalStrategy(visitedStates(iPeriod),:)
                IF (computeImpulseResponseToNash .EQ. -1) &
                    pPrime(iAgent) = pNash          ! The deviation to the Nash price is permanent
                !
            END DO
            !
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! End of Impulse Response analysis
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
        END DO                          ! End of loop aver shocking agent 
        !
    END DO        ! End of loop over games
    !$omp end parallel do
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Computing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Averages of prices and profits
    !
    avgPricesPre = avgPricesPre/DBLE(numGames)
    avgProfitsPre = avgProfitsPre/DBLE(numGames)
    avgPricesShock = avgPricesShock/DBLE(numGames)
    avgProfitsShock = avgProfitsShock/DBLE(numGames)
    avgPricesPost = avgPricesPost/DBLE(numGames)
    avgProfitsPost = avgProfitsPost/DBLE(numGames)
    avgPricesPercShock = avgPricesPercShock/DBLE(numGames)
    avgProfitsPercShock = avgProfitsPercShock/DBLE(numGames)
    avgPricesPercPost = avgPricesPercPost/DBLE(numGames)
    avgProfitsPercPost = avgProfitsPercPost/DBLE(numGames)
    avgPricesPreQ = avgPricesPreQ/DBLE(numGames)
    avgProfitsPreQ = avgProfitsPreQ/DBLE(numGames)
    avgPricesShockQ = avgPricesShockQ/DBLE(numGames)
    avgProfitsShockQ = avgProfitsShockQ/DBLE(numGames)
    avgPricesPostQ = avgPricesPostQ/DBLE(numGames)
    avgProfitsPostQ = avgProfitsPostQ/DBLE(numGames)
    avgPricesPercShockQ = avgPricesPercShockQ/DBLE(numGames)
    avgProfitsPercShockQ = avgProfitsPercShockQ/DBLE(numGames)
    avgPricesPercPostQ = avgPricesPercPostQ/DBLE(numGames)
    avgProfitsPercPostQ = avgProfitsPercPostQ/DBLE(numGames)
    !
    ! Computing aggregate (deviating and non-deviating) averages of prices and profits
    !
    AggrPricesPre = SUM(avgPricesPre)/DBLE(numAgents)
    AggrProfitsPre = SUM(avgProfitsPre)/DBLE(numAgents)
    AggrPricesPreQ = SUM(avgPricesPreQ)/DBLE(numAgents)
    AggrProfitsPreQ = SUM(avgProfitsPreQ)/DBLE(numAgents)
    !
    AggrDevPricesShock = 0.d0
    AggrDevProfitsShock = 0.d0
    AggrDevPricesShockQ = 0.d0
    AggrDevProfitsShockQ = 0.d0
    AggrDevPricesPercShock = 0.d0
    AggrDevProfitsPercShock = 0.d0
    AggrDevPricesPercShockQ = 0.d0
    AggrDevProfitsPercShockQ = 0.d0
    DO iPeriod = 1, numShockPeriodsPrint
        !
        DO iAgent = 1, numAgents
            !
            AggrDevPricesShock(iPeriod) = AggrDevPricesShock(iPeriod)+avgPricesShock(iPeriod,iAgent,iAgent)
            AggrDevProfitsShock(iPeriod) = AggrDevProfitsShock(iPeriod)+avgProfitsShock(iPeriod,iAgent,iAgent)
            AggrDevPricesShockQ(iPeriod) = AggrDevPricesShockQ(iPeriod)+avgPricesShockQ(iPeriod,iAgent,iAgent)
            AggrDevProfitsShockQ(iPeriod) = AggrDevProfitsShockQ(iPeriod)+avgProfitsShockQ(iPeriod,iAgent,iAgent)
            !
            AggrDevPricesPercShock(iPeriod) = AggrDevPricesPercShock(iPeriod)+avgPricesPercShock(iPeriod,iAgent,iAgent)
            AggrDevProfitsPercShock(iPeriod) = AggrDevProfitsPercShock(iPeriod)+avgProfitsPercShock(iPeriod,iAgent,iAgent)
            AggrDevPricesPercShockQ(iPeriod) = AggrDevPricesPercShockQ(iPeriod)+avgPricesPercShockQ(iPeriod,iAgent,iAgent)
            AggrDevProfitsPercShockQ(iPeriod) = AggrDevProfitsPercShockQ(iPeriod)+avgProfitsPercShockQ(iPeriod,iAgent,iAgent)
            !
        END DO
        AggrNonDevPricesShock(iPeriod) = (SUM(avgPricesShock(iPeriod,:,:))-AggrDevPricesShock(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesShock(iPeriod) = AggrDevPricesShock(iPeriod)/DBLE(numAgents)
        AggrNonDevProfitsShock(iPeriod) = (SUM(avgProfitsShock(iPeriod,:,:))-AggrDevProfitsShock(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsShock(iPeriod) = AggrDevProfitsShock(iPeriod)/DBLE(numAgents)
        AggrNonDevPricesShockQ(iPeriod) = (SUM(avgPricesShockQ(iPeriod,:,:))-AggrDevPricesShockQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesShockQ(iPeriod) = AggrDevPricesShockQ(iPeriod)/DBLE(numAgents)
        AggrNonDevProfitsShockQ(iPeriod) = (SUM(avgProfitsShockQ(iPeriod,:,:))-AggrDevProfitsShockQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsShockQ(iPeriod) = AggrDevProfitsShockQ(iPeriod)/DBLE(numAgents)
        !
        AggrNonDevPricesPercShock(iPeriod) = (SUM(avgPricesPercShock(iPeriod,:,:))-AggrDevPricesPercShock(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesPercShock(iPeriod) = AggrDevPricesPercShock(iPeriod)/DBLE(numAgents)
        AggrNonDevProfitsPercShock(iPeriod) = (SUM(avgProfitsPercShock(iPeriod,:,:))-AggrDevProfitsPercShock(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsPercShock(iPeriod) = AggrDevProfitsPercShock(iPeriod)/DBLE(numAgents)
        AggrNonDevPricesPercShockQ(iPeriod) = (SUM(avgPricesPercShockQ(iPeriod,:,:))-AggrDevPricesPercShockQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesPercShockQ(iPeriod) = AggrDevPricesPercShockQ(iPeriod)/DBLE(numAgents)
        AggrNonDevProfitsPercShockQ(iPeriod) = (SUM(avgProfitsPercShockQ(iPeriod,:,:))-AggrDevProfitsPercShockQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsPercShockQ(iPeriod) = AggrDevProfitsPercShockQ(iPeriod)/DBLE(numAgents)
        !
    END DO
    !
    AggrDevPricesPost = 0.d0
    AggrDevProfitsPost = 0.d0
    AggrDevPricesPostQ = 0.d0
    AggrDevProfitsPostQ = 0.d0
    AggrDevPricesPercPost = 0.d0
    AggrDevProfitsPercPost = 0.d0
    AggrDevPricesPercPostQ = 0.d0
    AggrDevProfitsPercPostQ = 0.d0
    DO iAgent = 1, numAgents
        !
        AggrDevPricesPost = AggrDevPricesPost+avgPricesPost(iAgent,iAgent)
        AggrDevProfitsPost = AggrDevProfitsPost+avgProfitsPost(iAgent,iAgent)
        AggrDevPricesPostQ = AggrDevPricesPostQ+avgPricesPostQ(iAgent,iAgent)
        AggrDevProfitsPostQ = AggrDevProfitsPostQ+avgProfitsPostQ(iAgent,iAgent)
        !
        AggrDevPricesPercPost = AggrDevPricesPercPost+avgPricesPercPost(iAgent,iAgent)
        AggrDevProfitsPercPost = AggrDevProfitsPercPost+avgProfitsPercPost(iAgent,iAgent)
        AggrDevPricesPercPostQ = AggrDevPricesPercPostQ+avgPricesPercPostQ(iAgent,iAgent)
        AggrDevProfitsPercPostQ = AggrDevProfitsPercPostQ+avgProfitsPercPostQ(iAgent,iAgent)
        !
    END DO
    AggrNonDevPricesPost = (SUM(avgPricesPost(:,:))-AggrDevPricesPost)/DBLE(numAgents*(numAgents-1))
    AggrDevPricesPost = AggrDevPricesPost/DBLE(numAgents)
    AggrNonDevProfitsPost = (SUM(avgProfitsPost(:,:))-AggrDevProfitsPost)/DBLE(numAgents*(numAgents-1))
    AggrDevProfitsPost = AggrDevProfitsPost/DBLE(numAgents)
    AggrNonDevPricesPostQ = (SUM(avgPricesPostQ(:,:))-AggrDevPricesPostQ)/DBLE(numAgents*(numAgents-1))
    AggrDevPricesPostQ = AggrDevPricesPostQ/DBLE(numAgents)
    AggrNonDevProfitsPostQ = (SUM(avgProfitsPostQ(:,:))-AggrDevProfitsPostQ)/DBLE(numAgents*(numAgents-1))
    AggrDevProfitsPostQ = AggrDevProfitsPostQ/DBLE(numAgents)
    !
    AggrNonDevPricesPercPost = (SUM(avgPricesPercPost(:,:))-AggrDevPricesPercPost)/DBLE(numAgents*(numAgents-1))
    AggrDevPricesPercPost = AggrDevPricesPercPost/DBLE(numAgents)
    AggrNonDevProfitsPercPost = (SUM(avgProfitsPercPost(:,:))-AggrDevProfitsPercPost)/DBLE(numAgents*(numAgents-1))
    AggrDevProfitsPercPost = AggrDevProfitsPercPost/DBLE(numAgents)
    AggrNonDevPricesPercPostQ = (SUM(avgPricesPercPostQ(:,:))-AggrDevPricesPercPostQ)/DBLE(numAgents*(numAgents-1))
    AggrDevPricesPercPostQ = AggrDevPricesPercPostQ/DBLE(numAgents)
    AggrNonDevProfitsPercPostQ = (SUM(avgProfitsPercPostQ(:,:))-AggrDevProfitsPercPostQ)/DBLE(numAgents*(numAgents-1))
    AggrDevProfitsPercPostQ = AggrDevProfitsPercPostQ/DBLE(numAgents)
    !
    ! Computing standard errors
    !
    avgPricesPreQ = SQRT((avgPricesPreQ-avgPricesPre**2)/DBLE(numGames))
    avgProfitsPreQ = SQRT((avgProfitsPreQ-avgProfitsPre**2)/DBLE(numGames))
    avgPricesShockQ = SQRT((avgPricesShockQ-avgPricesShock**2)/DBLE(numGames))
    avgProfitsShockQ = SQRT((avgProfitsShockQ-avgProfitsShock**2)/DBLE(numGames))
    avgPricesPostQ = SQRT((avgPricesPostQ-avgPricesPost**2)/DBLE(numGames))
    avgProfitsPostQ = SQRT((avgProfitsPostQ-avgProfitsPost**2)/DBLE(numGames))
    !
    AggrPricesPreQ = SQRT((AggrPricesPreQ-AggrPricesPre**2)/DBLE(numGames))
    AggrProfitsPreQ = SQRT((AggrProfitsPreQ-AggrProfitsPre**2)/DBLE(numGames))
    DO iPeriod = 1, numShockPeriodsPrint
        !
        AggrNonDevPricesShockQ(iPeriod) = SQRT((AggrNonDevPricesShockQ(iPeriod)-AggrNonDevPricesShock(iPeriod)**2)/DBLE(numGames))
        AggrDevPricesShockQ(iPeriod) = SQRT((AggrDevPricesShockQ(iPeriod)-AggrDevPricesShock(iPeriod)**2)/DBLE(numGames))
        AggrNonDevProfitsShockQ(iPeriod) = SQRT((AggrNonDevProfitsShockQ(iPeriod)-AggrNonDevProfitsShock(iPeriod)**2)/DBLE(numGames))
        AggrDevProfitsShockQ(iPeriod) = SQRT((AggrDevProfitsShockQ(iPeriod)-AggrDevProfitsShock(iPeriod)**2)/DBLE(numGames))
        !
        AggrNonDevPricesPercShockQ(iPeriod) = SQRT((AggrNonDevPricesPercShockQ(iPeriod)-AggrNonDevPricesPercShock(iPeriod)**2)/DBLE(numGames))
        AggrDevPricesPercShockQ(iPeriod) = SQRT((AggrDevPricesPercShockQ(iPeriod)-AggrDevPricesPercShock(iPeriod)**2)/DBLE(numGames))
        AggrNonDevProfitsPercShockQ(iPeriod) = SQRT((AggrNonDevProfitsPercShockQ(iPeriod)-AggrNonDevProfitsPercShock(iPeriod)**2)/DBLE(numGames))
        AggrDevProfitsPercShockQ(iPeriod) = SQRT((AggrDevProfitsPercShockQ(iPeriod)-AggrDevProfitsPercShock(iPeriod)**2)/DBLE(numGames))
        !
    END DO
    AggrNonDevPricesPostQ = SQRT((AggrNonDevPricesPostQ-AggrNonDevPricesPost**2)/DBLE(numGames))
    AggrDevPricesPostQ = SQRT((AggrDevPricesPostQ-AggrDevPricesPost**2)/DBLE(numGames))
    AggrNonDevProfitsPostQ = SQRT((AggrNonDevProfitsPostQ-AggrNonDevProfitsPost**2)/DBLE(numGames))
    AggrDevProfitsPostQ = SQRT((AggrDevProfitsPostQ-AggrDevProfitsPost**2)/DBLE(numGames))
    !
    AggrNonDevPricesPercPostQ = SQRT((AggrNonDevPricesPercPostQ-AggrNonDevPricesPercPost**2)/DBLE(numGames))
    AggrDevPricesPercPostQ = SQRT((AggrDevPricesPercPostQ-AggrDevPricesPercPost**2)/DBLE(numGames))
    AggrNonDevProfitsPercPostQ = SQRT((AggrNonDevProfitsPercPostQ-AggrNonDevProfitsPercPost**2)/DBLE(numGames))
    AggrDevProfitsPercPostQ = SQRT((AggrDevProfitsPercPostQ-AggrDevProfitsPercPost**2)/DBLE(numGames))    
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(100031,1) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (iPeriod, iPeriod = 1, numShockPeriodsPrint), (iPeriod, iPeriod = 1, numShockPeriodsPrint), &
            (ThresPeriodsLength(i), i = 1, numThresPeriodsLength), &
            ((iAgent, ThresPeriodsLength(i), i = 1, numThresPeriodsLength), iAgent = 1, numAgents), &
            ((iAgent, ThresPeriodsLength(i), i = 1, numThresPeriodsLength), iAgent = 1, numAgents), &
            ((iAgent, ThresPeriodsLength0(i), i = 1, numThresPeriodsLength0), iAgent = 1, numAgents), &
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), &  ! avgPrices
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), &  ! sePrices
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), &  ! avgProfits
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents)     ! seProfits        
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>(' MExplPar', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I2.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            'AggrPricePre ', <numShockPeriodsPrint>('AggrDevPriceShockPer', I3.3, ' '), 'AggrDevPricePost ', &
                <numShockPeriodsPrint>('AggrNonDevPriceShockPer', I3.3, ' '), 'AggrNonDevPricePost ', &
            'seAggrPricePre ', <numShockPeriodsPrint>('seAggrDevPriceShockPer', I3.3, ' '), 'seAggrDevPricePost ', &
                <numShockPeriodsPrint>('seAggrNonDevPriceShockPer', I3.3, ' '), 'seAggrNonDevPricePost ', &
            <numShockPeriodsPrint>('AggrDevPricePercShockPer', I3.3, ' '), 'AggrDevPricePercPost ', &
                <numShockPeriodsPrint>('AggrNonDevPricePercShockPer', I3.3, ' '), 'AggrNonDevPricePercPost ', &
            <numShockPeriodsPrint>('seAggrDevPricePercShockPer', I3.3, ' '), 'seAggrDevPricePercPost ', &
                <numShockPeriodsPrint>('seAggrNonDevPricePercShockPer', I3.3, ' '), 'seAggrNonDevPricePercPost ', &
            'AggrProfitPre ', <numShockPeriodsPrint>('AggrDevProfitShockPer', I3.3, ' '), 'AggrDevProfitPost ', &
                <numShockPeriodsPrint>('AggrNonDevProfitShockPer', I3.3, ' '), 'AggrNonDevProfitPost ', &
            'seAggrProfitPre ', <numShockPeriodsPrint>('seAggrDevProfitShockPer', I3.3, ' '), 'seAggrDevProfitPost ', &
                <numShockPeriodsPrint>('seAggrNonDevProfitShockPer', I3.3, ' '), 'seAggrNonDevProfitPost ', &
            <numShockPeriodsPrint>('AggrDevProfitPercShockPer', I3.3, ' '), 'AggrDevProfitPercPost ', &
                <numShockPeriodsPrint>('AggrNonDevProfitPercShockPer', I3.3, ' '), 'AggrNonDevProfitPercPost ', &
            <numShockPeriodsPrint>('seAggrDevProfitPercShockPer', I3.3, ' '), 'seAggrDevProfitPercPost ', &
                <numShockPeriodsPrint>('seAggrNonDevProfitPercShockPer', I3.3, ' '), 'seAggrNonDevProfitPercPost ', &
            <numThresPeriodsLength>('#PerLenPre=', I2.2, ' '), &
            <numAgents>(<numThresPeriodsLength>('#PerLenPostAg', I1, '=', I2.2, ' ')), &
            <numAgents>(<numThresPeriodsLength>('#PerLenShockAg', I1, '=', I2.2, ' ')), &
            <numAgents>(<numThresPeriodsLength+1>('#PunishStratAg', I1, '=', I2.2, ' ')), &
            <numAgents>(<numAgents>('Ag', I1, 'avgPricePre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'avgPriceShockAg', I1, 'Per', I3.3, ' '), &
                'Ag', I1, 'avgPricePostAg', I1, ' ')), &
            <numAgents>(<numAgents>('seAg', I1, 'avgPricePre', ' ', &
                <numShockPeriodsPrint>('seAg', I1, 'avgPriceShockAg', I1, 'Per', I3.3, ' '), &
                'seAg', I1, 'avgPricePostAg', I1, ' ')), &
            <numAgents>(<numAgents>('Ag', I1, 'avgProfitPre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'avgProfitShockAg', I1, 'Per', I3.3, ' '), &
                'Ag', I1, 'avgProfitPostAg', I1, ' ')), &
            <numAgents>(<numAgents>('seAg', I1, 'avgProfitPre', ' ', &
                <numShockPeriodsPrint>('seAg', I1, 'avgProfitShockAg', I1, 'Per', I3.3, ' '), &
                'seAg', I1, 'avgProfitPostAg', I1, ' ')) &
            )
        !
    END IF
    !
    WRITE(100031,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        AggrPricesPre, AggrDevPricesShock, AggrDevPricesPost, AggrNonDevPricesShock, AggrNonDevPricesPost, &
        AggrPricesPreQ, AggrDevPricesShockQ, AggrDevPricesPostQ, AggrNonDevPricesShockQ, AggrNonDevPricesPostQ, &
        AggrDevPricesPercShock, AggrDevPricesPercPost, AggrNonDevPricesPercShock, AggrNonDevPricesPercPost, &
        AggrDevPricesPercShockQ, AggrDevPricesPercPostQ, AggrNonDevPricesPercShockQ, AggrNonDevPricesPercPostQ, &
        AggrProfitsPre, AggrDevProfitsShock, AggrDevProfitsPost, AggrNonDevProfitsShock, AggrNonDevProfitsPost, &
        AggrProfitsPreQ, AggrDevProfitsShockQ, AggrDevProfitsPostQ, AggrNonDevProfitsShockQ, AggrNonDevProfitsPostQ, &
        AggrDevProfitsPercShock, AggrDevProfitsPercPost, AggrNonDevProfitsPercShock, AggrNonDevProfitsPercPost, &
        AggrDevProfitsPercShockQ, AggrDevProfitsPercPostQ, AggrNonDevProfitsPercShockQ, AggrNonDevProfitsPercPostQ, &
        FreqPeriodLengthPre, &
        ((FreqPeriodLengthPost(iAgent,i), i = 1, numThresPeriodsLength), iAgent = 1, numAgents), &
        ((FreqPeriodLengthShock(iAgent,i), i = 1, numThresPeriodsLength), iAgent = 1, numAgents), &
        ((FreqPunishmentStrategy(iAgent,i), i = 0, numThresPeriodsLength), iAgent = 1, numAgents), &
        ((avgPricesPre(jAgent), (avgPricesShock(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgPricesPost(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), &
        ((avgPricesPreQ(jAgent), (avgPricesShockQ(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgPricesPostQ(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), &
        ((avgProfitsPre(jAgent), (avgProfitsShock(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgProfitsPost(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), &
        ((avgProfitsPreQ(jAgent), (avgProfitsShockQ(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgProfitsPostQ(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents)
2   FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.7, 1X), &
        <6*numAgents>(F10.7, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        F12.7, 1X, <numShockPeriodsPrint>(F23.7,1X), F16.7, 1X, <numShockPeriodsPrint>(F26.7,1X), F19.7, 1X, &
        F14.7, 1X, <numShockPeriodsPrint>(F25.7,1X), F18.7, 1X, <numShockPeriodsPrint>(F28.7,1X), F21.7, 1X, &
        <numShockPeriodsPrint>(F27.7,1X), F20.7, 1X, <numShockPeriodsPrint>(F30.7,1X), F23.7, 1X, &
        <numShockPeriodsPrint>(F29.7,1X), F22.7, 1X, <numShockPeriodsPrint>(F32.7,1X), F25.7, 1X, &
        F13.7, 1X, <numShockPeriodsPrint>(F24.7,1X), F17.7, 1X, <numShockPeriodsPrint>(F27.7,1X), F20.7, 1X, &
        F15.7, 1X, <numShockPeriodsPrint>(F26.7,1X), F19.7, 1X, <numShockPeriodsPrint>(F29.7,1X), F22.7, 1X, &
        <numShockPeriodsPrint>(F28.7,1X), F21.7, 1X, <numShockPeriodsPrint>(F31.7,1X), F24.7, 1X, &
        <numShockPeriodsPrint>(F30.7,1X), F23.7, 1X, <numShockPeriodsPrint>(F33.7,1X), F26.7, 1X, &
        <numThresPeriodsLength>(I13, 1X), &
        <numAgents>(<numThresPeriodsLength>(I17, 1X)), &
        <numAgents>(<numThresPeriodsLength>(I18, 1X)), &
        <numAgents>(<numThresPeriodsLength0>(I18, 1X)), &
        <numAgents>(<numAgents>(F14.7, 1X, <numShockPeriodsPrint>(F25.7, 1X), F18.7, 1X)), &
        <numAgents>(<numAgents>(F16.7, 1X, <numShockPeriodsPrint>(F27.7, 1X), F20.7, 1X)), &
        <numAgents>(<numAgents>(F15.7, 1X, <numShockPeriodsPrint>(F26.7, 1X), F19.7, 1X)), &
        <numAgents>(<numAgents>(F17.7, 1X, <numShockPeriodsPrint>(F28.7, 1X), F21.7, 1X)) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeIRToNashAnalysis
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ImpulseResponseToNash
