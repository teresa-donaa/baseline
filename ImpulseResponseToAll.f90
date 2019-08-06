MODULE ImpulseResponseToAll
!
USE globals
USE QL_routines
!
! Computes analysis of the Impulse Response to a one-period deviation to all prices
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeIRToAllAnalysis ( iModel )
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
    INTEGER, PARAMETER :: numShockPeriodsPrint = 25
    INTEGER :: PeriodsLengthPre, PeriodsLengthPost, &
        visitedStatesPre(numStates+1), visitedStates(MAX(numShockPeriodsPrint,numStates+1)), &
        p(DepthState,numAgents), pPrime(numAgents), numPeriodsShockTmp(numShockPeriodsPrint,numAgents), &
        iStatePre, iGame, iAgent, iPrice, iPeriod, jAgent, &
        optimalStrategy(numStates,numAgents), LastObservedPrices(DepthState,numAgents), &
        indexShockState(LengthStates), iThres, i, j
    INTEGER, DIMENSION(numStates+1,numAgents) :: indexPricesPre
    REAL(8) :: nn
    REAL(8), DIMENSION(numStates+1,numAgents) :: visitedPrices, visitedProfits, PricesPre, ProfitsPre
    REAL(8), DIMENSION(numAgents) :: avgPricesPre, avgProfitsPre, avgPricesPreQ, avgProfitsPreQ
    REAL(8), DIMENSION(numPrices,numShockPeriodsPrint,numAgents,numAgents) :: &
        avgPricesShock, avgProfitsShock, avgPricesShockQ, avgProfitsShockQ, &
        avgPricesPercShock, avgProfitsPercShock, avgPricesPercShockQ, avgProfitsPercShockQ
    REAL(8), DIMENSION(numPrices,numAgents,numAgents) :: &
        avgPricesPost, avgProfitsPost, avgPricesPostQ, avgProfitsPostQ, &
        avgPricesPercPost, avgProfitsPercPost, avgPricesPercPostQ, avgProfitsPercPostQ
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents) :: &
        avgPricesShockTmp, avgProfitsShockTmp, avgPricesPercShockTmp, avgProfitsPercShockTmp
    LOGICAL :: FlagReturnedToState
	INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    REAL(8) :: AggrPricesPre, AggrProfitsPre
    REAL(8), DIMENSION(numPrices) :: AggrDevPricesPost, AggrNonDevPricesPost, AggrDevProfitsPost, AggrNonDevProfitsPost, &
        AggrDevPricesPercPost, AggrNonDevPricesPercPost, AggrDevProfitsPercPost, AggrNonDevProfitsPercPost
    REAL(8), DIMENSION(numPrices,numShockPeriodsPrint) :: &
        AggrDevPricesShock, AggrNonDevPricesShock, AggrDevProfitsShock, AggrNonDevProfitsShock, &
        AggrDevPricesPercShock, AggrNonDevPricesPercShock, AggrDevProfitsPercShock, AggrNonDevProfitsPercShock
    REAL(8) :: AggrPricesPreQ, AggrProfitsPreQ
    REAL(8), DIMENSION(numPrices) :: AggrDevPricesPostQ, AggrNonDevPricesPostQ, AggrDevProfitsPostQ, AggrNonDevProfitsPostQ, &
        AggrDevPricesPercPostQ, AggrNonDevPricesPercPostQ, AggrDevProfitsPercPostQ, AggrNonDevProfitsPercPostQ
    REAL(8), DIMENSION(numPrices,numShockPeriodsPrint) :: &
        AggrDevPricesShockQ, AggrNonDevPricesShockQ, AggrDevProfitsShockQ, AggrNonDevProfitsShockQ, &
        AggrDevPricesPercShockQ, AggrNonDevPricesPercShockQ, AggrDevProfitsPercShockQ, AggrNonDevProfitsPercShockQ
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Impulse Response functions to all prices'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !
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
    !$omp   visitedProfits,p,pPrime,iPeriod,iAgent,OptimalStrategyVec,LastStateVec, &
    !$omp   visitedStates,flagReturnedToState,jAgent,indexShockState,indexPricesPre,PricesPre,ProfitsPre, &
    !$omp   PeriodsLengthPre,iStatePre,PeriodsLengthPost, &
    !$omp   avgPricesShockTmp,avgProfitsShockTmp,avgPricesPercShockTmp,avgProfitsPercShockTmp, &
    !$omp   numPeriodsShockTmp,nn,iPrice) &
    !$omp firstprivate(PI,PricesGrids) &
    !$omp reduction(+ : avgPricesPre,avgProfitsPre,avgPricesShock,avgProfitsShock,avgPricesPost,avgProfitsPost, &
    !$omp   avgPricesPreQ,avgProfitsPreQ,avgPricesShockQ,avgProfitsShockQ,avgPricesPostQ,avgProfitsPostQ, &
    !$omp   avgPricesPercShock,avgProfitsPercShock,avgPricesPercPost,avgProfitsPercPost, &
    !$omp   avgPricesPercShockQ,avgProfitsPercShockQ,avgPricesPercPostQ,avgProfitsPercPostQ)
    DO iGame = 1, numGames        ! Start of loop over games
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
        indexPricesPre = 0
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
                indexPricesPre(iPeriod,iAgent) = pPrime(iAgent)
                PricesPre(iPeriod,iAgent) = PricesGrids(pPrime(iAgent),iAgent)
                ProfitsPre(iPeriod,iAgent) = PI(computeActionNumber(pPrime),iAgent)
                !
            END DO
            !
            ! Check if the state has already been visited
            !
            IF ((iPeriod .GE. 2) .AND. (ANY(visitedStatesPre(:iPeriod-1) .EQ. visitedStatesPre(iPeriod)))) EXIT
            !
            ! Update pPrime and iterate
            !
            pPrime = optimalStrategy(visitedStatesPre(iPeriod),:)
            !
        END DO
        !
        PeriodsLengthPre = &
            iPeriod-MINVAL(MINLOC((visitedStatesPre(:iPeriod-1)-visitedStatesPre(iPeriod))**2))
        !
        visitedStatesPre(:PeriodsLengthPre) = visitedStatesPre(iPeriod-PeriodsLengthPre+1:iPeriod)
        visitedStatesPre(PeriodsLengthPre+1:) = 0
        indexPricesPre(:PeriodsLengthPre,:) = indexPricesPre(iPeriod-PeriodsLengthPre+1:iPeriod,:)
        indexPricesPre(PeriodsLengthPre+1:,:) = 0
        PricesPre(:PeriodsLengthPre,:) = PricesPre(iPeriod-PeriodsLengthPre+1:iPeriod,:)
        PricesPre(PeriodsLengthPre+1:,:) = 0.d0
        ProfitsPre(:PeriodsLengthPre,:) = ProfitsPre(iPeriod-PeriodsLengthPre+1:iPeriod,:)
        ProfitsPre(PeriodsLengthPre+1:,:) = 0.d0
        !
        avgPricesPre = avgPricesPre+ &
            SUM(PricesPre(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
        avgPricesPreQ = avgPricesPreQ+ &
            (SUM(PricesPre(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre))**2
        avgProfitsPre = avgProfitsPre+ &
            SUM(ProfitsPre(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
        avgProfitsPreQ = avgProfitsPreQ+ &
            (SUM(ProfitsPre(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre))**2
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iAgent = 1, numAgents        ! Start of loop over shocking agent 
            !
            DO iPrice = 1, numPrices    ! Start of loop over possible deviations
                !
                avgPricesShockTmp = 0.d0
                avgProfitsShockTmp = 0.d0
                avgPricesPercShockTmp = 0.d0
                avgProfitsPercShockTmp = 0.d0
                numPeriodsShockTmp = 0
                !
                DO iStatePre = 1, PeriodsLengthPre      ! Start of loop over pre-shock cycle states
                    !
                    visitedStates = 0
                    pPrime = indexPricesPre(iStatePre,:)
                    !
                    ! Price selection in shock period:
                    ! Agent "iAgent" selects the each price in turn,
                    ! The other agents stick to the strategy at convergence
                    !
                    pPrime(iAgent) = iPrice
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
                                numPeriodsShockTmp(iPeriod,jAgent) = numPeriodsShockTmp(iPeriod,jAgent)+1
                                nn = DBLE(numPeriodsShockTmp(iPeriod,jAgent))
                                avgPricesShockTmp(iPeriod,jAgent) = &
                                    (nn-1.d0)/nn*avgPricesShockTmp(iPeriod,jAgent)+ &
                                    PricesGrids(pPrime(jAgent),jAgent)/nn
                                avgPricesPercShockTmp(iPeriod,jAgent) = &
                                    (nn-1.d0)/nn*avgPricesPercShockTmp(iPeriod,jAgent)+ &
                                    (PricesGrids(pPrime(jAgent),jAgent)/PricesPre(iStatePre,jAgent))/nn
                                avgProfitsShockTmp(iPeriod,jAgent) = &
                                    (nn-1.d0)/nn*avgProfitsShockTmp(iPeriod,jAgent)+ &
                                    PI(computeActionNumber(pPrime),jAgent)/nn
                                avgProfitsPercShockTmp(iPeriod,jAgent) = &
                                    (nn-1.d0)/nn*avgProfitsShockTmp(iPeriod,jAgent)+ &
                                    (PI(computeActionNumber(pPrime),jAgent)/ProfitsPre(iStatePre,jAgent))/nn
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
                            indexShockState = RESHAPE(p,(/ LengthStates /))
                            flagReturnedToState = .TRUE.
                            !
                        END IF
                        pPrime = optimalStrategy(visitedStates(iPeriod),:)
                        !
                    END DO
                    !
                    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Post-shock period analysis
                    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !
                    visitedPrices = 0.d0
                    visitedProfits = 0.d0
                    visitedStates = 0
                    p = RESHAPE(indexShockState, (/ DepthState,numAgents /) )
                    pPrime = optimalStrategy(computeStateNumber(p),:)
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
                        IF ((iPeriod .GE. 2) .AND. (ANY(visitedStates(:iPeriod-1) .EQ. visitedStates(iPeriod)))) EXIT
                        !
                        ! Update pPrime and iterate
                        !
                        pPrime = optimalStrategy(visitedStates(iPeriod),:)
                        !
                    END DO
                    !
                    PeriodsLengthPost = &
                        iPeriod-MINVAL(MINLOC((visitedStates(:iPeriod-1)-visitedStates(iPeriod))**2))
                    !
                    avgPricesPost(iPrice,iAgent,:) = avgPricesPost(iPrice,iAgent,:)+ &
                        SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)/DBLE(PeriodsLengthPre)
                    avgPricesPostQ(iPrice,iAgent,:) = avgPricesPostQ(iPrice,iAgent,:)+ &
                        (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost))**2/DBLE(PeriodsLengthPre)
                    avgPricesPercPost(iPrice,iAgent,:) = avgPricesPercPost(iPrice,iAgent,:)+ &
                        SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)/PricesPre(iStatePre,:)/DBLE(PeriodsLengthPre)
                    avgPricesPercPostQ(iPrice,iAgent,:) = avgPricesPercPostQ(iPrice,iAgent,:)+ &
                        (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)/PricesPre(iStatePre,:))**2/DBLE(PeriodsLengthPre)
                    !
                    avgProfitsPost(iPrice,iAgent,:) = avgProfitsPost(iPrice,iAgent,:)+ &
                        SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPost)/DBLE(PeriodsLengthPre)
                    avgProfitsPostQ(iPrice,iAgent,:) = avgProfitsPostQ(iPrice,iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPost))**2/DBLE(PeriodsLengthPre)
                    avgProfitsPercPost(iPrice,iAgent,:) = avgProfitsPercPost(iPrice,iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPost)/ProfitsPre(iStatePre,:)/DBLE(PeriodsLengthPre))
                    avgProfitsPercPostQ(iPrice,iAgent,:) = avgProfitsPercPostQ(iPrice,iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPost)/ProfitsPre(iStatePre,:))**2/DBLE(PeriodsLengthPre)
                    !
                END DO
                !
                ! Compute average prices and profits over pre-shock cycle states
                !
                avgPricesShock(iPrice,:,iAgent,:) = avgPricesShock(iPrice,:,iAgent,:)+avgPricesShockTmp
                avgPricesShockQ(iPrice,:,iAgent,:) = avgPricesShockQ(iPrice,:,iAgent,:)+avgPricesShockTmp**2
                avgProfitsShock(iPrice,:,iAgent,:) = avgProfitsShock(iPrice,:,iAgent,:)+avgProfitsShockTmp
                avgProfitsShockQ(iPrice,:,iAgent,:) = avgProfitsShockQ(iPrice,:,iAgent,:)+avgProfitsShockTmp**2
                avgPricesPercShock(iPrice,:,iAgent,:) = avgPricesPercShock(iPrice,:,iAgent,:)+avgPricesPercShockTmp
                avgPricesPercShockQ(iPrice,:,iAgent,:) = avgPricesPercShockQ(iPrice,:,iAgent,:)+avgPricesPercShockTmp**2
                avgProfitsPercShock(iPrice,:,iAgent,:) = avgProfitsPercShock(iPrice,:,iAgent,:)+avgProfitsPercShockTmp
                avgProfitsPercShockQ(iPrice,:,iAgent,:) = avgProfitsPercShockQ(iPrice,:,iAgent,:)+avgProfitsPercShockTmp**2
                !
                ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! End of Impulse Response analysis
                ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !
            END DO                      ! End of loop over deviations
            !
        END DO                          ! End of loop over shocking agent 
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
    DO iPrice = 1, numPrices
        !
        DO iPeriod = 1, numShockPeriodsPrint
            !
            DO iAgent = 1, numAgents
                !
                AggrDevPricesShock(iPrice,iPeriod) = AggrDevPricesShock(iPrice,iPeriod)+avgPricesShock(iPrice,iPeriod,iAgent,iAgent)
                AggrDevProfitsShock(iPrice,iPeriod) = AggrDevProfitsShock(iPrice,iPeriod)+avgProfitsShock(iPrice,iPeriod,iAgent,iAgent)
                AggrDevPricesShockQ(iPrice,iPeriod) = AggrDevPricesShockQ(iPrice,iPeriod)+avgPricesShockQ(iPrice,iPeriod,iAgent,iAgent)
                AggrDevProfitsShockQ(iPrice,iPeriod) = AggrDevProfitsShockQ(iPrice,iPeriod)+avgProfitsShockQ(iPrice,iPeriod,iAgent,iAgent)
                !
                AggrDevPricesPercShock(iPrice,iPeriod) = AggrDevPricesPercShock(iPrice,iPeriod)+avgPricesPercShock(iPrice,iPeriod,iAgent,iAgent)
                AggrDevProfitsPercShock(iPrice,iPeriod) = AggrDevProfitsPercShock(iPrice,iPeriod)+avgProfitsPercShock(iPrice,iPeriod,iAgent,iAgent)
                AggrDevPricesPercShockQ(iPrice,iPeriod) = AggrDevPricesPercShockQ(iPrice,iPeriod)+avgPricesPercShockQ(iPrice,iPeriod,iAgent,iAgent)
                AggrDevProfitsPercShockQ(iPrice,iPeriod) = AggrDevProfitsPercShockQ(iPrice,iPeriod)+avgProfitsPercShockQ(iPrice,iPeriod,iAgent,iAgent)
                !
            END DO
            AggrNonDevPricesShock(iPrice,iPeriod) = (SUM(avgPricesShock(iPrice,iPeriod,:,:))-AggrDevPricesShock(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevPricesShock(iPrice,iPeriod) = AggrDevPricesShock(iPrice,iPeriod)/DBLE(numAgents)
            AggrNonDevProfitsShock(iPrice,iPeriod) = (SUM(avgProfitsShock(iPrice,iPeriod,:,:))-AggrDevProfitsShock(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevProfitsShock(iPrice,iPeriod) = AggrDevProfitsShock(iPrice,iPeriod)/DBLE(numAgents)
            AggrNonDevPricesShockQ(iPrice,iPeriod) = (SUM(avgPricesShockQ(iPrice,iPeriod,:,:))-AggrDevPricesShockQ(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevPricesShockQ(iPrice,iPeriod) = AggrDevPricesShockQ(iPrice,iPeriod)/DBLE(numAgents)
            AggrNonDevProfitsShockQ(iPrice,iPeriod) = (SUM(avgProfitsShockQ(iPrice,iPeriod,:,:))-AggrDevProfitsShockQ(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevProfitsShockQ(iPrice,iPeriod) = AggrDevProfitsShockQ(iPrice,iPeriod)/DBLE(numAgents)
            !
            AggrNonDevPricesPercShock(iPrice,iPeriod) = (SUM(avgPricesPercShock(iPrice,iPeriod,:,:))-AggrDevPricesPercShock(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevPricesPercShock(iPrice,iPeriod) = AggrDevPricesPercShock(iPrice,iPeriod)/DBLE(numAgents)
            AggrNonDevProfitsPercShock(iPrice,iPeriod) = (SUM(avgProfitsPercShock(iPrice,iPeriod,:,:))-AggrDevProfitsPercShock(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevProfitsPercShock(iPrice,iPeriod) = AggrDevProfitsPercShock(iPrice,iPeriod)/DBLE(numAgents)
            AggrNonDevPricesPercShockQ(iPrice,iPeriod) = (SUM(avgPricesPercShockQ(iPrice,iPeriod,:,:))-AggrDevPricesPercShockQ(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevPricesPercShockQ(iPrice,iPeriod) = AggrDevPricesPercShockQ(iPrice,iPeriod)/DBLE(numAgents)
            AggrNonDevProfitsPercShockQ(iPrice,iPeriod) = (SUM(avgProfitsPercShockQ(iPrice,iPeriod,:,:))-AggrDevProfitsPercShockQ(iPrice,iPeriod))/DBLE(numAgents*(numAgents-1))
            AggrDevProfitsPercShockQ(iPrice,iPeriod) = AggrDevProfitsPercShockQ(iPrice,iPeriod)/DBLE(numAgents)
            !
        END DO
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
    DO iPrice = 1, numPrices
        !
        DO iAgent = 1, numAgents
            !
            AggrDevPricesPost(iPrice) = AggrDevPricesPost(iPrice)+avgPricesPost(iPrice,iAgent,iAgent)
            AggrDevProfitsPost(iPrice) = AggrDevProfitsPost(iPrice)+avgProfitsPost(iPrice,iAgent,iAgent)
            AggrDevPricesPostQ(iPrice) = AggrDevPricesPostQ(iPrice)+avgPricesPostQ(iPrice,iAgent,iAgent)
            AggrDevProfitsPostQ(iPrice) = AggrDevProfitsPostQ(iPrice)+avgProfitsPostQ(iPrice,iAgent,iAgent)
            !
            AggrDevPricesPercPost(iPrice) = AggrDevPricesPercPost(iPrice)+avgPricesPercPost(iPrice,iAgent,iAgent)
            AggrDevProfitsPercPost(iPrice) = AggrDevProfitsPercPost(iPrice)+avgProfitsPercPost(iPrice,iAgent,iAgent)
            AggrDevPricesPercPostQ(iPrice) = AggrDevPricesPercPostQ(iPrice)+avgPricesPercPostQ(iPrice,iAgent,iAgent)
            AggrDevProfitsPercPostQ(iPrice) = AggrDevProfitsPercPostQ(iPrice)+avgProfitsPercPostQ(iPrice,iAgent,iAgent)
            !
        END DO
        AggrNonDevPricesPost(iPrice) = (SUM(avgPricesPost(iPrice,:,:))-AggrDevPricesPost(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesPost(iPrice) = AggrDevPricesPost(iPrice)/DBLE(numAgents)
        AggrNonDevProfitsPost(iPrice) = (SUM(avgProfitsPost(iPrice,:,:))-AggrDevProfitsPost(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsPost(iPrice) = AggrDevProfitsPost(iPrice)/DBLE(numAgents)
        AggrNonDevPricesPostQ(iPrice) = (SUM(avgPricesPostQ(iPrice,:,:))-AggrDevPricesPostQ(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesPostQ(iPrice) = AggrDevPricesPostQ(iPrice)/DBLE(numAgents)
        AggrNonDevProfitsPostQ(iPrice) = (SUM(avgProfitsPostQ(iPrice,:,:))-AggrDevProfitsPostQ(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsPostQ(iPrice) = AggrDevProfitsPostQ(iPrice)/DBLE(numAgents)
        !
        AggrNonDevPricesPercPost(iPrice) = (SUM(avgPricesPercPost(iPrice,:,:))-AggrDevPricesPercPost(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesPercPost(iPrice) = AggrDevPricesPercPost(iPrice)/DBLE(numAgents)
        AggrNonDevProfitsPercPost(iPrice) = (SUM(avgProfitsPercPost(iPrice,:,:))-AggrDevProfitsPercPost(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsPercPost(iPrice) = AggrDevProfitsPercPost(iPrice)/DBLE(numAgents)
        AggrNonDevPricesPercPostQ(iPrice) = (SUM(avgPricesPercPostQ(iPrice,:,:))-AggrDevPricesPercPostQ(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevPricesPercPostQ(iPrice) = AggrDevPricesPercPostQ(iPrice)/DBLE(numAgents)
        AggrNonDevProfitsPercPostQ(iPrice) = (SUM(avgProfitsPercPostQ(iPrice,:,:))-AggrDevProfitsPercPostQ(iPrice))/DBLE(numAgents*(numAgents-1))
        AggrDevProfitsPercPostQ(iPrice) = AggrDevProfitsPercPostQ(iPrice)/DBLE(numAgents)
        !
    END DO
    !
    ! Computing standard errors
    !
    avgPricesPreQ = SQRT(ABS((avgPricesPreQ-avgPricesPre**2)/DBLE(numGames)))
    avgProfitsPreQ = SQRT(ABS((avgProfitsPreQ-avgProfitsPre**2)/DBLE(numGames)))
    avgPricesShockQ = SQRT(ABS((avgPricesShockQ-avgPricesShock**2)/DBLE(numGames)))
    avgProfitsShockQ = SQRT(ABS((avgProfitsShockQ-avgProfitsShock**2)/DBLE(numGames)))
    avgPricesPostQ = SQRT(ABS((avgPricesPostQ-avgPricesPost**2)/DBLE(numGames)))
    avgProfitsPostQ = SQRT(ABS((avgProfitsPostQ-avgProfitsPost**2)/DBLE(numGames)))
    !
    AggrPricesPreQ = SQRT(ABS((AggrPricesPreQ-AggrPricesPre**2)/DBLE(numGames)))
    AggrProfitsPreQ = SQRT(ABS((AggrProfitsPreQ-AggrProfitsPre**2)/DBLE(numGames)))
    DO iPrice = 1, numPrices
        !
        DO iPeriod = 1, numShockPeriodsPrint
            !
            AggrNonDevPricesShockQ(iPrice,iPeriod) = SQRT(ABS((AggrNonDevPricesShockQ(iPrice,iPeriod)-AggrNonDevPricesShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            AggrDevPricesShockQ(iPrice,iPeriod) = SQRT(ABS((AggrDevPricesShockQ(iPrice,iPeriod)-AggrDevPricesShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            AggrNonDevProfitsShockQ(iPrice,iPeriod) = SQRT(ABS((AggrNonDevProfitsShockQ(iPrice,iPeriod)-AggrNonDevProfitsShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            AggrDevProfitsShockQ(iPrice,iPeriod) = SQRT(ABS((AggrDevProfitsShockQ(iPrice,iPeriod)-AggrDevProfitsShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            !
            AggrNonDevPricesPercShockQ(iPrice,iPeriod) = SQRT(ABS((AggrNonDevPricesPercShockQ(iPrice,iPeriod)-AggrNonDevPricesPercShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            AggrDevPricesPercShockQ(iPrice,iPeriod) = SQRT(ABS((AggrDevPricesPercShockQ(iPrice,iPeriod)-AggrDevPricesPercShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            AggrNonDevProfitsPercShockQ(iPrice,iPeriod) = SQRT(ABS((AggrNonDevProfitsPercShockQ(iPrice,iPeriod)-AggrNonDevProfitsPercShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            AggrDevProfitsPercShockQ(iPrice,iPeriod) = SQRT(ABS((AggrDevProfitsPercShockQ(iPrice,iPeriod)-AggrDevProfitsPercShock(iPrice,iPeriod)**2)/DBLE(numGames)))
            !
        END DO
        AggrNonDevPricesPostQ(iPrice) = SQRT(ABS((AggrNonDevPricesPostQ(iPrice)-AggrNonDevPricesPost(iPrice)**2)/DBLE(numGames)))
        AggrDevPricesPostQ(iPrice) = SQRT(ABS((AggrDevPricesPostQ(iPrice)-AggrDevPricesPost(iPrice)**2)/DBLE(numGames)))
        AggrNonDevProfitsPostQ(iPrice) = SQRT(ABS((AggrNonDevProfitsPostQ(iPrice)-AggrNonDevProfitsPost(iPrice)**2)/DBLE(numGames)))
        AggrDevProfitsPostQ(iPrice) = SQRT(ABS((AggrDevProfitsPostQ(iPrice)-AggrDevProfitsPost(iPrice)**2)/DBLE(numGames)))
        !
        AggrNonDevPricesPercPostQ(iPrice) = SQRT(ABS((AggrNonDevPricesPercPostQ(iPrice)-AggrNonDevPricesPercPost(iPrice)**2)/DBLE(numGames)))
        AggrDevPricesPercPostQ(iPrice) = SQRT(ABS((AggrDevPricesPercPostQ(iPrice)-AggrDevPricesPercPost(iPrice)**2)/DBLE(numGames)))
        AggrNonDevProfitsPercPostQ(iPrice) = SQRT(ABS((AggrNonDevProfitsPercPostQ(iPrice)-AggrNonDevProfitsPercPost(iPrice)**2)/DBLE(numGames)))
        AggrDevProfitsPercPostQ(iPrice) = SQRT(ABS((AggrDevProfitsPercPostQ(iPrice)-AggrDevProfitsPercPost(iPrice)**2)/DBLE(numGames)))    
        !
    END DO
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(100032,1) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            ((iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, &
                (iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), iPrice, iPrice = 1, numPrices), &
            (((jAgent, (jAgent, iAgent, iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, iPrice, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices), &
            (((jAgent, (jAgent, iAgent, iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, iPrice, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices)
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I2.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            <numPrices>('AggrPricePre ', <numShockPeriodsPrint>('AggrDevPriceShockPr', I0.2, 'Per', I3.3, ' '), 'AggrDevPricePostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('AggrNonDevPriceShockPr', I0.2, 'Per', I3.3, ' '), 'AggrNonDevPricePostPr', I0.2, ' '), &
            <numPrices>('seAggrPricePre ', <numShockPeriodsPrint>('seAggrDevPriceShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrDevPricePostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('seAggrNonDevPriceShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrNonDevPricePostPr', I0.2, ' '), &
            <numPrices>(<numShockPeriodsPrint>('AggrDevPricePercShockPr', I0.2, 'Per', I3.3, ' '), 'AggrDevPricePercPostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('AggrNonDevPricePercShockPr', I0.2, 'Per', I3.3, ' '), 'AggrNonDevPricePercPostPr', I0.2, ' '), &
            <numPrices>(<numShockPeriodsPrint>('seAggrDevPricePercShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrDevPricePercPostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('seAggrNonDevPricePercShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrNonDevPricePercPostPr', I0.2, ' '), &
            <numPrices>('AggrProfitPre ', <numShockPeriodsPrint>('AggrDevProfitShockPr', I0.2, 'Per', I3.3, ' '), 'AggrDevProfitPostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('AggrNonDevProfitShockPr', I0.2, 'Per', I3.3, ' '), 'AggrNonDevProfitPostPr', I0.2, ' '), &
            <numPrices>('seAggrProfitPre ', <numShockPeriodsPrint>('seAggrDevProfitShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrDevProfitPostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('seAggrNonDevProfitShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrNonDevProfitPostPr', I0.2, ' '), &
            <numPrices>(<numShockPeriodsPrint>('AggrDevProfitPercShockPr', I0.2, 'Per', I3.3, ' '), 'AggrDevProfitPercPostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('AggrNonDevProfitPercShockPr', I0.2, 'Per', I3.3, ' '), 'AggrNonDevProfitPercPostPr', I0.2, ' '), &
            <numPrices>(<numShockPeriodsPrint>('seAggrDevProfitPercShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrDevProfitPercPostPr', I0.2, ' ', &
                <numShockPeriodsPrint>('seAggrNonDevProfitPercShockPr', I0.2, 'Per', I3.3, ' '), 'seAggrNonDevProfitPercPostPr', I0.2, ' '), &
            <numPrices>(<numAgents>(<numAgents>('Ag', I1, 'avgPricePre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'avgPriceShockAg', I1, 'Pr', I0.2, 'Per', I3.3, ' '), &
                'Ag', I1, 'avgPricePostAg', I1, 'Pr', I0.2, ' '))), &
            <numPrices>(<numAgents>(<numAgents>('Ag', I1, 'avgProfitPre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'avgProfitShockAg', I1, 'Pr', I0.2, 'Per', I3.3, ' '), &
                'Ag', I1, 'avgProfitPostAg', I1, 'Pr', I0.2, ' '))) &
            )
        !
    END IF
    !
    WRITE(100032,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        (AggrPricesPre, AggrDevPricesShock(iPrice,:), AggrDevPricesPost(iPrice), AggrNonDevPricesShock(iPrice,:), AggrNonDevPricesPost(iPrice), iPrice = 1, numPrices), &
        (AggrPricesPreQ, AggrDevPricesShockQ(iPrice,:), AggrDevPricesPostQ(iPrice), AggrNonDevPricesShockQ(iPrice,:), AggrNonDevPricesPostQ(iPrice), iPrice = 1, numPrices), &
        (AggrDevPricesPercShock(iPrice,:), AggrDevPricesPercPost(iPrice), AggrNonDevPricesPercShock(iPrice,:), AggrNonDevPricesPercPost(iPrice), iPrice = 1, numPrices), &
        (AggrDevPricesPercShockQ(iPrice,:), AggrDevPricesPercPostQ(iPrice), AggrNonDevPricesPercShockQ(iPrice,:), AggrNonDevPricesPercPostQ(iPrice), iPrice = 1, numPrices), &
        (AggrProfitsPre, AggrDevProfitsShock(iPrice,:), AggrDevProfitsPost(iPrice), AggrNonDevProfitsShock(iPrice,:), AggrNonDevProfitsPost(iPrice), iPrice = 1, numPrices), &
        (AggrProfitsPreQ, AggrDevProfitsShockQ(iPrice,:), AggrDevProfitsPostQ(iPrice), AggrNonDevProfitsShockQ(iPrice,:), AggrNonDevProfitsPostQ(iPrice), iPrice = 1, numPrices), &
        (AggrDevProfitsPercShock(iPrice,:), AggrDevProfitsPercPost(iPrice), AggrNonDevProfitsPercShock(iPrice,:), AggrNonDevProfitsPercPost(iPrice), iPrice = 1, numPrices), &
        (AggrDevProfitsPercShockQ(iPrice,:), AggrDevProfitsPercPostQ(iPrice), AggrNonDevProfitsPercShockQ(iPrice,:), AggrNonDevProfitsPercPostQ(iPrice), iPrice = 1, numPrices), &
        (((avgPricesPre(jAgent), (avgPricesShock(iPrice,iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgPricesPost(iPrice,iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices), &
        (((avgProfitsPre(jAgent), (avgProfitsShock(iPrice,iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgProfitsPost(iPrice,iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices)
2       FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        <numPrices>(F12.7, 1X, <numShockPeriodsPrint>(F27.7,1X), F20.7, 1X, <numShockPeriodsPrint>(F30.7,1X), F23.7, 1X), &
        <numPrices>(F14.7, 1X, <numShockPeriodsPrint>(F29.7,1X), F22.7, 1X, <numShockPeriodsPrint>(F32.7,1X), F25.7, 1X), &
        <numPrices>(<numShockPeriodsPrint>(F31.7,1X), F24.7, 1X, <numShockPeriodsPrint>(F34.7,1X), F27.7, 1X), &
        <numPrices>(<numShockPeriodsPrint>(F33.7,1X), F26.7, 1X, <numShockPeriodsPrint>(F36.7,1X), F29.7, 1X), &
        <numPrices>(F13.7, 1X, <numShockPeriodsPrint>(F28.7,1X), F21.7, 1X, <numShockPeriodsPrint>(F31.7,1X), F24.7, 1X), &
        <numPrices>(F15.7, 1X, <numShockPeriodsPrint>(F30.7,1X), F23.7, 1X, <numShockPeriodsPrint>(F33.7,1X), F26.7, 1X), &
        <numPrices>(<numShockPeriodsPrint>(F32.7,1X), F25.7, 1X, <numShockPeriodsPrint>(F35.7,1X), F28.7, 1X), &
        <numPrices>(<numShockPeriodsPrint>(F34.7,1X), F27.7, 1X, <numShockPeriodsPrint>(F37.7,1X), F30.7, 1X), &
        <numPrices>(<numAgents>(<numAgents>(F14.7, 1X, <numShockPeriodsPrint>(F29.7, 1X), F22.7, 1X))), &
        <numPrices>(<numAgents>(<numAgents>(F15.7, 1X, <numShockPeriodsPrint>(F30.7, 1X), F23.7, 1X))) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeIRToAllAnalysis
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ImpulseResponseToAll
