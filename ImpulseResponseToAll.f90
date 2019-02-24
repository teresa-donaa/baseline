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
        p(DepthState,numAgents), pPrime(numAgents), &
        iGame, iAgent, iPrice, iPeriod, jAgent, &
        optimalStrategy(numStates,numAgents), LastObservedPrices(DepthState,numAgents), &
        indexShockState(LengthStates), numPeriods, iThres, i, j
    REAL(8), DIMENSION(numStates+1,numAgents) :: visitedPrices, visitedProfits
    REAL(8), DIMENSION(numAgents) :: avgPricesPre, avgProfitsPre, avgPricesPreQ, avgProfitsPreQ
    REAL(8), DIMENSION(numPrices,numShockPeriodsPrint,numAgents,numAgents) :: &
        avgPricesShock, avgProfitsShock, avgPricesShockQ, avgProfitsShockQ
    REAL(8), DIMENSION(numPrices,numAgents,numAgents) :: &
        avgPricesPost, avgProfitsPost, avgPricesPostQ, avgProfitsPostQ
    LOGICAL :: FlagReturnedToState
	INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Impulse Response functions to all prices'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !
    numPeriods = numStates+1        ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
    avgPricesPre = 0.d0
    avgProfitsPre = 0.d0
    avgPricesShock = 0.d0
    avgProfitsShock = 0.d0
    avgPricesPost = 0.d0
    avgProfitsPost = 0.d0
    avgPricesPreQ = 0.d0
    avgProfitsPreQ = 0.d0
    avgPricesShockQ = 0.d0
    avgProfitsShockQ = 0.d0
    avgPricesPostQ = 0.d0
    avgProfitsPostQ = 0.d0
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
    !$omp private(OptimalStrategy,LastObservedPrices,visitedStatesPre,visitedPrices,visitedProfits, &
    !$omp   visitedStates,p,pPrime,iPrice,iAgent,iPeriod,jAgent,indexShockState, &
    !$omp   flagReturnedToState,PeriodsLengthPre,PeriodsLengthPost,OptimalStrategyVec,LastStateVec,i) &
    !$omp firstprivate(PI,PricesGrids) &
    !$omp reduction(+ : avgPricesPre,avgProfitsPre,avgPricesShock,avgProfitsShock,avgPricesPost,avgProfitsPost, &
    !$omp   avgPricesPreQ,avgProfitsPreQ,avgPricesShockQ,avgProfitsShockQ,avgPricesPostQ,avgProfitsPostQ)
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
                avgPricesPre = avgPricesPre+ &
                    SUM(visitedPrices(iPeriod-PeriodsLengthPre+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPre)
                avgPricesPreQ = avgPricesPreQ+ &
                    (SUM(visitedPrices(iPeriod-PeriodsLengthPre+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPre))**2
                avgProfitsPre = avgProfitsPre+ &
                    SUM(visitedProfits(iPeriod-PeriodsLengthPre+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPre)
                avgProfitsPreQ = avgProfitsPreQ+ &
                    (SUM(visitedProfits(iPeriod-PeriodsLengthPre+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPre))**2
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
        DO iAgent = 1, numAgents        ! Start of loop over shocking agent 
            !
            DO iPrice = 1, numPrices    ! Start of loop over possible deviations
                !
                visitedStates = 0
                p = LastObservedPrices
                !
                ! Price selection in shock period:
                ! Agent "iAgent" selects the each price in turn,
                ! The other agents stick to the strategy at convergence
                !
                pPrime = optimalStrategy(computeStateNumber(p),:)
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
                            avgPricesShock(iPrice,iPeriod,iAgent,jAgent) = & 
                                avgPricesShock(iPrice,iPeriod,iAgent,jAgent)+PricesGrids(pPrime(jAgent),jAgent)
                            avgPricesShockQ(iPrice,iPeriod,iAgent,jAgent) = & 
                                avgPricesShockQ(iPrice,iPeriod,iAgent,jAgent)+PricesGrids(pPrime(jAgent),jAgent)**2
                            avgProfitsShock(iPrice,iPeriod,iAgent,jAgent) = & 
                                avgProfitsShock(iPrice,iPeriod,iAgent,jAgent)+PI(computeActionNumber(pPrime),jAgent)
                            avgProfitsShockQ(iPrice,iPeriod,iAgent,jAgent) = & 
                                avgProfitsShockQ(iPrice,iPeriod,iAgent,jAgent)+PI(computeActionNumber(pPrime),jAgent)**2
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
                    IF ((iPeriod .GE. 2) .AND. (ANY(visitedStates(:iPeriod-1) .EQ. visitedStates(iPeriod)))) THEN
                        !
                        PeriodsLengthPost = &
                            iPeriod-MINVAL(MINLOC((visitedStates(:iPeriod-1)-visitedStates(iPeriod))**2))
                        avgPricesPost(iPrice,iAgent,:) = avgPricesPost(iPrice,iAgent,:)+ &
                            SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                                DBLE(PeriodsLengthPost)
                        avgPricesPostQ(iPrice,iAgent,:) = avgPricesPostQ(iPrice,iAgent,:)+ &
                            (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                                DBLE(PeriodsLengthPost))**2
                        avgProfitsPost(iPrice,iAgent,:) = avgProfitsPost(iPrice,iAgent,:)+ &
                            SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                                DBLE(PeriodsLengthPost)
                        avgProfitsPostQ(iPrice,iAgent,:) = avgProfitsPostQ(iPrice,iAgent,:)+ &
                            (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                                DBLE(PeriodsLengthPost))**2
                        EXIT
                        !
                    END IF
                    pPrime = optimalStrategy(visitedStates(iPeriod),:)
                    !
                END DO
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
    avgPricesPreQ = SQRT((avgPricesPreQ/DBLE(numGames)-avgPricesPre**2)/DBLE(numGames))
    avgProfitsPreQ = SQRT((avgProfitsPreQ/DBLE(numGames)-avgProfitsPre**2)/DBLE(numGames))
    avgPricesShockQ = SQRT((avgPricesShockQ/DBLE(numGames)-avgPricesShock**2)/DBLE(numGames))
    avgProfitsShockQ = SQRT((avgProfitsShockQ/DBLE(numGames)-avgProfitsShock**2)/DBLE(numGames))
    avgPricesPostQ = SQRT((avgPricesPostQ/DBLE(numGames)-avgPricesPost**2)/DBLE(numGames))
    avgProfitsPostQ = SQRT((avgProfitsPostQ/DBLE(numGames)-avgProfitsPost**2)/DBLE(numGames))
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
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (((jAgent, (jAgent, iAgent, iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, iPrice, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices), &
            (((jAgent, (jAgent, iAgent, iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, iPrice, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices), &
            (((jAgent, (jAgent, iAgent, iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, iPrice, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices), &
            (((jAgent, (jAgent, iAgent, iPrice, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, iPrice, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices)
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>(' MExplPar', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I2.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            <numPrices>(<numAgents>(<numAgents>('Ag', I1, 'avgPricePre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'avgPriceShockAg', I1, 'Pr', I0.2, 'Per', I3.3, ' '), &
                'Ag', I1, 'avgPricePostAg', I1, 'Pr', I0.2, ' '))), &
            <numPrices>(<numAgents>(<numAgents>('seAg', I1, 'avgPricePre', ' ', &
                <numShockPeriodsPrint>('seAg', I1, 'avgPriceShockAg', I1, 'Pr', I0.2, 'Per', I3.3, ' '), &
                'seAg', I1, 'avgPricePostAg', I1, 'Pr', I0.2, ' '))), &
            <numPrices>(<numAgents>(<numAgents>('Ag', I1, 'avgProfitPre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'avgProfitShockAg', I1, 'Pr', I0.2, 'Per', I3.3, ' '), &
                'Ag', I1, 'avgProfitPostAg', I1, 'Pr', I0.2, ' '))), &
            <numPrices>(<numAgents>(<numAgents>('seAg', I1, 'avgProfitPre', ' ', &
                <numShockPeriodsPrint>('seAg', I1, 'avgProfitShockAg', I1, 'Pr', I0.2, 'Per', I3.3, ' '), &
                'seAg', I1, 'avgProfitPostAg', I1, 'Pr', I0.2, ' '))) &
            )
        !
    END IF
    !
    WRITE(100032,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        (((avgPricesPre(jAgent), (avgPricesShock(iPrice,iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgPricesPost(iPrice,iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices), &
        (((avgProfitsPre(jAgent), (avgProfitsShock(iPrice,iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), avgProfitsPost(iPrice,iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), iPrice = 1, numPrices)
2       FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.3, 1X), &
        <6*numAgents>(F10.7, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        <numPrices>(<numAgents>(<numAgents>(F14.7, 1X, <numShockPeriodsPrint>(F29.7, 1X), F22.7, 1X))), &
        <numPrices>(<numAgents>(<numAgents>(F16.7, 1X, <numShockPeriodsPrint>(F31.7, 1X), F24.7, 1X))), &
        <numPrices>(<numAgents>(<numAgents>(F15.7, 1X, <numShockPeriodsPrint>(F30.7, 1X), F23.7, 1X))), &
        <numPrices>(<numAgents>(<numAgents>(F17.7, 1X, <numShockPeriodsPrint>(F32.7, 1X), F25.7, 1X))) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeIRToAllAnalysis
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ImpulseResponseToAll
