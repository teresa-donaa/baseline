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
        visitedStatesPre(numStates), visitedStates(numStates), &
        p(depthState,numAgents), pPrime(numAgents), &
        iPeriod, iAgent, jAgent, iPrice, tmp1(numAgents), &
        iGame, optimalStrategy(numStates,numAgents), lastObservedState(depthState,numAgents), &
        indexShockState(lengthStates), numPeriods, iThres, i, j
    INTEGER :: FreqPeriodLengthPre(numThresPeriodsLength)
    INTEGER, DIMENSION(numAgents,numThresPeriodsLength) :: FreqPeriodLengthShock, FreqPeriodLengthPost
    INTEGER :: FreqPunishmentStrategy(numAgents,0:numThresPeriodsLength)
    REAL(8) :: visitedPrices(numStates,numAgents), visitedProfits(numStates,numAgents), &
        pNash
    REAL(8), DIMENSION(numAgents) :: avgPricesPre, avgProfitsPre, avgPricesPreQ, avgProfitsPreQ
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents,numAgents) :: &
        avgPricesShock, avgProfitsShock, avgPricesShockQ, avgProfitsShockQ
    REAL(8), DIMENSION(numAgents,numAgents) :: avgPricesPost, avgProfitsPost, avgPricesPostQ, avgProfitsPostQ
    LOGICAL :: FlagReturnedToState
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(lengthStates)
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
    numPeriods = numStates          ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
    FreqPeriodLengthPre = 0
    FreqPeriodLengthShock = 0
    FreqPeriodLengthPost = 0
    FreqPunishmentStrategy = 0
    avgPricesPre = 0.d0
    avgPricesShock = 0.d0
    avgPricesPost = 0.d0
    avgProfitsPre = 0.d0
    avgProfitsShock = 0.d0
    avgProfitsPost = 0.d0
    avgPricesPreQ = 0.d0
    avgPricesShockQ = 0.d0
    avgPricesPostQ = 0.d0
    avgProfitsPreQ = 0.d0
    avgProfitsShockQ = 0.d0
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
    22  FORMAT(<lengthStates>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    PRINT*, 'Read indexLastState'
    CLOSE(UNIT = 999)                   ! Close indexLastState file
    !
    ! Beginning loop over games
    !
    !$omp parallel do &
    !$omp private(OptimalStrategy,lastObservedState,visitedStatesPre,visitedPrices, &
    !$omp   visitedProfits,p,pPrime,iPeriod,iAgent,pNash,OptimalStrategyVec,LastStateVec, &
    !$omp   visitedStates,iPrice,tmp1,flagReturnedToState,jAgent,indexShockState) &
    !$omp firstprivate(PI,PricesGrids,i) &
    !$omp reduction(+ : FreqPeriodLengthPre, &
    !$omp   avgPricesPre,avgProfitsPre,avgPricesShock,avgProfitsShock,avgPricesPost,avgProfitsPost, &
    !$omp   avgPricesPreQ,avgProfitsPreQ,avgPricesShockQ,avgProfitsShockQ,avgPricesPostQ,avgProfitsPostQ, &
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
        lastObservedState = RESHAPE(LastStateVec, (/ depthState,numAgents /) )
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Pre-shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        visitedStatesPre = 0
        visitedPrices = 0.d0
        visitedProfits = 0.d0
        p = lastObservedState
        pPrime = optimalStrategy(computeStateNumber(p),:)
        DO iPeriod = 1, numPeriods
            !
            IF (depthState .GT. 1) p(2:depthState,:) = p(1:depthState-1,:)
            p(1,:) = pPrime
            lastObservedState = p
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
        DO iAgent = 1, numAgents        ! Start of loop aver shocking agent 
            !
            visitedStates = 0
            p = lastObservedState
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
            DO iPeriod = 1, numPeriods
                !
                IF (depthState .GT. 1) p(2:depthState,:) = p(1:depthState-1,:)
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
                        avgProfitsShock(iPeriod,iAgent,jAgent) = & 
                            avgProfitsShock(iPeriod,iAgent,jAgent)+PI(computeActionNumber(pPrime),jAgent)
                        avgProfitsShockQ(iPeriod,iAgent,jAgent) = & 
                            avgProfitsShockQ(iPeriod,iAgent,jAgent)+PI(computeActionNumber(pPrime),jAgent)**2
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
                    indexShockState = RESHAPE(p,(/ lengthStates /))
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
                    indexShockState = RESHAPE(p,(/ lengthStates /))
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
            p = RESHAPE(indexShockState, (/ depthState,numAgents /) )
            pPrime = optimalStrategy(computeStateNumber(p),:)
            IF (computeImpulseResponseToNash .EQ. -1) &            
                pPrime(iAgent) = pNash          ! The deviation to the Nash price is permanent!
            !            
            DO iPeriod = 1, numPeriods
                !
                IF (depthState .GT. 1) p(2:depthState,:) = p(1:depthState-1,:)
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
                    avgPricesPost(iAgent,:) = avgPricesPost(iAgent,:)+ &
                        SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)
                    avgPricesPostQ(iAgent,:) = avgPricesPostQ(iAgent,:)+ &
                        (SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost))**2
                    avgProfitsPost(iAgent,:) = avgProfitsPost(iAgent,:)+ &
                        SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)
                    avgProfitsPostQ(iAgent,:) = avgProfitsPostQ(iAgent,:)+ &
                        (SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost))**2
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
        WRITE(100031,1) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
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
