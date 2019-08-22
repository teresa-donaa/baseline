MODULE DetailedImpulseResponseToAll
!
USE globals
USE QL_routines
USE ImpulseResponseToBR
!
! Computes disaggregated Impulse Responses to a one-period deviation to all prices
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeDetailedIRToAll ( iModel )
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
    INTEGER :: PeriodsLengthPre, PunishmentLength, PeriodsLengthPost, &
        VisitedStatesPre(numPeriods), VisitedStates(MAX(numShockPeriodsPrint,numPeriods)), &
        VisitedStatesTMP(numPeriods), SameCyclePrePost, &
        p(DepthState,numAgents), pPrime(numAgents), &
        iStatePre, iGame, iAgent, jAgent, iPrice, iPeriod, jPeriod, &
        OptimalStrategy(numStates,numAgents), LastObservedPrices(DepthState,numAgents), &
        indexShockState(LengthStates), i, j, PreCycleLength, CycleLength, IndexDynamicBR(numAgents)
	INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    INTEGER, DIMENSION(numPeriods,numAgents) :: IndPricesPre
    INTEGER, DIMENSION(numShockPeriodsPrint,numAgents) :: IndPricesShock, &
        StaticBRIndPrices, DynamicBRIndPrices
    REAL(8), DIMENSION(numPeriods,numAgents) :: visitedPrices, VisitedProfits, PricesPre, ProfitsPre
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents) :: PricesShock, ProfitsShock, &
        StaticBRPrices, DynamicBRPrices, OptStratQ, DynamicBRQ
    REAL(8), DIMENSION(numAgents) :: DeviationQ
    REAL(8), DIMENSION(numAgents) :: avgPricesPre, avgProfitsPre, avgPricesPost, avgProfitsPost
    LOGICAL :: FlagReturnedToState
    CHARACTER(len = 30) :: filename
    CHARACTER(len = 3) :: iPriceChar
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Detailed Impulse Responses to all prices'
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
    ! Writing header line in global output file
    !
    WRITE(100033,11) &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numAgents), (i, i = 1, numAgents)
11      FORMAT('    Game  DevTo_Price DevTo_IndPrice ', &
        <numAgents>(' NashProfit', I1, ' '), <numAgents>(' CoopProfit', I1, ' '), &
        'PreShock_CycleLength ', &
        <numAgents>('Avg_Pre_Price', I1, ' '), <numAgents>('Avg_Pre_Profit', I1, ' '), &
        'PreShock_NumInCycle ', &
        <numAgents>('PreShock_IndPrice', I1, ' '), <numAgents>('PreShock_Price', I1, ' '), <numAgents>('PreShock_Profit', I1, ' '), &
        'Shock_Agent Obs_Agent  Deviation_Q PunishmentLength SameCyclePrePost ', &
        <numShockPeriodsPrint>('Shock_IndPrice', I0.3, ' '), <numShockPeriodsPrint>('Shock_Price', I0.3, ' '), &
        <numShockPeriodsPrint>('Shock_Profit', I0.3, ' '), &
        <numShockPeriodsPrint>('StaticBR_IndPrice', I0.3, ' '), <numShockPeriodsPrint>('StaticBR_Price', I0.3, ' '), &
        <numShockPeriodsPrint>('DynamicBR_IndPrice', I0.3, ' '), <numShockPeriodsPrint>('DynamicBR_Price', I0.3, ' '), &
        <numShockPeriodsPrint>('OptStrat_Q', I0.3, ' '), <numShockPeriodsPrint>('DynamicBR_Q', I0.3, ' '), &
        'PostShock_CycleLength ', <numAgents>('Avg_Post_Price', I1, ' '), <numAgents>('Avg_Post_Profit', I1, ' '))
    !
    ! Beginning loop over deviation prices
    !
    DO iPrice = 1, numPrices    ! Start of loop over possible deviations
        !
        PRINT*, 'iPrice = ', iPrice
        !
        ! Opening output file
        !
        WRITE(iPriceChar,'(I0.3)') iPrice
        FileName = "IR_DevToPrice_" // iPriceChar // ".txt"
        OPEN(UNIT = 100,FILE = FileName)
        !
        ! Writing header line in separate deviation prices output files
        !
        WRITE(100,1) &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
            (i, i = 1, numShockPeriodsPrint), &
            (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
            (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
            (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
            (i, i = 1, numAgents), (i, i = 1, numAgents)
1       FORMAT('    Game ', &
            <numAgents>(' NashProfit', I1, ' '), <numAgents>(' CoopProfit', I1, ' '), &
            'PreShock_CycleLength ', &
            <numAgents>('Avg_Pre_Price', I1, ' '), <numAgents>('Avg_Pre_Profit', I1, ' '), &
            'PreShock_NumInCycle ', &
            <numAgents>('PreShock_IndPrice', I1, ' '), <numAgents>('PreShock_Price', I1, ' '), <numAgents>('PreShock_Profit', I1, ' '), &
            'Shock_Agent Obs_Agent  Deviation_Q PunishmentLength SameCyclePrePost ', &
            <numShockPeriodsPrint>('Shock_IndPrice', I0.3, ' '), <numShockPeriodsPrint>('Shock_Price', I0.3, ' '), &
            <numShockPeriodsPrint>('Shock_Profit', I0.3, ' '), &
            <numShockPeriodsPrint>('StaticBR_IndPrice', I0.3, ' '), <numShockPeriodsPrint>('StaticBR_Price', I0.3, ' '), &
            <numShockPeriodsPrint>('DynamicBR_IndPrice', I0.3, ' '), <numShockPeriodsPrint>('DynamicBR_Price', I0.3, ' '), &
            <numShockPeriodsPrint>('OptStrat_Q', I0.3, ' '), <numShockPeriodsPrint>('DynamicBR_Q', I0.3, ' '), &
            'PostShock_CycleLength ', <numAgents>('Avg_Post_Price', I1, ' '), <numAgents>('Avg_Post_Profit', I1, ' '))
        !
        ! Beginning loop over games
        !
        DO iGame = 1, numGames        ! Start of loop over games
            !
            OptimalStrategyVec = indexStrategies(:,iGame)
            LastStateVec = indexLastState(:,iGame)
            !
            OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
            IF (DepthState0 .EQ. 0) THEN
                !
                LastObservedPrices = OptimalStrategy
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
            VisitedStatesPre = 0
            PricesPre = 0.d0
            ProfitsPre = 0.d0
            avgPricesPre = 0.d0
            avgProfitsPre = 0.d0
            p = LastObservedPrices
            pPrime = OptimalStrategy(computeStateNumber(p),:)
            DO iPeriod = 1, numPeriods
                !
                IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                p(1,:) = pPrime
                LastObservedPrices = p
                VisitedStatesPre(iPeriod) = computeStateNumber(p)
                DO iAgent = 1, numAgents
                    !
                    PricesPre(iPeriod,iAgent) = PricesGrids(pPrime(iAgent),iAgent)
                    ProfitsPre(iPeriod,iAgent) = PI(computeActionNumber(pPrime),iAgent)
                    IndPricesPre(iPeriod,iAgent) = pPrime(iAgent)
                    !
                END DO
                !
                ! Check if the state has already been visited
                !
                IF ((iPeriod .GE. 2) .AND. (ANY(VisitedStatesPre(:iPeriod-1) .EQ. VisitedStatesPre(iPeriod)))) EXIT
                !
                ! Update pPrime and iterate
                !
                pPrime = OptimalStrategy(VisitedStatesPre(iPeriod),:)
                !
            END DO
            !
            PeriodsLengthPre = &
                iPeriod-MINVAL(MINLOC((VisitedStatesPre(:iPeriod-1)-VisitedStatesPre(iPeriod))**2))
            !
            VisitedStatesPre(:PeriodsLengthPre) = VisitedStatesPre(iPeriod-PeriodsLengthPre+1:iPeriod)
            VisitedStatesPre(PeriodsLengthPre+1:) = 0
            PricesPre(:PeriodsLengthPre,:) = PricesPre(iPeriod-PeriodsLengthPre+1:iPeriod,:)
            PricesPre(PeriodsLengthPre+1:,:) = 0.d0
            ProfitsPre(:PeriodsLengthPre,:) = ProfitsPre(iPeriod-PeriodsLengthPre+1:iPeriod,:)
            ProfitsPre(PeriodsLengthPre+1:,:) = 0.d0
            IndPricesPre(:PeriodsLengthPre,:) = IndPricesPre(iPeriod-PeriodsLengthPre+1:iPeriod,:)
            IndPricesPre(PeriodsLengthPre+1:,:) = 0.d0
            avgPricesPre = SUM(PricesPre(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
            avgProfitsPre = SUM(ProfitsPre(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
            !
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Shock period analysis
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            DO iStatePre = 1, PeriodsLengthPre      ! Start of loop over pre-shock cycle states
                !
                DO iAgent = 1, numAgents            ! Start of loop over deviating agent 
                    !
                    PricesShock = 0.d0
                    ProfitsShock = 0.d0
                    StaticBRPrices = 0.d0
                    DynamicBRPrices = 0.d0
                    IndPricesShock = 0
                    StaticBRIndPrices = 0
                    DynamicBRIndPrices = 0
                    OptStratQ = 0.d0
                    DynamicBRQ = 0.d0
                    DeviationQ = 0.d0
                    avgPricesPost = 0.d0
                    avgProfitsPost = 0.d0
                    VisitedStates = 0
                    !
                    p = RESHAPE(convertNumberBase(VisitedStatesPre(iStatePre)-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
                    !
                    ! pPrime and Q computed from Dynamic Best Reply
                    !
                    DO jAgent = 1, numAgents
                        !
                        CALL ComputeDynamicBestResponse(OptimalStrategy,VisitedStatesPre(iStatePre),jAgent,delta, &
                            IndexDynamicBR(jAgent),DynamicBRQ(1,jAgent))
                        DynamicBRIndPrices(1,jAgent) = IndexDynamicBR(jAgent)
                        !
                    END DO
                    !
                    ! pPrime and Q computed from Optimal Strategy
                    !
                    pPrime = OptimalStrategy(VisitedStatesPre(iStatePre),:)
                    DO jAgent = 1, numAgents
                        !
                        CALL computeQcell(OptimalStrategy,VisitedStatesPre(iStatePre),pPrime(jAgent),jAgent,delta, &
                            OptStratQ(1,jAgent),VisitedStatesTMP,PreCycleLength,CycleLength,iPeriod)                                !
                        !
                    END DO
                    !
                    ! Actual price selection in the first shock period:
                    ! Agent "iAgent" selects each price in turn,
                    ! The other agents stick to the strategy at convergence
                    !
                    pPrime(iAgent) = iPrice
                    DO jAgent = 1, numAgents
                        !
                        CALL computeQcell(OptimalStrategy,VisitedStatesPre(iStatePre),pPrime(jAgent),jAgent,delta, &
                            DeviationQ(jAgent),VisitedStatesTMP,PreCycleLength,CycleLength,iPeriod)                                !
                        !
                    END DO
                    !
                    ! Start of loop over shock periods
                    !
                    flagReturnedToState = .FALSE.
                    DO iPeriod = 1, MAX(numShockPeriodsPrint,numPeriods)
                        !
                        IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                        p(1,:) = pPrime
                        DO jAgent = 1, numAgents
                            !
                            IF (iPeriod .LE. numShockPeriodsPrint) THEN
                                !
                                indPricesShock(iPeriod,jAgent) = pPrime(jAgent)
                                PricesShock(iPeriod,jAgent) = PricesGrids(pPrime(jAgent),jAgent)
                                ProfitsShock(iPeriod,jAgent) = PI(computeActionNumber(pPrime),jAgent)
                                StaticBRIndPrices(iPeriod,jAgent) = ComputeStaticBestResponse(jAgent,p,PI(:,jAgent))
                                StaticBRPrices(iPeriod,jAgent) = PricesGrids(StaticBRIndPrices(iPeriod,jAgent),jAgent)
                                IF (iPeriod .GE. 2) THEN
                                    !
                                    ! pPrime and Q computed from Dynamic Best Reply
                                    !
                                    CALL ComputeDynamicBestResponse(OptimalStrategy,VisitedStates(iPeriod-1),jAgent,delta, &
                                        IndexDynamicBR(jAgent),DynamicBRQ(iPeriod,jAgent))
                                    !
                                    ! pPrime and Q computed from Optimal Strategy
                                    !
                                     CALL computeQcell(OptimalStrategy,VisitedStates(iPeriod-1),pPrime(jAgent),jAgent,delta, &
                                        OptStratQ(iPeriod,jAgent),VisitedStatesTMP,PreCycleLength,CycleLength,jPeriod)                                !
                                    !
                                END IF
                                DynamicBRIndPrices(iPeriod,jAgent) = IndexDynamicBR(jAgent)
                                DynamicBRPrices(iPeriod,jAgent) = PricesGrids(IndexDynamicBR(jAgent),jAgent)
                                !
                            END IF
                            !
                        END DO
                        VisitedStates(iPeriod) = computeStateNumber(p)
                        !
                        ! Check if the state has already been visited
                        ! Case 1: the state retuns to one of the states in the pre-shock cycle
                        !
                        IF ((.NOT.(flagReturnedToState)) .AND. &
                            (ANY(VisitedStatesPre(:PeriodsLengthPre) .EQ. VisitedStates(iPeriod)))) THEN
                            !
                            indexShockState = RESHAPE(p,(/ LengthStates /))
                            PunishmentLength = MAX(0,iPeriod-2)
                            SameCyclePrePost = 1
                            flagReturnedToState = .TRUE.
                            !
                        END IF
                        !
                        ! Case 2: after some time, the state starts cycling among a new set of states
                        !
                        IF ((iPeriod .GE. 2) .AND. (.NOT.(flagReturnedToState)) .AND. &
                            (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) THEN
                            !
                            indexShockState = RESHAPE(p,(/ LengthStates /))
                            PunishmentLength = &
                                iPeriod-MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
                            SameCyclePrePost = 0
                            flagReturnedToState = .TRUE.
                            !
                        END IF
                        pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
                        !
                    END DO
                    !
                    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Post-shock period analysis
                    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    !
                    visitedPrices = 0.d0
                    VisitedProfits = 0.d0
                    VisitedStates = 0
                    p = RESHAPE(indexShockState, (/ DepthState,numAgents /) )
                    pPrime = OptimalStrategy(computeStateNumber(p),:)
                    !
                    DO iPeriod = 1, numPeriods
                        !
                        IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                        p(1,:) = pPrime
                        VisitedStates(iPeriod) = computeStateNumber(p)
                        DO jAgent = 1, numAgents
                            !
                            visitedPrices(iPeriod,jAgent) = PricesGrids(pPrime(jAgent),jAgent)
                            VisitedProfits(iPeriod,jAgent) = PI(computeActionNumber(pPrime),jAgent)
                            !
                        END DO
                        !
                        ! Check if the state has already been visited
                        !
                        IF ((iPeriod .GE. 2) .AND. (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) EXIT
                        !
                        ! Update pPrime and iterate
                        !
                        pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
                        !
                    END DO
                    !
                    PeriodsLengthPost = &
                        iPeriod-MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
                    avgPricesPost = SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)
                    avgProfitsPost = SUM(VisitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPost)
                    !
                    ! Printing results to output files
                    !
                    DO jAgent = 1, numAgents    ! Start of loop over observed agent
                        !
                        WRITE(100033,12) iGame, PricesGrids(iPrice,iAgent), iPrice, &
                            NashProfits, CoopProfits, &
                            PeriodsLengthPre, &
                            avgPricesPre, avgProfitsPre, &
                            iStatePre, &
                            IndPricesPre(iStatePre,:), PricesPre(iStatePre,:), ProfitsPre(iStatePre,:), &
                            iAgent, jAgent, DeviationQ(jAgent), PunishmentLength, SameCyclePrePost, &
                            IndPricesShock(:,jAgent), PricesShock(:,jAgent), &
                            ProfitsShock(:,jAgent), &
                            StaticBRIndPrices(:,jAgent), StaticBRPrices(:,jAgent), &
                            DynamicBRIndPrices(:,jAgent), DynamicBRPrices(:,jAgent), &
                            OptStratQ(:,jAgent), DynamicBRQ(:,jAgent), &
                            PeriodsLengthPost, avgPricesPost, avgProfitsPost
12                      FORMAT(I8, 1X, F12.5, 1X, I14, 1X, &
                            <2*numAgents>(F12.5, 1X), &
                            I20, 1X, &
                            <numAgents>(F14.5, 1X), <numAgents>(F15.5, 1X), &
                            I19, 1X, &
                            <numAgents>(I18, 1X), <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X), &
                            I11, 1X, I9, 1X, F12.5, 1X, I16, 1X, I16, 1X, &
                            <numShockPeriodsPrint>(I17, 1X), <numShockPeriodsPrint>(F14.5, 1X), &
                            <numShockPeriodsPrint>(F15.5, 1X), &
                            <numShockPeriodsPrint>(I20, 1X), <numShockPeriodsPrint>(F17.5, 1X), &
                            <numShockPeriodsPrint>(I21, 1X), <numShockPeriodsPrint>(F18.5, 1X), &
                            <numShockPeriodsPrint>(F13.5, 1X), <numShockPeriodsPrint>(F14.5, 1X), &
                            I21, 1X, <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X))
                        !
                        WRITE(100,2) iGame, &
                            NashProfits, CoopProfits, &
                            PeriodsLengthPre, &
                            avgPricesPre, avgProfitsPre, &
                            iStatePre, &
                            IndPricesPre(iStatePre,:), PricesPre(iStatePre,:), ProfitsPre(iStatePre,:), &
                            iAgent, jAgent, DeviationQ(jAgent), PunishmentLength, SameCyclePrePost, &
                            IndPricesShock(:,jAgent), PricesShock(:,jAgent), &
                            ProfitsShock(:,jAgent), &
                            StaticBRIndPrices(:,jAgent), StaticBRPrices(:,jAgent), &
                            DynamicBRIndPrices(:,jAgent), DynamicBRPrices(:,jAgent), &
                            OptStratQ(:,jAgent), DynamicBRQ(:,jAgent), &
                            PeriodsLengthPost, avgPricesPost, avgProfitsPost
2                       FORMAT(I8, 1X, &
                            <2*numAgents>(F12.5, 1X), &
                            I20, 1X, &
                            <numAgents>(F14.5, 1X), <numAgents>(F15.5, 1X), &
                            I19, 1X, &
                            <numAgents>(I18, 1X), <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X), &
                            I11, 1X, I9, 1X, F12.5, 1X, I16, 1X, I16, 1X, &
                            <numShockPeriodsPrint>(I17, 1X), <numShockPeriodsPrint>(F14.5, 1X), &
                            <numShockPeriodsPrint>(F15.5, 1X), &
                            <numShockPeriodsPrint>(I20, 1X), <numShockPeriodsPrint>(F17.5, 1X), &
                            <numShockPeriodsPrint>(I21, 1X), <numShockPeriodsPrint>(F18.5, 1X), &
                            <numShockPeriodsPrint>(F13.5, 1X), <numShockPeriodsPrint>(F14.5, 1X), &
                            I21, 1X, <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X))
                        !
                    END DO                  ! End of loop over observed agent
                    !
                END DO                      ! End of loop over deviating agent
                !
            END DO                          ! End of loop over pre-shock cycle states
            !
        END DO                              ! End of loop over games
        !
        ! Closing output file
        !
        CLOSE(UNIT = 100)
        !
    END DO                                  ! End of loop over deviation prices
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeDetailedIRToAll
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE DetailedImpulseResponseToAll
