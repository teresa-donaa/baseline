MODULE DetailedImpulseResponseToAll
!
USE globals
USE QL_routines
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
    INTEGER :: PeriodsLengthPre, PeriodsLengthPost, &
        visitedStatesPre(numStates+1), visitedStates(MAX(numShockPeriodsPrint,numStates+1)), &
        p(DepthState,numAgents), pPrime(numAgents), &
        iStatePre, iGame, iAgent, iPrice, iPeriod, jAgent, &
        optimalStrategy(numStates,numAgents), LastObservedPrices(DepthState,numAgents), &
        indexShockState(LengthStates), numPeriods, i, j
	INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    INTEGER, DIMENSION(numStates+1,numAgents) :: indexPricesPre
    REAL(8) :: IC_Condition, avgIC_Condition
    REAL(8), DIMENSION(numStates+1,numAgents) :: visitedPrices, visitedProfits, PricesPre, ProfitsPre
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents) :: PricesShock, ProfitsShock
    REAL(8), DIMENSION(numAgents) :: avgPricesPre, avgProfitsPre, avgPricesPost, avgProfitsPost
    LOGICAL :: FlagReturnedToState
    CHARACTER(len = 30) :: filename
    CHARACTER(len = 3) :: iPriceChar
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Detailed Impulse Responses to all prices'
    !
    ! Initializing variables
    !
    numPeriods = numStates+1        ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
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
        ! Writing header line in output file
        !
        WRITE(100,1) (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numShockPeriodsPrint), (i, i = 1, numShockPeriodsPrint), &
            (i, i = 1, numAgents), (i, i = 1, numAgents)
1       FORMAT('    Game ', &
            <numAgents>(' NashProfit', I1, ' '), <numAgents>(' CoopProfit', I1, ' '), &
            'PreShock_CycleLength ', <numAgents>('Avg_Pre_Price', I1, ' '), <numAgents>('Avg_Pre_Profit', I1, ' '), 'PreShock_NumInCycle ', &
            <numAgents>('PreShock_Price', I1, ' '), <numAgents>('PreShock_Profit', I1, ' '), &
            'Shock_Agent Obs_Agent ', &
            <numShockPeriodsPrint>('Shock_Price', I0.3, ' '), <numShockPeriodsPrint>('Shock_Profit', I0.3, ' '), &
            'PostShock_CycleLength ', <numAgents>('Avg_Post_Price', I1, ' '), <numAgents>('Avg_Post_Profit', I1, ' '), &
            'IC_Condition avgIC_Condition ')
        !
        ! Beginning loop over games
        !
        DO iGame = 1, numGames        ! Start of loop over games
            !
            OptimalStrategyVec = indexStrategies(:,iGame)
            LastStateVec = indexLastState(:,iGame)
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
            PricesPre = 0.d0
            ProfitsPre = 0.d0
            avgPricesPre = 0.d0
            avgProfitsPre = 0.d0
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
                    avgPricesPost = 0.d0
                    avgProfitsPost = 0.d0
                    visitedStates = 0
                    pPrime = indexPricesPre(iStatePre,:)
                    !
                    ! Price selection in shock period:
                    ! Agent "iAgent" selects each price in turn,
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
                                PricesShock(iPeriod,jAgent) = PricesGrids(pPrime(jAgent),jAgent)
                                ProfitsShock(iPeriod,jAgent) = PI(computeActionNumber(pPrime),jAgent)
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
                    avgPricesPost = SUM(visitedPrices(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                            DBLE(PeriodsLengthPost)
                    avgProfitsPost = SUM(visitedProfits(iPeriod-PeriodsLengthPost+1:iPeriod,:),DIM = 1)/ &
                        DBLE(PeriodsLengthPost)
                    !
                    ! Printing results in output file
                    !
                    DO jAgent = 1, numAgents
                        !
                        IC_Condition = SUM(ProfitsShock(:,jAgent))-numShockPeriodsPrint*ProfitsPre(iStatePre,jAgent)
                        avgIC_Condition = SUM(ProfitsShock(:,jAgent))-numShockPeriodsPrint*avgProfitsPre(jAgent)
                        WRITE(100,2) iGame, NashProfits, CoopProfits, &
                            PeriodsLengthPre, avgPricesPre, avgProfitsPre, iStatePre, PricesPre(iStatePre,:), ProfitsPre(iStatePre,:), &
                            iAgent, jAgent, &
                            PricesShock(:,jAgent), ProfitsShock(:,jAgent), &
                            PeriodsLengthPost, avgPricesPost, avgProfitsPost, IC_Condition, avgIC_Condition
2                       FORMAT(I8, 1X, <2*numAgents>(F12.5, 1X), &
                            I20, 1X, <numAgents>(F14.5, 1X), <numAgents>(F15.5, 1X), I19, 1X, <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X), &
                            I11, 1X, I9, 1X, &
                            <numShockPeriodsPrint>(F14.5, 1X), <numShockPeriodsPrint>(F15.5, 1X), &
                            I21, 1X, <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X), F12.5, 1X, F15.5, 1X)
                        !
                    END DO
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
