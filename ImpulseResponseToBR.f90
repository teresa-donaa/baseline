MODULE ImpulseResponseToBR
!
USE globals
USE QL_routines
USE EquilibriumCheck
!
! Computes Impulse Response analysis to a one-period deviation to best response
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeIRAnalysisToBR ( iModel )
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
    !
    INTEGER :: i, j, iGame, iCycle, iAgent, jAgent, iPeriod, iStatePre, &
        PosThres, PeriodsLengthPre, DevPrice
    INTEGER, DIMENSION(numGames) :: CycleLength
    INTEGER, DIMENSION(numPeriods) :: VisitedStatesPre
    INTEGER, DIMENSION(numThresPeriodsLength) :: FreqPreLength
    INTEGER, DIMENSION(numPeriods,numGames) :: CycleStates
    INTEGER, DIMENSION(numAgents,numPeriods,numGames) :: CyclePrices
    INTEGER :: OptimalStrategyVec(lengthStrategies), OptimalStrategy(numStates,numAgents)
    INTEGER :: ShockLength, PostLength, PunishmentStrategy
    INTEGER, DIMENSION(numAgents,numThresPeriodsLength) :: FreqShockLength, FreqPostLength
    INTEGER :: FreqPunishmentStrategy(numAgents,0:numThresPeriodsLength)
    INTEGER :: ThresPeriodsLength(numThresPeriodsLength), ThresPeriodsLength0(numThresPeriodsLength0)
    !
    REAL(8) :: r_tmp, nn
    REAL(8), DIMENSION(numAgents,numPeriods,numGames) :: CycleProfits
    REAL(8), DIMENSION(numPeriods,numAgents) :: PrePrices, PreProfits
    REAL(8), DIMENSION(numAgents) :: AvgPrePrices, AvgPreProfits, AvgPrePricesQ, AvgPreProfitsQ
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents) :: ShockPrices, ShockProfits, &
        AvgShockPricesTmp, AvgShockProfitsTmp, AvgShockPricesPercTmp, AvgShockProfitsPercTmp
    REAL(8), DIMENSION(numAgents) :: PostPrices, PostProfits, &
        AvgPostPricesTmp, AvgPostProfitsTmp, AvgPostPricesPercTmp, AvgPostProfitsPercTmp
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents,numAgents) :: &
        AvgShockPrices, AvgShockProfits, AvgShockPricesQ, AvgShockProfitsQ, &
        AvgShockPricesPerc, AvgShockProfitsPerc, AvgShockPricesPercQ, AvgShockProfitsPercQ
    REAL(8), DIMENSION(numAgents,numAgents) :: &
        AvgPostPrices, AvgPostProfits, AvgPostPricesQ, AvgPostProfitsQ, &
        AvgPostPricesPerc, AvgPostProfitsPerc, AvgPostPricesPercQ, AvgPostProfitsPercQ
    REAL(8) :: AggrPrePrices, AggrDevPostPrices, AggrNonDevPostPrices, AggrPreProfits, AggrDevPostProfits, AggrNonDevPostProfits, &
        AggrDevPostPricesPerc, AggrNonDevPostPricesPerc, AggrDevPostProfitsPerc, AggrNonDevPostProfitsPerc
    REAL(8), DIMENSION(numShockPeriodsPrint) :: &
        AggrDevShockPrices, AggrNonDevShockPrices, AggrDevShockProfits, AggrNonDevShockProfits, &
        AggrDevShockPricesPerc, AggrNonDevShockPricesPerc, AggrDevShockProfitsPerc, AggrNonDevShockProfitsPerc
    REAL(8) :: AggrPrePricesQ, AggrDevPostPricesQ, AggrNonDevPostPricesQ, AggrPreProfitsQ, AggrDevPostProfitsQ, AggrNonDevPostProfitsQ, &
        AggrDevPostPricesPercQ, AggrNonDevPostPricesPercQ, AggrDevPostProfitsPercQ, AggrNonDevPostProfitsPercQ
    REAL(8), DIMENSION(numShockPeriodsPrint) :: &
        AggrDevShockPricesQ, AggrNonDevShockPricesQ, AggrDevShockProfitsQ, AggrNonDevShockProfitsQ, &
        AggrDevShockPricesPercQ, AggrNonDevShockPricesPercQ, AggrDevShockProfitsPercQ, AggrNonDevShockProfitsPercQ
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Impulse Response functions to BR'
    !
    ! Reading strategies and states at convergence from file
    !
    OPEN(UNIT = 998,FILE = FileNameIndexStrategies,STATUS = "OLD")    ! Open indexStrategies txt file
    DO i = 1, lengthStrategies
        !
        IF (MOD(i,10000) .EQ. 0) PRINT*, 'Read ', i, ' lines of indexStrategies'
        READ(998,21) (indexStrategies(i,iGame), iGame = 1, numGames)
    21  FORMAT(<numGames>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    CLOSE(UNIT = 998)                   ! Close indexStrategies txt file
    !
    OPEN(UNIT = 999,FILE = FileNamePriceCycles,STATUS = "OLD")     ! Open priceCycles txt file
    DO iGame = 1, numGames
        !
        READ(999,22) CycleLength(iGame)
    22  FORMAT(I8)    
        !
    END DO
    PRINT*, 'Read Cycles Length'
    CLOSE(UNIT = 999)                   ! Close priceCycles txt file
    !
    CycleStates = 0
    CyclePrices = 0
    CycleProfits = 0.d0
    OPEN(UNIT = 999,FILE = FileNamePriceCycles,STATUS = "OLD")     ! Re-Open priceCycles txt file
    DO iGame = 1, numGames
        !
        READ(999,23) CycleStates(:CycleLength(iGame),iGame), &
            ((CyclePrices(iAgent,iCycle,iGame), iCycle = 1, CycleLength(iGame)), iAgent = 1, numAgents), &
            ((CycleProfits(iAgent,iCycle,iGame), iCycle = 1, CycleLength(iGame)), iAgent = 1, numAgents)
23      FORMAT(9X, <CycleLength(iGame)>(I<lengthStatesPrint>, 1X), &
            <numAgents>(<CycleLength(iGame)>(I<lengthFormatActionPrint>, 1X)), &
            <numAgents>(<CycleLength(iGame)>(F8.5, 1X)))
        !
    END DO
    PRINT*, 'Read Cycles States'
    CLOSE(UNIT = 999)                   ! Re-Close priceCycles txt file
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !
    ThresPeriodsLength = (/ ( i, i = 1, numThresPeriodsLength ) /)
    ThresPeriodsLength0 = (/ 0, ThresPeriodsLength /)
    !
    FreqPreLength = 0
    AvgPrePrices = 0.d0
    AvgPreProfits = 0.d0
    AvgPrePricesQ = 0.d0
    AvgPreProfitsQ = 0.d0
    !
    FreqShockLength = 0
    FreqPunishmentStrategy = 0
    AvgShockPrices = 0.d0
    AvgShockPricesQ = 0.d0
    AvgShockProfits = 0.d0
    AvgShockProfitsQ = 0.d0
    AvgShockPricesPerc = 0.d0
    AvgShockPricesPercQ = 0.d0
    AvgShockProfitsPerc = 0.d0
    AvgShockProfitsPercQ = 0.d0
    !
    FreqPostLength = 0
    AvgPostPrices = 0.d0
    AvgPostPricesQ = 0.d0
    AvgPostProfits = 0.d0
    AvgPostProfitsQ = 0.d0
    AvgPostPricesPerc = 0.d0
    AvgPostPricesPercQ = 0.d0
    AvgPostProfitsPerc = 0.d0
    AvgPostProfitsPercQ = 0.d0
    !
    ! Beginning loop over games
    !
!    !$omp parallel do &
!    !$omp private(OptimalStrategy,LastObservedPrices,VisitedStatesPre,visitedPrices, &
!    !$omp   VisitedProfits,p,pPrime,iPeriod,iAgent,OptimalStrategyVec,LastStateVec, &
!    !$omp   VisitedStates,flagReturnedToState,jAgent,indexShockState,PrePrices,PreProfits, &
!    !$omp   PeriodsLengthPre,iStatePre,PeriodsLengthShock,PunishmentStrategy,PeriodsLengthPost, &
!    !$omp   AvgShockPricesTmp,AvgShockProfitsTmp,AvgShockPricesPercTmp,AvgShockProfitsPercTmp, &
!    !$omp   numPeriodsShockTmp,nn) &
!    !$omp firstprivate(PI,PricesGrids) &
!    !$omp reduction(+ : FreqPreLength, &
!    !$omp   AvgPrePrices,AvgPreProfits,AvgShockPrices,AvgShockProfits,AvgPostPrices,AvgPostProfits, &
!    !$omp   AvgPrePricesQ,AvgPreProfitsQ,AvgShockPricesQ,AvgShockProfitsQ,AvgPostPricesQ,AvgPostProfitsQ, &
!    !$omp   AvgShockPricesPerc,AvgShockProfitsPerc,AvgPostPricesPerc,AvgPostProfitsPerc, &
!    !$omp   AvgShockPricesPercQ,AvgShockProfitsPercQ,AvgPostPricesPercQ,AvgPostProfitsPercQ, &
!    !$omp   FreqPeriodLengthShock,FreqPunishmentStrategy,FreqPeriodLengthPost)
    DO iGame = 1, numGames        ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        !$omp critical
        OptimalStrategyVec = indexStrategies(:,iGame)
        !$omp end critical
        !
        OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Pre-shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        PeriodsLengthPre = CycleLength(iGame)
        PosThres = MIN(numThresPeriodsLength,PeriodsLengthPre)
        FreqPreLength(PosThres) = FreqPreLength(PosThres)+1
        VisitedStatesPre = 0
        PrePrices = 0.d0
        PreProfits = 0.d0
        DO iPeriod = 1, PeriodsLengthPre
            !
            VisitedStatesPre(iPeriod) = CycleStates(iPeriod,iGame)
            DO iAgent = 1, numAgents
                !
                PrePrices(iPeriod,iAgent) = PricesGrids(CyclePrices(iAgent,iPeriod,iGame),iAgent)
                PreProfits(iPeriod,iAgent) = CycleProfits(iAgent,iPeriod,iGame)
                !
            END DO
            !
        END DO
        !
        AvgPrePrices = AvgPrePrices+SUM(PrePrices(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
        AvgPrePricesQ = AvgPrePricesQ+(SUM(PrePrices(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre))**2
        AvgPreProfits = AvgPreProfits+SUM(PreProfits(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
        AvgPreProfitsQ = AvgPreProfitsQ+(SUM(PreProfits(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre))**2
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iAgent = 1, numAgents        ! Start of loop aver shocking agent 
            !
            AvgShockPricesTmp = 0.d0
            AvgShockProfitsTmp = 0.d0
            AvgShockPricesPercTmp = 0.d0
            AvgShockProfitsPercTmp = 0.d0
            AvgPostPricesTmp = 0.d0
            AvgPostProfitsTmp = 0.d0
            AvgPostPricesPercTmp = 0.d0
            AvgPostProfitsPercTmp = 0.d0
            !
            DO iStatePre = 1, PeriodsLengthPre      ! Start of loop over pre-shock cycle states
                !
                ! Agent "iAgent" selects the static best deviation price,
                ! The other agents stick to the strategy at convergence
                !
                CALL ComputeStaticBestResponse(OptimalStrategy,iStatePre,iAgent,PI(:,iAgent),DevPrice,r_tmp)
                CALL computeIndividualIR(OptimalStrategy,VisitedStatesPre(iStatePre),iAgent,DevPrice,numShockPeriodsPrint, &
                    PeriodsLengthPre,VisitedStatesPre(:PeriodsLengthPre), &
                    ShockPrices,ShockProfits,PostPrices,PostProfits, &
                    ShockLength,PunishmentStrategy,PostLength)
                !
                nn = DBLE(iStatePre)
                AvgShockPricesTmp = (nn-1.d0)/nn*AvgShockPricesTmp+ShockPrices/nn
                AvgShockPricesPercTmp = (nn-1.d0)/nn*AvgShockPricesPercTmp+ &
                    ShockPrices/SPREAD(PrePrices(iStatePre,:),DIM = 1,NCOPIES = numShockPeriodsPrint)/nn
                AvgShockProfitsTmp = (nn-1.d0)/nn*AvgShockProfitsTmp+ShockProfits/nn
                AvgShockProfitsPercTmp = (nn-1.d0)/nn*AvgShockProfitsPercTmp+ &
                    ShockProfits/SPREAD(PreProfits(iStatePre,:),DIM = 1,NCOPIES = numShockPeriodsPrint)/nn
                !
                AvgPostPricesTmp = (nn-1.d0)/nn*AvgPostPricesTmp+PostPrices/nn
                AvgPostPricesPercTmp = (nn-1.d0)/nn*AvgPostPricesPercTmp+PostPrices/PrePrices(iStatePre,:)/nn
                AvgPostProfitsTmp = (nn-1.d0)/nn*AvgPostProfitsTmp+PostProfits/nn
                AvgPostProfitsPercTmp = (nn-1.d0)/nn*AvgPostProfitsPercTmp+PostProfits/PreProfits(iStatePre,:)/nn
                !
                FreqShockLength(iAgent,MIN(numThresPeriodsLength,ShockLength)) = &
                    FreqShockLength(iAgent,MIN(numThresPeriodsLength,ShockLength))+1
                FreqPunishmentStrategy(iAgent,MIN(numThresPeriodsLength,PunishmentStrategy)) = &
                    FreqPunishmentStrategy(iAgent,MIN(numThresPeriodsLength,PunishmentStrategy))+1
                FreqPostLength(iAgent,MIN(numThresPeriodsLength,PostLength)) = &
                    FreqPostLength(iAgent,MIN(numThresPeriodsLength,PostLength))+1
                !
            END DO                          ! End of loop over pre-shock cycle states
            !
            ! Compute average prices and profits over pre-shock cycle states
            !
            AvgShockPrices(:,iAgent,:) = AvgShockPrices(:,iAgent,:)+AvgShockPricesTmp
            AvgShockPricesQ(:,iAgent,:) = AvgShockPricesQ(:,iAgent,:)+AvgShockPricesTmp**2
            AvgShockProfits(:,iAgent,:) = AvgShockProfits(:,iAgent,:)+AvgShockProfitsTmp
            AvgShockProfitsQ(:,iAgent,:) = AvgShockProfitsQ(:,iAgent,:)+AvgShockProfitsTmp**2
            AvgShockPricesPerc(:,iAgent,:) = AvgShockPricesPerc(:,iAgent,:)+AvgShockPricesPercTmp
            AvgShockPricesPercQ(:,iAgent,:) = AvgShockPricesPercQ(:,iAgent,:)+AvgShockPricesPercTmp**2
            AvgShockProfitsPerc(:,iAgent,:) = AvgShockProfitsPerc(:,iAgent,:)+AvgShockProfitsPercTmp
            AvgShockProfitsPercQ(:,iAgent,:) = AvgShockProfitsPercQ(:,iAgent,:)+AvgShockProfitsPercTmp**2
            AvgPostPrices(iAgent,:) = AvgPostPrices(iAgent,:)+AvgPostPricesTmp
            AvgPostPricesQ(iAgent,:) = AvgPostPricesQ(iAgent,:)+AvgPostPricesTmp**2
            AvgPostProfits(iAgent,:) = AvgPostProfits(iAgent,:)+AvgPostProfitsTmp
            AvgPostProfitsQ(iAgent,:) = AvgPostProfitsQ(iAgent,:)+AvgPostProfitsTmp**2
            AvgPostPricesPerc(iAgent,:) = AvgPostPricesPerc(iAgent,:)+AvgPostPricesPercTmp
            AvgPostPricesPercQ(iAgent,:) = AvgPostPricesPercQ(iAgent,:)+AvgPostPricesPercTmp**2
            AvgPostProfitsPerc(iAgent,:) = AvgPostProfitsPerc(iAgent,:)+AvgPostProfitsPercTmp
            AvgPostProfitsPercQ(iAgent,:) = AvgPostProfitsPercQ(iAgent,:)+AvgPostProfitsPercTmp**2
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
    AvgPrePrices = AvgPrePrices/DBLE(numGames)
    AvgPreProfits = AvgPreProfits/DBLE(numGames)
    AvgShockPrices = AvgShockPrices/DBLE(numGames)
    AvgShockProfits = AvgShockProfits/DBLE(numGames)
    AvgPostPrices = AvgPostPrices/DBLE(numGames)
    AvgPostProfits = AvgPostProfits/DBLE(numGames)
    AvgShockPricesPerc = AvgShockPricesPerc/DBLE(numGames)
    AvgShockProfitsPerc = AvgShockProfitsPerc/DBLE(numGames)
    AvgPostPricesPerc = AvgPostPricesPerc/DBLE(numGames)
    AvgPostProfitsPerc = AvgPostProfitsPerc/DBLE(numGames)
    AvgPrePricesQ = AvgPrePricesQ/DBLE(numGames)
    AvgPreProfitsQ = AvgPreProfitsQ/DBLE(numGames)
    AvgShockPricesQ = AvgShockPricesQ/DBLE(numGames)
    AvgShockProfitsQ = AvgShockProfitsQ/DBLE(numGames)
    AvgPostPricesQ = AvgPostPricesQ/DBLE(numGames)
    AvgPostProfitsQ = AvgPostProfitsQ/DBLE(numGames)
    AvgShockPricesPercQ = AvgShockPricesPercQ/DBLE(numGames)
    AvgShockProfitsPercQ = AvgShockProfitsPercQ/DBLE(numGames)
    AvgPostPricesPercQ = AvgPostPricesPercQ/DBLE(numGames)
    AvgPostProfitsPercQ = AvgPostProfitsPercQ/DBLE(numGames)
    !
    ! Computing aggregate (deviating and non-deviating) averages of prices and profits
    !
    AggrPrePrices = SUM(AvgPrePrices)/DBLE(numAgents)
    AggrPreProfits = SUM(AvgPreProfits)/DBLE(numAgents)
    AggrPrePricesQ = SUM(AvgPrePricesQ)/DBLE(numAgents)
    AggrPreProfitsQ = SUM(AvgPreProfitsQ)/DBLE(numAgents)
    !
    AggrDevShockPrices = 0.d0
    AggrDevShockProfits = 0.d0
    AggrDevShockPricesQ = 0.d0
    AggrDevShockProfitsQ = 0.d0
    AggrDevShockPricesPerc = 0.d0
    AggrDevShockProfitsPerc = 0.d0
    AggrDevShockPricesPercQ = 0.d0
    AggrDevShockProfitsPercQ = 0.d0
    DO iPeriod = 1, numShockPeriodsPrint
        !
        DO iAgent = 1, numAgents
            !
            AggrDevShockPrices(iPeriod) = AggrDevShockPrices(iPeriod)+AvgShockPrices(iPeriod,iAgent,iAgent)
            AggrDevShockProfits(iPeriod) = AggrDevShockProfits(iPeriod)+AvgShockProfits(iPeriod,iAgent,iAgent)
            AggrDevShockPricesQ(iPeriod) = AggrDevShockPricesQ(iPeriod)+AvgShockPricesQ(iPeriod,iAgent,iAgent)
            AggrDevShockProfitsQ(iPeriod) = AggrDevShockProfitsQ(iPeriod)+AvgShockProfitsQ(iPeriod,iAgent,iAgent)
            !
            AggrDevShockPricesPerc(iPeriod) = AggrDevShockPricesPerc(iPeriod)+AvgShockPricesPerc(iPeriod,iAgent,iAgent)
            AggrDevShockProfitsPerc(iPeriod) = AggrDevShockProfitsPerc(iPeriod)+AvgShockProfitsPerc(iPeriod,iAgent,iAgent)
            AggrDevShockPricesPercQ(iPeriod) = AggrDevShockPricesPercQ(iPeriod)+AvgShockPricesPercQ(iPeriod,iAgent,iAgent)
            AggrDevShockProfitsPercQ(iPeriod) = AggrDevShockProfitsPercQ(iPeriod)+AvgShockProfitsPercQ(iPeriod,iAgent,iAgent)
            !
        END DO
        AggrNonDevShockPrices(iPeriod) = (SUM(AvgShockPrices(iPeriod,:,:))-AggrDevShockPrices(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockPrices(iPeriod) = AggrDevShockPrices(iPeriod)/DBLE(numAgents)
        AggrNonDevShockProfits(iPeriod) = (SUM(AvgShockProfits(iPeriod,:,:))-AggrDevShockProfits(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockProfits(iPeriod) = AggrDevShockProfits(iPeriod)/DBLE(numAgents)
        AggrNonDevShockPricesQ(iPeriod) = (SUM(AvgShockPricesQ(iPeriod,:,:))-AggrDevShockPricesQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockPricesQ(iPeriod) = AggrDevShockPricesQ(iPeriod)/DBLE(numAgents)
        AggrNonDevShockProfitsQ(iPeriod) = (SUM(AvgShockProfitsQ(iPeriod,:,:))-AggrDevShockProfitsQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockProfitsQ(iPeriod) = AggrDevShockProfitsQ(iPeriod)/DBLE(numAgents)
        !
        AggrNonDevShockPricesPerc(iPeriod) = (SUM(AvgShockPricesPerc(iPeriod,:,:))-AggrDevShockPricesPerc(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockPricesPerc(iPeriod) = AggrDevShockPricesPerc(iPeriod)/DBLE(numAgents)
        AggrNonDevShockProfitsPerc(iPeriod) = (SUM(AvgShockProfitsPerc(iPeriod,:,:))-AggrDevShockProfitsPerc(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockProfitsPerc(iPeriod) = AggrDevShockProfitsPerc(iPeriod)/DBLE(numAgents)
        AggrNonDevShockPricesPercQ(iPeriod) = (SUM(AvgShockPricesPercQ(iPeriod,:,:))-AggrDevShockPricesPercQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockPricesPercQ(iPeriod) = AggrDevShockPricesPercQ(iPeriod)/DBLE(numAgents)
        AggrNonDevShockProfitsPercQ(iPeriod) = (SUM(AvgShockProfitsPercQ(iPeriod,:,:))-AggrDevShockProfitsPercQ(iPeriod))/DBLE(numAgents*(numAgents-1))
        AggrDevShockProfitsPercQ(iPeriod) = AggrDevShockProfitsPercQ(iPeriod)/DBLE(numAgents)
        !
    END DO
    !
    AggrDevPostPrices = 0.d0
    AggrDevPostProfits = 0.d0
    AggrDevPostPricesQ = 0.d0
    AggrDevPostProfitsQ = 0.d0
    AggrDevPostPricesPerc = 0.d0
    AggrDevPostProfitsPerc = 0.d0
    AggrDevPostPricesPercQ = 0.d0
    AggrDevPostProfitsPercQ = 0.d0
    DO iAgent = 1, numAgents
        !
        AggrDevPostPrices = AggrDevPostPrices+AvgPostPrices(iAgent,iAgent)
        AggrDevPostProfits = AggrDevPostProfits+AvgPostProfits(iAgent,iAgent)
        AggrDevPostPricesQ = AggrDevPostPricesQ+AvgPostPricesQ(iAgent,iAgent)
        AggrDevPostProfitsQ = AggrDevPostProfitsQ+AvgPostProfitsQ(iAgent,iAgent)
        !
        AggrDevPostPricesPerc = AggrDevPostPricesPerc+AvgPostPricesPerc(iAgent,iAgent)
        AggrDevPostProfitsPerc = AggrDevPostProfitsPerc+AvgPostProfitsPerc(iAgent,iAgent)
        AggrDevPostPricesPercQ = AggrDevPostPricesPercQ+AvgPostPricesPercQ(iAgent,iAgent)
        AggrDevPostProfitsPercQ = AggrDevPostProfitsPercQ+AvgPostProfitsPercQ(iAgent,iAgent)
        !
    END DO
    AggrNonDevPostPrices = (SUM(AvgPostPrices(:,:))-AggrDevPostPrices)/DBLE(numAgents*(numAgents-1))
    AggrDevPostPrices = AggrDevPostPrices/DBLE(numAgents)
    AggrNonDevPostProfits = (SUM(AvgPostProfits(:,:))-AggrDevPostProfits)/DBLE(numAgents*(numAgents-1))
    AggrDevPostProfits = AggrDevPostProfits/DBLE(numAgents)
    AggrNonDevPostPricesQ = (SUM(AvgPostPricesQ(:,:))-AggrDevPostPricesQ)/DBLE(numAgents*(numAgents-1))
    AggrDevPostPricesQ = AggrDevPostPricesQ/DBLE(numAgents)
    AggrNonDevPostProfitsQ = (SUM(AvgPostProfitsQ(:,:))-AggrDevPostProfitsQ)/DBLE(numAgents*(numAgents-1))
    AggrDevPostProfitsQ = AggrDevPostProfitsQ/DBLE(numAgents)
    !
    AggrNonDevPostPricesPerc = (SUM(AvgPostPricesPerc(:,:))-AggrDevPostPricesPerc)/DBLE(numAgents*(numAgents-1))
    AggrDevPostPricesPerc = AggrDevPostPricesPerc/DBLE(numAgents)
    AggrNonDevPostProfitsPerc = (SUM(AvgPostProfitsPerc(:,:))-AggrDevPostProfitsPerc)/DBLE(numAgents*(numAgents-1))
    AggrDevPostProfitsPerc = AggrDevPostProfitsPerc/DBLE(numAgents)
    AggrNonDevPostPricesPercQ = (SUM(AvgPostPricesPercQ(:,:))-AggrDevPostPricesPercQ)/DBLE(numAgents*(numAgents-1))
    AggrDevPostPricesPercQ = AggrDevPostPricesPercQ/DBLE(numAgents)
    AggrNonDevPostProfitsPercQ = (SUM(AvgPostProfitsPercQ(:,:))-AggrDevPostProfitsPercQ)/DBLE(numAgents*(numAgents-1))
    AggrDevPostProfitsPercQ = AggrDevPostProfitsPercQ/DBLE(numAgents)
    !
    ! Computing standard errors
    !
    AvgPrePricesQ = SQRT(ABS((AvgPrePricesQ-AvgPrePrices**2)/DBLE(numGames)))
    AvgPreProfitsQ = SQRT(ABS((AvgPreProfitsQ-AvgPreProfits**2)/DBLE(numGames)))
    AvgShockPricesQ = SQRT(ABS((AvgShockPricesQ-AvgShockPrices**2)/DBLE(numGames)))
    AvgShockProfitsQ = SQRT(ABS((AvgShockProfitsQ-AvgShockProfits**2)/DBLE(numGames)))
    AvgPostPricesQ = SQRT(ABS((AvgPostPricesQ-AvgPostPrices**2)/DBLE(numGames)))
    AvgPostProfitsQ = SQRT(ABS((AvgPostProfitsQ-AvgPostProfits**2)/DBLE(numGames)))
    !
    AggrPrePricesQ = SQRT(ABS((AggrPrePricesQ-AggrPrePrices**2)/DBLE(numGames)))
    AggrPreProfitsQ = SQRT(ABS((AggrPreProfitsQ-AggrPreProfits**2)/DBLE(numGames)))
    DO iPeriod = 1, numShockPeriodsPrint
        !
        AggrNonDevShockPricesQ(iPeriod) = SQRT(ABS((AggrNonDevShockPricesQ(iPeriod)-AggrNonDevShockPrices(iPeriod)**2)/DBLE(numGames)))
        AggrDevShockPricesQ(iPeriod) = SQRT(ABS((AggrDevShockPricesQ(iPeriod)-AggrDevShockPrices(iPeriod)**2)/DBLE(numGames)))
        AggrNonDevShockProfitsQ(iPeriod) = SQRT(ABS((AggrNonDevShockProfitsQ(iPeriod)-AggrNonDevShockProfits(iPeriod)**2)/DBLE(numGames)))
        AggrDevShockProfitsQ(iPeriod) = SQRT(ABS((AggrDevShockProfitsQ(iPeriod)-AggrDevShockProfits(iPeriod)**2)/DBLE(numGames)))
        !
        AggrNonDevShockPricesPercQ(iPeriod) = SQRT(ABS((AggrNonDevShockPricesPercQ(iPeriod)-AggrNonDevShockPricesPerc(iPeriod)**2)/DBLE(numGames)))
        AggrDevShockPricesPercQ(iPeriod) = SQRT(ABS((AggrDevShockPricesPercQ(iPeriod)-AggrDevShockPricesPerc(iPeriod)**2)/DBLE(numGames)))
        AggrNonDevShockProfitsPercQ(iPeriod) = SQRT(ABS((AggrNonDevShockProfitsPercQ(iPeriod)-AggrNonDevShockProfitsPerc(iPeriod)**2)/DBLE(numGames)))
        AggrDevShockProfitsPercQ(iPeriod) = SQRT(ABS((AggrDevShockProfitsPercQ(iPeriod)-AggrDevShockProfitsPerc(iPeriod)**2)/DBLE(numGames)))
        !
    END DO
    AggrNonDevPostPricesQ = SQRT(ABS((AggrNonDevPostPricesQ-AggrNonDevPostPrices**2)/DBLE(numGames)))
    AggrDevPostPricesQ = SQRT(ABS((AggrDevPostPricesQ-AggrDevPostPrices**2)/DBLE(numGames)))
    AggrNonDevPostProfitsQ = SQRT(ABS((AggrNonDevPostProfitsQ-AggrNonDevPostProfits**2)/DBLE(numGames)))
    AggrDevPostProfitsQ = SQRT(ABS((AggrDevPostProfitsQ-AggrDevPostProfits**2)/DBLE(numGames)))
    !
    AggrNonDevPostPricesPercQ = SQRT(ABS((AggrNonDevPostPricesPercQ-AggrNonDevPostPricesPerc**2)/DBLE(numGames)))
    AggrDevPostPricesPercQ = SQRT(ABS((AggrDevPostPricesPercQ-AggrDevPostPricesPerc**2)/DBLE(numGames)))
    AggrNonDevPostProfitsPercQ = SQRT(ABS((AggrNonDevPostProfitsPercQ-AggrNonDevPostProfitsPerc**2)/DBLE(numGames)))
    AggrDevPostProfitsPercQ = SQRT(ABS((AggrDevPostProfitsPercQ-AggrDevPostProfitsPerc**2)/DBLE(numGames)))   
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10003,1) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
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
                jAgent = 1, numAgents), iAgent = 1, numAgents), &  ! AvgPrices
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), &  ! sePrices
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents), &  ! AvgProfits
            ((jAgent, (jAgent, iAgent, iPeriod, iPeriod = 1, numShockPeriodsPrint), jAgent, iAgent, &
                jAgent = 1, numAgents), iAgent = 1, numAgents)     ! seProfits        
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), &
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
            <numAgents>(<numAgents>('Ag', I1, 'AvgPricePre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'AvgPriceShockAg', I1, 'Per', I3.3, ' '), &
                'Ag', I1, 'AvgPricePostAg', I1, ' ')), &
            <numAgents>(<numAgents>('seAg', I1, 'AvgPricePre', ' ', &
                <numShockPeriodsPrint>('seAg', I1, 'AvgPriceShockAg', I1, 'Per', I3.3, ' '), &
                'seAg', I1, 'AvgPricePostAg', I1, ' ')), &
            <numAgents>(<numAgents>('Ag', I1, 'AvgProfitPre', ' ', &
                <numShockPeriodsPrint>('Ag', I1, 'AvgProfitShockAg', I1, 'Per', I3.3, ' '), &
                'Ag', I1, 'AvgProfitPostAg', I1, ' ')), &
            <numAgents>(<numAgents>('seAg', I1, 'AvgProfitPre', ' ', &
                <numShockPeriodsPrint>('seAg', I1, 'AvgProfitShockAg', I1, 'Per', I3.3, ' '), &
                'seAg', I1, 'AvgProfitPostAg', I1, ' ')) &
            )
        !
    END IF
    !
    WRITE(10003,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        AggrPrePrices, AggrDevShockPrices, AggrDevPostPrices, AggrNonDevShockPrices, AggrNonDevPostPrices, &
        AggrPrePricesQ, AggrDevShockPricesQ, AggrDevPostPricesQ, AggrNonDevShockPricesQ, AggrNonDevPostPricesQ, &
        AggrDevShockPricesPerc, AggrDevPostPricesPerc, AggrNonDevShockPricesPerc, AggrNonDevPostPricesPerc, &
        AggrDevShockPricesPercQ, AggrDevPostPricesPercQ, AggrNonDevShockPricesPercQ, AggrNonDevPostPricesPercQ, &
        AggrPreProfits, AggrDevShockProfits, AggrDevPostProfits, AggrNonDevShockProfits, AggrNonDevPostProfits, &
        AggrPreProfitsQ, AggrDevShockProfitsQ, AggrDevPostProfitsQ, AggrNonDevShockProfitsQ, AggrNonDevPostProfitsQ, &
        AggrDevShockProfitsPerc, AggrDevPostProfitsPerc, AggrNonDevShockProfitsPerc, AggrNonDevPostProfitsPerc, &
        AggrDevShockProfitsPercQ, AggrDevPostProfitsPercQ, AggrNonDevShockProfitsPercQ, AggrNonDevPostProfitsPercQ, &
        FreqPreLength, &
        ((FreqPostLength(iAgent,i), i = 1, numThresPeriodsLength), iAgent = 1, numAgents), &
        ((FreqShockLength(iAgent,i), i = 1, numThresPeriodsLength), iAgent = 1, numAgents), &
        ((FreqPunishmentStrategy(iAgent,i), i = 0, numThresPeriodsLength), iAgent = 1, numAgents), &
        ((AvgPrePrices(jAgent), (AvgShockPrices(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), AvgPostPrices(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), &
        ((AvgPrePricesQ(jAgent), (AvgShockPricesQ(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), AvgPostPricesQ(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), &
        ((AvgPreProfits(jAgent), (AvgShockProfits(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), AvgPostProfits(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents), &
        ((AvgPreProfitsQ(jAgent), (AvgShockProfitsQ(iPeriod,iAgent,jAgent), iPeriod = 1, numShockPeriodsPrint), AvgPostProfitsQ(iAgent,jAgent), &
            jAgent = 1, numAgents), iAgent = 1, numAgents)
2   FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
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
    END SUBROUTINE computeIRAnalysisToBR
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ComputeStaticBestResponse ( OptimalStrategy, iState, iAgent, PIiAgent, IndexStaticBR, PIStaticBR )
    !
    ! Computes static best response of iAgent given all agents' strategies 
    ! 'Best' means that the selected price maximizes iAgent's profits assuming 
    ! that rivals play according to their strategies
    !
    IMPLICIT NONE
    !
    ! Declare dummy variables
    !
    INTEGER, INTENT(IN) :: OptimalStrategy(numStates,numAgents)
    INTEGER, INTENT(IN) :: iState, iAgent
    REAL(8), INTENT(IN) :: PIiAgent(numActions)
    INTEGER, INTENT(OUT) :: IndexStaticBR
    REAL(8), INTENT(OUT) :: PIStaticBR
    !
    ! Declare local variables
    !
    INTEGER :: iPrice
    INTEGER, DIMENSION(numAgents) :: pPrime
    REAL(8), DIMENSION(numPrices) :: selProfits
    !
    ! Beginning execution
    !
    pPrime = OptimalStrategy(iState,:)
    selProfits = 0.d0
    DO iPrice = 1, numPrices
        !
        pPrime(iAgent) = iPrice
        selProfits(iPrice) = PIiAgent(computeActionNumber(pPrime))
        !
    END DO
    IndexStaticBR = MINVAL(MAXLOC(selProfits))
    PIStaticBR = MAXVAL(selProfits)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ComputeStaticBestResponse    
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ComputeDynamicBestResponse ( OptimalStrategy, iState, iAgent, delta, IndexDynamicBR, QDynamicBR )
    !
    ! Computes dynamic best response of one agent given all agents' strategies 
    ! 'Best' means that the selected price maximizes Q given the state and assuming 
    ! that opponents play according to their strategies
    !
    IMPLICIT NONE
    !
    ! Declare dummy variables
    !
    INTEGER, INTENT(IN) :: OptimalStrategy(numStates,numAgents)
    INTEGER, INTENT(IN) :: iState
    INTEGER, INTENT(IN) :: iAgent
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    INTEGER, INTENT(OUT) :: IndexDynamicBR
    REAL(8), INTENT(OUT) :: QDynamicBR
    !
    ! Declare local variables
    !
    INTEGER :: iPrice, PreCycleLength, CycleLength, iPeriod
    INTEGER, DIMENSION(numPeriods) :: VisitedStates
    REAL(8), DIMENSION(numPrices) :: selQ
    !
    ! Beginning execution
    !
    selQ = 0.d0
    DO iPrice = 1, numPrices
        !
        CALL computeQcell(OptimalStrategy,iState,iPrice,iAgent,delta, &
            selQ(iPrice),VisitedStates,PreCycleLength,CycleLength,iPeriod)
        !
    END DO
    IndexDynamicBR = MINVAL(MAXLOC(selQ))
    QDynamicBR = MAXVAL(selQ)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ComputeDynamicBestResponse    
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeIndividualIR ( OptimalStrategy, InitialState, DevAgent, DevPrice, ShockObsLength, &
        PreCycleLength, PreCycleStates, &
        ShockPrices, ShockProfits, AvgPostPrices, AvgPostProfits, &
        ShockLength, PunishmentStrategy, PostLength )
    !
    ! INPUT:
    !
    ! - OptimalStrategy    : a strategy for all agents
    ! - InitialState       : an initial state
    ! - DevAgent           : the deviating agent index
    ! - DevPrice           : the deviation price index
    ! - ShockObsLength     : the length of the observation interval of the deviation period
    ! - PreCycleLength     : the length of the pre-deviation cycle
    ! - PreCycleStates     : the pre-deviation cycle states
    !
    ! OUTPUT:
    !
    ! - ShockPrices        : the trajectory of all agents' prices in the deviation interval
    ! - ShockProfits       : the trajectory of all agents' profits in the deviation interval
    ! - AvgPostPrices      : the average of all agents' prices in the post-deviation cycle
    ! - AvgPostProfits     : the average of all agents' profits in the post-deviation cycle
    ! - ShockLength        : the length of the non-cyclic deviation interval
    ! - PunishmentStrategy : indicator. After the deviation:
    !                        = 0: the system returns to a cycle different from the pre-deviation cycle
    !                        > 0: the system returns to the pre-deviation cycle after PunishmentStrategy periods
    ! - PostLength         : the length of the post-deviation cycle
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(numStates,numAgents), INTENT(IN) :: OptimalStrategy
    INTEGER, INTENT(IN) :: InitialState, DevAgent, DevPrice, ShockObsLength, PreCycleLength
    INTEGER, DIMENSION(PreCycleLength), INTENT(IN) :: PreCycleStates
    REAL(8), DIMENSION(ShockObsLength,numAgents), INTENT(OUT) :: ShockPrices, ShockProfits
    REAL(8), DIMENSION(numAgents), INTENT(OUT) :: AvgPostPrices, AvgPostProfits
    INTEGER, INTENT(OUT) :: ShockLength, PunishmentStrategy, PostLength
    !
    ! Declaring local variables
    !
    INTEGER :: iPeriod, jAgent
    INTEGER :: p(DepthState,numAgents), pPrime(numAgents) 
    INTEGER :: VisitedStates(MAX(ShockObsLength,numPeriods))
    INTEGER :: indexShockState(LengthStates)
    !
    REAL(8), DIMENSION(numPeriods,numAgents) :: visitedPrices, VisitedProfits
    !
    LOGICAL :: FlagReturnedToState
    !
    ! Beginning execution
    !
    p = RESHAPE(convertNumberBase(InitialState-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
    pPrime = OptimalStrategy(InitialState,:)
    !
    ! Agent "DevAgent" selects the best deviation price,
    ! the other agents stick to the strategy at convergence
    !
    pPrime(DevAgent) = DevPrice
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over deviation period
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    VisitedStates = 0
    flagReturnedToState = .FALSE.
    DO iPeriod = 1, MAX(ShockObsLength,numPeriods)
        !
        IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
        p(1,:) = pPrime
        VisitedStates(iPeriod) = computeStateNumber(p)
        DO jAgent = 1, numAgents
            !
            IF (iPeriod .LE. ShockObsLength) THEN
                !
                ShockPrices(iPeriod,jAgent) = PricesGrids(pPrime(jAgent),jAgent)
                ShockProfits(iPeriod,jAgent) = PI(computeActionNumber(pPrime),jAgent)
                !
            END IF
            !
        END DO
        !
        ! Check if the state has already been visited
        ! Case 1: the state retuns to one of the states in the pre-shock cycle
        !
        IF ((.NOT.(flagReturnedToState)) .AND. (ANY(PreCycleStates .EQ. VisitedStates(iPeriod)))) THEN
            !
            ShockLength = iPeriod
            PunishmentStrategy = iPeriod
            indexShockState = RESHAPE(p,(/ LengthStates /))
            flagReturnedToState = .TRUE.
            !
        END IF
        !
        ! Case 2: after some time, the state starts cycling among a new set of states
        !
        IF ((iPeriod .GE. 2) .AND. (.NOT.(flagReturnedToState)) .AND. &
            (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) THEN
            !
            ShockLength = MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
            PunishmentStrategy = 0
            indexShockState = RESHAPE(p,(/ LengthStates /))
            flagReturnedToState = .TRUE.
            !
        END IF
        pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
        !
    END DO
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Post-shock period 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    VisitedStates = 0
    VisitedPrices = 0.d0
    VisitedProfits = 0.d0
    p = RESHAPE(indexShockState, (/ DepthState,numAgents /) )
    pPrime = OptimalStrategy(computeStateNumber(p),:)
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
    PostLength = iPeriod-MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
    !
    AvgPostPrices = SUM(visitedPrices(iPeriod-PostLength+1:iPeriod,:),DIM = 1)/DBLE(PostLength)
    AvgPostProfits = SUM(visitedProfits(iPeriod-PostLength+1:iPeriod,:),DIM = 1)/DBLE(PostLength)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeIndividualIR    
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ImpulseResponseToBR
