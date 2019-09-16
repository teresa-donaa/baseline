MODULE globals
!
USE generic_routines
!
! Declare global parameters and variables 
!
IMPLICIT NONE
!
! Parameters
!
INTEGER, PARAMETER :: numShockPeriodsPrint = 25
INTEGER, PARAMETER :: numThresCycleLength = 10
!
REAL(8), PARAMETER :: SlackOnPath = 0.05d0
REAL(8), PARAMETER :: SlackOffPath = 0.05d0
REAL(8), PARAMETER :: twoOverPi = 0.636619772367581343075535053490d0    ! 2/Pi
REAL(8), PARAMETER :: piOverTwo = 1.57079632679489661923132169164d0     ! Pi/2
!
! Variables
!
INTEGER :: numModels, numCores, numGames, itersPerYear, maxNumYears, maxIters, &
    itersInPerfMeasPeriod, printQ, printP, codModel, PerfMeasPeriodTime, numPrices, lengthFormatActionPrint, &
    typeExplorationMechanism, DepthState0, DepthState, LengthStates, lengthStatesPrint, numStates, lengthStrategies, &
    typePayoffInput, numAgents, numActions, numDemandParameters, numPeriods, &
    numExplorationParameters, SwitchQLearningResults, SwitchConvergenceResults, &
    SwitchImpulseResponseToBR, SwitchImpulseResponseToNash, SwitchImpulseResponseToAll, &
    SwitchEquilibriumCheck, SwitchPIGapToMaximum, SwitchQGapToMaximum, SwitchDetailedAnalysis, SwitchRestart
REAL(8) :: PerfMeasPeriodLength, meanNashProfit, meanCoopProfit, gammaSinghVives
CHARACTER(len = 50) :: ModelNumber, FileNameIndexStrategies, FileNameIndexLastState, FileNamePriceCycles
!
INTEGER, ALLOCATABLE :: converged(:), indexActions(:,:), indexLastState(:,:), indexStrategies(:,:), &
    cStates(:), cActions(:), priceCycles(:,:), sampledIndexStrategies(:,:), sampledPriceCycles(:,:), &
    indexStates(:,:), indexEquivalentStates(:,:), indexNashPrices(:), indexCoopPrices(:), &
    SwitchMixedStrategies(:), QMatrixInitializationT(:), QMatrixInitializationF(:)
REAL(8), ALLOCATABLE :: timeToConvergence(:), NashProfits(:), CoopProfits(:), &
    maxValQ(:,:), NashPrices(:), CoopPrices(:), &
    PI(:,:), PIQ(:,:), avgPI(:), avgPIQ(:), alpha(:), delta(:), DiscountFactors(:,:), & 
    meanProfit(:), seProfit(:), meanProfitGain(:), seProfitGain(:), DemandParameters(:), &
    NashMarketShares(:), CoopMarketShares(:), PricesGrids(:,:), MExpl(:), ExplorationParameters(:), &
    QMatrixInitializationR(:,:), QMatrixInitializationU(:)
CHARACTER(len = :), ALLOCATABLE :: labelStates(:)
CHARACTER(len = :), ALLOCATABLE :: QFileFolderName(:)
CHARACTER(len = 1), ALLOCATABLE :: typeQInitialization(:)

!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readBatchVariables ( unitNumber )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unitNumber
    !
    ! Declaring local variables
    !
    INTEGER :: i, iAgent, iPrices, jPrices, iState, iAction, jState
    INTEGER, ALLOCATABLE :: switchedState(:), indA1(:), indA2(:)
    !
    ! Beginning execution
    !
    READ(unitNumber,'(1X)') 
    READ(unitNumber,*) numModels
    READ(unitNumber,'(1X)') 
    READ(unitNumber,*) numCores
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numGames
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) itersPerYear
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) maxNumYears
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) PerfMeasPeriodTime
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) PerfMeasPeriodLength
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numAgents
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) DepthState0
    DepthState = MAX(1,DepthState0)          ! Accomodates the DepthState = 0 case
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numPrices
    !
    ! Global variables
    !
    maxIters = maxNumYears*itersPerYear
    itersInPerfMeasPeriod = INT(PerfMeasPeriodLength*itersPerYear)
    LengthStates = MAX(1,numAgents*DepthState0)
    lengthStatesPrint = LengthStates*(1+FLOOR(LOG10(DBLE(numPrices))))+LengthStates-1
    numStates = numPrices**(numAgents*DepthState0)
    numPeriods = numStates+1
    numActions = numPrices**numAgents   ! Actions contain combinations of prices;
                                        ! they coincide with states when DepthState == 1
    lengthStrategies = numAgents*numStates
    lengthFormatActionPrint = FLOOR(LOG10(DBLE(numPrices)))+1
    ALLOCATE(SwitchMixedStrategies(numAgents))
    READ(unitNumber,'(1X)')
    !
    ! Read type of exploration mechanism
    !
    READ(unitNumber,*) typeExplorationMechanism
    IF (typeExplorationMechanism .EQ. 1) numExplorationParameters = 1*numAgents           ! Constant 
    IF (typeExplorationMechanism .EQ. 2) numExplorationParameters = 1*numAgents           ! Exponentially decreasing (beta)
    IF (typeExplorationMechanism .EQ. 3) numExplorationParameters = 1*numAgents           ! Exponentially decreasing (m)
    READ(unitNumber,'(1X)')
    !
    ! Read type of payoff input
    !
    READ(unitNumber,*) typePayoffInput
    IF (typePayoffInput .EQ. 0) numDemandParameters = numActions                ! Pi1 matrix
    IF (typePayoffInput .EQ. 1) numDemandParameters = 1+2                       ! gamma, extend
    IF (typePayoffInput .EQ. 2) numDemandParameters = 2*numAgents+4             ! a0, ai, ci, sigma, extend
    IF (typePayoffInput .EQ. 3) numDemandParameters = 2*numAgents+4             ! a0, ai, ci, sigma = 0, extend
    READ(unitNumber,'(1X)')
    !
    ! Read type of Q matrix initialization
    !
    ALLOCATE(typeQInitialization(numAgents),QMatrixInitializationR(numAgents,2), &
        QMatrixInitializationU(numAgents),QMatrixInitializationT(numAgents), &
        QMatrixInitializationF(numAgents))
    READ(unitNumber,*) typeQInitialization
    DO iAgent = 1, numAgents
        !
        IF (typeQInitialization(iAgent) .EQ. 'R') READ(unitNumber,*) QMatrixInitializationR(iAgent,:)
        IF (typeQInitialization(iAgent) .EQ. 'U') READ(unitNumber,*) QMatrixInitializationU(iAgent)
        IF (typeQInitialization(iAgent) .EQ. 'T') READ(unitNumber,*) QMatrixInitializationT(iAgent)
        IF (typeQInitialization(iAgent) .EQ. 'F') READ(unitNumber,*) QMatrixInitializationF(iAgent)
        !
    END DO
    READ(unitNumber,'(1X)')
    !
    ! Continue reading input settings
    !
    READ(unitNumber,*) SwitchQLearningResults
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchConvergenceResults
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchImpulseResponseToBR
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchImpulseResponseToNash
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchImpulseResponseToAll
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchEquilibriumCheck
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchQGapToMaximum
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchPIGapToMaximum
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchDetailedAnalysis
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchRestart
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchMixedStrategies
    READ(unitNumber,'(1X)')
    !
    ! Allocating matrices and vectors
    !
    ALLOCATE(indexActions(numActions,numAgents), indexStates(numStates,LengthStates), &
        indexEquivalentStates(numStates,numAgents), timeToConvergence(numGames), &
        converged(numGames),cStates(LengthStates),cActions(numAgents),DiscountFactors(0:numStates,numAgents), &
        maxValQ(numStates,numAgents), DemandParameters(numDemandParameters), &
        ExplorationParameters(numExplorationParameters), MExpl(numExplorationParameters), &
        alpha(numAgents),delta(numAgents),NashProfits(numAgents),CoopProfits(numAgents), &
        meanProfit(numAgents),seProfit(numAgents),meanProfitGain(numAgents),seProfitGain(numAgents), &
        PI(numActions,numAgents),PIQ(numActions,numAgents),avgPI(numActions),avgPIQ(numActions), &
        indexNashPrices(numAgents),indexCoopPrices(numAgents), &
        NashPrices(numAgents),CoopPrices(numAgents), &
        NashMarketShares(numAgents),CoopMarketShares(numAgents),PricesGrids(numPrices,numAgents))
    ALLOCATE(CHARACTER(len = 3+lengthStatesPrint) :: labelStates(numStates))
    ALLOCATE(CHARACTER(len = 200) :: QFileFolderName(numAgents))
    !
    cStates = (/ (numPrices**i, i = LengthStates-1, 0, -1) /)
    cActions = (/ (numPrices**i, i = numAgents-1, 0, -1) /)
    !
    ! Actions contain the most recent prices of all agents. Actions and States
    ! coincide when DepthState == 1
    !
    DO iAction = 1, numActions
        !
        indexActions(iAction,:) = convertNumberBase(iAction-1,numPrices,numAgents)
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readBatchVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE closeBatch ( )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    DEALLOCATE(indexActions,timeToConvergence,converged,cStates,cActions,maxValQ, &
        labelStates,indexStates,NashProfits,CoopProfits,QFileFolderName,DiscountFactors, &
        alpha,MExpl,ExplorationParameters,delta,indexEquivalentStates, &
        meanProfit,seProfit,meanProfitGain,seProfitGain,DemandParameters,PI,PIQ,avgPI,avgPIQ, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids, &
        SwitchMixedStrategies,typeQInitialization,QMatrixInitializationR,QMatrixInitializationU, &
        QMatrixInitializationT,QMatrixInitializationF)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE closeBatch
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readModelVariables ( unitNumber )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unitNumber
    !
    ! Declaring local variables
    !
    INTEGER :: i
    !
    ! Beginning execution
    !
    READ(unitNumber,*) codModel, printQ, printP, alpha, MExpl, delta, &
        DemandParameters, NashPrices, CoopPrices
    IF (typeExplorationMechanism .EQ. 2) THEN
        !
        ExplorationParameters = MExpl
        !
    ELSE IF (typeExplorationMechanism .EQ. 3) THEN
        !
        ExplorationParameters = -DBLE(itersPerYear)/DBLE(numAgents+1)* &
            LOG(1.d0-(DBLE(numPrices-1)/DBLE(numPrices))**numAgents/(DBLE(numStates*numPrices)*MExpl))
        !
    END IF
    DiscountFactors = TRANSPOSE(RESHAPE((/ (delta**i, i = 0, numPeriods-1) /),(/ numAgents,numPeriods /)))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readModelVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE globals