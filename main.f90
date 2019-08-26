PROGRAM main
!
USE globals
USE LearningSimulation
USE LearningSimulationRestart
USE ConvergenceResults
USE ImpulseResponse
USE DetailedImpulseResponseToAll
USE EquilibriumCheck
USE DetailedEquilibriumCheck
USE QGapToMaximum
USE PIGapToMaximum
USE MixedStrategiesResults
USE QL_routines
USE PI_routines
USE generic_routines
!
IMPLICIT NONE
!
! Declaring variables and parameters
!
INTEGER :: iModel, iter, i, j, h, iGame, numGamesConverged, iAgent, errcode, nlines
REAL(8) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
REAL(8) :: meanAvgProfit, seAvgProfit
REAL(8) :: meanAvgProfitGain, seAvgProfitGain
REAL(8) :: herfFreq, entropyFreq, giniFreq, freqSymmetricStrategies
CHARACTER(len = 5) :: iChar
CHARACTER(len = 50) :: ModelName, FileName, FileNameMSR
REAL(8), ALLOCATABLE :: alpha_tmp(:), beta_tmp(:), delta_tmp(:)
!
! Beginning execution
!
! Opening files
!
ModelName = "detailed_IR_analysis.txt"
FileName = "mod_" // ModelName
!
OPEN(UNIT = 10001,FILE = FileName)
CALL readBatchVariables(10001)
DO iAgent = 1, numAgents
    !
    IF (typeQInitialization(iAgent) .EQ. 'T') THEN
        !
        WRITE(iChar,'(I0.5)') QMatrixInitializationT(iAgent)
        QFileFolderName(iAgent) = &
            'C:/Users/sergio.pastorello/Documents/jobs/dynamic pricing/qlearning/baseline/mixed_Q_analysis/Q_' // &
                iChar // '/'
        !
    END IF
    !
END DO
!
IF (computeQLearningResults .EQ. 1) THEN
    !
    FileName = "res_" // ModelName
    OPEN(UNIT = 10002,FILE = FileName)
    !
END IF
IF (computeConvergenceResults .EQ. 1) THEN
    !
    FileName = "ConvResults_" // ModelName
    OPEN(UNIT = 100022,FILE = FileName)
    !
END IF
IF (computeImpulseResponseToBR .EQ. 1) THEN
    !
    FileName = "irToBR_" // ModelName
    OPEN(UNIT = 10003,FILE = FileName)
    !
END IF
IF (computeImpulseResponseToNash .GE. 1) THEN
    !
    FileName = "irToNash_" // ModelName
    OPEN(UNIT = 100031,FILE = FileName)
    !
END IF
IF (computeImpulseResponseToAll .EQ. 1) THEn
    !
    FileName = "irToAll_" // ModelName
    OPEN(UNIT = 100032,FILE = FileName)
    !
END IF
IF (computeEquilibriumCheck .EQ. 1) THEN
    !
    FileName = "ec_" // ModelName
    OPEN(UNIT = 10004,FILE = FileName)
    !
END IF
IF (computeQGapToMaximum .EQ. 1) THEN
    !
    FileName = "qg_" // ModelName
    OPEN(UNIT = 10006,FILE = FileName)
    !
END IF
IF (computePIGapToMaximum .EQ. 1) THEN
    !
    FileName = "pg_" // ModelName
    OPEN(UNIT = 10007,FILE = FileName)
    !
END IF
IF (computeDetailedImpulseResponseToAll .EQ. 1) THEN
    !
    FileName = "det_irToAll_" // ModelName
    OPEN(UNIT = 100033,FILE = FileName)
    !
END IF
IF (computeDetailedEquilibriumCheck .EQ. 1) THEN
    !
    FileName = "det_ec_" // ModelName
    OPEN(UNIT = 100043,FILE = FileName)
    !
END IF
labelStates = computeStatesCodePrint()
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Loop over models
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
IF (computeMixedStrategies(1) .EQ. 0) THEN
    !
    DO iModel = 1, numModels
        !
        ! Read model parameters
        !
        CALL readModelVariables(10001)
        !
        ! Creating the PI matrix
        !
        ALLOCATE(indexStrategies(lengthStrategies,numGames),indexLastState(LengthStates,numGames))
        IF (typePayoffInput .EQ. 0) CALL computePIMatricesGiven(DemandParameters,NashPrices,CoopPrices,&
            PI,NashProfits,CoopProfits, &
            indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
        IF (typePayoffInput .EQ. 1) CALL computePIMatricesSinghVives(DemandParameters,NashPrices,CoopPrices,&
            PI,NashProfits,CoopProfits, &
            indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
        IF (typePayoffInput .EQ. 2) CALL computePIMatricesLogit(DemandParameters,NashPrices,CoopPrices,&
            PI,NashProfits,CoopProfits, &
            indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
        IF (typePayoffInput .EQ. 3) CALL computePIMatricesLogitSigma0(DemandParameters,NashPrices,CoopPrices,&
            PI,NashProfits,CoopProfits, &
            indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
        PIQ = PI**2
        avgPI = SUM(PI,DIM = 2)/numAgents
        avgPIQ = avgPI**2
        !
        ! Creating I/O filenames
        !
        i = 1+INT(LOG10(DBLE(numModels)))
        WRITE(ModelNumber, "(I0.<i>, A4)") codModel, ".txt"
        FileNameIndexStrategies = "indexStrategiesTransposed_" // ModelNumber
        FileNameIndexLastState = "indexLastState_" // ModelNumber
        FileNamePriceCycles = "priceCycles_" // ModelNumber
        !
        ! Print message
        !
        WRITE(*,11) iModel, numModels, numCores
    11  FORMAT('model = ', I6, ' / numModels = ', I6, ' / numCores = ', I6)  
        !
        ! Compute QL strategy 
        !
        IF (computeQLearningResults .EQ. 1) THEN
            !
            IF (computeRestart .EQ. 0) THEN
                !
                CALL computeModel(iModel,codModel,alpha,ExplorationParameters,delta, &
                    converged,indexLastState,timeToConvergence)
                !
            ELSE IF (computeRestart .EQ. 1) THEN
                !
                CALL computeModelRestart(iModel,codModel,alpha,ExplorationParameters,delta, &
                    converged,indexLastState,timeToConvergence)
                !
            END IF
            !
            ! Output of convergence results
            !
            CALL computeIndicators(iModel,converged,timeToConvergence) 
            !
        END IF
        !
        ! Results at convergence
        ! 
        IF (computeConvergenceResults .EQ. 1) CALL ComputeConvResults(iModel)
        !
        ! Impulse Response analysis to one-period deviation to static best response
        ! 
        IF (computeImpulseResponseToBR .EQ. 1) CALL computeIRAnalysis(iModel,10003,-1)
        !
        ! Impulse Response to a permanent or transitory deviation to Nash prices
        !
        IF (computeImpulseResponseToNash .GE. 1) CALL computeIRAnalysis(iModel,100031,computeImpulseResponseToNash)
        !
        ! Impulse Response analysis to one-period deviation to all prices
        !
        IF (computeImpulseResponseToAll .EQ. 1) THEN
            !
            DO i = 1, numPrices
                !
                CALL computeIRAnalysis(iModel,100032,-i)
                !
            END DO
            !
        END IF
        !
        ! Equilibrium Check
        !
        IF (computeEquilibriumCheck .EQ. 1) CALL computeEqCheck(iModel)
        !
        ! Q and Average PI Gap w.r.t. Maximum
        !
        IF (computeQGapToMaximum .EQ. 1) CALL computeQGapToMax(iModel)
        !
        ! Q and Average PI Gap w.r.t. Maximum
        !
        IF (computePIGapToMaximum .EQ. 1) CALL computeAvgPIGapToMax(iModel)
        !
        ! Detailed Impulse Response analysis to one-period deviation to all prices
        !
        IF (computeDetailedImpulseResponseToAll .EQ. 1) CALL computeDetailedIRToAll(iModel)
        !
        ! Detailed EquilibriumCheck analysis 
        !
        IF (computeDetailedEquilibriumCheck .EQ. 1) CALL computeDetailedEC(iModel)
        !
        ! Deallocate arrays
        !
        DEALLOCATE(indexStrategies,indexLastState)
        !    
        ! End of loop over models
        !
    END DO
    !
END IF
!
! Mixed strategies results
!
IF (computeMixedStrategies(1) .GT. 0) THEN
    !
    ALLOCATE(indexStrategies(lengthStrategies,numGames),priceCycles(numPeriods*numAgents,numGames), &
        sampledIndexStrategies(lengthStrategies,numGames),sampledPriceCycles(numStates+2,numAgents*numGames), &
        alpha_tmp(numAgents),beta_tmp(numAgents),delta_tmp(numAgents))
    !
    ! Read model parameters
    !
    FileName = "ConvResults_" // ModelName
    FileNameMSR = 'MSRes'
!    i = 1+INT(LOG10(DBLE(numModels)))
    OPEN(UNIT = 100022,FILE = FileName,READONLY,IOSTAT = errcode)
    DO iAgent = 1, numAgents
        !
        REWIND(UNIT = 100022)
        IF (computeMixedStrategies(iAgent) .GT. 1) THEN
            !
            READ(100022,2534) 
2534        FORMAT(<computeMixedStrategies(iAgent)-1+1>(/))
            BACKSPACE(UNIT = 100022)
            !
        END IF
        READ(100022,*) i, alpha_tmp, beta_tmp, delta_tmp, &
            DemandParameters, NashPrices, CoopPrices
        IF (typeExplorationMechanism .EQ. 3) beta_tmp = -DBLE(itersPerYear)/DBLE(numAgents+1)* &
                LOG(1.d0-(DBLE(numPrices-1)/DBLE(numPrices))**numAgents/(DBLE(numStates*numPrices)*beta_tmp))
        alpha(iAgent) = alpha_tmp(iAgent)
        MExpl(iAgent) = beta_tmp(iAgent)
        delta(iAgent) = delta_tmp(iAgent)
        !
        WRITE(ModelNumber, "(I0.<1+INT(LOG10(DBLE(numModels)))>, A4)") computeMixedStrategies(iAgent)
        FileNameMSR = TRIM(FileNameMSR) // '_'
        FileNameMSR = TRIM(FileNameMSR) // ModelNumber
        !
    END DO
    CLOSE(UNIT = 100022)
    FileNameMSR = TRIM(FileNameMSR) // '.txt'
    !
    ! Creating the PI matrix
    !
    IF (typePayoffInput .EQ. 0) CALL computePIMatricesGiven(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 1) CALL computePIMatricesSinghVives(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 2) CALL computePIMatricesLogit(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 3) CALL computePIMatricesLogitSigma0(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    !
    ! Compute results
    !
    OPEN(UNIT = 111,FILE = FileNameMSR)
    CALL ComputeMixStratResults ( )
    CLOSE(UNIT = 111)
    !
    ! Deallocate variables
    !
    DEALLOCATE(indexStrategies,priceCycles,sampledIndexStrategies,sampledPriceCycles, &
        alpha_tmp,beta_tmp,delta_tmp)
    !
END IF
!
! Deallocating arrays
!
CALL closeBatch()
!
! Closing output files
!
CLOSE(UNIT = 10001)
IF (computeQLearningResults .EQ. 1) CLOSE(UNIT = 10002)
IF (computeConvergenceResults .EQ. 1) CLOSE(UNIT = 100022)
IF (computeImpulseResponseToBR .EQ. 1) CLOSE(UNIT = 10003)
IF (computeImpulseResponseToNash .GE. 1) CLOSE(UNIT = 100031)
IF (computeImpulseResponseToAll .EQ. 1) CLOSE(UNIT = 100032)
IF (computeEquilibriumCheck .EQ. 1) CLOSE(UNIT = 10004)
IF (computeQGapToMaximum .EQ. 1) CLOSE(UNIT = 10006)
IF (computePIGapToMaximum .EQ. 1) CLOSE(UNIT = 10007)
IF (computeDetailedImpulseResponseToAll .EQ. 1) CLOSE(UNIT = 100033)
IF (computeDetailedEquilibriumCheck .EQ. 1) CLOSE(UNIT = 100043)
!
! End of execution
!
END PROGRAM main