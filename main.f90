PROGRAM main
!
USE globals
USE LearningSimulation
USE LearningSimulationRestart
USE ConvergenceResults
USE PolicyImprovement
USE ImpulseResponse
USE EquilibriumCheck
USE QGapToMaximum
USE PIGapToMaximum
USE DetailedAnalysis
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
ModelName = "mixedQ_5920_5920.txt"
FileName = "A_mod_" // ModelName
!
OPEN(UNIT = 10001,FILE = FileName)
CALL readBatchVariables(10001)
!
IF (SwitchQLearningResults .EQ. 1) THEN
    !
    FileName = "A_res_" // ModelName
    OPEN(UNIT = 10002,FILE = FileName)
    !
END IF
IF (SwitchConvergenceResults .EQ. 1) THEN
    !
    FileName = "A_convResults_" // ModelName
    OPEN(UNIT = 100022,FILE = FileName)
    !
END IF
IF (SwitchImpulseResponseToBR .EQ. 1) THEN
    !
    FileName = "A_irToBR_" // ModelName
    OPEN(UNIT = 10003,FILE = FileName)
    !
END IF
IF (SwitchImpulseResponseToNash .GE. 1) THEN
    !
    FileName = "A_irToNash_" // ModelName
    OPEN(UNIT = 100031,FILE = FileName)
    !
END IF
IF (SwitchImpulseResponseToAll .EQ. 1) THEn
    !
    FileName = "A_irToAll_" // ModelName
    OPEN(UNIT = 100032,FILE = FileName)
    !
END IF
IF (SwitchEquilibriumCheck .EQ. 1) THEN
    !
    FileName = "A_ec_" // ModelName
    OPEN(UNIT = 10004,FILE = FileName)
    !
END IF
IF (SwitchQGapToMaximum .EQ. 1) THEN
    !
    FileName = "A_qg_" // ModelName
    OPEN(UNIT = 10006,FILE = FileName)
    !
END IF
IF (SwitchPIGapToMaximum .EQ. 1) THEN
    !
    FileName = "A_pg_" // ModelName
    OPEN(UNIT = 10007,FILE = FileName)
    !
END IF
IF (SwitchDetailedAnalysis .EQ. 1) THEN
    !
    FileName = "A_det_" // ModelName
    OPEN(UNIT = 100033,FILE = FileName)
    !
END IF
labelStates = computeStatesCodePrint()
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Loop over models
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
DO iModel = 1, numModels
    !
    ! Read model parameters
    !
    CALL readModelVariables(10001)
    !
    ! When using trained Q matrices, store the name of the directory containing them
    !
    DO iAgent = 1, numAgents
        !
        IF (typeQInitialization(iAgent) .EQ. 'T') THEN
            !
            WRITE(iChar,'(I0.5)') NINT(parQInitialization(iAgent,1))
            QFileFolderName(iAgent) = &
                'C:/Users/sergio.pastorello/Documents/jobs/dynamic pricing/qlearning/baseline/mixed_Q_analysis/Q_' // &
                    iChar // '/'
            !
        END IF
        !
    END DO
    !
    ! Creating the PI matrix
    !
    ALLOCATE(indexStrategies(lengthStrategies,numGames),indexLastState(lengthStates,numGames))
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
    i = 1+INT(LOG10(DBLE(totModels)))
    WRITE(ModelNumber, "(I0.<i>, A4)") codModel, ".txt"
    FileNameInfoModel = "InfoModel_" // ModelNumber
    !
    ! Print message
    !
    WRITE(*,11) iModel, numModels, numCores
11  FORMAT('model = ', I6, ' / numModels = ', I6, ' / numCores = ', I6)  
    !
    ! Compute QL strategy 
    !
    IF (SwitchQLearningResults .EQ. 1) THEN
        !
        IF (SwitchRestart .EQ. 0) THEN
            !
            CALL computeModel(iModel,codModel,alpha,ExplorationParameters,delta)
            !
        ELSE IF (SwitchRestart .EQ. 1) THEN
            !
            CALL computeModelRestart(iModel,codModel,alpha,ExplorationParameters,delta)
            !
        END IF
        !
    END IF
    !
    ! Results at convergence
    ! 
    IF (SwitchConvergenceResults .EQ. 1) CALL ComputeConvResults(iModel)
!@SP
!CALL ImprovePolicy ( )
!@SP
    !
    ! Impulse Response analysis to one-period deviation to static best response
    ! NB: The last argument in computeIRAnalysis is "IRType", and it's crucial:
    ! IRType < 0 : One-period deviation to the price IRType
    ! IRType = 0 : One-period deviation to static BR
    ! IRType > 0 : IRType-period deviation to Nash 
    ! 
    IF (SwitchImpulseResponseToBR .EQ. 1) CALL computeIRAnalysis(iModel,10003,0)
    !
    ! Impulse Response to a permanent or transitory deviation to Nash prices
    !
    IF (SwitchImpulseResponseToNash .GE. 1) CALL computeIRAnalysis(iModel,100031,SwitchImpulseResponseToNash)
    !
    ! Impulse Response analysis to one-period deviation to all prices
    !
    IF (SwitchImpulseResponseToAll .EQ. 1) THEN
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
    IF (SwitchEquilibriumCheck .EQ. 1) CALL computeEqCheck(iModel)
    !
    ! Q and Average PI Gap w.r.t. Maximum
    !
    IF (SwitchQGapToMaximum .EQ. 1) CALL computeQGapToMax(iModel)
    !
    ! Q and Average PI Gap w.r.t. Maximum
    !
    IF (SwitchPIGapToMaximum .EQ. 1) CALL computeAvgPIGapToMax(iModel)
    !
    ! Detailed Impulse Response analysis to one-period deviation to all prices
    !
    IF (SwitchDetailedAnalysis .EQ. 1) CALL ComputeDetailedAnalysis(iModel)
    !
    ! Deallocate arrays
    !
    DEALLOCATE(indexStrategies,indexLastState)
    !    
    ! End of loop over models
    !
END DO
!
! Deallocating arrays
!
CALL closeBatch()
!
! Closing output files
!
CLOSE(UNIT = 10001)
IF (SwitchQLearningResults .EQ. 1) CLOSE(UNIT = 10002)
IF (SwitchConvergenceResults .EQ. 1) CLOSE(UNIT = 100022)
IF (SwitchImpulseResponseToBR .EQ. 1) CLOSE(UNIT = 10003)
IF (SwitchImpulseResponseToNash .GE. 1) CLOSE(UNIT = 100031)
IF (SwitchImpulseResponseToAll .EQ. 1) CLOSE(UNIT = 100032)
IF (SwitchEquilibriumCheck .EQ. 1) CLOSE(UNIT = 10004)
IF (SwitchQGapToMaximum .EQ. 1) CLOSE(UNIT = 10006)
IF (SwitchPIGapToMaximum .EQ. 1) CLOSE(UNIT = 10007)
IF (SwitchDetailedAnalysis .EQ. 1) CLOSE(UNIT = 100033)
!
! End of execution
!
END PROGRAM main