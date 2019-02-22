MODULE IO_routines
!    
USE globals
!
! Various input/output routines 
!
IMPLICIT NONE
!
CONTAINS
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStrategyLabelPrint ( )
    !
    ! Generates the stratebies label in the output file
    !
    IMPLICIT NONE
    !
    ! Declaring function's type
    !
    CHARACTER(len = 9), DIMENSION(numPrintStrategies) :: computeStrategyLabelPrint
    !
    ! Declaring local variables
    !
    INTEGER :: i
    !
    ! Beginning execution
    !
    DO i = 1, numPrintStrategies
        !
        IF ((useNashStrategies .EQ. 1) .AND. (i .LE. numNashStrategies)) THEN
            !
            WRITE(computeStrategyLabelPrint(i),"(A8,I1)") ' EqStrat', i
            !
        ELSE
            !
            WRITE(computeStrategyLabelPrint(i),"(A8,I1)") 'OthStrat', i-useNashStrategies*numNashStrategies
            !
        END IF
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStrategyLabelPrint
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStrategyCodePrint ( indexStrategy )
    !
    ! Transforms the strategy index vector code in a printable string (with '.' and a '-')
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(lengthStrategies), INTENT(IN) :: indexStrategy
    !
    ! Declaring function's type
    !
    CHARACTER(len = lengthStrategiesPrint) :: computeStrategyCodePrint
    !
    ! Declaring local variables
    !
    INTEGER :: i
    CHARACTER(len = lengthFormatActionPrint) :: tmp
    !
    ! Beginning execution
    !
    DO i = 1, lengthStrategies
        !
        WRITE(tmp,'(I<lengthFormatActionPrint>)') indexStrategy(i)
        IF (i .EQ. 1) THEN 
            !
            computeStrategyCodePrint = TRIM(tmp)   
            !
        ELSE IF (MOD(i,numStates) .NE. 1) THEN
            !
            computeStrategyCodePrint = TRIM(computeStrategyCodePrint) // '.' // TRIM(tmp)   
            !
        ELSE
            !
            computeStrategyCodePrint = TRIM(computeStrategyCodePrint) // '-' // TRIM(tmp)   
            !
        END IF
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStrategyCodePrint
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE PrintResFile ( iModel, numGamesConverged, &
        meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence ) 
    !
    ! Prints the res*.txt output file
    ! This file contains indicators computed by computeIndicators
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel, numGamesConverged
    REAL(8), INTENT(IN) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, h
    !
    ! Beginning execution
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10002,891) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents)
891         FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>(' MExplPar', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I0.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I0.2, ' ')), &
            '   numConv ', &
            '    avgTTC      seTTC     medTTC ')
        !
    END IF
    !
    WRITE(10002,991) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        numGamesConverged, &
        meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
991 FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.3, 1X), &
        <6*numAgents>(F10.7, 1X), &
        <numPrices*numAgents>(F10.7, 1X), &
        I10, 1X, &
        <3>(F10.2, 1X))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE PrintResFile
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE PrintExtraResFile ( iModel, labelStates, &
        converged, indexStrategies, &
        timeToConvergence, vecProfit, freqStates )
    !
    ! Prints the res*.txt output file
    ! This file contains indicators computed by computeIndicators
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel
    CHARACTER(len = 3+lengthStatesPrint), DIMENSION(numStates), INTENT(IN) :: labelStates
    INTEGER, DIMENSION(numGames), INTENT(IN) :: converged
    INTEGER, DIMENSION(lengthStrategies,numGames), INTENT(IN) :: indexStrategies
    REAL(8), DIMENSION(numGames), INTENT(IN) :: timeToConvergence
    REAL(8), DIMENSION(numGames,numAgents), INTENT(IN) :: vecProfit
    REAL(8), DIMENSION(numStates,numGames), INTENT(IN) :: freqStates
    !
    ! Declaring local variables
    !
    CHARACTER(len = 3) :: iModelsChar
    CHARACTER(len = 20) :: printExpFileName
    INTEGER :: j, iGames
    !
    ! Beginning execution
    !
    IF (printExp .EQ. 1) THEN
        !
        WRITE(iModelsChar,'(I0.3)') iModel
        printExpFileName = 'results' // iModelsChar // '.txt'
        OPEN(UNIT = 102,FILE = printExpFileName)
        WRITE(102,1020) (labelStates(j), j = 1, numStates)
1020    FORMAT('  Game Converged ', <MAX(2,lengthStrategiesPrint-8)>X, 'Strategy ', &
            ' YrsToConv    Profit1    Profit2 ', <numStates>(A<MAX(10,3+lengthStatesPrint)>, ' ')) 
        WRITE(102,1021) (iGames, converged(iGames), ADJUSTR(computeStrategyCodePrint(indexStrategies(:,iGames))), &
            timeToConvergence(iGames), vecProfit(iGames,:), freqStates(:,iGames), &
            iGames = 1, numGames)
1021    FORMAT(<numGames>(I6, 1X, I9, 1X, A<lengthStrategiesPrint>, 1X, &
            <3+numStates>(F10.3, 1X), /))
        CLOSE(UNIT = 102)
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE PrintExtraResFile
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! End of execution
!
END MODULE IO_routines