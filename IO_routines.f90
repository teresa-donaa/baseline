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
        <3*numAgents+numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
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
! End of execution
!
END MODULE IO_routines