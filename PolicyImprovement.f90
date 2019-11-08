MODULE PolicyImprovement
!
USE globals
USE EquilibriumCheck
USE QL_routines
!
! Applies Generalized Policy Improvement on detailed strategies at convergence
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ImprovePolicy ( )
    !
    ! Computes equilibrium check for an individual replication
    !
    ! INPUT:
    !
    ! - Strategy         : initial strategy for all agents
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: i, iGame, iAgent, iState, StrategyPrice, iPrice, VisitedStates(numPeriods), &
        PreCycleLength, GameCycleLength, iPeriod, numImprovedPrices, ImprovedPrices(numStates), &
        numStatesBRAll(numAgents), numStatesBROnPath(numAgents), numStatesBROffPath(numAgents), &
        numStatesEQAll, numStatesEQOnPath, numStatesEQOffPath, numImprovements
	INTEGER :: OptimalStrategyVec(lengthStrategies), OptimalStrategy(numStates,numAgents)
    INTEGER, DIMENSION(numStates,numAgents) :: OldStrategy, NewStrategy
    REAL(8) :: StateValueFunction(numPrices), MaxStateValueFunction, TestDiff, TestCrit
    !
    ! Beginning execution
    !
    ! Reading strategies and states at convergence from file
    !
    CALL ReadInfoModel()
    !
    ! Beginning loop over games
    !
    OPEN(UNIT = 996,FILE = "PolicyImprovementStrategies.txt", STATUS = "REPLACE")
    DO iGame = 1, numGames              ! Start of loop over games
        !
        OptimalStrategyVec = indexStrategies(:,iGame)
        OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        NewStrategy = OptimalStrategy
        DO i = 1, 200                       ! Start of contraction loop
            !
            OldStrategy = NewStrategy
            !
            ! For each agent A and each state S, find the optimal price
            !
            DO iState = 1, numStates            ! Start of loop over states
                !
                ! Compute state value function for OptimalStrategy in iState, for all prices and agents
                !
                DO iAgent = 1, numAgents            ! Start of loop over agents
                    !
                    StateValueFunction = 0.d0
                    DO iPrice = 1, numPrices            ! Start of loop over prices to compute a row of Q
                        !
                        CALL computeQCell(OldStrategy,iState,iPrice,iAgent,delta, &
                            StateValueFunction(iPrice),VisitedStates,PreCycleLength,GameCycleLength)
                        !
                    END DO                              ! End of loop over prices to compute a row of Q
                    !
                    MaxStateValueFunction = MAXVAL(StateValueFunction)
                    numImprovedPrices = 0
                    ImprovedPrices = 0
                    DO iPrice = 1, numPrices            ! Start of loop over prices to find optimal price(s)
                        !
                        TestDiff = ABS(StateValueFunction(iPrice)-MaxStateValueFunction)
                        IF (TestDiff .LE. TestCrit) THEN
                            !
                            numImprovedPrices = numImprovedPrices+1                        
                            ImprovedPrices(numImprovedPrices) = iPrice
                            !
                        END IF
                        !
                    END DO                              ! End of loop over prices to find optimal price(s)
                    !
                    StrategyPrice = OldStrategy(iState,iAgent)
                    IF (ANY(MAXLOC(ImprovedPrices(:numImprovedPrices)) .EQ. StrategyPrice)) THEN
                        !
                        NewStrategy(iState,iAgent) = OldStrategy(iState,iAgent)
                        !
                    ELSE
                        !
                        NewStrategy(iState,iAgent) = MINVAL(MAXLOC(StateValueFunction))
                        !
                    END IF
                    !
                END DO                              ! End of loop over agents
                !
            END DO                              ! End of loop over states
            numImprovements = COUNT(NewStrategy .NE. OldStrategy)
!            PRINT*, 'iGame = ', iGame, ' ; i = ', i, ' ; numImprovements = ', numImprovements
            IF (numImprovements .EQ. 0) EXIT
            !
        END DO                              ! End of contraction loop
        PRINT*, 'iGame = ', iGame, ' ; i = ', i, ' ; numImprovements = ', numImprovements
        WRITE(996,100) numImprovements, i, (NewStrategy(:,iAgent), iAgent = 1, numAgents)
100     FORMAT(I5, 1X, I5, 1X, <numAgents*numStates>(I<lengthFormatActionPrint>,1X))    
        !
    END DO                              ! End of loop over iGames
    CLOSE(UNIT = 996)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ImprovePolicy
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE PolicyImprovement
