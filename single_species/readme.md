In this folder, the codes describe how the optimal harvesting policy for the single species case is generated among codes that that are somehow useful for the single species case.

In the HJB_WO_Discount.m the value function and optimal policy is found for the logistic growth system without discounting. It also plots the population dynamics using the generated policy.

HJB_With_Discounting.m generates the value function and optimal policy while considering discounting. It also plots the population dynamics using the generated policy.

Stability_test.m computes number of numerical solutions under different conditions to visually confirm the region where the numerical solution is stable.

msy_single.m plots the maximum sustainable yield state as well as the optimal harvest state.

ScalarStdWienerProcess.m is the noise generator for the system. This code snippet was provided in a course I attended in the spring of 2018. The course is called 02685 Scientific Computing for differential equations and was taught by Allan Peter Engsig-Karup, John Bagterp Jørgensen and Dimitri Boiroux at the Technical Yniversity of Denmark
