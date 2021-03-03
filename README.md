# DPTrafficFlowControl
These tools rely on dynamic programming for obtaining a closed-loop controller that aims to minimize the Total Travel Time (TTT) of the traffic flow in a three-link merge junction, an elemental component in transportation networks.

The main script is "MergeJunctionMDP.m".

The dynamics of the three-link merge junction are described in:
Samuel Coogan & Murat Arcak. A Benchmark Problem in Transportation Networks. ACM, 2016. 

Please refer to "DPTrafficFlowControl.pdf" for a detailed explanation of the one-step-lookahead value iteration method that it is implemented. You will also find some simulation results, technical notes, and a brief description of the Matlab functions in this repository.
