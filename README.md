# SU2-CoSim Interface
The repository presents several examples to demonstrate the possebilities of coupling [SU2](https://su2code.github.io/) with [KratosMultiphysics](https://github.com/KratosMultiphysics) by using [CoSimApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication) to solve FSI problems.


## Overview
- [CoSimulation Application Overview](#cosimulation-application-overview)
- [SU2 - Kratos Examples](#su2-kratos-examples)
- [Developer Guide](#developer-guide)
- [Links to Softwares](#links-to-softwares)

<a name="cosimulation-application-overview"></a>

## CoSimulation Application Overview

The CoSimulation Application contains the core developments in coupling black-box solvers and other software-tools within Kratos Multiphysics. The detailed documentation about the application can be found [here](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication#overview), Kratos CoSim examples can be found [here](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication#examples), OpenFoam - Kratos CoSim FSI examples [here](https://github.com/ashishdarekar/Kratos_OpenFOAM_adapter).

<a name="su2-kratos-examples"></a>

## SU2 Kratos Examples

The repository contains several examples:
- FSI simulation of the flexible Onera M6 wing
- Mixed fidelity models:
  - Low-Low fidelity model
  - Mid-Low fidelity model
  - High-Low fidelity model
  - High-Mid-Low fidelity model
- Optimization
- Kratos FSI example with Compresable Potential and Structural Mechanics applications

"Parallel" examples demonstrate the cases with mixed parallelizations: SU2 uses "mpi" and Kratos "OpenMP" parallel types.

<a name="developer-guide"></a>

## Developer Guide

### General structure of Remote Controlled Adapter for external solver
The demonstrated examples are implemented using **remote_controlled_solver_wrapper**. In this scenario, Kratos controlles the external solvers and send commands to it via additional solver adapter. The solver adapter directly controls the external solver and works as an interface to send and recive data. The main functions that has to be implemented are:

- **AdvanceInTime**: Move forward the solution in time.
- **InitializeSolutionStep**: Solver initialization before solving time step.
- **Predict**: Exchange residual information to accelarate coupling solution.
- **SolveSolutionStep**: Solve current state of the model.
- **FinalizeSolutionStep**: Solver finalization after solving time step.
- **OutputSolutionStep**: Output solution
- **Import/ExportData**: Collect data information from common interface, combine it as a data vector and send it to Kratos or vice versa.
- **Import/ExportMesh**: Read common interface from external solver and create ***Kratos.Model_Part*** to exchange mesh infromation with Kratos.

### Using Different parallelization across solvers

In several scenrious, it is suitable to use different partitions for different solvers. In this case, the communication is esteblished only on the root processor. Others are waiting for signal from root to execute required functions. One implementation example can be found in **Onera_FSI** folder.

### Alternative implementations

Alternatevily, one can implement the **Remote Controlled Adapter** directly as a ***remote_controlled_solver_wrapper*** in Kratos. Please, refer to **CoSimApp** [developer guide](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication#user-guide).

### Remark about Node Ids

While creating ***Kratos.ModelPart*** one should make sure to add nodes with sorted node ids. Otherwise, creating of the ***Kratos.ModelPart*** can take a long time because internaly the nodes are going to be sorted every time as a new node is added.

<a name="links-to-softwares"></a>

## Links to Softwares
The original link to softwares are: [SU2](https://su2code.github.io/), [KratosMultiphysics](https://github.com/KratosMultiphysics).