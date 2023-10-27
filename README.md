# SU2-CoSim Interface
The repository presents several examples to demonstrate the possebilities of coupling [SU2](https://su2code.github.io/) with [KratosMultiphysics](https://github.com/KratosMultiphysics) by using [CoSimApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication) to solve FSI problems.


## Overview
- [CoSimulation Application Overview](#cosimulation-application)
- [Developer Guide](#developer-guide)
    - [Structure of the Application](#structure-of-the-application)
    - [How to couple a new solver / software-tool?](#how-to-couple-a-new-solver--software-tool)
      - [Interface of SolverWrapper](#interface-of-solverwrapper)
      - [Remote controlled CoSimulation](#remote-controlled-cosimulation)
    - [Using a solver in MPI](#using-a-solver-in-mpi)
  - [References](#references)
- [SU2 - Kratos Examples](#SU2-KratosExamples)

<a name="cosimulation-application"></a>

## CoSimulation Application Overview

The CoSimulation Application contains the core developments in coupling black-box solvers and other software-tools within Kratos Multiphysics. The detailed documentation about the application can be found [here](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication#overview), Kratos CoSim examples can be found [here](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CoSimulationApplication#examples), OpenFoam - Kratos CoSim FSI examples [here](https://github.com/ashishdarekar/Kratos_OpenFOAM_adapter).

<a name="developer-guide"></a>


## Developer Guide


<a name="developer-guide_structure-of-the-application"></a>


The original link to softwares are: [SU2](https://su2code.github.io/), [KratosMultiphysics](https://github.com/KratosMultiphysics).