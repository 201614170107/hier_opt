# Constrained multilevel hierarchical optimization for energy markets 
Distributed control strategy for the design of an energy market. The method relies on a hierarchical structure of aggregators (Queens) for the coordination of prosumers (Woekers). The hierarchy reflects the voltage level separations of the electrical grid and allows aggregating prosumers in pools, while taking into account the grid operational constraints. To reach optimal coordination, the prosumers communicate their forecasted power profile to the upper level of the hierarchy. Each time the information crosses upwards a level of the hierarchy, it is first aggregated, both to strongly reduce the data flow and to preserve the privacy. The mathematical derivation of the algorithm is described in [1].

Requirements
=====
The script uses the [YALMIP](https://yalmip.github.io/) matlab package, which is freely available. The script uses the [GUROBI](http://www.gurobi.com/) solver by default, for which free academic licences are available. If not installed, the script reverts to `quadprog`.    

Usage
=====
The `main.m` script generate a communication hierarchy, defined by the `grid_scheleton` matlab structure. An instance of the class `hieragg` is then built and used to solve the first time step of the cooridnation problem decomposed among the Workers. Note that, depending on the uncontrolled power profiles of the Workers and on the initial state of charge of the batteries, the grid constraints cannot always be satisfied. In this case the convergence criterion is on the difference between primal error and constraints violations. 
Feasibility can be restored adding slack batteries in each branch of the structure. This can be done specifying 
```matlab
slack_b = true;
```
in the `main.m` script. In this case the proximal parameter `t_prox` should be changed to avoid bumpy behavior in the convergence. 

Licence
=====
The code is distributed under MIT licence.

References
=====

 [[1]](https://arxiv.org/abs/1803.03560) *L.Nespoli, V.Medici, Constrained hierarchical networked optimization for energy markets, 2018*
