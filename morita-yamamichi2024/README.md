# Character displacement or priority effects: Immigration timing can affect community assembly with rapid evolution

These codes are to simulate and draw figures of our paper, "Character displacement or priority effects: Immigration timing can affect community assembly with rapid evolution." Our paper is an updated version of Morita and Yamamichi (2023) Popul. Ecol. by introducing evolution of a quantitative trait (character displacement and local adaptation) to a population dynamics model of two competitor species assuming a trade-off between fecundity and interspecific resource competition.

## For recording results

First of all, you must put a directory, `morita-yamamichi2024`, under the home directory. All the results of simulation are put into an existing directory, `result` under `morita-yamamichi2024`.

# 1. Source File

The following table shows that source files correspond to figures. Please source `.R` file and evaluate `.nb` files below.

+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Figure                         | file             | function                                                                                                                                                                                                                          |
+:==============================:+:================:+:=================================================================================================================================================================================================================================:+
| 2a, 2b, 2c, S10a.b.d.e, S16a-d | `figure02abc.R`  | Draw population and trait dynamics of two species in the (discrete-time) Leslie-Gower model.                                                                                                                                      |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2d, S9a-g, S10c.f              | `figure02d.R`    | Draw trajectories of simulation results on a nullcline plot in the (discrete-time) Leslie-Gower model.                                                                                                                            |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3a, S4a-e, S16e                | `figure03a.R`    | Draw eco-evolutionary outcomes following immigration timings of species 2.                                                                                                                                                        |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3b, S4f, S7, S8h, S9h          | `figure03b.R`    | Calculate analytical solutions of parameter thresholds on which coexistence is always possible and extinction always occurs, and draw the solutions and simulation results.                                                       |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S2a, S2b, S2c                  | `figureS02abc.R` | Draw population and trait dynamics of two species in the (continuous-time) Lotka-Volterra competition model.                                                                                                                      |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S2d                            | `figureS02d.R`   | Draw trajectories of simulation results on a nullcline plot in the (continuous-time) Lotka-Volterra competition model.                                                                                                            |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S3                             | `figureS03.nb`   | Draw simulation results that were consistent with analytical solutions in the coexistence and extinction borders in the Leslie-Gower model.                                                                                       |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S5                             | `figureS05.R`    | Calculate analytical solutions of parameter thresholds on which coexistence is always possible and extinction always occurs, and draw the solutions and simulation results.                                                       |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S6                             | `figureS06.R`    | Draw the difference of eco-evolutionary outcomes between at early and late immigration of species 2.                                                                                                                              |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S8a-g                          | `figureS08.nb`   | Draw nullclines in various parameter values in the Leslie-Gower model. `figureS6.nb` gives a figure of S6 (d).                                                                                                                    |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S11                            | `figureS11.nb`   | Draw conditions for ecological priority effects in the consumer-resource model.                                                                                                                                                   |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S12                            | `figureS12.R`    | Draw simulation results of population and trait dynamics in the continuous-time consumer-resource model.                                                                                                                          |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S13                            | `figureS13.R`    | `figureS11.R` draws a figure of a phase diagram of coexistence, competitive exclusion, and ecological priority effects in the consumer-resource model. `figureS11.nb` gives equilibria.                                           |
|                                |                  |                                                                                                                                                                                                                                   |
|                                | `figureS13.nb`   |                                                                                                                                                                                                                                   |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S14                            | `figureS14.R`    | Draw a figure of a phase diagram of niche and competitive ability differences in the Chessonian framework of the Leslie-Gower model.                                                                                              |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S15                            | `figureS15.nb`   | Draw a figure of the change in equilibria of the frequency of species 1 following maximum fecundity and simulation results of the dynamics of populations of species 1 and 2, trait of species 1, and the frequency of species 1. |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| S17, 18                        | `figureS17.nb`   | Draw a figure showing the simulation results of clonal models in the discrete-time Leslie-Gower model and continuous-time Lotka-Volterra model.                                                                                   |
+--------------------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

# 2. file list

## a. R language (function)

-   `figure02abc.R` (Source file for setting parameters & initial values)

-   `figure02d.R` (Source file for setting parameters & initial values)

    -   `research_dynamics_optEtrait.R` (Simulate the dynamics of population and traits)
    -   `make_color_map.R` (Draw a heat map of figure 2d)

-   `figureS05.R` (Source file for setting parameters & initial values)

    -   `extinctBoundary.R` (Calculate parameter conditions determining whether "alternative stable states occurs" or "coexistence is always possible and extinction of an immigrating species occurs")

-   `figureS02d.R` (Source file for setting parameters & initial values)

    -   `calc_model02.R` (Simulate the dynamics of population and traits)

-   `figureS12.R` (Source file for setting parameters & initial values)

    -   `calc_S12.R` (Simulate the dynamics of population and traits)
    -   `graphics_S12.R` (Draw a figure)

## b. Mathematica

-   `figureS03.nb` (Draw a figure of analytical results)
-   `figureS08.nb` (Draw a nullcline plot)
-   `figureS11.nb` (Draw a phase diagram)
-   `figureS13.nb` (Calculate equilibria)

## 3. Versions of Code languages

### a. R langiage

```         
> version
               _                           
platform       x86_64-apple-darwin15.6.0   
arch           x86_64                      
os             darwin15.6.0                
system         x86_64, darwin15.6.0        
status                                     
major          3                           
minor          6.3                         
year           2020                        
month          02                          
day            29                          
svn rev        77875                       
language       R                           
version.string R version 3.6.3 (2020-02-29)
nickname       Holding the Windsock  
>
```

### b. Mathematica

```         
12.1 (Japanese)
```
