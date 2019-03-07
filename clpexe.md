The CLP Executable
==================

Quick Start
===========

The result of `make unitTest` (executed in the build directory) is an executable
`clp` as well as the CLP and COIN-OR libraries. The executable can be used
to perform various unit tests, but can also be used as a standalone
solver. As the executable has a very simple solution file format, the
user may wish to modify `ClpMain.cpp`, which contains the
source of the executable (modifications could even be offered as a
contribution to CLP).

The `clp` executable operates in command line mode or prompted mode.
Entering `clp` will invoke the prompted mode, while `clp <filename>`
will import a problem in MPS format from `filename`, solve it using the
dual simplex method and exit. The command
`clp <filename> -primalsimplex` instructs the executable tp import a
file and solve using the primal simplex method. An additional solitary
dash ("`-`") starts the prompt mode once the execution of the initial
command has been completed. The "`-`" is necessary as part of the
command; invoking prompt mode as a separate command will result in the
loss of problem information related to the initial command. So, the
following sequences of commands are equivalent in the sense that both
maximize a problem using the dual simplex method and write a solution to
file: `solfile`:
```
$ clp filename -maximize -dualsimplex -solution solfile
```

```
$ clp filename -maximize -
Clp:duals
Clp:solution solfile
Clp:quit
```

The executable is at a very early stage of development. Comments and
suggestions would be appreciated.

Online Help and Basic Usage
===========================

The executable has some command-completion functionality as well as some
online help. Below is a table with some examples which summarize these
capabilities.

| Command      | Result
|--------------|-------------------------------------------------------------------------------------------------------------------------------------
| `?`          | Gives a list of all commands |
| `p?`         | Gives a list of all commands which begin with `p`. |
| `p??`        | Gives a list of all commands which begin with `p`, with a short explanation for each. |
| `primals??`  | If is this is enough to uniquely determine a command (in this example, `primalS`, for primal simplex), a long explanation is given. |

In addition, matching a name without a `?` will either execute the command
or give the value of the corresponding parameter as follows: `primalw`
will give the current value of the `primalWeight` parameter while
`primalw 1.0e7` will change it to `1.0e7`.

A Sample Session
================

Below is a sample CLP executable prompt-mode session. A small problem is
loaded and solved under various conditions with the primal and dual
simplex methods. Note the use of the `allslack` command; it sets the
basis to all slacks and resets the solution.

```
$ clp
Coin LP version 0.99.9, build Sep 14 2004
Clp takes input from arguments ( - switches to stdin)
Enter ? for list of commands or help
Clp:import../Data/Sample/p0033.mps
At line 15 NAME  P0033
At line 16 ROWS
At line 34 COLUMNS
At line 109 RHS
At line 118 BOUNDS
At line 152 ENDATA
Problem P0033 has 16 rows, 33 columns and 98 elements
Model was imported from ./../Data/Sample/p0033.mps in 0 seconds
Clp:primals
Presolve 15 (-1) rows, 32 (-1) columns and 97 (-1) elements
0  Obj 0 Primal inf 27.2175 (10) Dual inf 6.42094e+11 (32)
32  Obj 2520.57
Optimal - objective value 2520.57
After Postsolve, objective 2520.57, infeasibilities - dual 0 (0), primal 0 (0)
Optimal objective 2520.571739 - 32 iterations time 0.012, Presolve 0.01
Clp:max
Clp:primals
Presolve 11 (-5) rows, 25 (-8) columns and 84 (-14) elements
0  Obj 4807.92 Dual inf 1700.71 (15)
End of values pass after 2 iterations
2  Obj 4921.7 Dual inf 580.637 (5)
9  Obj 5299.7
Optimal - objective value 5299.7
After Postsolve, objective 5299.7, infeasibilities - dual 643.608 (9), primal 27.0826 (10)
Presolved model was optimal, full model needs cleaning up
0  Obj 5299.7
0  Obj 5299.7
Optimal - objective value 5299.7
Optimal objective 5299.698868 - 9 iterations time 0.022, Presolve 0.02
Clp:allslack
Clp:duals
Presolve 11 (-5) rows, 25 (-8) columns and 84 (-14) elements
0  Obj 2752 Primal inf 24.4867 (6) Dual inf 4280.55 (25)
8  Obj 5299.7
Optimal - objective value 5299.7
After Postsolve, objective 5299.7, infeasibilities - dual 704.58 (8), primal 27.0792 (10)
Presolved model was optimal, full model needs cleaning up
0  Obj 5299.7
0  Obj 5299.7
Optimal - objective value 5299.7
Optimal objective 5299.698868 - 8 iterations time 0.032, Presolve 0.01
Clp:min
Clp:duals
Presolve 15 (-1) rows, 32 (-1) columns and 97 (-1) elements
0  Obj 5299.7 Dual inf 4632.26 (28)
16  Obj 2520.57
Optimal - objective value 2520.57
After Postsolve, objective 2520.57, infeasibilities - dual 2052.5 (13), primal 27.1143 (10)
Presolved model was optimal, full model needs cleaning up
0  Obj 2520.57
0  Obj 2520.57
Optimal - objective value 2520.57
Optimal objective 2520.571739 - 16 iterations time 0.012, Presolve 0.01
Clp:allslack
Clp:presolve off
Clp:primals
0  Obj 0 Primal inf 27.2175 (10) Dual inf 6.39167e+11 (32)
32  Obj 2520.57
Optimal - objective value 2520.57
Optimal objective 2520.571739 - 32 iterations time 0.002
Clp:allslack
Clp:maxIt 10
maxIterations was changed from 99999999 to 10
Clp:primals
0  Obj 0 Primal inf 27.2175 (10) Dual inf 6.39167e+11 (32)
Stopped - objective value 4.24664e+10
Stopped objective 4.246637759e+10 - 10 iterations time 0.002
Clp:quit
```
