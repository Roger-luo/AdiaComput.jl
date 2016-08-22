# AdiaComput.jl

[Quantum Adiabatic Computing](https://en.wikipedia.org/wiki/Adiabatic_quantum_computation) simulator

[![Build Status](https://travis-ci.org/Roger-luo/AdiaComput.jl.svg?branch=master)](https://travis-ci.org/Roger-luo/AdiaComput.jl)

AdiaComput.jl is a package for simulating quantum adiabatic computing implemented in Julia. The design of AdiaComput.jl is aimed to implementing an efficient, pluggable simulator which allows users create their own adiabatic computing routine. Current AdiaComput.jl only supports CPU computing in the simulation. GPU branch of the simulator is old and needs to be further developed.

# Installation

## Dependencies

- [QuComputStates.jl](https://github.com/Roger-luo/QuComputStates.jl)
- [QuDynamics.jl](https://github.com/JuliaQuantum/QuDynamics.jl)
- [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl)
- [QuSAT.jl](https://github.com/Roger-luo/QuSAT.jl)
- [QuBase.jl](https://github.com/JuliaQuantum/QuBase.jl)
- [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl.git)

Use the following commands in Julia to install dependencies.

```julia
Pkg.clone("https://github.com/Roger-luo/QuComputStates.jl")
Pkg.clone("https://github.com/JuliaQuantum/QuDynamics.jl")
Pkg.add("JuMP")
Pkg.clone("https://github.com/Roger-luo/QuSAT.jl")
Pkg.clone("https://github.com/JuliaQuantum/QuBase.jl")
Pkg.add("Ipopt")
```

The package has not been submitted into METADATA.jl, therefore, if you want to use this package, you need to use `Pkg.clone`

```julia
Pkg.clone("https://github.com/Roger-luo/AdiaComput.jl")
```

# User Guide

For adiabatic computing, you need to generate your own problem Hamiltonian, denoted as `HP` in following codes.

For example, we could use the build-in `QuSAT.jl` interface to generate a 4-qubit Hamiltonian for [exact cover problem](https://en.wikipedia.org/wiki/Exact_cover).

```julia
using QuDynamics
using AdiaComput
using QuSAT
using QuBase
using QuComputStates

ins,ans = generate(4)
pH = QuArray(pHamiltonian(ins,4),(comput_basis(4),comput_basis(4)))
```

Then you will need to create an AQC obeject, there are few key words that you could tweak at the beginning.

- `maxtime` max evolution time.
- `method` method for Shrodinger equation solver. All the method in `QuDynamics.jl` is allowed.
    - ODE solvers `QuODE45`,`QuODE23s`,`QuODE78`
    - Expm solvers `QuExpokit`,`QuExpmV`

```julia
aqc = AQC(pH,4;maxtime=100)
```

use `|>` operator to insert routines into your simulated quantum adiabatic computer. In this example, we simply simulate the normal routine of an adiabatic computing.

```julia
(0.0:1e-2:100,QuODE45())|>aqc
for (t,psi) in aqc end
```

The tuple on the left of `|>` consist of two parts

- time list for your propagator (evolution)
- your propagator/evolution routine

The AQC's states is stored in `aqc.state`, you could check the success probability of your evolution.

```
p = norm(aqc.state[ans[1]+1])^2
```

# Custom routines

To custome the routines, a subtype of `QuDynamics.QuPropagatorMethod` will need to be defined. And you will need to define your own `propagate` function, which returns an `QuArray` as next state.
