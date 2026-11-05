# OpenRXN
A free, open-source tool for modeling chemical reaction networks in Python.  Offers a single platform for pharmacokinetics, pharmacodynamics, and reaction-diffusion equations. 

## Features
* definition of compartments, reactions and species as Python objects
* capable of stochastic (e.g. Gillespie) or deterministic modeling (e.g. ODEs)
* uses Pint for units throughout
* easy definition of 1D, 2D or 3D compartment arrays
* can visualize compartment connectivity as graphs (NetworkX)
* ligand rebinding near a target (`Vicinity`, ODE only)

## Dependencies
* NetworkX
* Numpy
* Scipy
* Pint
* pandas
* matplotlib

## Installation with pip

First, clone the source from github:

```
git clone https://github.com/ADicksonLab/OpenRXN.git
```

Change into that directory and install using pip:

```
cd OpenRXN
pip install .
```

Check the install by running one of the examples.  E.g.:
```
python examples/1D_diffusion.py
```

## Quick start

```python
from openrxn import unit
from openrxn.reactions import Reaction, Species
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem

D, R, DR = Species('D'), Species('R'), Species('DR')
binding = Reaction('binding', [D, R], [DR], [1, 1], [1],
                   kf=1e6/(unit.molar*unit.sec), kr=0.1/unit.sec)

well = Compartment('well', volume=1e-4*unit.L)
well.add_rxn_to_compartment(binding)
system = ODESystem(Model(compartments=[well]).flatten())   # always flatten first

idx = system.state.index['well']                 # species -> position in the state vector
system.set_q([idx['D']], 10e-9*unit.molar)       # concentrations are multiplied by the volume
system.set_q([idx['R']], 1e-9*unit.molar)
result = system.run(600)                         # seconds; result.t, result.y (molecule counts)
R_tot = result.y[idx['R'], -1] + result.y[idx['DR'], -1]
print('occupancy at 600 s:', result.y[idx['DR'], -1] / R_tot)   # ~0.09 = 10 nM / (10 nM + Kd)
```

Use `from openrxn import unit` for every quantity: Pint quantities from another
unit registry cannot be combined with OpenRXN's.

## Running the tests

```
pip install pytest
python -m pytest tests
```

The suite takes about 2 minutes; `examples/membrane_slab.py` alone takes about 50 s.

## Conventions

* `set_q` accepts a number of molecules, an amount (`mol`, `mmol`, ...) or a
  concentration (`mol/L`, `unit.molar`, `unit.nanomolar`, ...). A concentration needs a
  compartment with a 3D volume.
* `system.state.index[compartment][species]` gives the position of a species in the
  state vector. Species are ordered as they first appear in the compartment's
  connections and reactions, so the order is the same in every Python session; still,
  look positions up by name rather than hard-coding them.
* `run(T)` returns the whole trajectory from 0 to T. Reporters record the state at
  t = 0 and then every `freq`.
* A `ResConnection` to a `Reservoir` holds the reservoir concentration at the shared
  face, so the distance is half the cell width. Rates are twice those of versions
  before November 2026, which used the full width.
* A link from A to B with rates `(k_out, k_in)` needs the link from B to A to carry
  `(k_in, k_out)`. `stack` and `join3D` set this up (k_out is from `self` into the other
  array); links inside one array must have equal rates both ways. `Model.flatten()`
  warns about any link that does not mirror its reverse, since such a link does not
  conserve mass.

## Limitations

* The ODE Jacobian is dense, so models with thousands of states are slow and use every
  BLAS thread; cap them with e.g. `OPENBLAS_NUM_THREADS=8`.
* `GillespieSystem` supports reservoirs with constant concentrations only
  (`Reservoir(concs=...)`), does not support `Vicinity`, and draws from NumPy's global
  random state (seed it with `np.random.seed`).
