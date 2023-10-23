# ffit

Quick and dirty code to fit emperical forcefield nonbonded parameters to ab initio single points.

## Building

```
git clone https://github.com/aehogan/ffit.git
cd ffit
mkdir build
cd build
cmake ..
make
```

## Testing

```
cd ffit/examples
../build/ffit he.inp
```

## Data files

See examples/co2.data or examples/he.data for the format of ab initio input data. Briefly, the data files are modified XYZ format files with the ab initio energy in K/kB on the comment line and a molecule id number added before the XYZ data and an atomic charge added after.

## Input commands

- Lines starting with ! or # are treated as comments
- steps [int]
  - number of fitting steps, can be zero to only print out information of the current fit
- output_freq [int]
  - after how many steps to print fitting information
- max_energy [float]
  - the energy scaling parameter to start reducing the important of very repulsive configurations
- data_file [str]
  - filename for an ab initio data file
- es [on|off]
  - permanent electrostatics on / off
- pol [on|off]
  - induced dipole electrostatics on / off
- lj [on|off]
  - Lennard-Jones potential on / off
- phahst [on|off]
  - modified tang toennies (exponential 6-8-10) on / off
- type_map [str] [int]
  - atom name to atom id
- lj_sig [int] [float]
  - atom id and initial Lennard-Jones sigma (or exponential rho) in angstrom
- lj_eps [int] [float]
  - atom id and initial Lennard-Jones epsilon (or exponential beta) in K/kB (or angstrom^-1)
- c6 [int] [float]
  - atom id and initial c6 in a.u.
- c8 [int] [float]
  - atom id and initial c8 in a.u.
- c10 [int] [float]
  - atom id and initial c10 in a.u.
- alpha [int] [float]
  - atom id and initial polarizability in angstrom^3
