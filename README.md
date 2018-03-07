# BetaScatt GEANT4 simulation App for scattered electrons

BetaScatt is a GEANT4 simulation application that is used
to estimate the amount of scattered betas that reach the
detector in several the emission channeling setups from
IKS, CTN and ISOLDE.


## Installation

Before installation the user needs to have ROOT and GEANT4
installed. Also make sure environment variables for ROOT
and GEANT4 are correctly set.

The software has been tested with ROOT 6.08
and GEANT4 10.4.

To download clone the repository

```console
$ git clone https://github.com/eric-presbitero/BetaScatt.git
```

To install create a build directory next to the source.

```console
$ cd BetaScatt
$ cd build
$ cmake -DGeant4_DIR=(...)/geant4.10.04-install/lib/Geant4-10.4.0 ../G4BetaScat
$ make
```

## Configuring a Run

Run configuration is, for coherence, saved in a file ending in
 "\_r?.data", where ? can be any number.
Several examples can be found in folder "BetaScatt/Test"

Configuration values for the parameters MeanDepth and SigmaDepth
need to be taken from TRIM simulation with the correct incidence angle

## Running a simulation

To run a simulation locally do

```console
$ ./BetaScatt file_r.data "$PWD"
```

The first argument is the run configuration file.

The second argument is the path to the output folder.
In the example it is set to the current directory.


## Contributors
- Bart de Vries
- Ligia Amorim
- Angelo Costa
- Eric David Bosne


## License

The software is under an MIT license which can be found
in the LICENSE file
