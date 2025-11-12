# Master Thesis Project on 22Al

This is my masters thesis project, which concerns itself with data of the beta decay of 22Al taken at FRIB in June-July 2023.

The analysis consists of an analysis of the beta-delayed single proton emission and of an analysis of the beta-delayed two-proton emission of the 22Al nucleus.
Inspiration is drawn from previous work on the same nucleus, where the focus was the beta-delayed two-proton emission.
This work can be found in Erik's git repo: https://gitlab.au.dk/ausa/erik/e21010

## Prerequisites
In order to run the scripts presented some prerequisites are needed:

### Packages needed for data analysis
* ROOT
* GSL
* [Ausalib](https://gitlab.au.dk/ausa/ausalib)
* [Telescope](https://gitlab.au.dk/ausa/erik/telescope) package from Erik
* [libconfig++](https://hyperrealm.github.io/libconfig/)

### Other folders and files needed

* Updated [setup](https://gitlab.au.dk/ausa/e21010) folder:
* Unpacked data or sorted data.
* Ausalib [Sorter](https://gitlab.au.dk/ausa/sorter) (unless in possession of already sorted data)

## Data Analysis

### Building the project 
From the root of the project folder, do

```shell
mkdir build
cd build
cmake ..
cmake --build .
```


### Running the analysis
Base C++ analysis is run using: 

```shell
./path-to-execution-script <config-file.cfg> <run-number>
```
Additionally one can utilize parallel to run the analysis more efficiently.
The following is an example of the run command of a specific analysis of 22Al that runs the specified analysis on all the run numbers from the 22Al data:

```shell
parallel -u ./analysis/genanalysis/mainanalyzer singleprotons.cfg ::: $(seq 40 40 ; seq 45 52 ; seq 55 55 ; seq 57 59 ; seq 61 61 ; seq 63 74 )
```

Additional Python analysis in run from within the specific notebooks. From the input folders in the python scripts one can see what the required C++ analysis is required for running the individual scripts.

## Thesis
The thesis can be found [here](https://wiki.kern.phys.au.dk/Unpublished%20reports%20and%20notes).


