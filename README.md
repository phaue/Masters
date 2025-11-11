This is my masters thesis project, which concerns itself with data of the beta decay of 22Al taken at FRIB in June-July 2023.

The analysis consists of an analysis of the beta-delayed single proton emission and of an analysis of the beta-delayed two-proton emission of the 22Al nucleus.
Inspiration is drawn from previous work on the same nucleus, where the focus was the beta-delayed two-proton emission.
This work can be found in Erik's git repo: https://gitlab.au.dk/ausa/erik/e21010

In order to run the scripts presented some prerequisites are needed:

### Packages needed for data analysis
ROOT
GSL
[Ausalib](https://gitlab.au.dk/ausa/ausalib)
[Telescope](https://gitlab.au.dk/ausa/erik/telescope) package from Erik
[libconfig++](https://hyperrealm.github.io/libconfig/)

### Other folders and files needed

Updated [setup](https://gitlab.au.dk/ausa/e21010) folder:
Unpacked data. Can be found on Aarhus subatomic physics group's erda.


Analysis is run using ./path-to-execution-script some-config.cfg some-run-number
Additionally one can utilize parallel to run more analysis more efficiently.
The following is an example of the run command of a specific analysis of 22Al that runs the specified analysis on all the specified run numbers:
parallel -u ./analysis/genanalysis/mainanalyzer singleprotons.cfg ::: $(seq 40 40 ; seq 45 52 ; seq 55 55 ; seq 57 59 ; seq 61 61 ; seq 63 74 )