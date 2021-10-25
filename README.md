# RESPFitting

The example here creates and runs the necessary RESP input files for AmberTools. 
The generation of ESP data is designed to be used with Q-Chem program, although ESP data can be generated with other softwares and imported into the fitting procedure.

A typical call to the program will generate ESP points using the Merz-Kollman procedure, call Q-Chem and evalue the ESP at each of those points, and then call AmberTools to generate the RESP charges.

At the moment, only single conformational fits can be performed, and geometry optimization can not be performed durring the fitting process.

# Usage
```
usage: resp_fitting.py [-h] [-ipt IPT] [-pdb PDB] [-out OUT] [-esp ESP]
                       [-chg CHG] -mol MOL

optional arguments:
  -h, --help  show this help message and exit
  -ipt IPT  (REQUIRED if running Q-Chem) input file to change program settings
  -pdb PDB  PDB file with topology information
  -out OUT  file to print program output to
  -esp ESP  ESP data file using Q-Chem output format
  -chg CHG  pre-computed ESP charges, for analysis only
  -mol MOL  (REQUIRED) molecule file, either .pdb or .xyz
  ```
```-mol``` option can either be a PDB or XYZ file. \
```-out``` output file is optional, and if not supplied, the program's output will be directoed to the sonsole instead.\
```-esp``` ESP data file that contains the electrostatic potential at varous points around the molecule. The file must be formatted using the Q-Chem output format:
```
# comment line 1
# comment line 2
# comment line 3
# comment line 4
<x_1> <y_1> <z_1> <esp_1>
<x_2> <y_2> <z_2> <esp_2>
...
<x_N> <y_N> <z_N> <esp_N>
```
All are assumed to be in atomic units. The same format is also used by Q-Chem when generating ESP data.

# input file
The input file, specified with the ```-ipt``` flag, contains the options for both the RESP python program and Q-Chem. 
The format is designed to mimic that of Q-Chem's input file format, and an example of what one might look like is as follows:
```
$resp
    name        adenine                 ! save results in the directory 'resp_adenine'
    charge      0                       ! total molecular charge of 0 e
    vdw_ratios  [1.4, 1.6, 1.8, 2.0]    ! default vdW radii scaling factors
    mk_density  20                      ! 20 Merz-Kollman ESP points per angstrom^2 on each sphere
$end

$rem
    jobtype         sp     ! single point energy calculation
    method          HF     ! Hartree-Fock exchange
    basis           6-31G* ! Popel basis set
    sym_ignore      true   ! do not translate or rotate coords
    resp_charges    true   ! Q-Chem RESP charges (for comparison only)
    chelpg          true   ! Q-Chem CHELPG charges (for comparison only)
$end
```
Using the ```$resp``` section specifies options for the python program.  
All available options and their default values are specified in the example file above. If an input file is not provided, the program will use the default values.
The ```$rem``` section specifies the options for the Q-Chem job.  
Any option that Q-Chem can use for it's input can also be used here, and only the Q-Chem executable itself will check the input file for valid formatting and options.
And additional Q-Chem secions that begin with ```$``` can also be specified in this file.
For example, the ```$cdft``` section can be specified to perform a constrained DFT calcualtion, and if electrostatic data is to be generated, it will do so using the resulting constrained QM calculation.  
**The "ESPGrid" Q-Chem option does not need to be specified here!** The Python program will automatically create an input file with the correct ESP input option specified. Doing so yourself may generate an error with Q-Chem, depending on the version used. 
# Examples
Examples of how to run the program can be found in the "examples" directory and contain aditional READMEs for how to run them.


