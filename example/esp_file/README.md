# Creating RESP charges from an existing ESP data file

The example here creates and runs the necessary RESP input files for AmberTools using an already exsting ESP data file. 
This can be extremly usefull if an expensive quantum chemistry calculation was already performed to generate the ESP data from another source, or if the resp_fittiny.py program was run and an error occured and you want to continue the computation procedure.
The data file (plot.esp) contains the electrostatic potential evaluated at various points around the molecule. 
The format for this file must follow that of Q-Chem's output:
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

To re-run example, delete or rename the 'mol_resp' directory and call the program with

```
python3 resp_fitting.py -mol mol.pdb -esp plot.esp -out <output_file>
```
The output file is optional. If not included, output will print to the colsole instead.
