# Generating RESP Charges with Q-Chem

This example generates RESP charges using Q-Chem to generate ESP data. The ESP data (stored in plot.esp) is then passed to AmberTools to generate RESP charges using the default method settings.

To re-run example, delete the 'adenine_resp' directory (or change the 'name' variable in the input file to something other than 'adenine'), and call the program with

```
python3 resp_fitting.py -mol adenine.pdb -ipt input -out <output_file>
```
