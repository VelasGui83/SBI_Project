Basic procedure for structure minimization
==================================================

We will run a very simple minimization process. Keep in mind that it is still in progress, you may find crashes or bugs. 

It is not appliable to all the complexes. It will only work for basic structures made of common aminoacids. This procedure was tested succesfully with the complexes given in python class.
Please consider if your input will fit the MD requirement.

PDB Files preparation
======================

A pdb file is required. All the cystalographic waters and duplicated residues will be remove with the following script. Furthemore it also will add the hydrogens and SS bridges in the complex

    pdb4amber -i input.pdb -o sample.pdb --dry --reduce

Tleap will charge a general ForceField to the complex. Also it will automatically solvate and neutralize the sample. It writes the output files "parm: Topology" "rst: Coordinates"

    xleap leaprc

LEAPRC CONTENT. This argument file can be found in the git. Consider changing the path for the input for tleap in line 2.

    source leaprc.protein.ff14SB
    sample = loadPDB "/home/brian/Escritorio/test/sample.pdb"
    source leaprc.water.tip3p
    solvatebox sample TIP3PBOX 10.0
    addions sample Na+ 0
    addions sample Cl- 0
    saveAmberParm sample sample_xleap.parm sample_xleap.crd
    quit


Minimization process
=====================
2 files can be found in this git, one for the minimization process parameters an the other one for the heating process parameters. If you desire to change the parameters, take a look at this website tutorial:

http://ambermd.org/tutorials/basic/tutorial0/index.htm

Minimization: Avoid contacts and clashes. To change parameters, go to 01_Min.in file

    sander -O -i 01_Min.in -o 01_Min.out -p sample_xleap.parm -c sample_xleap.crd -r 01_Min.ncrst -inf 01_Min.mdinfo


Heating: Introduces Kinetic energy to the system. Tries to stabilize the system energy. To change parameters, go to 02_Heat.in file.

    sander -O -i 02_Heat.in -o 02_Heat.out -p sample_xleap.parm -c 01_Min.ncrst -r 02_Heat.ncrst -x 02_Heat.nc -inf 02_Heat.mdinfo


Retrieve PDB model
============
Once the inimization processes are done, one needs to retrieve a pdb file. To do so, we use the following comand, which uses the topology and coordinate files for generating the pdb model

    ambpdb -p sample_xleap.parm -c out_coord.mdcrd > final_min.pdb

Now it is time to clear the pdb file. We need to elimiate the solvent molecules, the counterions added, and the hydrogens in the chain. By executing the following script we end up with the final PDB minimized model

    pdbstrip -ysao /path/to/file.pdb

One can see how the energy of the system is minimized (looking at output files). Side chains are rlocated acording to a more realistic situation.


## __All this procedure was intended to be automated, but we ran out of time before it can be done.__