Theoretical background & scientific explanation
===
## Introduction
One of the forms of the quaternary structure are the protein complexes. Individual proteins can participate in the formation of a large variety of protein complexes, differences in the composition of these macro-complexes result in a difference in their function. That is the reason because of why is vital to know the structure of these macro-complexes.

That's where the idea of this project came from. In this project an algorithm to create a macro-complex is proposed, based on the interaction / position of complexes with two interacting proteins subunits. Provided a list of these complexes, the program will propose a single solution for the macro-complex structure.

In a general overview, the program recognizes and groups each chain of the PDB files provided. The chains are compared between themselves first by a pairwise alignment and second a structural alignment if certain conditions are met. By this time, the program knows which chains are identical between them, making them a perfect target-template system for the superposition. The program starts from a single complex and begins to build the structure superposing those complexes that share an identical chain provided that no clashes are met. The program ends when it can no longer place more complexes in the structure because of steric clashes (a given complex can be used more than once). The user will obtain the best macro-complex solution according to the algorithm.

Algorithmic aproach
===

## Data analysis

The workflow of the program starts with the proportioning of the input.
Once the user has provided it, the program converts and extracts the information from the PDB files classifying them into program classes. Then, information of individual chains and interacting complexes is saved separately.

## Alignment

First of all, pairwise alignments for each chain, between all of the provided, are carried out using a module named *Levenshtein*. This method is a string metric for measuring the difference between two sequences, in other words, calculates the minimum number of single character edits for the sequences to match, which can be translated into a similarity percentage. The unique purpose of the alignment is seeing which have **more than 90% sequence identity**.

Once the Levenshtein alignment method has been done, the program treats proteins and DNA/RNA separately for the next analysis. For DNA/RNA chains, another Levenshtein alignment will be performed with its sequence (this time the sequence similarity percentage is required to be 100%). However, for protein chains, an structural alignment will be performed.

The structural alignment is carried out using STAMP Software. The program will compare each chain with the ones  which shares more than 90% sequence similarity. The aim of using STAMP is to know which chains are structurally identical to each other (STAMP gives the maximum score to them), by this way we will find common chains shared by the complexes, useful for the superposition process. Furthermore, to avoid unnecessary STAMP calculations, if two chains are identical between them, they will share the same output data (it will be a waste to calculate the same structural alignment once more).
We are aware that the pairwise alignment is irrelevant for the structural alignment to work, but we decided to implement it for future improvements as shown later. 

At this point the program only is able to superpose identical chains (maximum STAMP score). It would become more useful for analyzing orthologs ( little sequence and structure differences). Nevertheless, the pairwise alignment is useful for reducing the chains analyzed by STAMP.

## Superposition

Once all the chains hev been analysed and we know which are identical among them, the program proceeds to the process of superposition. It **starts from a random complex**, (which has two chains to evaluate) from which all the macro-complex will be built. From the initial complex, a chain is selected. From that specific chain a list of possible candidates to superpose will appear, which are the complex that share an identical chain. The program will iterate over those candidates, superposing the identical chains and evaluating if there is any remaining chain from the candidate complex able to relocate without generating steric clashes (the steric clash system is detailed in the *“steric clashes”* section). If no clashes are detected, only the remaining candidate chain is added to the final model. In addition, the recently added chain will also become a new chain to evaluate (same process as the two first chains).

The process will continue iterating every time a new chain is added. As soon as a new chain is added to the final model, there are new possibilities for the candidate complex to fit in. Once the program has checked all the chains to evaluate (repeating the ones added in previous iterations) and there are no more new chains added, it will stop iterating. By this way, we make sure that the program explores all the possible space ending up with an unique model.

In order to apply the transformation coordinates for the superposition, a matrix is calculated. Without going further into details, for the proteins, the matrix is calculated using the CA atoms, however, for DNA/RNA all the atoms are used. As it compares identical chains, we will end up with the same matrix as if we were using all protein atoms but it will be more computational expensive.

It is important note that the algorithm can reuse complexes previously added to the final model, meaning that is able to reconstruct macro-molecular structures if a few PDB complexes are given. Keep in mind that those few PDB files must fulfill the minimum essential interactions between all chains in order to generate the final model.

## Steric clashes

In order to evaluate the steric clashes when superposing structures, 2 approaches were designed. However, only one of them is used on the code.
Given the candidate chain that is not superposing, it calculates the distance between its center of mass and the furthest atoms from the same chain. Then, the program draws an imaginary sphere with the previous distance, plus two angstroms, from the center of mass (those 2 angstroms are required for avoiding possible clashes near the furthest atoms). All the atoms from other chains included in the drawn sphere will be evaluated one by one versus the atoms from the candidate chain. If two atoms are in a distance closer than 1.1 angstroms, it will be considered as a clash. Once there has been ten clashes, the program will discard the candidate chain and will try another chain. 
The cut off distance was selected due to the fact that 1.09 angstroms is approximately the minimum distance on a common protein, between the C(sp3)-H. After testing a bunch of models were only one clash was avoiding the final model to be correct, we decided to accept a **maximum of 10 distances** to be less than 1.1 angstroms (It needs further analysis). It may require further testing.
The other hidden approach compares all the atoms from the candidate chain and the ones which are in the final model.

## Saving the results

By default, the output is in PDB format. However, files having more than 52 chains are saved in CIF format. 

Features
===

Here we present some of the features that macromaker has:
- **No Stoichiometry is required** for the program to run. As you may have deduced, the program will not stop until all the space is explored and no more chains could be added.  As we want the user to do the less work possible, as a team we decided that the stoichiometry could be annoying for the user (as a required input). 
- The user can **limit the number of chains** placed in the final model. By default, the maximum number of structures saved is 200, but it can be easily changed. If you desire to limit or expand the maximum chains in the model, please head to the “README”. This is to avoid infinite macro-complex structures to continue building (I.e. microtubules). 
- When building huge macro-molecular systems, the user may need to stop the job without losing all the computational time done, that is why a **Pause-Recover system** is included in the program. If you want to take advantage of this option, please go visit the manual “README”.
- If you desire to **remove the crystallographic waters, solvent residues, hydrogens, heteroatoms or alternative atom locations** from your input PDB files and the Output model, please head to the “README” usage. Keep in mind that weird residues positions may derive into clashes, deriving in a bad model.

## Working in progress (not implemented)

We decided to implement an option where the user is able to submit the output model to a **minimization dynamic**, useful for retrieving a more accurate model where the side chains and possible clashes are relaxed. Keep in mind that this part is on development and will not work for all models. 
- PDB output model will be used to perform a minimization. 
- It cleans the output PDB file, adds hydrogens and prepares the files for an **Amber MD** (tleap). By now it only applies a general **ForceField (SB14FF)**, so it will only work for proteins made of only common aminoacids. 
- The complex is submitted to a simple **Minimization-Heating-Production system**. When it finishes, all the solvent and counterions are removed and we end obtaining a PDB file with the minimized complex.

The procedure was only tested successfully in the example provided for the task. It probably will not work for other models because we are not able to predict the Force Fields needed to explain all the atoms. The aim is just trying to improve the quality of the model predicted, please do not understand it as a real MD because most of the procedure may need mayor changes.

## Improvements

- **Being able to superpose high similarity chains. The code is already suitable for it.**
- **Analyze the output model. Make it easier for the user to decide if the model is good or wrong.**
- **Add a good MD / minimization system, applicable to more models**
- **Tolerant to weird PDB formats**
- **Improve test coverages**
- **Further and deeper analysis of clash cut off and maximum number of clashes per chain in order to englove a larger set of macro-complexes**