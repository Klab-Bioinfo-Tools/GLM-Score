GLM-Score
========

A set of empirical scoring functions for predicting Receptor-Ligand binding affinities (in pKd units): protein-DNA/RNA, protein-small molecule, and protein-protein complexes.
Each scoring function was built using generalized linear models (GLMs) applied to specific and curated data sets.

#Requirements

- GCC (http://gcc.gnu.org).
- R with GLM package installed (http://cran.r-project.org).

#Download

(LINUX):

    wget https://github.com/Klab-Bioinfo-Tools/GLM-Score/archive/master.zip -O GLM-Score.zip

(MAC):

    curl https://github.com/Klab-Bioinfo-Tools/GLM-Score/archive/master.zip -o GLM-Score.zip

#Extract the files:

    unzip GLM-Score.zip
    
#Install
    
    cd GLM-Score-master
    g++ GLM-Score.cpp -o GLM-Score
    
#Run
    
This is an example of how to run GLM-Score:

    ./GLM-Score protein.pdb ligand.pdb ligand.mol2 DNA
    
There are 4 input parameters for GLM-Score:

1. Protein or receptor structure in PDB format (www.pdb.org);
2. Ligand structure in PDB format (the ligand can also be another protein);
3. The same ligand structure, but in MOL2 format. There are many stand alone and online tools available for converting PDB files into MOL2 format. Some examples are OpenBabel (http://openbabel.org) and the online molecular formats converter (http://www.webqc.org/molecularformatsconverter.php).
4. The type of ligand molecule (DNA, protein, or small_molecule). Dependending on the type of molecule you provide, the program will choose the best statistical model for predicting the binding affinity of your binding complex.


    

