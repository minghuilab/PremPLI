# PremPLI
## About
<font size=4> 
  
PremPLI predicts the effects of single mutations on protein-ligand interactions by calculating the binding affinity changes quantitatively. It can be used for guiding the design of ligand-binding proteins, identifying and understanding disease driver mutations, and finding potential resistance mutations for different drugs. The 3D structure of a protein-ligand complex is required for performing the prediction.
  
</font>

## Scoring mutations with PremPLI
<font size=4> 

We recommend that most users who just want to obtain PremPLI predictions use [PremPLI website](https://lilab.jysw.suda.edu.cn/research/PremPLI/) to obtain scores.

</font>

## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/PremPLI/releases) on github.

</font>

## Installation

#### I. PREREQUISITES

<font size=4>
 
PremPLI requires the following software and packages.

1. PROVEAN

   This is available at the PROVEAN website.

   http://provean.jcvi.org/index.php/

2. NCBI BLAST 2.4.0

   This is available at the NCBI ftp site.

   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/

3. DSSP

   This is available at the DSSP website.

   https://swift.cmbi.umcn.nl/gv/dssp/

4. FoldX

   This is available at the FoldX website.

   http://foldxsuite.crg.eu/

5. VMD

   This is available at the VMD website.

   https://www.ks.uiuc.edu/Research/vmd/

6. UCSF Chimera

   This is available at the Chimera website.

   http://www.cgl.ucsf.edu/chimera/

7. Arpeggio

   This is available at the Arpeggio website.

   https://bitbucket.org/harryjubb/arpeggio/

8. XLOGP3

   This is available at the XLOGP3 website.

   http://sioc-ccbg.ac.cn/skins/ccbgwebsite/software/xlogp3/

9. Python packages: pandas, biopython, sklearn and rpy2

   To install these packages you can use the following command:
</font>

<font size=4>

	$ conda install -c conda-forge pandas
	$ conda install -c conda-forge biopython
	$ conda install -c conda-forge scikit-learn
	$ conda install -c r rpy2

</font> 

<font size=4>

10. R packages: randomForest and stringr

</font>

<font size=4>

	$ install.packages('randomForest')
	$ install.packages('stringr')

</font> 

#### II. INSTALLATION INSTRUCTIONS

<font size=4>

1. Download and/or install prerequisites described above.

2. Download and unpack the distribution:

</font>

<font size=4>

	$ wget https://github.com/minghuilab/PremPLI/archive/v1.0.tar.gz
	$ tar -zxvf v1.0.tar.gz

</font> 

<font size=4>

3. Change to the source directory:

</font>

<font size=4>

	$ cd PremPLI-1.0/v1.0

</font> 

<font size=4>

4. Change the path parameters in PremPLI.py:

</font>

<font size=4>

	workdir = Your working directory
	pathvmd = path for running VMD software  # /usr/local/bin/vmd
	pathmkdssp = path for running DSSP software  # /usr/local/bin/mkdssp
	pathpsiblast = path for running PSI-BLAST software  # /usr/local/bin/blast/psiblast
	pathblastdb = path for blastdb  # /usr/local/bin/blastdb/nr
	pathprovean = path for PROVEAN software  # /usr/bin/provean.sh
	patharpeggio = path for Arpeggio software  # /usr/local/bin/arpeggio/arpeggio.py
	pathxlogp3 = path for running XLOGP3 software  # /usr/local/bin/
	
</font>

&nbsp; &nbsp; The FoldX software needs to be installed in the working directory.

#### III. RUNNING PremPLI

<font size=4>

	$ python PremPLI.py -i 2021031002563493343426069

</font> 

## Platform

<font size=4>

PremPLI is only intended to run on *linux* operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 2.7 and R 3.4.0 (or higher) installed.

</font>
