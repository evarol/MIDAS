MIDAS

  Section of Biomedical Image Analysis
  Department of Radiology
  University of Pennsylvania
  Richard Building
  3700 Hamilton Walk, 7th Floor
  Philadelphia, PA 19104

  Web:   https://www.med.upenn.edu/sbia/
  Email: sbia-software at uphs.upenn.edu

  Copyright (c) 2018 University of Pennsylvania. All rights reserved.
  See https://www.med.upenn.edu/sbia/software-agreement.html or COPYING file.

Author:
Erdem Varol
software@cbica.upenn.edu

===============
1. INTRODUCTION
===============
This software yield statistical maps for group comparisons or regressions. The statistical maps are based on regionally linear multivariate discriminative analysis. P-value maps are based on an analytic approximation of permutation testing.

===============
2. TESTING & INSTALLATION
===============

This software has been primarily implemented in MATLAB for Linux operating systems.

----------------
 Dependencies
----------------
- NIFTI Matlab toolbox (necessary files have been included in package)


----------------
 Installation
----------------

Midas can be run directly in a matlab environment without compilation.

OPTIONAL:

If user wants to run midas as a standalone executable, then it must be compiled as following (using the additionally obtained matlab compiler "mcc"):

Run the following command in a MATLAB environment:

   mcc -m midas.m

-----------------
 Test
-----------------
We provided a test sample in the test folder.

To test in matlab enviroment, use the command:

midas('-i','test.csv','-o','.','-r',15,'-p',200,'-c',0.1)

To test in command line using compiled executable, use the command:

midas -i test.csv -o . -r 15 -p 200 -c 0.1

This runs a MIDAS experiment which may take a few minutes. The test case contains a gray matter RAVENS maps of 46 subjects from the ADNI dataset. The output should yield statistical maps for diagnosis, age, and sex along with p-value maps in .nii.gz format. An accompanying MATLAB .mat file that stores these results is also output as MIDAS_results.mat in the output directory.

-----------------
 Test Verification
-----------------

Pre-computed statistical and p-value maps along with the .mat file have been included in directory "Pre_computed_test_results". The user may verify that their test results match the pre-computed results to confirm proper set-up.

==========
3. USAGE
==========


I. Input images

MIDAS requires input images in nifti gz format (.nii.gz file extension). The images should be registered to a common template space and share the same dimensionality.


II. Running "MIDAS":

Here is a brief introduction to running MIDAS. For a complete list of parameters, see --help option.

To run this software, you will need an input csv file, with the following mandatory fields in the following column order:
(Column 1) ID: ID for subject
(Column 2) filepath: path for input image for subject
(Column 3 and on) covariate_1,covariate_2,...: covariates to obtain statistical maps for (Needs to be numerical)

NOTE: column header names can be arbitrary, only the order matters.

An example input csv file looks as following:
    
ID,        filepath,            DIAG,    AGE,    SEX
subject_1,    ./data/subject_1.nii.gz,    1,    79.3,    -1
subject_2,    ./data/subject_2.nii.gz,    1,    71.4,    1
subject_3,    ./data/subject_3.nii.gz,    1,    82.7,    -1

    
If you install the package successfully, there will be two ways of running MIDAS:

1. Running MIDAS in a matlab enviroment, a simple example:

        midas('-i','test.csv','-o','.','-r',15,'-p',200,'-c',0.1)

2. Running matlab compiled MIDAS executables in the command line, a simple example:

    midas -i test.csv -o . -r 15 -p 200 -c 0.1


The software returns:
1. statistical maps for provided covariates along with p-value maps in .nii.gz format.
2. Same results in a matlab format in MIDAS_results.mat 


===========
4. REFERENCE
===========

If you find this software useful, please cite:

Varol, Erdem, Aristeidis Sotiras, Christos Davatzikos. "MIDAS: regionally linear multivariate discriminative statistical mapping." NeuroImage (2018)

===========
5. LICENSING
===========

  See https://www.med.upenn.edu/sbia/software-agreement.html or COPYING.txt file.

