rPANGEAHIV Installation
======================

To install this R package from the command line: 

 1. cd to the directory that contains 'rPANGEAHIVsim' 
 2. type 'R CMD build rPANGEAHIVsim'  
 3. type 'R CMD INSTALL rPANGEAHIVsim_1.0-0.tar.gz'

----------
 
rPANGEAHIV Help files
======================
Then run R and load the package with 'library(rPANGEAHIVsim)'. The package content that is currently documented is shown with 'library(help=rPANGEAHIVsim)'. You can type e.g.  '?prog.HPTN071.input.parser.v2' to see how the calling syntax for different functions. The help file for the virus tree simulator is '?cmd.VirusTreeSimulator'.

----------

rPANGEAHIV Simulations
======================
The central function is rPANGEAHIVsim.pipeline(), type '?rPANGEAHIVsim.pipeline' for help. This function creates a UNIX batch file.

> **PANGEAphy-01: simulated genome**
> The following parts of the genome are simulated:

> - gag: p17 start to pol PROT start; length 1440 nucleotides. This is shorter than HXB2-K03455-gag due to several deletions in HIV-1C. The the simulated gag gene does not include the last 14 amino acids of p6, due to the overlap with pol.

> - pol: PROT start to Integrase end; length 2844 nucleotides. 

> - env: CDS signal peptide start to gp41 end; length 2523 nucleotides. This shorter than HXB2-K03455-env due to several deletions in HIV-1C.