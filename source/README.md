
PANGEA.HIV.sim Installation
======================

To install this R package from the command line: 

 1. cd to the directory that contains 'PANGEA.HIV.sim' 
 2. type 'R CMD build PANGEA.HIV.sim'  
 3. type 'R CMD INSTALL PANGEA.HIV.sim_1.0-0.tar.gz'


----------

PANGEA.HIV.sim in action
======================
The central function is PANGEA.HIV.sim.pipeline(). This function creates a UNIX batch file that must be started manually. The user syntax (in R) for example input files is as follows:

    library(PANGEA.HIV.sim)
    indir		<- system.file(package="PANGEA.HIV.sim", "misc")
    infile.ind	<- '140716_RUN001_IND.csv'
    infile.trm	<- '140716_RUN001_TRM.csv'
    outdir	<- '/Users/Oliver/git/HPTN071sim/tmp140909'
    dir.create(outdir, showWarnings=FALSE)
    file		<- PANGEA.HIV.sim.pipeline(indir, infile.ind, infile.trm, outdir)

See also '?PANGEA.HIV.sim.pipeline' for help.


----------
 
PANGEA.HIV.sim Documentation
======================

**PANGEAphy-01: help with the software**
> The package documentation can be loaded with `library(help=PANGEA.HIV.sim)`. The main function is `PANGEA.HIV.sim.pipeline`, which glues various subprograms together. Type `?PANGEA.HIV.sim.pipeline` for help. You can also type, for example,  `?prog.HPTN071.input.parser.v2` to see how the different subprograms can be executed. 

----------

**PANGEAphy-01: simulated genome**
> The following parts of the genome are simulated:

> - gag: p17 start to pol PROT start; length 1440 nucleotides. This is shorter than HXB2-K03455-gag due to several deletions in HIV-1C. The the simulated gag gene does not include the last 14 amino acids of p6, due to the overlap with pol.

> - pol: PROT start to Integrase end; length 2844 nucleotides. 

> - env: CDS signal peptide start to gp41 end; length 2523 nucleotides. This shorter than HXB2-K03455-env due to several deletions in HIV-1C.

----------

**PANGEAphy-01: genome simulation**
> Concatenated gag+pol+env HIV sequences are generated along each simulated transmission chain (from either the DSPS or HPTN071 epidemic simulator). Briefly, the simulation proceeds along the following steps:
> 
> - Simulate a viral genealogy through sampled and unsampled individuals of a transmission chain that allows for within-host evolution of the virus. Branch lengths are in units of calendar time. The VirusTreeSimulator JAVA program is used for this purpose (Matthew Hall).
> - Sample a starting sequence for each transmission chain at the simulated root time. A pool of ~250k starting sequences was generated with BEAST ancestral state reconstruction methods from the available 390 HIV-1C full genome sequences (see above).
> - Sample within-host evolutionary rates to compute the expected number of substitutions along each branch. Rates along transmission edges are dampened to account for the rate discrepancy between within-host and between-host HIV evolution.
> - Simulate partial HIV sequences for each gene and each codon position with the SeqGen program (Andrew Rambaut) and concatenate the partial sequences into the full genome.

----------

**PANGEAphy-01: further details**
> XX TODO XX
>