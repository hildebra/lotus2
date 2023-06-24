# LotuS2
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lotus2/badges/downloads.svg)](https://anaconda.org/bioconda/lotus2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lotus2/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/lotus2)

LotuS2 is a amplicon sequencing pipeline, that is programmed to be lightweight, easy to use, fast without comprimising the quality of reconstructing microbial communitites. It supports 16S, 18S and ITS amplicons, and support for other amplicon targets is available. Currently five different sequence clustering algorithms (DADA2, uparse, unoise3, cd-hit, vsearch) are supported as well as multiple options to assigning taxonomic annotations. LotuS2 output can be imported directly into R or as text file into other programs.
Full documentation on http://lotus2.earlham.ac.uk/

### REQUIREMENTS
LotuS2 requires a perl installation and sdm requires a fairly recent C++ compiler (like gcc or clang) that supports C++11.
Lambda currently only works under linux, the option to use lambda is not available on mac os. Instead, Blast, vsearch or usearch can be used.

### INSTALL LotuS2
LotuS2 can be installed via conda https://anaconda.org/bioconda/lotus2 
```{sh}
conda install -c bioconda lotus2
```
If there should be problems with the conda solver, try:
```{sh}
conda create -c conda-forge -c bioconda --strict-channel-priority -n lotus2 lotus2
conda activate lotus2
```

Alternatively, often the github contains pre-release versions and can be installed via:
```{sh}
git clone https://github.com/hildebra/lotus2.git
perl autoInstall.pl
```

All required software will be downloaded and installed in this directory.

If you want to install software and databases to other locations, follow installation instructions on: 
http://lotus2.earlham.ac.uk/main.php?site=documentation

LotuS2 packs a static compiled linux sdm binary, this should be useable out of the box. However, to manually compile **sdm**  (the autoinstaller can also do this) go to the lotus subdirectory *sdm_src* and run 
**make** to compile the sdm binary. Next copy the binary into the lotus directory using **cp**. For using MACs and LotuS2, this would be a requirement as the static binary will only work Linux systems.
```{sh}
cd sdm_src
make
cp sdm ../sdm
```

###  UPDATE LotuS2
If you installed LotuS2 via "git clone", you can get the latest release via "git pull".
LotuS2 has a built in mechanism to upgrade LotuS2, so that properitary programs & databases (once installed) don't have to be downloaded again. To use this feature a) install LotuS2 for the first time using the autoinstaller. 
Once you know or want to check for new updates, simply excecute the autoinstaller again and you will be prompted if LotuS2 should be updated. In case no new updates are available, the autoinstaller will exit without making changes, so this function can be used frequently. New updates will be avaialble from the LotuS2 webpage (http://lotus2.earlham.ac.uk/) or from github (https://github.com/hildebra/lotus2).

## EXAMPLES & CONFIGURATION
To test your installation, run a minimal example using test files distributed with LotuS2
*Note: files can be found relative to the lotus2 executable. If you used a conda installation, use "which lotus2", the databases and examples will be present in the shared conda folder*
:
```{sh}
./lotus2 -i Example/ -m Example/miSeqMap.sm.txt -o myTestRun
```

LotuS2 will try to choose default options. However, note that you should have seen a warning that no PCR primers were provided.

In the next example, we will explicitly configure the read filtering by providing sdm_miSeq2.txt, explicitly defining this to be 16S data from an illumina miSeq machine, to remove PCR primers used in this experiment, to use DADA2 instead of UPARSE clustering algorithm, to use alginments of ASVs against SILVA reference database instead of RDPclassifier taxonomic annotations:
```{sh}
./lotus2 -i Example/ -m Example/miSeqMap.sm.txt -o myTestRun2 -s configs/sdm_miSeq2.txt -p miSeq -amplicon_type SSU -forwardPrimer GTGYCAGCMGCCGCGGTAA -reversePrimer GGACTACNVGGGTWTCTAAT -CL dada2 -refDB SLV -taxAligner lambda
```
Building the lambda formatted SILVA reference database will take a long time the first time you run this. Please ensure that during installation you selected that the SILVA database will be installed (otherwise this example will not work).

There are >60 flags with which you can further customize each LotuS2 run, but we try to optimize LotuS2 to work pretty well with just default options. Please run ./lotus2 to see these options.

### Coustom Reference Database
To use your own reference database for LotuS2, you can employ the flags -refDB and -tax4refDB to supply a fasta formated reference database and tab-delimited taxonomy file, respectively. The format of these is the same as in the database already installed with the LotuS2 autoinstall, e.g. have a look at DB/SLV_138_SSU.fasta and DB/SLV_138_LSU.tax for examples. In the \*.tax file, the levels are fixed to 7 levels (kingdom, phylum, class, order, family, genus, species). These are demarked by tags k__; p__ ; etc and separated by a ";" character. In case tax information is missing, use "?" to inset this information, for example:

FJ588878	k__Eukaryota; p__Phragmoplastophyta; c__?; o__?; f__?; g__?; s__Osyris wightiana

Let's simulate using SILVA138 (SLV) as custom database now, with vsearch as search algorithm and uparse OTU clusering, using the following example (this might take some time )
*Note: databases are installed relative to the lotus2 excecutable (see above)*:
```{sh}
./lotus2 -tax4refDB DB/SLV_138_SSU.tax -refDB DB/SLV_138_SSU.fasta -i Example/ -m Example/miSeqMap.sm.txt -o myTestRun3 -forwardPrimer GTGYCAGCMGCCGCGGTAA -reversePrimer GGACTACNVGGGTWTCTAAT -CL uparse -taxAligner vsearch 
```

We can make this even more complicated, by having both HitDB and SLV as complimentary databases searched, using either
```{sh}
./lotus2 -tax4refDB DB/SLV_138_SSU.tax,DB/HITdb/HITdb_taxonomy.txt -refDB DB/SLV_138_SSU.fasta,DB/HITdb/HITdb_sequences.fna -i Example/ -m Example/miSeqMap.sm.txt -o myTestRun3 -forwardPrimer GTGYCAGCMGCCGCGGTAA -reversePrimer GGACTACNVGGGTWTCTAAT -CL uparse -taxAligner vsearch 
```
or (this is a shortcut, possible because GG and SLV are in-built):

```{sh}
./lotus2 -refDB SLV,HITdb -i Example/ -m Example/miSeqMap.sm.txt -o myTestRun3 -forwardPrimer GTGYCAGCMGCCGCGGTAA -reversePrimer GGACTACNVGGGTWTCTAAT -CL uparse -refDB SLV -taxAligner vsearch 
```

Note that "-refDB GG,SLV" and "-refDB SLV,GG" would likely give a different result, as GG is the primary annotation source in the first, SLV is primary in the second case.

###  PacBio CCS amplicon sequence processing with LotuS2

When using PacBio CCS (HiFi) amplicons, set the ** -p PacBio ** flag first, this will be set important default options that in our hands work relatively well for PacBio. As clustering algortithm, CD-HIT might be a good choice because it's independent of read quality and read length differences, but in the end any of the clusterings will work.

However, you might want to adopt the sdm read quality filtering options further. Have a look at either **configs/sdm_PacBio_ITS.txt** or **configs/sdm_PacBio_LSSU.txt**, these are the default configuration for ITS or SSU/LSU PacBio amplicons, respectively. Parameters can be freely modified, but of specific importance are:
```
minSeqLength	700
maxSeqLength	2000
TruncateSequenceLength	-1
```

The first two describe the length distribution you expect of your amplicon (e.g. full length 16S would be ~1500 bp). Narrowing this further down can help avoid false positives, but beeing too narrow will also remove some natural biological variation. The range 700-1500bp in the LSSU file is chosen to be inclusive, but if you know you have full length 16S, 1300-1600 might be a better choice. For 18S amplicons, 1600-2000 might be a better choice. It's important to note that **TruncateSequenceLength -1** deactivates sequence truncation. This option is very important for illumina amplicons, because truncating amplicons is extremely important for DADA2 and UPARSE and UNOISE3 clusterings (please refer to [R Edgars excellent UPARSE paper](https://www.nature.com/articles/nmeth.2604)). However, since PacBio reads do not have a length dependent decay in quality scores and under the assumption that the complete amplicon is sequenced, truncating sequences for the clustering step is not necessary.

Also note the option 
```
ExtensivePrimerChecks	T
RejectSeqWithoutFwdPrim	T
RejectSeqWithoutRevPrim	T
```
This options are important to remove faulty amplicons that can result from longer PCRs and suboptimal primers, PCR and sample conditions and faulty CCS derivations. This doesn't need to be the case for your experiment, but sometimes they occur and hence we decided to use a very strict quality filtering for PacBio data. Note that you need to pass your amplicon primers to LotuS2, otherwise all reads will be removed during quality filtering (**RejectSeqWithoutFwdPrim	T** option).

So to summarize, maybe create your copy of the default **configs/sdm_PacBio_ITS.txt** or **configs/sdm_PacBio_LSSU.txt**, modify it to your needs, and run LotuS2 with a command similar to:

```{sh}
./lotus2 -i PacBioDir/ -p PacBio -id 0.97 -CL cdhit -s configs/sdm_PacBio_my_copy.txt -refDB SLV  -m PacBioDir/my_PacBio.map -o /my/PacBio/LotuS2 -forwardPrimer XYZ -reversePrimer XYZ ...
```


## Publications related to LotuS2
LotuS2: https://www.biorxiv.org/content/10.1101/2021.12.24.474111v1

offtarget removal: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01012-1

LotuS: http://www.microbiomejournal.com/content/2/1/30

## Acknowledgements 
LotuS2 was developed at Quadram Institute Bioscience (QIB) & Earlham Institaute (EI), Norwich, UK. Various members of the Hildebrand group contributed to the pipelines (see publication Özkurt et al., 2022).
(c) Falk.Hildebrand {at} gmail.com

**Please cite LotuS2 with:**

**Pipeline** - Özkurt E, Fritscher J, et al. (2022) LotuS2: An ultrafast and highly accurate tool for amplicon sequencing analysis. Microbiome 10:176  doi:10.1186/s40168-022-01365-1.

**offtarget removal** - Bedarf JR, Beraza N, Khazneh H, Özkurt E, et al (2021) Much ado about nothing? Off-target amplification can lead to false-positive bacterial brain microbiome detection in healthy and Parkinson’s disease individuals. Microbiome ;9:75.

We would like to acknowledge the following proprietary software, that is used in LotuS2. Please acknowledge these if your LotuS2 run was using them (listed in LotuSLogS/citations.txt for each LotuS2 run):

* **DADA2** - Callahan, B., McMurdie, P., Rosen, M. et al. 2016. DADA2: High resolution sample inference from Illumina amplicon data. Nat Methods, 13. 581–583 (2016).

* **UPARSE** - Edgar RC. 2013. UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nat Methods, 10, 996–998 (2013).

* **VSEARCH** - Rognes T, Flouri T, Nichols B, Quince C, Mahé F.2016. VSEARCH: a versatile open source tool for metagenomics PeerJ. vol. 4 e2584.

* **swarm** - Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. 2014. Swarm: robust and fast clustering method for amplicon-based studies. PeerJ 2: e593.

* **CD-HIT** - Fu L, Niu B, Zhu Z, Wu S, Li W. 2012. CD-HIT: Accelerated for clustering the next-generation sequencing data. Bioinformatics 28: 3150–3152.

* **uchime** - Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. 2011. UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27: 2194–200.

* **RDP classifier** - Wang Q, Garrity GM, Tiedje JM, Cole JR. 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Env Microbiol 73: 5261–5267; DOI: 10.1128/AEM.00062-07.

* **lambda aligner** - Hauswedell H, Singer J, Reinert K. 2014. Lambda: the local aligner for massive biological data. Bioinformatics 30: i349–i355. 

* **Blast+** - Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215: 403–10.

* **Clustal Omega** - Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, et al. 2011. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol 7: 539.

* **MAFFT** - Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013;30:772–80.

* **fasttree2** - Price MN, Dehal PS, Arkin AP. 2010. FastTree 2--approximately maximum-likelihood trees for large alignments. ed. A.F.Y. Poon. PLoS One 5: e9490.

* **IQ-TREE 2** - Nguyen L-T, Schmidt HA, von Haeseler A, Minh BQ. IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Mol Biol Evol. 2015;32:268–74.

## Databases

* **Greengenes** - McDonald D, Price MN, Goodrich J, Nawrocki EP, DeSantis TZ, Probst A, Andersen GL, Knight R, Hugenholtz P. 2012. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. ISME J 6: 610–8.

* **SILVA** - Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glockner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucleic Acid Res. 42:D643-D648 

* **HITdb** - Ritari J, Salojärvi J, Lahti L & de Vos WM. Improved taxonomic assignment of human intestinal 16S rRNA sequences by a dedicated reference database. BMC Genomics. 2015 Dec 12;16(1):1056.

* **beetax** - Jones, JC, Fruciano, C, Hildebrand, F, et al. Gut microbiota composition is associated with environmental landscape in honey bees. Ecol Evol. 2018; 8: 441– 451.

* **PR2** - Guillou L, Bachar D, Audic S, et al. The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. Nucleic Acids Res. 2013;41(Database issue):D597-D604

## ITS specific

* **UNITE ITS chimera DB** - Nilsson et al. 2015. A comprehensive, automatically updated fungal ITS sequence dataset for reference-based chimera control in environmental sequencing efforts. Microbes and Environments 

* **UNITE ITS taxonomical refDB** - Koljalg, Urmas, et al. "Towards a unified paradigm for sequence-based identification of fungi." Molecular Ecology 22.21 (2013): 5271-5277.

* **ITSx** - Bengtsson‐Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A., De Wit, P., Sánchez‐García, M., Ebersberger, I., de Sousa, F., Amend, A., Jumpponen, A., Unterseher, M., Kristiansson, E., Abarenkov, K., Bertrand, Y.J.K., Sanli, K., Eriksson, K.M., Vik, U., Veldre, V. and Nilsson, R.H. 2013. Improved software detection and extraction of ITS1 and ITS 2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods in Ecology and Evolution


## Mathematical models

* Puente-Sánchez 2016 - A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies. Nucleic Acids Res. 2016;44(4):e40. 

## C++ libraries

* **gzip libraries** - (gzstream.h) https://gist.github.com/piti118/1508048 and zlib library (http://www.zlib.net/)
* **robin_hood hash map libraries** - (robinhood.h) https://github.com/martinus/robin-hood-hashing
