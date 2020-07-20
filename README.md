               +?+                                        =+???+                
               I7I                   +II                ,?77++I7                
               I7I         ~==~      ?77     ,      ,   =77                     
               I7I       =I777777  =I77777  =77    ?77  =77I,                   
               I7I      +77    I7I   ?77    =77    ?77   +7777I=                
               I7I      I7,    =77   ?77    =77    ?77      I777I               
               I7I      I7=    =77   ?77    =77    ?77        ?77               
               I7I      ?77    ?7I   ?77    =I7=  =I77  =,    +77               
               I777777I  I77III7I    +777I   I77777I77  ?77II777                

----------------------------------
# LotuS
----------------------------------
(c) Falk.Hildebrand {at} gmail.com

http://lotus2.earlham.ac.uk/

http://www.microbiomejournal.com/content/2/1/30


### REQUIREMENTS
LotuS requires a perl installation and sdm requires a fairly recent C++ compiler (like gcc or clang) that supports C++11.
Lambda currently only works under linux, the option to use lambda is not available on mac os :( Instead, Blast can be installed.

### INSTALL LotuS
To install LotuS and required properitary software execute
```{sh}
perl autoInstall.pl
```
All required software will be downloaded and installed in this directory.
A video tutorial is available here: 
http://lotus2.earlham.ac.uk/documentation.html

If you want to install software and databases to other locations, follow installation instructions on: 
http://lotus2.earlham.ac.uk/documentation.html

Compiling **sdm** manually (the autoinstaller is alternatively doing this):
go to the lotus subdirectory *sdm_src* and run 
**make** to compile the sdm binary. Next copy the binary into the lotus directory using **cp**.
```{sh}
cd sdm_src
make
cp sdm ../sdm
```

###  UPDATE LotuS
LotuS has a built in mechanism to upgrade LotuS, so that properitary programs & databases (once installed) don't have to be downloaded again. To use this feature a) install LotuS for the first time using the autoinstaller. 
Once you know or want to check for new updates, simply excecute the autoinstaller again and you will be prompted if LotuS should be updated. In case no new updates are available, the autoinstaller will exit without making changes, so this function can be used frequently. New updates will always be announced on LotuS webpage (http://lotus2.earlham.ac.uk/).

### EXAMPLES
To test your installation, run the example test set:
```{sh}
./lotus.pl -i Example/ -m Example/miSeqMap.sm.txt -s sdm_miSeq.txt -p miSeq -o myTestRun
```


**Please cite LotuS with:**

Hildebrand F, Tadeo RY, Voigt AY, Bork P, Raes J. 2014. LotuS: an efficient and user-friendly OTU processing pipeline. Microbiome 2: 30. 


## Acknowledgements of the softwares

I would like to acknowledge the following software, that is used in LotuS, please also acknowledge these if you use them:

* **DADA2** - Callahan, B., McMurdie, P., Rosen, M. et al. 2016. DADA2: High resolution sample inference from Illumina amplicon data. Nat Methods, 13. 581–583 (2016).

* **UPARSE** - Edgar RC. 2013. UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nat Methods, 10, 996–998 (2013).

* **VSEARCH** - Rognes T, Flouri T, Nichols B, Quince C, Mahé F.2016. VSEARCH: a versatile open source tool for metagenomics PeerJ. vol. 4 e2584.

* **swarm** - Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. 2014. Swarm: robust and fast clustering method for amplicon-based studies. PeerJ 2: e593.

* **CD-HIT** - Fu L, Niu B, Zhu Z, Wu S, Li W. 2012. CD-HIT: Accelerated for clustering the next-generation sequencing data. Bioinformatics 28: 3150–3152.

* **DNACLUST** - Ghodsi, M., Liu, B., & Pop, M. (2011). DNACLUST: accurate and efficient clustering of phylogenetic marker genes. BMC Bioinformatics, 12, 271.

* **uchime** - Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. 2011. UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27: 2194–200.

* **flash** - Magoc T, Salzberg SL. 2011. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27: 2957–63.

* **RDP classifier** - Wang Q, Garrity GM, Tiedje JM, Cole JR. 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Env Microbiol 73: 5261–5267; DOI: 10.1128/AEM.00062-07.

* **lambda aligner** - Hauswedell H, Singer J, Reinert K. 2014. Lambda: the local aligner for massive biological data. Bioinformatics 30: i349–i355. 

* **Blast+** - Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215: 403–10.

* **Clustal Omega** - Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, et al. 2011. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol 7: 539.

* **fasttree2** - Price MN, Dehal PS, Arkin AP. 2010. FastTree 2--approximately maximum-likelihood trees for large alignments. ed. A.F.Y. Poon. PLoS One 5: e9490.

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

