#!/usr/bin/perl
# lOTUs2 - less OTU scripts
# Copyright (C) 2020 Falk Hildebrand

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

# contact
# ------
# Falk.Hildebrand [at] gmail.com
#

use strict;
use warnings;
use Getopt::Long qw( GetOptions );


#use threads;
use 5.012;
use FindBin qw($RealBin);
#my $LWPsimple = eval {require LWP::Simple;LWP::Simple->import();1;};

sub usage;
sub help;
sub frame;sub printL;
sub announce_options;
sub buildOTUs;
sub mergeUCs;sub delineateUCs;sub cutUCstring;sub uniq;
sub prepLtsOptions;

sub makeAbundTable2;
sub readMap;sub writeMap; sub addSeqRun;
sub assignTaxOnly; sub runRDP;

sub writeUTAXhiera;
sub doDBblasting;
sub numberOTUs;
sub getTaxForOTUfromRefBlast;
sub get16Sstrand;
sub splitBlastTax;
sub extractTaxoForRefs;
sub calcHighTax;
sub buildTree;
sub combine;
sub contamination_rem;
sub readOTUmat;
sub runPhyloObj;
sub readFastaHd;
sub replaceFastaHd;
sub splitFastas;
sub extractFastas;
sub writeFasta;
sub readFasta;
sub revComplFasta;
sub reverse_complement_IUPAC;
sub newUCSizes;
sub readLinkRefFasta;
sub readTaxIn;
sub biomFmt;
sub finWarn;
sub printWarnings;    #1 to add warnings, the other to print these
sub mergeRds;
sub clean_otu_mat;
sub systemL;
sub swarmClust;
sub swarm4us_size;
sub dnaClust2UC;
sub chimera_ref_rem; sub chimera_denovo;
sub checkLtsVer;
sub ITSxOTUs; sub vfxtractor;
sub lulu;
sub checkXtalk;
sub annotateFaProTax;
sub systemW;

#keep track of time
my $start = time;
my $duration;
$| = 1; #don't buffer ostream

#print qx/ps -o args $$/."\n";
my $cmdCall = qx/ps -o args $$/;

# --------------------
# Progams Pathways  -- get info from lotus.cfg
my $LCABin      = "";
my $sdmBin      = "./sdm";
my $usBin       = ""; #absolute path to usearch binary
my $dada2Scr    = ""; #absolute path to dada2 pipeline script
my $LULUscript  = ""; #script for LULU OTU matrix cleanup
my $phyloLnk    = ""; #R helper script to create phyloSeq objects
my $Rscript     = ""; #absolute path to Rscript -> also serves as check that Rscript exists on cluster..
my $defRscriptOpt  = " --vanilla ";
my $swarmBin    = "";
my $VSBin       = "";
my $VSBinOri    = "";
my $VSused      = 1;
my $cdhitBin    = "";
my $dnaclustBin = "";
my $mini2Bin    = "minimap2";
my $mjar = ""; #RDP multiclassifier java archieve (e.g. /YY/MultiClassifier.jar)
my $rdpjar = ""
  ; #alternatively leave this "" and set environmental variabel to RDP_JAR_PATH, as described in RDP documentary
my $blastBin     = "";  #Blastn+ binary (e.g. /YY/ncbi-blast-2.2.XX+/bin/blastn)
my $mkBldbBin    = "";
my $lambdaBin    = "";  #lambda ref DB search
my $lambdaIdxBin = "";
my $clustaloBin  = ""; #clustalo multiple alignment binary (e.g. /YY/clustaloBin-1.2.0-OS-x86_64)
my $mafftBin     = ""; #different aligner than clustalo, better performance
my $fasttreeBin  = ""; #FastTree binary - for speed improvements /YY/FastTReeMP recommended -- requires clastalo
my $iqTreeBin    = ""; #iqTree, better choice instead of fasttree
my $flashBin     = "";    #flash for merging paired reads
my $itsxBin      = "";    #identifies ITS regions
my $hmmsVX = "" ; my $vxtrBin = ""; #vxtractor

my $hmmsrchBin   = "";    #hmm required for itsx

# --------------------
#   databases
my $UCHIME_REFDB = ""; #reference database for ref chimera filtering - suggested is gold.fa or LSU93
my @TAX_REFDB = (); #greengenes or SILVA fasta ref database - requires $blastBin
my @TAX_RANKS = ();    #greengenes or SILVA taxonomic assignments for ref fasta database
my $CONT_REFDB_PHIX = "";           #ref Genome for PhiX contamination checks
my @refDBname       = ();
my $refDBwanted     = "";
my $refDBwantedTaxo = "";           #only used with custom ref DB files
my $ampliconType    = "SSU";
my $organism        = "bacteria";
my $FaProTax        = undef; #faprotax -> functional annotations from taxonomy~~

# --------------------

# --------------------
#general lOTUs parameters
my $selfID = "LotuS 2.01"; #release candidate: 2.0
my $citations = "$selfID: Hildebrand F, Tadeo RY, Voigt AY, Bork P, Raes J. 2014. LotuS: an efficient and user-friendly OTU processing pipeline. Microbiome 2: 30.\n";
my $noChimChk = 0; #deactivate all chimera checks 1=no nothing, 2=no denovo, 3=no ref; default = 0
my $mainLogFile = "";
my $osname      = $^O;
my $verbosity = 1;
my $versionOut=0;

#LotuS file related
my $rmOutDir = 0;
my $lotusCfg  = "$RealBin/lOTUs.cfg";

my $ClusterPipe_pre = "?"; #default is unknown..
my $ClusterPipe     = 1; #use UPARSE (1) or otupipe(0) or SWARM (2) or cd-hit(3), dnaclust (4), micca (5)
my $clusteringNameStr = "otupipe";
my $lotus_tempDir    = "";
my $sdmOpt           = "";    #default sdm (simple demultiplexer) options
my $damagedFQ        = 0;
my $BlastCores       = 12;    #number of cores to use for BLASTn
my $clQue            = "";    #"all.q";#"highmem";
my ${OTU_prefix}     = "OTU";    # Prefix for OTU label, default is OTU_ giving OTU_1, OTU_2 ...
my $chimera_prefix   = "CHIMERA_";# Prefix for chimera label, default is CHIM_ giving CHIM_1, OTU_2 ...
my $sep_smplID       = "___";    #separator of smpID & ori fasta id
my $extendedLogs     = 1; #write chimeric OTUs, exact blast/RDP hits to extra dir
my $keepTmpFiles     = 0; #more of a debugging option
my $checkForUpdates  = 1; #check online if a new lotus version is avaialble
my $maxReadOverlap   = 250;          #flash parameter
my $mergePreCluster  = 0;  #merge reads in sdm already (derep, demulti, cleaned reads)?
my $takeNonMerge     = 1; #also use single end that don't merge
my $maxHitOnly       = 0;
my $greengAnno       = 0;            #if 1, annotate OTUs to best greengenes hit
my $pseudoRefOTU     = 0;            #replace OTU ids with best hit (LCA)
my $numInput         = 2;
my $saveDemulti      = 0;            #save a copy of demultiplexed reads??
my $check_map        = "";
my $create_map       = "";
my $curSdmV          = 1.46;
my $platform         = "miSeq";        #454, miSeq, hiSeq, PacBio
my $keepUnclassified = 1;
my $keepOfftargets   = 0;
my $doITSx           = 1;            #run itsx in its mode?
my $doVXTR           = "0";   #run v xtractor?
my $doLULU           = 0;   #run LULU post matrix filter?
my $ITSpartial       = 0;            #itsx --partial parameter
my $finalWarnings    = "";
my $remFromEnd       = ""; #fix for strange behavior of flash, where overlaps of too short amplicons can include rev primer / adaptor
my $doPhiX			 = 1;
my $dada2Seed        = 0; #seed for dada2 to produce reproducible results
my $buildPhylo       = 1; #iqtree(2), fastree (1)


#my $combineSamples = 0; #controls if samples are combined
my $chimCnt = 0;    #should chimeric OTU counts be split up among parents?

#flow controls
my $onlyTaxRedo = 0;    #assumes that outdir already contains finished lotus run
my $TaxOnly     = "0"; #will override most functionality and skip directly ahead to tax assignment part
						#can also be used to assign tax directly to fa given as argument, e.g. -taxOnly /xx/test.fa

# --------------------
#similarity taxo search options
my $RDPCONF     = 0.8;    #RDP confidence threshhold
my $utaxConf    = 0.8;    #UTAX confidence threshhold
my $LCAfraction = 0.8;    #number of matching taxa to accept taxa level
my $LCAcover = 0.5;
my $minBit      = 120;
my $minEval     = 1e-14;    #blast filtering; should be optimized to different platforms
my @idThr = ( 97, 93, 93, 91, 88, 78, 0 );
my $lengthTolerance =
  0.85;     #length of hits to still consider valid, compared to longest match
my $linesPerFile = 4000000;
my $otuRefDB =  "denovo"; #OTU building strategy: "denovo" (default), "ref_closed", "ref_open"
my $doBlasting = -1;       #$doBlasting 2:lambda, 1:blast, 0:RDP, -1: ini value, 3:utax 4: vsearch 5: usearch
my $doRDPing          = -1;    #0: not, 1:do RDP
my $doBlasting_pre    = -1;
my $custContamCheckDB = "";
my $lowMemLambI       = 0;

# --------------------
#uparse / otupipe / cdhit / swarm options
my $useVsearch = 0; # 1(vsearch) or 0(usearch) for chim check are being used, searching, matching tasks
#my $preferVsearch   = 0;      #0=use usearch always; 1=use vesarch always
my $chimera_absskew = 2;      # Abundance skew for de novo chimera filtering # 2
my $id_OTU          = .97;    # Id threshold for OTU clustering  #Default: .97
my $swarmClus_d     = 1;
my $id_OTU_noise =
  .99;   # Id threshold for error correction (noise removal) round #Default: .99
my $uthreads                = 1;
my $dereplicate_minsize_def = 2
  ; # Discard clusters < dereplicate_minsize in error correction round #Default: 2 for UPARSE, 4 for otupipe
my $dereplicate_minsize = -1;
my $doXtalk             = -1;    #check for cross talk in OTU tables
my $usearchVer          = 7; my $usearchsubV = 0;
my $REFflag = 0;#reference based OTU clustering requested?

#flash control
my $flashCustom = "";

#### DEPRECEATED ###
my $truncLfwd    = 130;    #250 for fwd illumina miSeq & 454; 90 for hiSeq fwd
my $truncLrev    = 200;    #200 for rev illumina miSeq; 90 for hiSeq rev
my $truncQual    = 25;     #15 for illumina miSeq, 25 for 454
my $UPARSEfilter = 0
  ; #use additional quality filter from uparse (1) or use sdm only for filtering (0) # Default: 0
#### DEPRECEATED ###

my $input = "";
my $outdir = "";
my $inq         = "";
my $barcodefile = "";
my $mapFile;
my $exec       = 0;    #my $highmem = 0;
my $sdmDerepDo = 2;
if ( !@ARGV ) {
    usage();
    exit 1;
}

GetOptions(
    "help|?"                => \&help,
    "i=s"                   => \$input,
    "o=s"                   => \$outdir,
    "barcode|MID=s"         => \$barcodefile,
    "m|map=s"               => \$mapFile,
    "taxOnly|TaxOnly=s"     => \$TaxOnly,
    "redoTaxOnly=i"         => \$onlyTaxRedo,
    "check_map=s"           => \$check_map,
	"create_map=s"          => \$create_map,
    "q|qual=s"              => \$inq,
    "s|sdmopt=s"            => \$sdmOpt,
    "tmp|tmpDir=s"          => \$lotus_tempDir,
    "c|config=s"            => \$lotusCfg,
    "exe|executionMode=i"   => \$exec,
    "keepTmpFiles=i"        => \$keepTmpFiles,
    "extendedLogs=i"        => \$extendedLogs,
    "cl|CL|clustering|UP|UPARSE=s" => \$ClusterPipe_pre,
    "t|thr|threads=i"       => \$uthreads,
	"v|version"             => \$versionOut,
	"verbosity=i"           => \$verbosity,
    "highmem=i"             => \$sdmDerepDo,
    "useBestBlastHitOnly=i" => \$maxHitOnly,
    "pseudoRefOTUcalling=i" => \$pseudoRefOTU,
    "greengenesSpecies=i"   => \$greengAnno,
    "id=f"                  => \$id_OTU,
    "xtalk=i"               => \$doXtalk,
    "saveDemultiplex=i"     => \$saveDemulti,        #1=yes, 0=not, 2=yes,unfiltered
    "rdp_thr=f"             => \$RDPCONF,
    "ITSx|itsextraction=i"  => \$doITSx,
	"VXtr=s"                => \$doVXTR,
	"lulu=i"                => \$doLULU,
    "itsx_partial=i"        => \$ITSpartial,
    "utax_thr=f"            => \$utaxConf,
    "LCA_frac=f"            => \$LCAfraction,
	"LCA_cover=f"           => \$LCAcover,
    "chim_skew=f"           => \$chimera_absskew,		#miSeq,hiSeq,pacbio
    "p|platform=s"          => \$platform,
    "tolerateCorruptFq=i"   => \$damagedFQ,
    "keepUnclassified=i"    => \$keepUnclassified,
	"keepOfftargets=i"      => \$keepOfftargets,
    "derepMin=s"            => \$dereplicate_minsize,
    "doBlast|simBasedTaxo|taxAligner=s"=> \$doBlasting_pre,
    "refDB=s"               => \$refDBwanted,
    "tax4refDB=s"           => \$refDBwantedTaxo,
    "amplicon_type=s"       => \$ampliconType,       #SSU LSU ITS ITS1 ITS2
    "tax_group=s"           => \$organism,           #fungi bacteria, euakaryote
    "readOverlap=i"         => \$maxReadOverlap,
	"mergePreClusterReads=i"     => \$mergePreCluster,
    "endRem=s"              => \$remFromEnd,
    "swarm_distance=i"      => \$swarmClus_d,
	"dada2seed=i"           => \$dada2Seed,
    "OTUbuild=s"            => \$otuRefDB,
    "count_chimeras"        => \$chimCnt,            #T or F
    "offtargetDB=s"         => \$custContamCheckDB,
    "flash_param=s"         => \$flashCustom,
    "deactivateChimeraCheck=i" => \$noChimChk,
	#"VsearchChimera=i"		=> \$useVsearch,
	"useVsearch=i"		=> \$useVsearch,
	"removePhiX=i"			=> \$doPhiX,
	"buildPhylo=i"          => \$buildPhylo,

    # "flashAvgLength" => \$flashLength,
    #"flashAvgLengthSD" => \$flashSD,
) or usage();


#routines that are external of a proper 16S run
if ($versionOut){
	print "$selfID\n";
	exit(0);
}
if ( $check_map ne "" ) {
	$mapFile = $check_map;
	my ( $x11, $x22, $x33 ) = readMap();
	print "\n\nmapping file seems correct\n";
	exit(0);
}
if ($create_map ne ""){
	autoMap($input, $create_map);
	exit(0); 
}
#declare global vars
my $logDir       = $outdir . "/LotuSLogS/";
my $extendedLogD = $outdir . "/ExtraFiles/";
$mainLogFile = $logDir . "LotuS_run.log";
my $cmdLogFile = $logDir . "LotuS_cmds.log";
my $progOutPut = $logDir . "LotuS_progout.log";
#create dirs
if ($TaxOnly eq "0" ){
	system("rm -f -r $outdir") if ( $exec == 0 && $onlyTaxRedo == 0 );
	system "rm -f $progOutPut\n"; #from here log systemL
	
} else {
	unless (-d $outdir && -d $logDir){
		if (-f $TaxOnly){
			$outdir = $TaxOnly;$outdir =~ s/[^\/]+$//;
			print "TaxOnly option specified, but not an output dir. Assumming: $outdir\n" ;
		} else {
			die "No output dir defined, only possible if -taxOnly argument is a file\n";
		}
		
		$logDir       = $outdir . "/LotuSLogS/";
		$mainLogFile = $logDir . "LotuS_run.log";
		$extendedLogD = $outdir . "/ExtraFiles/" ;
		$cmdLogFile = $logDir . "LotuS_cmds.log";
		$progOutPut = $logDir . "LotuS_progout.log";
	}
}
#reset logfile
#create output folders
system("mkdir -p $outdir/primary/;") unless ( -d "$outdir/primary" );
system("mkdir -p $outdir") unless ( -d $outdir );
system("mkdir -p $logDir") unless ( -d $logDir );
system("mkdir -p $extendedLogD;") if ( !-d $extendedLogs && $extendedLogs);	
open LOG, ">", $mainLogFile or die "Can't open Logfile $mainLogFile\n";
open cmdLOG , ">$cmdLogFile" or die "Can't open cmd log $cmdLogFile\n";



#just setup options to default pars, check consistency of options etc
prepLtsOptions();



#-----------------
#all parameters set; pipeline starts from here
printL( frame( $selfID . "\n". $cmdCall ), 0 );

my @inputArray = split( /,/, $input );
my @inqArray   = split( /,/, $inq );
$numInput = scalar(@inputArray);
if ( scalar(@inqArray) > 0
    && $numInput != scalar(@inqArray)
    && -f $inqArray[0] ){
    printL("Error: fasta input file number does not correspond ot quality file number.\n", 1);
}

#unless ($platform eq "miSeq" || $platform eq "hiSeq") {#only support paired reads for hi/miSeq
#	$numInput = 1;
#}

#read map, also check map file format
my ( $mapHref, $combHref, $hasCombiSmpls ) = ( {}, {}, 0 );
my %mapH; my $cpMapFile = "$outdir/primary/in.map";
if ( $TaxOnly eq "0" ) {
    ( $mapHref, $combHref, $hasCombiSmpls ) = readMap();
	my ($hr,$hr2 )= addSeqRun($mapHref);#SequencingRun
	writeMap($hr2,$cpMapFile);
	%mapH = %{$hr};
	#systemL("cp $map $outdir/primary\n");
}

#die "$cpMapFile\n\n";

#die $exec."\n";
if ( $ClusterPipe == 0 ) {
    printL( "Warning: otupipe sequence clustering mode is depreceated.\n\n",0 );
}
if ( $ClusterPipe == 0 && $id_OTU_noise < $id_OTU ) {
    printL( "id_OTU must be bigger-or-equal than id_OTU_noise\n", 2 );
}
if ( $id_OTU > 1 || $id_OTU < 0 ) {
    printL( "\"-id\" set to value <0 or >1: $id_OTU\nHas to be between 0 and 1\n", 2 );
}

systemL("rm -f -r $lotus_tempDir;") if ( $exec == 0 );
systemL("mkdir -p $lotus_tempDir;") unless ( -d $lotus_tempDir );
if ( !-d $lotus_tempDir ) {
    die( "Failed to make tmp dir or doesn't exist: \"" . $lotus_tempDir . "\"\n" );
}

#CONSTANT file paths
my $highLvlDir   = $outdir . "/higherLvl/";
my $FunctOutDir = $outdir . "/derrivedFunctions/";
my $RDP_hierFile = "$outdir/hiera_RDP.txt";
my $SIM_hierFile = "$outdir/hiera_BLAST.txt";


#my $currdir=`pwd`;
my $clustMode = "de novo";
$clustMode = "reference closed" if ( $otuRefDB eq "ref_closed" );
$clustMode = "reference open"   if ( $otuRefDB eq "ref_open" );
announce_options($ClusterPipe);

#=========
# Pipeline

# ////////////////////////// sdm 1st time (demult,qual etc) /////////////////////////////////////////////
#  cmdArgs["-i_MID_fastq"]
my ($sdmIn,$derepOutHQ,$qualOffset,$filterOutAdd,$filterOut,
		$derepOutHQ2,$derepOutMap,$sdmDemultiDir)=sdmStep1();

#die();
#exit;
# ////////////////////////// OTU building ///////////////////////////////////////


my $A_UCfil;
my $tmpOTU = "$lotus_tempDir/tmp_otu.fa";
my $OTUfa  = "$outdir/otus.fa";



my $ucFinalFile = "";
$duration = time - $start;
if ( $ClusterPipe != 0 && $onlyTaxRedo == 0 && $TaxOnly eq "0" ) {
	#$ClusterPipe == 0 
	($A_UCfil) = buildOTUs($tmpOTU);
	my @tmp = @$A_UCfil;
	#die "@tmp\n";
	$ucFinalFile = $tmp[0];
}
elsif ( $onlyTaxRedo == 1 ) {
    printL "Clustering step was skipped\n", 0;
}
elsif ( $TaxOnly ne "0" ) {
    printL
"No qual filter, demultiplexing or clustering required, as taxonomy only requested\n",
      0;
}
else {
    printL"clustering method (-CL/-clustering) unknown, given argument was '$ClusterPipe'\n",5;
}

#/////////////////////////////////////////  SEED extension /////////////////////
#$derepOutMap,$derepOutHQ
my @preMergeSeedsFiles; my @mergeSeedsFilesSing;
my $OTUSEED    = "$lotus_tempDir/otu_seeds.fna";
my $OTUrefSEED = "$lotus_tempDir/otu_seeds.fna.ref";

#uc additions have to be in matrix creation (>0.96 through sdm)
#my $UCadditions = $ucFinalFile.".ADD";
my $OTUmFile = "$outdir/OTU.txt";
my $OTUmRefFile =   "$outdir/OTU_psRef.txt";    #diff to OTU.txt: collapse entries with same ref
$OTUmRefFile = "" unless ($pseudoRefOTU);
my $seedExtDone = 1;
my $refClusSDM  = "";
my $didMerge = 0; #were reads already merged before this step? -> so far always no


if ($REFflag){ #these need to be treated extra, as no new optimal ref seq needs to be identified..
    $refClusSDM .="-optimalRead2Cluster_ref $ucFinalFile.ref -OTU_fallback_refclust $tmpOTU.ref";
    $refClusSDM .=" -ucAdditionalCounts_refclust $ucFinalFile.ADDREF -ucAdditionalCounts_refclust1 $ucFinalFile.RESTREF";
    #systemL "cat $UCadditions"."REF"." >> $UCadditions"
}

my $sdmOut2 = "-o_fna $OTUSEED";
if ( $numInput == 2 ) {
    push( @preMergeSeedsFiles, "$lotus_tempDir/otu_seeds.1.fq" );
    push( @preMergeSeedsFiles, "$lotus_tempDir/otu_seeds.2.fq" );
	@mergeSeedsFilesSing = @preMergeSeedsFiles; 
	$mergeSeedsFilesSing[0] =~ s/\.1\.fq/\.1\.singl\.fq/;$mergeSeedsFilesSing[1] =~ s/\.2\.fq/\.2\.singl\.fq/;

    $sdmOut2    = "-o_fastq " . $preMergeSeedsFiles[0] . "," . $preMergeSeedsFiles[1];
    $OTUrefSEED = "$lotus_tempDir/otu_seeds.1.fq.ref";
    #TODO : $sdmIn
}
my $upVer = "";
$upVer = "-uparseVer $usearchVer " if ($ClusterPipe == 1 || $ClusterPipe == 7); #use same format for .uc file in uparse and dada2 mode
$upVer = "-uparseVer N11 "if ($ClusterPipe == 6);
#note that -mergedPairs $didMerge refers to old implementation with read merge done before sdm (not supported any longer, legacy option)
my $mergeOptions = "-merge_pairs_seed 1 ";
my $sdmcmd="";

#----------------  SEED extension -------------------
#sdm based Seed extension of OTU clusters
if ($sdmDerepDo) {
    $sdmIn = "-i_fastq $derepOutHQ";
    if ( $numInput == 2 ) { $sdmIn = "-i_fastq $derepOutHQ,$derepOutHQ2"; }
    $sdmcmd ="$sdmBin $sdmIn $sdmOut2 $upVer -optimalRead2Cluster $ucFinalFile -paired $numInput -sample_sep $sep_smplID -derep_map $derepOutMap -options $sdmOpt $mergeOptions $qualOffset -mergedPairs $didMerge -log $logDir/SeedExtensionStats.txt  -OTU_fallback $tmpOTU -ucAdditionalCounts $ucFinalFile.ADD -ucAdditionalCounts1 $ucFinalFile.REST -otu_matrix $OTUmFile -count_chimeras $chimCnt $refClusSDM";
} else {
    $sdmcmd = "$sdmBin $sdmIn $sdmOut2 $upVer -optimalRead2Cluster $ucFinalFile -paired $numInput -sample_sep $sep_smplID -map $cpMapFile -options $sdmOpt $qualOffset $mergeOptions -mergedPairs $didMerge -log $logDir/SeedExtensionStats.txt -OTU_fallback $tmpOTU -otu_matrix $OTUmFile -ucAdditionalCounts $ucFinalFile.ADD -ucAdditionalCounts1 $ucFinalFile.REST -count_chimeras $chimCnt $refClusSDM";
}
my $status = 0;

#die $sdmcmd."\n";
if ( $exec == 0 && $onlyTaxRedo == 0 && $TaxOnly eq "0" ) {
	my $mrgSentence = ""; $mrgSentence = "and merging pairs of " if ($mergeOptions ne "");
    printL (frame("Extending $mrgSentence${OTU_prefix} Seeds"), 0);
    $status = systemL($sdmcmd);
}
elsif ($onlyTaxRedo) { printL "Skipping Seed extension step\n", 0; }
if ($status) {
    printL "Failed $sdmcmd\n",                    0;
    printL "Fallback to ${OTU_prefix} median sequences.\n", 0;
    $seedExtDone = 0;
	$OTUSEED = $tmpOTU;
    #exit(11);
} else {
    $seedExtDone = 1;
}
undef $tmpOTU;



#/////////////////////////////////////////  paired read merging /////////////////////
mergeRds();
#die "$OTUSEED\n";
if ($TaxOnly ne "0" && -f $TaxOnly){
	$OTUfa = $TaxOnly;
	if ($TaxOnly =~ m/\.gz$/){$OTUfa =~ s/\.gz//; systemL "gunzip -c $TaxOnly >$OTUfa "; }
}
#die "$OTUSEED\n";
#
#merge not combined reads & check for reverse adaptors/primers
#do this also in case no pairs were found at all.. Singletons need to be fq->fna translated into $OTUSEED

#at this point $tmpOTU (== $OTUSEED) contains the OTU ref seqs

#add additional sequences to unfinal file
#my $lnkHref = numberOTUs($tmpOTU,$OTUfa,${OTU_prefix});
#new algo 0.97: don't need to rename reads any longer (and keep track of this), done by sdm - but no backlinging to uparse any longer
#/////////////////////////////////////////  check for contaminants /////////////////////
my $OTUrefDBlnk;
if ( $exec == 0 && $onlyTaxRedo == 0 && $TaxOnly eq "0" ) {
	#first remove denovo OTUs.. pretty much a no-brainer
	chimera_denovo($OTUSEED);
    #remove chimeras on longer merged reads
    my $refChims = chimera_ref_rem( $OTUSEED);
    #ITSx - ITS region fungi only 
    my $nonITShref = ITSxOTUs($OTUSEED);
	#V-Xtractor
	my $nonVXTRs = vfxtractor($OTUSEED);
    #phiX
	my $phiXhref = contamination_rem( $OTUSEED, $CONT_REFDB_PHIX, "phiX" );
    #custom DB for contamination - off-target
	my $xtraConthref = contamination_rem( $OTUSEED, $custContamCheckDB, "offTarget" );
	#merge contamination results..
	${$xtraConthref}{"phiX"} = ${$phiXhref}{"phiX"} ;
	${$xtraConthref}{"ITSx"} = $nonITShref;

    #link between OTU and refDB seq - replace with each other
    $OTUrefDBlnk = readLinkRefFasta( $OTUrefSEED . ".lnks" );
	
	#finalize OTU file;
	systemL "cp $OTUSEED $OTUfa";
	#die "$OTUfa\n";

    #but OTU matrix already written, need to remove these
    my $XtalkRef = checkXtalk( $OTUfa, $OTUmFile );
	${$xtraConthref}{"Xtalk"} = $XtalkRef;
	
	my $luluref= lulu($OTUfa,$OTUmFile);
	${$xtraConthref}{"LULU"} = $luluref;
	
	#integrate all these filters now
    $OTUrefDBlnk = clean_otu_mat($OTUfa,$OTUrefSEED, $OTUmFile,
        $OTUrefDBlnk,$xtraConthref) if ($seedExtDone);    #&& $otuRefDB ne "ref_closed" );
        # $OTUfa.ref contains reference Seqs only, needs to be merged later..
        #and last check for cross talk in remaining OTU match
		
	
		
} elsif ($onlyTaxRedo) {
    printL "Skipping removal of contaminated OTUs step\n", 0;
}
# ////////////////////////// TAXONOMY ////////////////////////////////////////////
my $RDPTAX = 0;
my $REFTAX = 0;

runRDP();

my $cmd="";

#some warnings to throw
if ( !$doBlasting && substr( $ampliconType, 0, 3 ) eq "ITS" ) {
    my $failedBlastITS = "ITS region was chosen as target; this requires a similarity based taxonomic annotation and excludes RDP tax annotation.\n";
    $failedBlastITS .=       "Blast similarity based annotation is not possible due to: ";
    if ( !$doBlasting ) {
        $failedBlastITS .= "Similarity search being deactivated.";
    }
    $failedBlastITS .= "\nTherefore LotuS had to abort..\n";
    printL $failedBlastITS, 87;
}
if ($TaxOnly ne "0") {
    my $hier = assignTaxOnly( $OTUfa, $outdir );
	if ($rmOutDir){
		systemL("rm -rf $outdir;"); #online in extreme cases, keep well defined & controlled
		printL frame("Taxonomy for \n$TaxOnly\n has been assigned to \n$hier"),0;
	} else {
		printL frame("Taxonomy has been assigned to $input, output in \n$outdir\n"), 0;
	}
    exit(0);
}

#pre 0.97
#my $lnkHref="";
#my ($OTUmatref,$failsR) = makeAbundTable($taxblastf,"$lotus_tempDir/RDPotus.tax",$A_UCfil,$OTUmFile,$lnkHref,\@avSmps);
my ( $OTUmatref, $avOTUsR ) = readOTUmat($OTUmFile);

#debug
#my %retMat = %{$OTUmatref};my @hdss = keys %retMat;my @fdfd = keys %{$retMat{bl21}};die "@hdss\n@fdfd\n$retMat{bl14}{OTU_1} $retMat{bl14}{OTU_2} $retMat{bl14}{OTU_3}\n";

#this subroutine also has blast/LCA algo inside
#die $OTUfa."\n";
my ($failsR) = makeAbundTable2( "$lotus_tempDir/RDPotus.tax", $OTUmatref );    #,\@avSmps);

if ( !-d $highLvlDir ) {
    if ( systemL("mkdir -p $highLvlDir;") ) {
        printL("Could not create Higher level abundance matrix directory $highLvlDir.", 23);
    }
}
if (0 && !-d $FunctOutDir ) { #currently not used, deactivate
    if ( systemL("mkdir -p $FunctOutDir;") ) {
        printL("Could not create Functional level abundance matrix directory $FunctOutDir.", 29);
    }
}
systemL("cp $OTUmFile $highLvlDir;");

#higher taxonomy & bioms
if ( $REFTAX || $RDPTAX ) {
	my $taxRefHR;
    my $table_dir = "$outdir/Tables";
    if ($REFTAX) {
        printL frame("Calculating Taxonomic Abundance Tables from @refDBname assignments"), 0;
        $taxRefHR = calcHighTax( $OTUmatref, $SIM_hierFile, $failsR, 1, $OTUmRefFile );
        biomFmt( $OTUmatref, $SIM_hierFile, "$outdir/OTU.biom", 1, {} );
        if ($pseudoRefOTU) {
            my ( $OTUmatref2, $avOTUsR2 ) = readOTUmat($OTUmRefFile);
            biomFmt( $OTUmatref2, $SIM_hierFile, "$outdir/OTU_psRef.biom", 1,
                $taxRefHR );
        }

    } elsif ($RDPTAX) {
        printL(frame("Calculating Taxonomic Abundance Tables from RDP \nclassifier assignments, Confidence $RDPCONF "),0);
        $taxRefHR = calcHighTax( $OTUmatref, $RDP_hierFile, $failsR, 0, "" );
        biomFmt( $OTUmatref, $RDP_hierFile, "$outdir/OTU.biom", 0, {} );
    }

	#   TODO 
	#annotate OTU's with functions
	#my $OTU2Funct = annotateFaProTax($taxRefHR,$FaProTax); #TODO
	#calcHighFunc($OTU2Funct,$FunctOutDir); #TODO
}

#merge in ref seq fastas
if ($REFflag) {
    systemL("cat $OTUfa.ref >> $OTUfa;");
    unlink "$OTUfa.ref";
}

#building tree & MSA
my $treeF = "";
$treeF = buildTree( $OTUfa, $outdir );

#citations file
open O, ">$logDir/citations.txt" or printL" Failed opening $logDir/citations.txt\n",0;
print O $citations;
close O;

my $phyloseqCreated=runPhyloObj($treeF);

systemL("rm -rf $lotus_tempDir;") if ($exec == 0 && !$keepTmpFiles);  #printL "Delete temp dir $lotus_tempDir\n", 0; }
systemL("rm -rf $outdir;") if ($rmOutDir); #online in extreme cases, keep well defined & controlled
printL(    frame("LotuS2 finished. Output:\n$outdir\n\- LotuSLogS/ contains run statistics (useful for describing data/amount of reads/quality\n- LotuSLogS/citations.txt: papers of programs used in this run\nNext steps: you can use the rtk program in this pipeline, to generate rarefaction curves and diversity estimates of your samples.\n"    ),0);
my $phyloHlp="";
$phyloHlp="- Phyloseq: $outdir/phyloseq.Rdata can be directly loaded in R\n" if ($phyloseqCreated);
printL(frame("          Next steps:          \n- Rarefaction analysis: can be done with rtk (avaialble in R or use bin/rtk)\n$phyloHlp- Phylogeny: ${OTU_prefix} phylogentic tree available in $outdir/tree.nwk\n- .biom: $outdir/OTU.biom contains biom formated output\n- tutorial: Visit http://lotus2.earlham.ac.uk for more tutorials in data analysis\n"));
printWarnings();

close LOG; close cmdLOG;











































#--------------------------############################----------------------------------###########################



sub sdmStep1{
	my $sdmcmd       = "";
	my $filterOut    = "$lotus_tempDir/demulti.fna";
	my $filterOutAdd = "$lotus_tempDir/demulti.add.fna";
	if ( $numInput > 1 ) {
		$filterOut    = "$lotus_tempDir/demulti.1.fna,$lotus_tempDir/demulti.2.fna";
		$filterOutAdd = "$lotus_tempDir/demulti.1.add.fna,$lotus_tempDir/demulti.2.add.fna";
	}
	my $filOutCmd = "-o_fna ";
	if ($UPARSEfilter) {
		$filterOut    = "$lotus_tempDir/demulti.fastq";
		$filOutCmd    = "-o_fastq ";
		$filterOutAdd = "$lotus_tempDir/demulti.add.fastq";
		if ( $numInput > 1 ) {
			$filterOut    = "$lotus_tempDir/demulti.1.fastq,$lotus_tempDir/demulti.2.fastq";
			$filterOutAdd = "$lotus_tempDir/demulti.1.add.fastq,$lotus_tempDir/demulti.2.add.fastq";
		}
	}
	my $qualOffset = "-o_qual_offset 33";       #33 for UPARSE
	my $sdmOut     = $filOutCmd . $filterOut;
	my $sdmIn      = "";
	my $paired     = $numInput;

	#for now: only use fwd pair
	#if ($paired != 1){$paired.= " -onlyPair 1";}
	if ( $barcodefile ne "" ) {
		$sdmIn .= " -i_MID_fastq $barcodefile";
	}
	my $derepCmd    = "";
	my $paired_sdm  = "";
	my $derepOutHQ  = "";
	my $derepOutHQ2 = "";
	my $derepOutMap = "";

	#sdm merge options: -merge_pairs_filter -merge_pairs_demulti -merge_pairs_derep
	my $sdmMergeOpt = "";
	if ($mergePreCluster){
		$sdmMergeOpt = "-merge_pairs_derep 1 ";
		if ($ClusterPipe == 7){#dada2 requires more output..
		$sdmMergeOpt .= "-merge_pairs_demulti 1 ";
		}
	}

	my $sdmDemultiDir = "";
	my $sdmOptStr     = "-options $sdmOpt ";
	if ( $saveDemulti == 2 || $saveDemulti == 1 || $ClusterPipe == 7 ) { #dada2 also requires filtered raw reads
		$sdmDemultiDir = "$outdir/demultiplexed/";
		if ( $saveDemulti == 1 ) {
			printL frame("Demultiplexed input files into single samples, no quality filtering done\nWill abort after demultiplexing (mode \"1\" is only for preparing publication ready, raw sequences)!\n"),0;
			$sdmOptStr = "";
		}
		if ($ClusterPipe == 7){ #dada2.. no rd pair info in head!
			$sdmOptStr .= "-pairedRD_HD_out 0 -pairedDemulti 1 -derep_format fq -derepPerSR 1 -DemultiBPperSR 1e8";
			$sdmDemultiDir = "$lotus_tempDir/demultiplexed/" if ($saveDemulti == 0);
		}
	}
	
	#die "$sdmOptStr\n";

	if ($sdmDerepDo) {
		$derepCmd = "-o_dereplicate $lotus_tempDir/derep.fas ";
		if ( 0 && $ClusterPipe == 2 ) {
			$derepCmd .= "-dere_size_fmt 1 ";
		}
		else { $derepCmd .= "-dere_size_fmt 0 "; }
		$derepCmd .= " -min_derep_copies $dereplicate_minsize ";
		$derepOutHQ  = "$lotus_tempDir/derep.1.hq.fq";
		$derepOutMap = "$lotus_tempDir/derep.map";
		$derepOutHQ2 = "$lotus_tempDir/derep.2.hq.fq";
		if ( $saveDemulti == 3 ) {   #temp deactivated, since I have  $sdmDemultiDir
			$derepCmd .= " -suppressOutput 0";
		}
		else {
			$derepCmd .= " -suppressOutput 1";
		}
	}
	else {
		if ( $paired != 1 ) { $paired_sdm .= " -onlyPair 1"; }
	}
	my $demultiSaveCmd = "";
	if ( $sdmDemultiDir ne "" ) {
		systemL "mkdir -p $sdmDemultiDir;" unless ( -d $sdmDemultiDir );
		$demultiSaveCmd .= " -o_demultiplex $sdmDemultiDir";
	}
	my $dmgCmd = "";
	if ($damagedFQ) { $dmgCmd = "-ignore_IO_errors 1"; }
	my $mainSDMlog = "$logDir/demulti.log";
	if ($TaxOnly ne "0"){
		printL("Skipping Quality Filtering & demultiplexing & dereplication step\n",0);
		return ($sdmIn,$derepOutHQ,$qualOffset,$filterOutAdd,$filterOut,$derepOutHQ2,$derepOutMap,$sdmDemultiDir);
	}
	
	if ( -d $input ) {
		$sdmIn = "-i_path $input ";
	} elsif ( -f $inputArray[0] && $inq ne "" && -f $inqArray[0] ) {
		if ( $paired == 1 ) {
			$sdmIn = "-i_fna $inputArray[0] -i_qual $inqArray[0] ";
		} elsif ( $paired == 2 ) {
			$sdmIn = "-i_fna $inputArray[0],$inputArray[1] -i_qual $inqArray[0],$inqArray[1] ";
		}
	} else {
		if ( $paired == 1 ) {
			$sdmIn = "-i_fastq $inputArray[0]";
		} elsif ( $paired == 2 ) {
			$sdmIn = "-i_fastq $inputArray[0],$inputArray[1]";
		}
	}
	my $mrgOpt = "";#"-merge_pairs 1";
	my $multCore = ""; $multCore = " -threads $uthreads";

	#primary sequence filtering + demultiplexing + dereplication
	$sdmcmd = "$sdmBin $sdmIn $sdmOut -sample_sep $sep_smplID  -log $mainSDMlog -map $cpMapFile $sdmMergeOpt $sdmOptStr $demultiSaveCmd $derepCmd $dmgCmd $qualOffset -paired $paired $paired_sdm -maxReadsPerOutput $linesPerFile -oneLineFastaFormat 1 $mrgOpt $multCore ";    #4000000
	#die $sdmcmd."\n";


	if ( $exec == 0 && $onlyTaxRedo == 0 && $TaxOnly eq "0" ) {
		#$duration = time - $start;
		printL( frame("Demultiplexing, filtering, dereplicating input files, this might take some time..\ncheck progress at $progOutPut\n",1,2),0 );
		systemL("cp $sdmOpt $outdir/primary");
		if ( systemL($sdmcmd) != 0 ) {
			printL "FAILED sdm demultiplexing step: " . $sdmcmd . "\n";
			exit(4);
		}
		my $fileNum = `ls -1 $mainSDMlog* | wc -l`;
		my @sdmOfiles = glob("${mainSDMlog}0*");

		if ( $fileNum > 0 && @sdmOfiles>0) {
			systemL "mkdir -p $logDir/SDMperFile/; mv $mainSDMlog". "0* $logDir/SDMperFile/";
		}
		if ( $fileNum > 10 ) {
			systemL "tar zcf $logDir/SDMperFile.tar.gz $logDir/SDMperFile/; rm -r $logDir/SDMperFile/;";
		}
		my $readSdmLog = `cat $mainSDMlog`;
		if ( $readSdmLog =~ m/binomial est\. errors/ ) {
			$citations .= "Poisson binomial model based read filtering: Fernando Puente-Sánchez, Jacobo Aguirre, Víctor Parro (2015).A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies. Nucleic Acids Res.(2015).\n";
		}
		
		#final report
		$readSdmLog =~ m/(Reads processed:.*)/;
		my $shrtRpt = $1;
		$readSdmLog =~ m/(Rejected:.*)/;
		$shrtRpt .= "\n".$1;
		$readSdmLog =~ m/(Accepted \(Mid\+High qual\):.*)/;
		$shrtRpt .= "\n".$1;
		printL( frame("Finished primary read processing with sdm:\n".$shrtRpt."\nFor an extensive report see $mainSDMlog\n",1,3),0 );

		#postprocessing of output files
		if ( $saveDemulti == 1 || $saveDemulti == 2 ) {    #gzip stuff
			printL "Zipping demultiplex output..\n";
			systemL "gzip $sdmDemultiDir/*.fq;";
		}
		if ( $saveDemulti == 1 ) {
			printL frame("Demultiplexed intput files with no quality filtering to:\n$sdmDemultiDir\nFinished task, if you want to have a complete LotuS run, change option \"-saveDemultiplex 0\".\n",1),0;
			exit(0);
		}
		if ( $saveDemulti == 3 && $exec == 0 )
		{    #&& $onlyTaxRedo==0){#gzip demultiplexed fastas
			systemL "mkdir -p $outdir/demultiplexed; ";
			my @allOuts = split /,/, $filterOut;
			foreach (@allOuts) {
				$_ =~ m/\/([^\/]+$)/;
				my $fn  = $1;
				my $fns = $fn;
				$fns =~ s/\.(f[^\.]+)$/\.singl\.$1/;

	#			systemL "gzip -c $_*.singl > $outdir/demultiplexed/$fn.singl.gz; rm -f $_*.singl";
	#			systemL "gzip -c $_* > $outdir/demultiplexed/$fn.gz; rm -f $_*";
				die "DEBUG new singl filenames\n";
				systemL "gzip -c $fns > $outdir/demultiplexed/$fn.singl.gz; rm -f $_*.singl;";
				systemL "gzip -c $_* > $outdir/demultiplexed/$fn.gz; rm -f $_*;";
			}
			@allOuts = split /,/, $filterOutAdd;
			foreach (@allOuts) {
				$_ =~ m/\/([^\/]+$)/;
				systemL "gzip -c $_* > $outdir/demultiplexed/$1.gz; rm -f $_*;";
			}
		}
		if ($ClusterPipe == 7){#derep by single file, combine map at end
			#$cmd = "cat $lotus_tempDir/derep.*.map > $lotus_tempDir/derep.map;";
			#$cmd = "cat $lotus_tempDir/derep.*.fas > $lotus_tempDir/derep.fas;";
			#systemL $cmd;
		}
		


	} elsif ($onlyTaxRedo) {
		printL("Skipping Quality Filtering & demultiplexing & dereplication step\n",0);
	}
	
	return ($sdmIn,$derepOutHQ,$qualOffset,$filterOutAdd,$filterOut,$derepOutHQ2,$derepOutMap,$sdmDemultiDir);
}


sub cntFastaEntrs($){
	my $tmp = `grep -c '^>' $_[0]`;chomp $tmp;
	return $tmp;
}

sub runPhyloObj{
	my ($treeF) = @_;
	return 0 if (!-f $phyloLnk || !-f $Rscript);
	die "Incorrect phyloLnk script defined $phyloLnk" unless (-f $phyloLnk);
	die "Incorrect R installation (can't find Rscript)" unless (-f $Rscript);
	my $cmd = "$Rscript $defRscriptOpt $phyloLnk $OTUmFile $SIM_hierFile $cpMapFile $treeF;";
	#die "$cmd\n";
	if ((systemL $cmd)){
		printL "Could not create phyloseq object\n","w";
		return 0;
	}
	return 1;
}

sub printWarnings() {
    if ( $finalWarnings eq "" ) { return; }
    printL "\nThe following WARNINGS occured:\n", 0;
    printL $finalWarnings. "\n", 0;
}

sub readLinkRefFasta() {
    my ($inF) = @_;
    my %ret;
    if ( !-e $inF ) { return \%ret; }
    open I, "<$inF";
    while ( my $line = <I> ) {
        chomp $line;
        my @spl = split( /\t/, $line );
        $ret{ $spl[0] } = $spl[1];
    }
    close I;

    #unlink $inF;
    return \%ret;
}

sub loadFaProTax($){
	my $opt_db = $_[0];
	#my $DB = retrieve($opt_db);
	#return $DB;
}

sub vfxtractor{
	my ($otusFA) = @_;
	my %ret;
	return \%ret unless ( substr( $ampliconType, 0, 3 ) eq "SSU" );
	return \%ret if ( $doVXTR eq "0" );
	if ( !-e $vxtrBin ) {
		printL "Did not find V-xtractor binary at $vxtrBin\nNo V extraction used\n","w";
		return \%ret;
	}
	if (!-d $hmmsVX){
		printL "Did not find V-xtractor HMMs at $hmmsVX\nNo V extraction used\n","w";
		return \%ret;
	}
	if ( !-e $hmmsrchBin ) {
		printL "Did not find hmmscan binary at $hmmsrchBin\nNo V extraction used\n","w";
		return \%ret;
	}
	my $outBFile = $otusFA . ".vxtr";
	my $outCFile = $otusFA . ".vxtr.csv";
	my $defReg = ".V1-V9.";
	my $region = $defReg;
	$region = $doVXTR unless ($doVXTR eq "1");

	my $cmd = "$vxtrBin -o $outBFile -r $region -c $outCFile -hmmdir $hmmsVX/".substr( $ampliconType, 0, 3 )."/bacteria/ -nc $uthreads -hmmsc $hmmsrchBin $otusFA;"; #-nc $uthreads
	$cmd .= "cp $outCFile $logDir/VXtractor.summary.txt\n";
	if ( systemL($cmd) != 0 ) { printL( "Failed command:\n$cmd\n", 1 ); }

#die "VXdie ".$cmd;
#if (-z "$outBFile.full.fasta"){printL "Could not find any valid ITS OTU's. Aborting run.\n Remaining OTUs can be found in $outBFile*\n",923;}
	my $hr;
	my %ITSo;
	
	$hr   = readFasta("$outBFile");
	%ITSo = %{$hr};

	#my $ITSfa = `grep -c '^>' $outBFile.full.fasta`;
	#my $orifa = `grep -c '^>' $otusFA`;chomp $orifa;    #chomp $ITSfa;
	$hr = readFasta($otusFA); my %FNA = %{$hr};
	
	if ( scalar( keys(%ITSo) ) == 0 ) {
		printL "Could not find any valid ${OTU_prefix}'s within variable regions $region (via V-Xtractor). Aborting run.\n Remaining OTUs can be found in $outBFile*\n", 923;
	}
	
	#check which OTUs existed before, but are not longer listed in ITSx output
	my @prevOTUs = keys %FNA;
	my %lookup;
	foreach my $k ( keys %ITSo ) {
		my $hdde = $k;# substr($k,1);#"${OTU_prefix}x_1";
		#print $hdde;
		$lookup{$hdde} = 1;
	}
	
	my $delOTUs=0;
	for my $otu ( @prevOTUs ) {
		if (!exists($lookup{ $otu })){
			$ret{$otu} = 1 ;
			$delOTUs++;
		}
	}
	

	printL frame( "V region extraction: Kept " . scalar( keys(%ITSo) ) . ", deleted $delOTUs ${OTU_prefix}'s identified as regions $doVXTR (of "  . scalar( keys(%FNA) ) . " ${OTU_prefix}'s).\n"),0;
	systemL "rm -f $outBFile*;";
	$citations .= "V extraction of ribosomal regions: V-Xtractor: an open-source, high-throughput software tool to identify and extract hypervariable regions of small subunit (16S/18S) ribosomal RNA gene sequences. Martin Hartmann 1, Charles G Howes, Kessy Abarenkov, William W Mohn, R Henrik Nilsson\n";
	return \%ret;


}

sub ITSxOTUs {
	my ($otusFA) = @_;
	my %ret;
	return \%ret unless ( substr( $ampliconType, 0, 3 ) eq "ITS" );
	return \%ret if ( !$doITSx );
	if ( !-e $itsxBin ) {
		printL "Did not find ITSx binary at $itsxBin\nNo ITS extraction used\n","w";
		return \%ret;
	}
	if ( !-e $hmmsrchBin ) {
		printL "Did not find hmmscan binary at $hmmsrchBin\nNo ITS extraction used\n","w";
		return \%ret;
	}

	my $outBFile = $otusFA . ".itsX";
	my $ITSxReg  = "ITS1,ITS2";
	if ( $ampliconType eq "ITS2" ) {
		$ITSxReg = "ITS2";
		printL "Setting to ITS2 region\n";
	}
	if ( $ampliconType eq "ITS1" ) {
		$ITSxReg = "ITS1";
		printL "Setting to ITS1 region\n";
	}
	my $itsxOrg = "all";
	$itsxOrg = "F" if ( lc($organism) eq "fungi" );
	my $xcmd = "";
	
	#die "$itsxOrg\n";
	#perl ITSx -h
	#ITSx -- Identifies ITS sequences and extracts the ITS region
	#by Johan Bengtsson-Palme et al., University of Gothenburg
	#Version: 1.0.11

	my $preChk = `perl $itsxBin -h 2>&1`;
	if ($preChk =~ m/Version: 1\.([\.0-9]+)/){
		if ($1 >= 1.3){$preChk = "--temp $lotus_tempDir/itsx/";systemL "mkdir -p $lotus_tempDir/itsx/";}
	}
	my $cmd = "perl $itsxBin -i $otusFA -o $outBFile --cpu $uthreads --multi_thread T --heuristics T -t $itsxOrg --silent T --fasta T --save_regions $ITSxReg --partial $ITSpartial --hmmBin $hmmsrchBin --preserve T;";
	$cmd .= "cp $outBFile.summary.txt $logDir/ITSx.summary.txt\n";
	if ( systemL($cmd) != 0 ) { printL( "Failed command:\n$cmd\n", 1 ); }

#die "ITSxdie ".$cmd;
#if (-z "$outBFile.full.fasta"){printL "Could not find any valid ITS OTU's. Aborting run.\n Remaining OTUs can be found in $outBFile*\n",923;}
	my $hr;
	my %ITSo;
	
	if ($ampliconType eq "ITS"){
		#printL "cat $outBFile.ITS1.fasta $outBFile.ITS2.fasta >> $outBFile.full.fasta\n\n",0;
		systemL "cat $outBFile.ITS1.fasta $outBFile.ITS2.fasta >> $outBFile.full.fasta";
		$hr   = readFasta("$outBFile.full.fasta");
		%ITSo = %{$hr};
	} else{
		if ( $ampliconType eq "ITS1" || $ampliconType eq "ITS" ) {
			$hr   = readFasta("$outBFile.ITS1.fasta");
			%ITSo = %{$hr};
		} 
		if ( $ampliconType eq "ITS2" || $ampliconType eq "ITS" ) {
			$hr   = readFasta("$outBFile.ITS2.fasta");
			%ITSo = ( %ITSo, %{$hr} );
		}
	}

	#my $ITSfa = `grep -c '^>' $outBFile.full.fasta`;
	#my $orifa = `grep -c '^>' $otusFA`;chomp $orifa;    #chomp $ITSfa;
	$hr = readFasta($otusFA); my %FNA = %{$hr};
	
	if ( scalar( keys(%ITSo) ) == 0 ) {
		printL "Could not find any valid ITS ${OTU_prefix}'s. Aborting run.\n Remaining OTUs can be found in $outBFile*\n", 923;
	}
	
	#check which OTUs existed before, but are not longer listed in ITSx output
	my @prevOTUs = keys %FNA;
	my %lookup;
	foreach my $k ( keys %ITSo ) {
		my $hdde = $k;# substr($k,1);#"${OTU_prefix}x_1";
		
		#known format, but remove the ItSx info
		if ($k =~ m/^(ASV|[z]?OTU_\d+)\|/){
			$hdde = $1;
		}
		#same again, different known format, but remove the ItSx info
		if ($k =~ m/^(ASV|[z]?OTU_\d+);size=\d+\|/){
			$hdde = $1;
		}
		#print $hdde;
		$lookup{$hdde} = 1;
	}
	
	my $delOTUs=0;
	for my $otu ( @prevOTUs ) {
		if (!exists($lookup{ $otu })){
			$ret{$otu} = 1 ;
			$delOTUs++;
		}
	}
	

	printL frame( "ITSx analysis: Kept " . scalar( keys(%ITSo) ) . ", deleted $delOTUs ${OTU_prefix}'s identified as $ITSxReg (of "  . scalar( keys(%FNA) ) . " ${OTU_prefix}'s).\n"),0;
	#die scalar( keys(%ret) )."\n";

	#systemL "cat $outBFile.full.fasta > $otusFA";
	#open O, ">$otusFA" or die "Can't open output $otusFA\n";
	#foreach my $k ( keys %ITSo ) {
	#	my $hdde = "${OTU_prefix}x_1";
	#	if ($k =~ m/^(ASV|[z]?OTU_\d+)\|/){
	#		$hdde = $1;
	#	}
	#	print O ">$hdde\n$ITSo{$k}\n";
	#}
	#close O;
	#die "rm -f $outBFile*;";
	systemL "rm -f $outBFile*;";
	$citations .= "ITSx removal of non ITS OTUs: ITSx: Johan Bengtsson-Palm et al. (2013) Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for use in environmental sequencing. Methods in Ecology and Evolution, 4: 914-919, 2013\n";
	return \%ret;
	#die $cmd."\n";
}

sub contamination_rem($ $ $ ) {
    my ( $otusFA, $refDB1, $nameRDB1 ) = @_;
	my @refDBs = split /,/,$refDB1;
	my $contRemStr=0;
	my $hr = readFasta($otusFA);
	my %OTUs = %{$hr};
	my %totHits;
	my %ret;
	
	if (!$doPhiX && $nameRDB1 =~ m/^phiX/){
		return \%ret;
	}
	
	my $Lreport="";#text for frame at end of routine
	my $anycontRem = 0; #total count of contaminants across databases

	foreach (my $i=0;$i<@refDBs;$i++){
		my $refDB = $refDBs[$i];
		my $nameRDB = $nameRDB1.$i;

		my $outContaminated = "$logDir/out.cont.$nameRDB.fna";
		my $required        = 1;
		if ( $refDB eq "" ) { $required = 0; }
		systemL "rm -f $outContaminated;" if ( -e $outContaminated );
		my $contRem = 0;
		my $matchAlgo="";

		#printL frame "Searching for contaminant OTUs with $nameRDB ref DB";
		if ( $refDB ne "" && -f $refDB && -s $otusFA ) {
			#die "hex seq to 50kb pieces\n";
			my $hexDB = $refDB;
			$hexDB .= ".lts.fna";
			my $hitsFile = $otusFA . ".$nameRDB.cont_hit.uc";
			
			if ($nameRDB ne "PhiX"){ #minimap2, now default
				$hitsFile= "$otusFA.$nameRDB.cont_hit.paf";
				$cmd = "$mini2Bin  -N 1 -t $uthreads -o $hitsFile $otusFA $refDB ;";
				$matchAlgo = "minimap2";
				#die "$cmd\n";
			} elsif ( $VSused == 0 ) { #deprecated
				$cmd = "$usBin -usearch_local $otusFA -db $refDB -uc $hitsFile  -query_cov .8 -log $logDir/$nameRDB"
				  . "_contami_align.log ";
				$cmd .= "-id .9 -threads $uthreads -strand both ;";
				$matchAlgo = "usearch";
			}
			else {
				$cmd = "$VSBin -usearch_local $otusFA -db $refDB --maxseqlength 99999999999 -uc $hitsFile --query_cov .8 -log $logDir/$nameRDB" . "_contami_align.log ";
				$cmd .= "-maxhits 1 -top_hits_only -strand both -id .9 -threads $uthreads --dbmask none --qmask none; ";    #.95
				$matchAlgo = "vsearch";
			}

			#die $cmd."\n";
			if ( systemL($cmd) != 0 ) { printL( "Failed command:\n$cmd\n", 1 ); }
			#create tmp
			my %hits;
			if ($hitsFile =~ m/\.paf/){
				open I, "<", $hitsFile or die "Can't open search result file $hitsFile";
				while (<I>) {my @spl = split(/\t/); 
					#look for matches at 60% overlap
					if ( $spl[9] > $spl[6] * 0.5 || (($spl[9] / $spl[10]) > 0.9) ) { $hits{$spl[5]}=1; $contRem++; }
				}
				close I;
			} else {
				open I, "<", $hitsFile or die "Can't open search result file $hitsFile";
				while (<I>) {my @spl = split(/\t/);
					if ( $spl[0] eq "H" ) { $hits{$spl[8]}=1; $contRem++; }
				}
				close I;
			}
			my @hitA = keys %hits;
			my $storeContamIDsF = "$logDir/$nameRDB1.$nameRDB.otus";
			my $Xtxt ="\n";
			#$Xtxt .= "\nIDs stored in $storeContamIDsF" if ($contRem);
			if ($i==(@refDBs-1) && $anycontRem){
				$Xtxt = "";
				if ($keepOfftargets){$Xtxt .= ", will be retained in ${OTU_prefix} matrix.";} else {$Xtxt .= ", will be removed from ${OTU_prefix} matrix.";}
				$Xtxt .= "\nEnsure to use option -keepUnclassified to retain all these in output matrix." if (!$keepUnclassified && $contRem );
				$Xtxt .= "\n";
			}
			$Lreport .= "Found ". scalar(@hitA). " ${OTU_prefix}'s using $matchAlgo ($nameRDB: $refDB)$Xtxt";

			#die ("@hits\n");
			$anycontRem += $contRem;
			$totHits{$nameRDB} = \%hits;
			next if (scalar(@hitA) == 0);# { return(0);}
			systemL ("gzip $hitsFile\nmv $hitsFile.gz $logDir;");
			#open LLOG,">$storeContamIDsF";
			#print LLOG "$nameRDB\t$refDB\n@hitA\n";
			#close LLOG;
			#report
			#(add - deleted before) contaminated Fastas to file
			#open Ox, ">$outContaminated" or printL "Can't open contaminated OTUs file $outContaminated\n", 39;
			#foreach my $hi (@hitA) {
			#	print Ox ">" . $hi."_".$nameRDB . "\n" . $OTUs{$hi} . "\n";
			#}
			#close Ox;
			#systemL "gzip $outContaminated";
			#if (($nonChim-$emptyOTUcnt)==0){
			#	printL "Empty OTU matrix.. aborting LotuS run!\n",87;
			#}
			#print "\n\n\n\n\n\nWARN DEBUG TODO .95 .45 $outContaminated\n";
		} elsif ($required) {
			my $warnStr = "Could not check for contaminated OTUs, because ";
			unless ( $refDB ne "" && -f $refDB ) {
				$warnStr .= "\"$nameRDB\" reference database \"$refDB\"did not exist.\n";
			} else {  $warnStr .= "${OTU_prefix} fasta file was empty\n";
			}
			printL $warnStr,"w";
			$contRem=0;
			#systemL("cp $otusFA $outfile");
			#$outfile = "$lotus_tempDir/uparse.fa";
		}
		#print "CR $contRem\n";
		$contRemStr += $contRem;
	}
	#nothing found? don't bother writing file again..
	printL (frame( $Lreport), 0) if ($Lreport ne "");
	my @kh = keys %totHits;
	#print @kh. "   @kh\n";
	return \%totHits; #if ($contRemStr == 0);
	
	die; #should no longer go here..
	$contRemStr = scalar(keys(%totHits));
	#print "\n\n\"". scalar(keys %totHits) ."\n";
	#print remaining OTUs
	open Oa, ">$otusFA" or printL "Can't open ${OTU_prefix} file $otusFA\n", 39;
	my $cnt=0;
	foreach my $hi ( keys %OTUs ) {
		next if (exists($totHits{$hi}));
		print Oa ">" . $hi . "\n" . $OTUs{$hi} . "\n";
		$cnt++;
	}
	close Oa;
	printL "Writing $cnt/" . scalar(keys(%OTUs)) ." ${OTU_prefix} seeds after $nameRDB1 removal\n";
    return $contRemStr;
}
sub getOTUsize($){
	my ($otuM) = @_;
	my %ret;
	open I,"<$otuM" or die "can't open matrix $otuM\n";
	my $cnt=0;
	while(<I>){
		chomp;
		$cnt++; next if ($cnt == 1);
		my @spl = split /\t/;
		my $key = shift @spl;
		my $oc=0;
		foreach my $occ (@spl){$oc += $occ;}
		$ret{$key} = $oc;
	}
	close I;
	return \%ret;
}
sub addSize2OTUfna($ $){
	my ($otuFA,$hr) = @_;
	my %siz = %{$hr};
	
	replaceFastaHd($hr,$otuFA);
	
}

sub chimera_denovo($){
	my ($OTUfa) = @_;
	#uparse, unoise, dada2 have their own chimera checks
    if ($ClusterPipe == 1 ||$ClusterPipe == 7||$ClusterPipe == 6 || !-e $OTUfa || $noChimChk == 1 || $noChimChk == 2 ){
		return $OTUfa;
	}
	my $iniEntries = cntFastaEntrs($OTUfa);
	my $chimOut = "$lotus_tempDir/chimeras_denovo.fa";
	if ($extendedLogs) {
		$chimOut = "$extendedLogD/chimeras_denovo.fa";
	}
	
	my $progUsed = "vsearch uchime";
	$cmd = "$VSBin -uchime_denovo $OTUfa -chimeras $chimOut -nonchimeras $lotus_tempDir/tmp1.fa -abskew $chimera_absskew -log $logDir/chimera_dn.log;";

	#die "\n\n$usearchVer\n";
	if (!$useVsearch && $usearchVer >= 10 && !$VSused ) {
		if ( $usearchVer == 10.0 && $usearchsubV <= 240 ) {
			#really dirty hack..
			$cmd ="$VSBinOri -uchime_denovo $OTUfa -chimeras $chimOut -nonchimeras $lotus_tempDir/tmp1.fa -abskew $chimera_absskew -log $logDir/chimera_dn.log;";
			printL "Can't do de novo chimer filter, since usearch 10.0.240 currently has a bug with this\nUsing vsearch chimera detection instead\n","w";
		}
		else {
			#needs to have OTUsizes attached to OTUs
			my $hr = getOTUsize($OTUmFile);
			addSize2OTUfna($OTUfa,$hr);
			$cmd = "$usBin -sortbysize $OTUfa -fastaout $OTUfa.srt; \n";
			$cmd .= "$usBin -uchime3_denovo $OTUfa.srt -chimeras $chimOut -nonchimeras $lotus_tempDir/tmp1.fa -log $logDir/chimera_dn.log;";
			$cmd .= "rm $OTUfa.srt\n";
			$progUsed = "usearch uchime3";
			#die "$cmd\n";
		}

		#replace until bug is fixed
	}
	elsif ( $usearchVer >= 9 && !$VSused ) {
		$cmd = "$usBin -uchime2_denovo $OTUfa -abskew 16 -chimeras $chimOut -nonchimeras $lotus_tempDir/tmp1.fa -log $logDir/chimera_dn.log;";
		$progUsed = "usearch uchime2";
	}

	#die $cmd."\n";
	$cmd .= "\nrm $OTUfa\nmv -f $lotus_tempDir/tmp1.fa $OTUfa";
	if ( systemL($cmd) != 0 ) {
		printL( "uchime de novo failed// aborting\n", 1 );
	}
	if ( $usearchVer >= 9 ) {
		$citations .= "uchime2 chimera detection deNovo: Edgar, R.C. (2016), UCHIME2: Improved chimera detection for amplicon sequences, http://dx.doi.org/10.1101/074252..\n";
	}
	elsif ($VSused) {
		$citations .= "Vsearch chimera detection deNovo: [VSEARCH paper]\n";
	}
	else {
		$citations .= "uchime chimera detection deNovo: Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. 2011. UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27: 2194–200.\n";
	}
	
	my $endEntries = cntFastaEntrs($OTUfa);
	printL( frame("De novo chimera filter using $progUsed\nTotal removed ${OTU_prefix}s: (". ($iniEntries-$endEntries) ."/$iniEntries)"), 0 );


}

sub chimera_ref_rem($) {
    my ( $otusFA ) = @_;
    my $chimOut = "";
	my $iniEntries = cntFastaEntrs($otusFA);
	my $progUsed = "";
	my $outfile = $otusFA.".tmp.fa";
    if ($ClusterPipe == 1 ||$ClusterPipe == 7 || $UCHIME_REFDB eq ""  || !-f $UCHIME_REFDB ||  $noChimChk == 1 ){
		printL "No ref based chimera detection\n";
		return $chimOut;
	}
	#start ref chimera fasta
	if ($extendedLogs) {
		$chimOut = "$extendedLogD/chimeras_ref.fa";
	} else {$chimOut = "$lotus_tempDir/chimeras_ref.fa";}
	printL "Could not find fasta otu file $otusFA. Aborting..\n", 33 unless ( -s $otusFA );
	$cmd = "$VSBin -uchime_ref  $otusFA -db $UCHIME_REFDB -strand plus -chimeras $chimOut -nonchimeras $outfile -threads $uthreads -log $logDir/uchime_refdb.log;";
	$progUsed = "vsearch uchime_ref";
	if (!$useVsearch && $usearchVer >= 9 && !$VSused ) {
		$cmd = "$usBin -uchime2_ref  $otusFA -db $UCHIME_REFDB -mode balanced -strand plus -chimeras $chimOut -notmatched $outfile -threads $uthreads -log $logDir/uchime_refdb.log;";
		$progUsed = "usearch uchime2_ref";
	}
	#die $cmd."\n";
	if ( systemL($cmd) != 0 ) { printL( "Failed command:\n$cmd\n", 1 ); }
	if ( $usearchVer >= 9 ) {
		$citations .= "uchime2 chimera detection deNovo: Edgar, R.C. (2016), UCHIME2: Improved chimera detection for amplicon sequences, http://dx.doi.org/10.1101/074252..\n";
	}
	elsif ( $VSused == 0 ) {
		$citations .= "uchime reference based chimera detection: Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. 2011. UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27: 2194–200.\n";
	}
	else {
		$citations .= "Vsearch reference based chimera detection: \n";
	}
	systemL "rm $otusFA; mv $outfile $otusFA";
	#print $outfile."\n";

	
	my $endEntries = cntFastaEntrs($otusFA);
	printL( frame("Ref chimera filter using $progUsed\nTotal removed ${OTU_prefix}s: (". ($iniEntries-$endEntries) ."/$iniEntries)"), 0 );

	return $chimOut;
}

sub checkXtalk($ $) {
    my ( $otuFA, $otuM ) = @_;
	my %ret;
    if ( !$doXtalk ) { return \%ret; }
    if ( $usearchVer < 11 ) {
        printL "cannot check for cross-talk, as only implemented in usearch version > 11\n",83;
    }
    my $otuM1 = $otuM . ".noXref";
    #systemL "cp $otuM $otuM1;";
    my $cmd = "$usBin -otutab_xtalk $otuM -otutabout $otuM1 -report $logDir/crossTalk_analysis.txt;";
	my $uncrCnt = 0;
    #get the OTUs that are not in both tables
    systemL($cmd,"crosstalk",0);
    if ( -z "$logDir/crossTalk_analysis.txt" ) {
        printL "Cross talk unsuccessful, continue without\n","w";
    } else {
		my $tmp = `cut -f1 $otuM`; my @prevOTUs = split /\n/,$tmp;
		$tmp = `cut -f1 $otuM1`; my @newOTUs = split /\n/,$tmp;
		my %lookup = map { $_ => 1 } @newOTUs;
		for my $thing ( @prevOTUs ) {
			if (!exists($lookup{ $thing })){
				$ret{$thing} = 1 ;
				$uncrCnt++;
			}
		}
        $citations .= "CrossTalk ${OTU_prefix} removal: UNCROSS2: identification of cross-talk in 16S rRNA OTU tables. Robert C. Edgar . Bioarxiv (https://www.biorxiv.org/content/biorxiv/early/2018/08/27/400762.full.pdf)\n";
    }
	systemL "rm -f $otuM1;";
	printL(frame( "$uncrCnt ${OTU_prefix}'s removed with UNCROSS2\n"),0);
	return \%ret;
}

sub lulu{
	my %ret;
	return \%ret if (!$doLULU);
	my ($otuFas,$otuM) = @_;
	my $xtr="";
	if ($VSused){$xtr = "--iddef 1";
	} else {$xtr = "--strand plus";}
	my $idMin = 0.84;
	my $cmd = "$VSBin --usearch_global $otuFas --db $otuFas --self --id $idMin $xtr --userout $lotus_tempDir/lulu_match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10;";
	systemL $cmd,"vsearch for LULU",1;
	$cmd = "$Rscript $defRscriptOpt $LULUscript $lotus_tempDir/lulu_match_list.txt $otuM $logDir;";
	
	#die $cmd."\n";
	
	
	systemL $cmd,"LULU Rscript",1;
	my $rmLst = `cat $lotus_tempDir/lulu_match_list.txt.rm`;
	chomp $rmLst;
	my $LULUcnt=0;
	foreach my $x (split /\n/,$rmLst){
		$ret{$x} = 1;
		$LULUcnt ++;
	}
	
	printL(frame( "$LULUcnt ${OTU_prefix}'s removed with LULU\n"),0);

	$citations .= "LULU: Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188.\n";
	systemL "rm $lotus_tempDir/lulu_match_list.txt*";
	return \%ret;
}


#remove entries from OTU matrix
sub clean_otu_mat($ $ $ $) {
	my ( $OTUfa, $OTUrefFa, $OTUmFile, $OTUrefDBlnk, $ContaHR ) = @_;#$phiXCnt, $xtraContCnt
	my %Contam = %{$ContaHR};
	my @contaK = keys %Contam;
	my $phiXCnt = 0;my $xtraContCnt = 0;
	my %contT; my %contaCnt; my %contaCntRds; my %contRN; #new system of tracking contaminant..
	foreach my $conta (@contaK){
		if ($conta eq "phiX"){
			$phiXCnt += scalar(keys(%{$Contam{"phiX"}}));
		} else {
			$xtraContCnt += scalar(keys(%{$Contam{$conta}}));
		}
		foreach my $otus (keys %{$Contam{$conta}}){
			$contT{$otus} .= "$conta.";
		}

	}
	#search for local matchs of chimeras in clean reads
	#get headers of OTUs
	my $hr  = readFasta($OTUfa);
	my %hds = %{$hr};
	my $cnt = -1;

	#my @kkk =  keys %hds; die scalar @kkk ."\n";
	my %refHds;
	if (-e $OTUrefFa){
		$hr = readFasta($OTUrefFa);
		my %refHds = %{$hr};
	}
	my %ORDL   = %{$OTUrefDBlnk};
	my %ORDL2;
	my $chimRm      = 0; my $chimRdCnt   = 0;
	my %OTUcnt;
	my %OTUmat;
	
	my $Lreport = "";
	

	#die "$OTUmFile\n$OTUfa\n";
	open I, "<$OTUmFile";
	if ($extendedLogs) {
		open O2x, ">$extendedLogD/otu_mat.chim.txt" or printL "Failed to open $extendedLogD/otu_mat.chim.txt", 33;
	}
	while ( my $line = <I> ) {
		$cnt++;
		chomp $line;
		if ( $cnt == 0 ) {
			$OTUmat{head} = $line;
			if ($extendedLogs) { print O2x $line . "\n"; }
			next;
		}
		my @spl = split( /\t/, $line );
		my $ot  = shift @spl;
		#specific OTU read count
		if ( exists( $hds{$ot} ) || exists( $refHds{$ot} ) ){    #exists in non-chimeric set or ref DB set
			my $rdCnt =0;$rdCnt += $_ for @spl;
			#die "$rdCnt @spl\n";
			$OTUcnt{$ot} = $rdCnt;
			$OTUmat{$ot} = join( "\t", @spl );
		} else {
			#print $ot."\n";
			#die $line."\n";
			$chimRm++;    #don't include
			$contRN{$ot} = "Chim.".$chimRm;
			$chimRdCnt += $_ for @spl;
			if ($extendedLogs) { print O2x $contRN{$ot} ."\t". join( "\t", @spl ) . "\n"; }
		}
	}
	close I;

	if ($extendedLogs) { close O2x; }

	#print resorted OTU matrix & contaminant matrix..
	open OF,  ">$OTUfa";
	open OOFF,  ">$logDir/OTU.contaminants.fa";
	open OFR, ">$OTUfa.ref" if ( $REFflag );
	open O,   ">$OTUmFile";
	open O3x, ">$extendedLogD/otu_mat.contam.txt" or printL "Failed to open $extendedLogD/otu_mat.contam.txt", 33 if ($extendedLogs);
	my $emptyOTUcnt = 0;my $nonChim   = 0; my $nonChimRds = 0;
	print O $OTUmat{head} . "\n";
	print O3x $OTUmat{head} . "\n" if ($extendedLogs); 
	my @sorted_otus = ( sort { $OTUcnt{$b} <=> $OTUcnt{$a} } keys %OTUcnt);    #sort(keys(%OTUcnt));
	my $OTUcntd = 1;    # my $maxOTUdig = length (keys %OTUcnt)
	my $intoMatLast = "";
	foreach my $ot (@sorted_otus) {
		next unless (exists($OTUcnt{$ot})); #chimera
		if (exists($contT{$ot}) ){
			$contaCnt{$contT{$ot}} ++;
			$contaCntRds{$contT{$ot}} += $OTUcnt{$ot};
			$contRN{$ot} = ${OTU_prefix}.".".$contT{$ot}.$contaCnt{$contT{$ot}};
			if ($extendedLogs) { print O3x $contRN{$ot}  ."\t". $OTUmat{$ot} . "\n"; }
			print OOFF ">".$contRN{$ot}."\n$hds{$ot}\n"; #write fasta
			#retain them? Then include in final matrix, otu fasta etc -> will also mean that we get taxo annotation/tree etc for these
			if ($keepOfftargets){
				$intoMatLast .= $contRN{$ot}  ."\t". $OTUmat{$ot} . "\n"; 
				print OF ">".$contRN{$ot}."\n$hds{$ot}\n";
			}
			next;
		} 
		if ($OTUcnt{$ot} == 0){
			$emptyOTUcnt++;
			next;
		}
		#remaining OTUs that passed various checks are renamed & written here
		$nonChim++;$nonChimRds+=$OTUcnt{$ot};
		#$newOname = sprintf("%08d", $OTUcnt);
		my $newOname = ${OTU_prefix} . $OTUcntd;
		$OTUcntd++;
		print O $newOname . "\t" . $OTUmat{$ot} . "\n";
		if ( exists( $hds{$ot} ) ) {
			print OF ">" . $newOname . "\n$hds{$ot}\n";
		}
		elsif ( exists( $refHds{$ot} ) ) {
			print OFR ">" . $newOname . "\n$refHds{$ot}\n" if ( $REFflag );
			$ORDL2{$newOname} = $ORDL{$ot};

			#die ("new:".$ORDL2{$newOname}."  old: ".$ORDL{$ot}."\n");
		} else {
			printL "Fatal error, cannot find ${OTU_prefix} $ot\n", 87;
		}
	}
	print O $intoMatLast;
	close O;
	close OF;
	close OOFF;
	close OFR if ( $REFflag );
	close O3x if ($extendedLogs) ;
	#if ( $REFflag ) { unlink "$OTUfa.ref"; }


	#writeFasta(\%newOTUs,$OTUfa);
	$Lreport .= "Postfilter:\n"; 
	my $chimTag = "chimeric";#$chimTag .= "/ITSx" if ($doITSx);
	$Lreport .= "Chimeras: $chimRm $chimTag found ($chimRdCnt reads)\n" if ($chimRm);
	my $strTmp = "Contaminants: ";
	my $lcnt =0;
	foreach my $k (keys %contaCnt){
		$strTmp .=  "$k: $contaCnt{$k} ${OTU_prefix}'s removed ($contaCntRds{$k} reads); ";
		$lcnt += $contaCntRds{$k}+$contaCnt{$k};
	}
	$Lreport .= "Extended logs active, contaminant and chimeric matrix will be created.\n" if ($extendedLogs);
	$strTmp.="\n";
	$strTmp ="" if ($lcnt==0);
	$Lreport .= $strTmp;
	$Lreport.= "Removed mismapped OTUs ($emptyOTUcnt) with 0 matrix counts..\n" if ( $emptyOTUcnt > 0 );
	$Lreport .= "After filtering $nonChim $OTU_prefix ($nonChimRds reads) remaining in matrix.\n";

	if ( ( $nonChim - $emptyOTUcnt ) == 0 ) {
		#print "$nonChim - $emptyOTUcnt\n";
		printL "Empty ${OTU_prefix} matrix.. aborting LotuS run!\n", 87;
	}

	printL frame($Lreport),0;
	return \%ORDL2;
}

sub forceMerge_fq2fna($ $ $ $ $) {
    my ( $ifq1, $ifq2, $mfq, $sdms, $out ) = @_; #not merged1, not m2, merged, sdm sing, SEEDFNA
	my $seedCnt = 0;
    #print $mfq."\n";
    open O, ">", $out or die "Can't open seedFNA $out\n";
    my $hd    = "";
    my $lnCnt = 0;

    #check for rev adaptors
    my $check4endStr = 0;
    my $endRmvs      = "";
    if ( $remFromEnd ne "" ) {
        $check4endStr = 1;
        my @spl = split( /,/, $remFromEnd );
        $endRmvs .= $spl[0] . '.*$';
        for ( my $i = 1 ; $i < @spl ; $i++ ) {
            $endRmvs .= "|" . $spl[$i] . '.*$';
        }

        #print $remFromEnd."  C ".$endRmvs."\n";
        #$rawFileSrchStr1 = '.*1\.f[^\.]*q\.gz$';
    }
#main file with merged reads
	if (-e $mfq){
		open I, "<", $mfq or die "Couldn't open merge file $mfq\n";
		while ( my $line = <I> ) {
			chomp $line;
			if ( $line =~ m/^@/ && $lnCnt == 0 ) {
				$line =~ s/^@/>/;
				$line =~ s/.\d$// if ($line =~ m/[OTUASV]+_\d+\.\d+/);
				print O $line . "\n"; $seedCnt++;
				#print $line . "\n"; 
			}
			if ( $lnCnt == 1 ) {
				$line =~ s/$endRmvs// if ($check4endStr);
				print O $line . "\n";
			}
			$lnCnt++;
			$lnCnt = 0 if ( $lnCnt == 4 );
		}
		close I;
	}

    #merge two unmerged reads
    $lnCnt = 0;
	#die "XX  $mfq\n";

    #print $ifq1."\n";
	if (-s $ifq1){
		#die"XXXAS\n";
		open I,  "<", $ifq1;
		#open I2, "<", $ifq2;
		while ( my $line = <I> ) {
			#my $line2 = <I2>;
			chomp $line;
			#chomp $line2;
			if ( $line =~ m/^@/ && $lnCnt == 0 ) {
				$line =~ s/^@/>/;
				$line =~ s/.\d$//;

				#print $line."\n";
				print O $line . "\n"; $seedCnt++;
			}
			if ( $lnCnt == 1 ) {
				$line =~ s/$endRmvs// if ($check4endStr);
				print O $line . "\n";    #."NNNN".$line2."\n";
			}
			$lnCnt++;
			$lnCnt = 0 if ( $lnCnt == 4 );
		}
		close I;
		#close I2;
	}
    if ( -f $sdms ) { 
	#print "sdms $sdms\n";
		open I, "<", $sdms;
		$lnCnt = 0;
		while ( my $line = <I> ) {
			chomp $line;
			if ( $line =~ m/^@/ && $lnCnt == 0 ) {
				$line =~ s/$endRmvs// if ($check4endStr);
				$line =~ s/^@/>/;
				$line =~ s/.\d$//;
				print O $line . "\n";$seedCnt++;
			}
			if ( $lnCnt == 1 ) {
				print O $line . "\n";
			}
			$lnCnt++;
			$lnCnt = 0 if ( $lnCnt == 4 );
		}
		close I;
	}
    close O;
	printL frame("Found $seedCnt fasta seed sequences based on seed extension and read merging\n"),0;
}

sub mergeRds{
	return if ($numInput != 2 );
	my $single1  = "";my $single2  = "";
	#print "mergeOptions  $mergeOptions\n";
	my $mergeFile = "";
	if($mergeOptions eq ""){#non-sdm merge, use external programs
		my $key = "merged";
		#$mergeFile = "$lotus_tempDir/$key";
		$mergeFile   = "$lotus_tempDir/$key.extendedFrags.fastq";
		$single1 = "$lotus_tempDir/$key.notCombined_1.fastq";
		$single2 = "$lotus_tempDir/$key.notCombined_2.fastq";
		if (@preMergeSeedsFiles > 0 && -s $preMergeSeedsFiles[0] > 0 #check that file even exists..
				&& (-f $VSBin || -f $flashBin) #as well as the alignment programs
				&& $otuRefDB ne "ref_closed" && $onlyTaxRedo == 0  && $TaxOnly eq "0" ) #and otherwise also wanted step..
		{
			my $mergCmd = "";
			
			#switch to vsearch for lotus2
			if (-f $VSBin){
				$mergCmd = "$VSBin  -fastq_mergepairs $preMergeSeedsFiles[0]  -reverse  $preMergeSeedsFiles[1] --fastqout $mergeFile --fastqout_notmerged_fwd $single1 --fastqout_notmerged_rev $single2 --threads $BlastCores;";
				if ($VSused == 1){
					$citations .= "VSEARCH read pair merging: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584\n";
				}else{
					$citations .= "USEARCH read pair merging: R.C. Edgar (2010), Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19) 2460-2461\n";
				}
				#die "$mergCmd\n"; #find out what the unmerged reads are..
			} elsif ( -f $flashBin ) {
				my $flOvOption = "-m 10 -M $maxReadOverlap ";
				if ( $flashCustom ne "" ) {

					#$flOvOption = "-r $flashLength -s $flashSD";
					$flashCustom =~ s/\"//g;
					$flOvOption = " $flashCustom ";
				}
				$mergCmd = "$flashBin $flOvOption -o $key -d $lotus_tempDir -t $BlastCores "
				  . $preMergeSeedsFiles[0] . " " . $preMergeSeedsFiles[1];
				$mergCmd .= "cp $lotus_tempDir/$key.hist $logDir/FlashPairedSeedsMerges.hist";
				$citations .="Flash read pair merging: Magoc T, Salzberg SL. 2011. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27: 2957-63\n";
			} 

			#print $mergCmd."\n";
			if ( $exec == 0 ) {
				#die $mergCmd."\n".$flashCustom."\n";
				printL (frame("Merging ${OTU_prefix} seed paired reads"), 0);
				if ( !systemL($mergCmd) == 0 ) {
					printL( "Merge command failed: \n$mergCmd\n", 3 );
				}
			}
			$didMerge = 1;
		} elsif ( $onlyTaxRedo && $numInput == 2 ) {
			printL "Skipping paired end merging step\n", 0;
		}
	} elsif ($numInput == 2) { #so this is from the sdm merge
#		$single1  = $preMergeSeedsFiles[0]; $single2  = $preMergeSeedsFiles[1];
		$single1  =""; $single2  = "";
		$mergeFile = $preMergeSeedsFiles[0];
		$mergeFile =~ s/\.1\.fq$/\.merg\.fq/;
	}
	
	#needs to be run in any case for paired input, just to make sure all reads are correctly listed in $OTUSEED
	#die "$single1, $single2, $input,$mergeSeedsFilesSing[0], $OTUSEED\n";
	#TODO: add option to sdm to give out unmerged r1
	forceMerge_fq2fna( $single1, $single2, $mergeFile , $mergeSeedsFilesSing[0], $OTUSEED );
}

sub readTaxIn($ $ $ $ ) {
    my ( $inT, $LCAtax, $biomFmt, $calcHit2DB ) = @_;
    my %Taxo;    #
    my @lvls =
      ( "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species" );
    my @lvls2 = ( "k__", "p__", "c__", "o__", "f__", "g__", "s__" );
    my %hit2DB;
    my %hit2DBtax;
    my $taxStart   = 0;    #RDP case
    my $taxEnd     = 7;
    my $tarEntr    = 7;
    my $hit2DBEntr = 8;

    if ( $LCAtax == 1 ) {    #LCA case
        $taxStart = 1;
        $taxEnd   = 8;
        $tarEntr  = 0;
        unshift @lvls,  "OTU";
        unshift @lvls2, "";
    }
    open I, "<", $inT or die("$inT not found: $!\n");
    my $markLvl = -1;
    my $cnt     = -1;
    while ( my $line = <I> ) {
        chomp($line);
        $cnt++;
        my @spl = split( /\t/, $line );
        if ( $cnt == 1 && $markLvl == -1 ) {
            $markLvl = 0;
            $markLvl = 1 unless ( $spl[0] =~ m/[kpc]__/ );
        }
        for ( my $i = 0 ; $i < @spl ; $i++ ) {
            $spl[$i] = "?" if ( $spl[$i] eq "" );
        }
        next
          if ( $cnt == 0 )
          ; #|| $spl[0] eq "domain" || $spl[0] eq "Domain" || $spl[0] eq "OTU");
            #phylum unknown: ignore entry
            #if ($spl[1] eq "?" && $spl[2] eq "?"){next;}
        if ( @spl < 6 ) { die $line; }
        if ( $markLvl == 1 && $biomFmt == 1 ) {    #add k__ etc tags
            for ( my $i = 0 ; $i < @lvls2 ; $i++ ) {
                chomp $spl[$i];
                $spl[$i] = $lvls2[$i] . $spl[$i];
            }
        }
        if ( $biomFmt == 1 ) {                     #simple joining of levels
                #print ("@spl $taxStart .. $taxEnd\n");
            $Taxo{ $spl[$tarEntr] } =
              join( "\", \"", @spl[ $taxStart .. ( $taxEnd - 1 ) ] );
        }
        else {    #complex add up of levels
            my $foreRun = $spl[$taxStart];
            for ( my $i = $taxStart + 1 ; $i < $taxEnd ; $i++ ) {
                my $tartax = $foreRun . ";" . $spl[$i];
                if ( exists( $Taxo{ $lvls[$i] }{$tartax} ) ) {
                    push( @{ $Taxo{ $lvls[$i] }{$tartax} }, $spl[$tarEntr] );
                }
                else {
                    my @tmp = ( $spl[$tarEntr] );
                    $Taxo{ $lvls[$i] }{$tartax} = \@tmp;
                }
                $foreRun = $tartax;
            }
        }

        #prepare tag for greengenes hit
        if ($calcHit2DB) {
            if ( $spl[$hit2DBEntr] ne "?" ) {

                #previous one link.. no longer needed
                #$hit2DB{$spl[$tarEntr]} = $spl[$hit2DBEntr];}
                if ( exists( $hit2DB{ $spl[$hit2DBEntr] } ) ) {
                    push( @{ $hit2DB{ $spl[$hit2DBEntr] } }, $spl[$tarEntr] );
                }
                else {
                    my @tmp = ( $spl[$tarEntr] );
                    $hit2DB{ $spl[$hit2DBEntr] } = \@tmp;
                }

                #for bioim format, need to change the lvls
                if ( $markLvl == 1 ) {
                    for ( my $i = 0 ; $i < @lvls2 ; $i++ ) {
                        chomp $spl[$i];
                        $spl[$i] = $lvls2[$i] . $spl[$i];
                    }
                }
                $hit2DBtax{ $spl[$hit2DBEntr] } =
                  join( "\", \"", @spl[ $taxStart .. ( $taxEnd - 1 ) ] )
                  ;    #in biom format

            }
            else { $hit2DB{ $spl[$tarEntr] } = [ $spl[$tarEntr] ]; }
        }
    }
    close I;

    #create high  lvl matrix
    if ( $LCAtax == 1 ) {    #LCA case
        shift @lvls;
        shift @lvls2;
    }
    return ( \%Taxo, \@lvls, \%hit2DB, \%hit2DBtax );
}

sub biomFmt($ $ $ $ $) {
    my ( $otutab, $inT, $bioOut, $LCAtax, $bioTaxXHR ) = @_;

    if (0) {                 #$combineSamples==1){
        open O, ">$bioOut.not";
        print O "Samples are being combined; biom format is not supported by LotuS in this case\n";
        close O;
        printL "Samples are being combined; biom format is not supported by LotuS in this case\n";
        return;
    }
    my %bioTaxX = %{$bioTaxXHR};
    my %mapH    = %{$mapHref};
    my %combH   = %{$combHref};
    if ($hasCombiSmpls) {
        printL "Combined samples in lotus run.. attempting merge of metadata in .biom file\n","w";
    }
    my $otutab2 = $lotus_tempDir . "/OTUpTax_tmp.txt";
    my %OTUmat  = %{$otutab};
    my @avSmps  = sort( keys %OTUmat );
    my @avOTUs  = sort( keys %{ $OTUmat{ $avSmps[0] } } );

    #die "@avSmps\n";
    my @colNms = @{ $mapH{'#SampleID'} };

    #my @lvls = ("Domain","Phylum","Class","Order","Family","Genus","Species");

   # read tax normal in every case, hit2DB would be the new OTU IDs, if required
    my ( $hr1, $ar1, $hr2, $hr3 ) =
      readTaxIn( $inT, $LCAtax, 1, 0 );    #hr3 not used here
    my %Tax    = %{$hr1};
    my %hit2db = %{$hr2};                  #my @lvls = @{$ar1};

    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =  gmtime();
    my $biomo =
"{\n \"id\":null,\n \"format\": \"Biological Observation Matrix 0.9.1-dev\",
	 \"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",
	 \"type\": \"OTU table\",
	 \"generated_by\": \"$selfID\",
	 \"date\": \"" . sprintf( "%04d-%02d-%02d", $year + 1900, $mon, $mday ) . "T";

    #$biomo .= sprintf("%04d-%02d-%02d", $year,$mon,$mday);
    $biomo .= sprintf( "%02d:%02d:%02d", $hour, $min, $sec );
    $biomo .= "\",\n \"rows\":[\n";

#    {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}},

    my $cnt = 0;
    foreach my $id (@avOTUs) {
        $cnt++;

        #print $id."\n";
        my $taxstr;
        my $id2 = $id;

        if ( exists( $Tax{$id} ) ) {
            $taxstr = $Tax{$id};
        }
        elsif ( exists( $bioTaxX{$id} ) ) {

            #print "Heureka! $bioTaxX{$id}\n $id\n";
            $taxstr = $bioTaxX{$id};

            #if ($replOTUids4Hit){$id2 = $hit2db{$id};}
        }
        else {
            $taxstr =
"k__?\", \"p__?\", \"c__?\", \"o__?\", \"f__?\", \"g__?\", \"s__?";
        }
        if ( $cnt > 1 ) {
            $biomo .=
                ",\n            {\"id\":\""
              . $id2
              . "\", \"metadata\":{\"taxonomy\":[\""
              . $taxstr . "\"]}}";
        }
        else {
            $biomo .=
                "            {\"id\":\""
              . $id2
              . "\", \"metadata\":{\"taxonomy\":[\""
              . $taxstr . "\"]}}";
        }
    }
    $biomo .= "\n],\n \"columns\":[\n";

   #             {"id":"Sample1", "metadata":null},
   #   {"id":"Sample1", "metadata":{
   #                             "BarcodeSequence":"CGCTTATCGAGA",
   #                             "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
   #                             "BODY_SITE":"gut",
   #                             "Description":"human gut"}},
   #    {"id":"Sample2", "metadata":{

    $cnt = 0;

    #print $avSmps[0]."\n";
    foreach my $smpl (@avSmps) {
        my $prepMeta = "{\n                             ";    #"null";
                                                              #print $smpl."\n";
        my @curMeta  = @{ $mapH{ $combH{$smpl} } };

        my $cnt2 = -1;
        foreach my $cn (@colNms) {
            $cnt2++;
            $prepMeta .= ",\n                             " if ( $cnt2 > 0 );
            $prepMeta .= "\"" . $cn . "\":\"" . $curMeta[$cnt2] . "\"";
        }
        $prepMeta .= "}";
        $cnt++;
        if ( $cnt > 1 ) {
            $biomo .=
                ",\n            {\"id\":\""
              . $smpl
              . "\", \"metadata\":$prepMeta}";
        }
        else {    #cnt==1
            $biomo .=
              "            {\"id\":\"" . $smpl . "\", \"metadata\":$prepMeta}";
        }
    }
    my $rnum = scalar(@avOTUs);
    my $cnum = scalar(@avSmps);

    $biomo .=
"],\n\"matrix_type\": \"dense\",\n    \"matrix_element_type\": \"int\",\n    \"shape\": ["
      . $rnum . ","
      . $cnum . "],
		\"data\":  [";
    $cnt = 0;
    foreach my $otu (@avOTUs) {
        $cnt++;
        if ( $cnt == 1 ) {
            $biomo .= "[";
        }
        else {
            $biomo .= ",\n[";
        }
        my $cnt2 = 0;
        foreach my $smpl (@avSmps) {
            $cnt2++;
            if ( exists( $OTUmat{$smpl}{$otu} ) ) {
                if ( $cnt2 == 1 ) {
                    $biomo .= $OTUmat{$smpl}{$otu};
                }
                else {
                    $biomo .= "," . $OTUmat{$smpl}{$otu};
                }
            }
            else {
                #print "$smpl  $otu\n";
                if   ( $cnt2 == 1 ) { $biomo .= "0"; }
                else                { $biomo .= ",0"; }
            }
        }
        $biomo .= "]";
    }

    $biomo .= "]\n}";
    open O, ">", $bioOut;
    print O $biomo;
    close O;
    #printL frame("biom file created: $bioOut"), 0;

#die;
#my $cmd = "biom convert -i $otutab2 -o $bioOut --table-type \"otu table\" --process-obs-metadata taxonomy";
}

sub truePath($){
	my ($tmp) = @_;
	my $oriTmp = $tmp;
	return "$RealBin/$tmp" if (-f "$RealBin/$tmp");#relative path thing..
	$tmp = `which $tmp 2>/dev/null`;chomp($tmp);
	$tmp="" if ($tmp =~ /^which: no/);
	return $tmp;
}

sub readPaths_aligners($) {
    my ($inF) = @_;
    die("$inF does not point to a valid lotus configuration\n")
      unless ( -f $inF );
	  #die "$inF\n";
	$Rscript = truePath("Rscript");
    open I, "<", $inF;
    while ( my $line = <I> ) {
		chomp $line;
		next if ( $line =~ m/^#/ );
		next if ( length($line) < 5 );    #skip empty lines
		$line =~ s/\"//g;
		if ( $line =~ m/^usearch\s(\S+)/ ) {
			$usBin = truePath($1);
		} elsif ( $line =~ m/^dada2R\s+(\S+)/ ) {
			$dada2Scr = truePath $1;
		} elsif ( $line =~ m/^phyloLnk\s+(\S+)/ ) {
			$phyloLnk = truePath $1;
		} elsif ( $line =~ m/^LULUR\s+(\S+)/ ) {
			$LULUscript = truePath $1;
		} elsif ( $line =~ m/^vsearch\s+(\S+)/ ) {
			$VSBin = truePath($1);
			$VSBinOri = $VSBin;  #deactivate default vsearch again.. prob with chimera finder.
		}
		elsif ( $line =~ m/^LCA\s+(\S+)/ ) {
			$LCABin = truePath($1);
		}
		elsif ( $line =~ m/^multiRDPjar\s+(\S+)/ ) {
			$mjar = truePath($1);
		}
		elsif ( $line =~ m/^RDPjar\s+(\S+)/ ) {
			$rdpjar = truePath($1);
		}
		elsif ( $line =~ m/^blastn\s+(\S+)/ ) {
			$blastBin = truePath($1);
		}
		elsif ( $line =~ m/^makeBlastDB\s+(\S+)/ ) {
			$mkBldbBin = truePath($1);
		}
		elsif ( $line =~ m/^mafft\s+(\S+)/ ) {
			$mafftBin = truePath($1);
		}
		elsif ( $line =~ m/^clustalo\s+(\S+)/ ) {
			$clustaloBin = truePath($1);
		}
		elsif ( $line =~ m/^iqtree\s+(\S+)/ ) {
			$iqTreeBin = truePath($1);
		}
		elsif ( $line =~ m/^fasttree\s+(\S+)/ ) {
			$fasttreeBin = truePath($1);
		}
		elsif ( $line =~ m/^sdm\s+(\S+)/ ) {
			$sdmBin = truePath($1);
		}
		elsif ( $line =~ m/^flashBin\s+(\S+)/ ) {
			$flashBin = truePath($1);
		}
		elsif ( $line =~ m/^swarm\s+(\S+)/ ) {
			$swarmBin = truePath($1);
		}
		elsif ( $line =~ m/^cd-hit\s+(\S+)/ ) {
			$cdhitBin = truePath($1);
		}
		elsif ( $line =~ m/^dnaclust\s+(\S+)/ ) {
			$dnaclustBin = truePath($1);
		}
		elsif ( $line =~ m/^minimap2\s+(\S+)/ ) {
			$mini2Bin = truePath($1);
		}
		elsif ( $line =~ m/^lambda\s+(\S+)/ ) {
			$lambdaBin = truePath($1);
		}
		elsif ( $line =~ m/^lambda_index\s+(\S+)/ ) {
			$lambdaIdxBin = truePath($1);
		} elsif ( $line =~ m/^vxtractor\s+(\S+)/ ) {$vxtrBin = truePath($1);
		} elsif ( $line =~ m/^vxtractorHMMs\s+(\S+)/ ) {$hmmsVX = $1;
		} elsif ( $line =~ m/^itsx\s+(\S+)/ ) {$itsxBin = truePath($1);
		} elsif ( $line =~ m/^hmmsearch\s+(\S+)/ ) { $hmmsrchBin = truePath($1);
		} elsif ( $line =~ m/^CheckForUpdates\s+(\S+)/ ) {$checkForUpdates = $1;
		}
	}

    #check that usearch is execuatble
    if ( $doBlasting == 3 && !-f $usBin ) {
        printL "UTAX tax classification requested, but no usearch binary found at $usBin\nAborting..\n", 93;
    }
    if ( -f $usBin && !-X $usBin ) {
        printL "It seems like your usearch binary is not executable, attempting to change this (needs sufficient user rights)\n",0;
        printL "Failed \"chmod +x $usBin\"\n"  if ( systemL "chmod +x $usBin;" ), 33;
    }
	if ($LCABin eq ""){
		printL "Essential LCA program not found, please check that it is in the locations given in lOTUs.cfg\n",28;
	}
	if ($sdmBin eq ""){
		printL "Essential sdm program not found, please check that it is in the locations given in lOTUs.cfg\n",28;
	}
	#die;
}

#read database paths on hdd
sub readPaths {    #read tax databases and setup correct usage
    my ($inF,$defDBset) = @_;
    die("$inF does not point to a valid lotus configuration\n")
      unless ( -f $inF );
    open I, "<", $inF;
    my $TAX_REFDB_GG        = "";
    my $GGinfile            = "";
    my $TAX_REFDB_SLV       = "";
    my $SLVinfile           = "";
    my $TAX_REFDB_SLV_LSU   = "";
    my $SLVinfile_LSU       = "";
    my $UCHIME_REFssu       = "";
    my $UCHIME_REFlsu       = "";
    my $UCHIME_REFits       = "";
    my $UCHIME_REFits2      = "";
    my $UCHIME_REFits1      = "";
    my $TAX_RANK_ITS_UNITE  = "";
    my $TAX_REFDB_ITS_UNITE = "";
    my $TAX_UTAX_ITS        = "";
    my $TAX_UTAX_16S        = "";
    my $TAX_REFDB_16S_HITdb = "";
    my $TAX_RANK_16S_HITdb  = "";
    my $TAX_REFDB_16S_PR2   = "";
    my $TAX_RANK_16S_PR2    = "";
    my $TAX_REFDB_BT        = "";
    my $TAX_RANK_16S_BT     = "";

    while ( my $line = <I> ) {
        chomp $line;
        next if ( $line =~ m/^#/ );
        next if ( length($line) < 5 );    #skip empty lines
        $line =~ s/\"//g;

        #databases
        #$UCHIME_REFits $TAX_RANK_ITS_UNITE $TAX_REFDB_ITS_UNITE
        if ( $line =~ m/^UCHIME_REFDB\s+(\S+)/ ) {
            $UCHIME_REFssu = $1;
        }
        elsif ( $line =~ m/^UCHIME_REFDB_LSU\s+(\S+)/ ) {
            $UCHIME_REFlsu = $1;
        }
        elsif ( $line =~ m/^FAPROTAXDB\s+(\S+)/ ) {
            my $FaProTaxDBfile = $1;
			$FaProTax = loadFaProTax($FaProTaxDBfile);
        }
        elsif ( $line =~ m/^UCHIME_REFDB_ITS\s+(\S+)/ ) {
            $UCHIME_REFits = $1;
        }
        elsif ( $line =~ m/^UCHIME_REFDB_ITS1\s+(\S+)/ ) {
            $UCHIME_REFits1 = $1;
        }
        elsif ( $line =~ m/^UCHIME_REFDB_ITS2\s+(\S+)/ ) {
            $UCHIME_REFits2 = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_GG\s+(\S+)/ ) {
            $TAX_REFDB_GG = $1;
        }
        elsif ( $line =~ m/^TAX_RANK_GG\s+(\S+)/ ) {
            $GGinfile = $1;
        }
        elsif ($line =~ m/^TAX_REFDB_SSU_SLV\s+(\S+)/
            || $line =~ m/^TAX_REFDB_SLV\s+(\S+)/ )
        {
            $TAX_REFDB_SLV = $1;
        }
        elsif ($line =~ m/^TAX_RANK_SSU_SLV\s+(\S+)/
            || $line =~ m/^TAX_RANK_SLV\s+(\S+)/ )
        {
            $SLVinfile = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_LSU_SLV\s+(\S+)/ ) {
            $TAX_REFDB_SLV_LSU = $1;
        }
        elsif ( $line =~ m/^TAX_RANK_LSU_SLV\s+(\S+)/ ) {
            $SLVinfile_LSU = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_ITS_UNITE\s+(\S+)/ ) {
            $TAX_REFDB_ITS_UNITE = $1;
        }
        elsif ( $line =~ m/^TAX_RANK_ITS_UNITE\s+(\S+)/ ) {
            $TAX_RANK_ITS_UNITE = $1;

            #HITdb
        }
        elsif ( $line =~ m/^TAX_RANK_HITdb\s+(\S+)/ ) {
            $TAX_RANK_16S_HITdb = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_HITdb\s+(\S+)/ ) {
            $TAX_REFDB_16S_HITdb = $1;

            #PR2
        }
        elsif ( $line =~ m/^TAX_RANK_PR2\s+(\S+)/ ) {
            $TAX_RANK_16S_PR2 = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_PR2\s+(\S+)/ ) {
            $TAX_REFDB_16S_PR2 = $1;

            #bee tax
        }
        elsif ( $line =~ m/^TAX_RANK_BEE\s+(\S+)/ ) {
            $TAX_RANK_16S_BT = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_BEE\s+(\S+)/ ) {
            $TAX_REFDB_BT = $1;

            #UTAX
        }
        elsif ( $line =~ m/^TAX_REFDB_ITS_UTAX\s+(\S+)/ ) {
            $TAX_UTAX_ITS = $1;
        }
        elsif ( $line =~ m/^TAX_REFDB_SSU_UTAX\s+(\S+)/ ) {
            $TAX_UTAX_16S = $1;

            #ref DB for PhiX conatminations
        }
        elsif ( $line =~ m/^REFDB_PHIX\s+(\S+)/ ) {
            $CONT_REFDB_PHIX = $1;
        }
        elsif ( $line =~ m/^IdentityThresholds\s+(\S+)/ ) {
            my @spl = split( /,/, $1 );
            if ( scalar(@spl) == 6 ) {
                push( @spl, 0 );
            }
            if ( scalar(@spl) == 7 ) {
                for ( my $i = 0 ; $i < @spl ; $i++ ) {
                    if ( $spl[$i] > 100 || $spl[$i] < 0 ) {
                        print "Error reading \"IdentityThresholds\" from configuration file: every entry has to be <100 and >0. Found " . $spl[$i] . ".\n";
                        exit(22);
                    }
                    $idThr[$i] = $spl[$i];
                }
            }
            else {
                printL "Error reading \"IdentityThresholds\" from configuration file: has to contain 6 levels seperated by \",\".\n",23;
            }

            #die "@idThr\n";
        }
    }
    close I;
    if ( !-e $usBin && !-e $rdpjar && !-e $sdmBin ) {
        printL
"\nWARNING:: Several essential auxiliary programs are missing: you can always install these and configure into your lotus installation by excecuting ./autoinstall.pl\n\n",
          0;
    }

    $UCHIME_REFDB = $UCHIME_REFssu;
    if ( $ampliconType eq "LSU" ) {
        $UCHIME_REFDB = $UCHIME_REFlsu;
    }
    elsif ( $ampliconType eq "ITS" ) {
        $UCHIME_REFDB = $UCHIME_REFits;
    }
    elsif ( $ampliconType eq "ITS2" ) {
        $UCHIME_REFDB = $UCHIME_REFits2;
    }
    elsif ( $ampliconType eq "ITS1" ) {
        $UCHIME_REFDB = $UCHIME_REFits1;
    }
    if ( !-e $UCHIME_REFDB ) {
        printL "Requested DB for uchime ref at\n$UCHIME_REFDB\ndoes not exist; LotuS will run without reference based ${OTU_prefix} chimera checking.\n","w";
    }

    if (   ( $refDBwanted !~ m/SLV/ && $refDBwanted !~ m/PR2/ )
        && $ampliconType eq "LSU"
        && !-f $refDBwanted )
    {
        printL "-refDB \"GG\" does not contain taxonomic infomation for amplicon type $ampliconType \n switchung to \"SLV\".\n",
          0;
        $refDBwanted = "SLV";
    }

    #die $refDBwanted."$doBlasting\n";
    #printL "XX\nXX\n",0;
    if ( $doBlasting > 0 ) {
        #set up default paths
        if ( $doBlasting == 3 ) {
            $refDBwanted = "UTAX";
            if ( $usearchVer < 8.1 ) {
                printL "Your usearch version ($usearchVer) is not compatible with utax, please upgrade: http://www.drive5.com/usearch/\n",63;
            }
            if ( $usearchVer >= 9 ) {
                printL "Please use usearch ver 8 for now, if you want to use utax\n", 22;
            }
        }
        elsif ($defDBset) {
            if ( !-f $TAX_REFDB_GG && $refDBwanted eq "GG" ) {
                $refDBwanted = "SLV";
            }
        }
        for my $subrdbw ( split /,/, $refDBwanted ) {

            #print $subrdbw."\n";
            if ( -f $subrdbw ) {    #custom DB
                push @refDBname, "custom";
                push( @TAX_REFDB, $subrdbw );
                if ( $refDBwantedTaxo eq "" ) {
                    printL "Please provide a taxonomy file for custom ref DB @TAX_REFDB\n",45;
                }
                elsif ( !-f $refDBwantedTaxo ) {
                    printL "Taxonomy file for custom ref DB $refDBwantedTaxo does not exist\n",45;
                }
                push @TAX_RANKS, $refDBwantedTaxo;
                next;
            }
            $subrdbw = uc $subrdbw;
            if ( $subrdbw =~ m/^UTAX$/ ) {    #utax algo
                push @refDBname, "UTAX";
                if ( substr( $ampliconType, 0, 3 ) eq "ITS" ) {
                    push( @TAX_REFDB, $TAX_UTAX_ITS );
                }
                else {
                    push( @TAX_REFDB, $TAX_UTAX_16S );
                }
            }
            if ( $subrdbw =~ m/^UNITE$/
                || substr( $ampliconType, 0, 3 ) eq "ITS" )
            {                                 #overrides everything else
                #printL "Using UNITE ITS ref seq database.\n", 0;
                $ampliconType = "ITS";    #unless (substr($ampliconType,0,3) eq "ITS");
                push @TAX_REFDB, $TAX_REFDB_ITS_UNITE;
                push @TAX_RANKS, $TAX_RANK_ITS_UNITE;
                if ( !-f $TAX_REFDB[-1] || !-f $TAX_RANKS[-1] ) {
                    printL "Could not find UNITE ITS DB files at \n$TAX_REFDB[-1]\n$TAX_RANKS[-1]\nPlease check that these files exist.\n",55;
                }
                push @refDBname, "UNITE";
                $subrdbw = "";
                $citations .=
"UNITE ITS chimera DB - Nilsson et al. 2015. A comprehensive, automatically updated fungal ITS sequence dataset for reference-based chimera control in environmental sequencing efforts. Microbes and Environments. \n";
                $citations .=
"UNITE ITS taxonomical refDB - Koljalg, Urmas, et al. Towards a unified paradigm for sequence-based identification of fungi. Molecular Ecology 22.21 (2013): 5271-5277.\n";
            }
            if ( $subrdbw =~ m/^HITDB$/ ) {
                push @TAX_REFDB, $TAX_REFDB_16S_HITdb;
                push @TAX_RANKS, $TAX_RANK_16S_HITdb;
                push @refDBname, "HITdb";
                if ( !-f $TAX_REFDB[-1] || !-f $TAX_RANKS[-1] ) {
                    printL "Could not find HITdb files at \n$TAX_REFDB[-1]\n$TAX_RANKS[-1]\nPlease check that these files exist.\n",55;
                }
                $citations .= "HITdb gut specific tax database - J. Ritari, J. Salojärvi, L. Lahti, W. M. de Vos, Improved taxonomic assignment of human intestinal 16S rRNA sequences by a dedicated reference database. BMC Genomics. 16, 1056 (2015). \n";
            }
            if ( $subrdbw =~ m/^PR2$/ ) {
                push @TAX_REFDB, $TAX_REFDB_16S_PR2;
                push @TAX_RANKS, $TAX_RANK_16S_PR2;
                if ( !-f $TAX_REFDB[-1] || !-f $TAX_RANKS[-1] ) {
                    printL "Could not find PR2 files at \n$TAX_REFDB[-1]\n$TAX_RANKS[-1]\nPlease check that these files exist.\n",55;
                }
                $citations .= "PR2 LSU specific tax database - The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Laure Guillou, Dipankar Bachar, Stephane Audic, David Bass, Cedric Berney, Lucie Bittner, Christophe Boutte, Gaetan Burgaud, Colomban de Vargas, Johan Decelle, Javier del Campo, et al. Nucleic Acids Res (2013) Volume 41 Issue D1: Pp. D597-D604. \n";
                push @refDBname, "PR2";
            }
            if ( $subrdbw =~ m/^GG$/ ) {
                if ( !-f $TAX_REFDB_GG || !-f $GGinfile ) {
                    printL "Could not find greengenes DB files at \n$TAX_REFDB_GG\n$GGinfile\nPlease check that these files exist.\n",55;
                }

#if (-f $TAX_REFDB_SLV && -f $SLVinfile){				#printL "Found both SILVA and greengenes databases. Using Greengenes for this run.\n",0;			}
                push @refDBname, "greengenes";
                push @TAX_REFDB, $TAX_REFDB_GG;
                push @TAX_RANKS, $GGinfile;
                $citations .=
"greengenes 16S database - McDonald D, Price MN, Goodrich J, Nawrocki EP, DeSantis TZ, Probst A, Andersen GL, Knight R, Hugenholtz P. 2012. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. ISME J 6: 610–8.\n";
            }

            if ( $subrdbw =~ m/^BEETAX$/ ) {
                if ( !-f $TAX_REFDB_GG || !-f $GGinfile ) {
                    printL"Could not find greengenes DB files at \n$TAX_REFDB_GG\n$GGinfile\nPlease check that these files exist.\n",55;
                }

#if (-f $TAX_REFDB_SLV && -f $SLVinfile){				#printL "Found both SILVA and greengenes databases. Using Greengenes for this run.\n",0;			}
                push @refDBname, "BeeTax";
                push @TAX_REFDB, $TAX_REFDB_BT;
                push @TAX_RANKS, $TAX_RANK_16S_BT;
                $citations .= "Bee specific reference database: Jones, Julia C et al. 2017. “Gut Microbiota Composition Is Associated with Environmental Landscape in Honey Bees.” Ecology and Evolution (October): 1–11.\n";
            }
            if ( $subrdbw =~ m/^SLV$/ ) {
                push @refDBname, "SILVA";
                if ( $ampliconType eq "SSU" ) {
                    printL "Using Silva SSU ref seq database.\n", 0;
                    push @TAX_REFDB, $TAX_REFDB_SLV;
                    push @TAX_RANKS, $SLVinfile;
                }
                elsif ( $ampliconType eq "LSU" ) {
                    printL "Using Silva LSU ref seq database.\n", 0;
                    push @TAX_REFDB, $TAX_REFDB_SLV_LSU;
                    push @TAX_RANKS, $SLVinfile_LSU;
                }
                if ( !-f $TAX_REFDB[-1] || !-f $TAX_RANKS[-1] ) {
                    printL "Could not find Silva DB files at \n$TAX_REFDB[-1]\n$TAX_RANKS[-1]\nPlease check that these files exist.\n",55;
                }
                $citations .= "SILVA 16S/18S database - Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glockner FO (2014) The SILVA and \"All-species Living Tree Project (LTP)\" taxonomic frameworks. Nucleic Acid Res. 42:D643-D648 \n";

            }
            if ( @refDBname == 0 ) {
                printL"Found no valid path to a ref database \"$subrdbw\" or could not identify term. No reference based ${OTU_prefix} picking.Aborting LotuS\n",3;
            }

        }
    }

    #die "@TAX_REFDB\n";
    #die $ampliconType.substr($ampliconType,0,3)."\n";
    if ( substr( $ampliconType, 0, 3 ) eq "ITS" ) { $doRDPing = 0; }

    if ( !-f $sdmBin ) {
        printL "Could not find sdm binary at\n\"$sdmBin\"\n.Aborting..\n", 3;
    }
    if ( !-f $usBin ) {
        printL "Could not find usearch binary at\n\"$usBin\"\n.Aborting..\n", 3;
    }
    if ( !-f $VSBin ||  $useVsearch == 0) {
        $VSBin  = $usBin;
        $VSused = 0;
        printL ("Could not find vsearch binaries at\n$VSBin\n, switching to usearch binaries instead\n","w") if ( !-f $VSBin );
    }
    else {
        $citations .= "VSEARCH 1.13 (chimera de novo / ref; ${OTU_prefix} alignments): Rognes T (2015) https://github.com/torognes/vsearch\n";
    }
    if ( !-f $mjar && !-f $rdpjar && $doRDPing > 0 ) {
        printL "Could not find rdp jar at\n\"$mjar\"\nor\n\"$rdpjar\"\n.Aborting..\n", 3;
    }

}


sub announce_options{
	if ( !$onlyTaxRedo && $TaxOnly eq "0" ) {
	printL "$clustMode sequence clustering with $clusteringNameStr into ${OTU_prefix}'s\n",0;
		if ( $ClusterPipe == 7 ) {
			die "Incorrect dada2 script defined $dada2Scr" unless (-f $dada2Scr);
			die "Incorrect R installation (can't find Rscript)" unless (-f $Rscript);
		}
		#if   ($sdmDerepDo) { printL ("Running fast LotuS mode..\n",0); }else{ printL ("Running low-mem LotuS mode..\n",0); }
	} elsif ($TaxOnly eq "1") {    #prep for redoing tax, save previous tax somehwere
		printL "Re-Running only tax assignments, no de novo clustering\n", 0;
		my $k       = 0;
		my $newLDir = "$outdir/prevLtsTax_$k";
		while ( -d $newLDir ) { $k++; $newLDir = "$outdir/prevLtsTax_$k"; }
		printL frame("Saving previous Tax to $newLDir"), 0;
		systemL "mkdir -p $newLDir/LotuSLogS/;";
		systemL "mv $highLvlDir $FunctOutDir $RDP_hierFile $SIM_hierFile $outdir/hierachy_cnt.tax $outdir/cnadjusted_hierachy_cnt.tax $newLDir/;";
		systemL "mv $outdir/LotuSLogS/* $newLDir/LotuSLogS/;";
		systemL("mkdir -p $logDir;") unless ( -d $logDir );

		if ($extendedLogs) {
			systemL("mkdir -p $extendedLogD;") unless ( -d $extendedLogs );
		}

	}
	printL "------------ I/O configuration --------------\n", 0;
	printL( "Input=   $input\nOutput=  $outdir\n", 0 );    #InputFileNum=$numInput\n
	if ( $barcodefile ne "" ) {
		printL("Barcodes= $barcodefile\n",0);
	}
	printL "TempDir= $lotus_tempDir\n",                                     0;
	printL "------------ Configuration LotuS --------------\n", 0;
	printL( "Sequencing platform=$platform\nAmpliconType=$ampliconType\n", 0 );
	if ( $ClusterPipe == 2 ) {
		printL "Swarm inter-cluster distance=$swarmClus_d\n", 0;
	}
	else {
		printL "OTU id=$id_OTU\n", 0;
	}
	printL "min unique read abundance=" . ($dereplicate_minsize) . "\n", 0;
	if ( $noChimChk == 1 || $noChimChk == 3 ) {
		printL "No RefDB based Chimera checking\n", 0;
	}
	elsif ( $noChimChk == 0 || $noChimChk == 2 ) {
		printL "UCHIME_REFDB, ABSKEW=$UCHIME_REFDB, $chimera_absskew\nOTU, Chimera prefix=${OTU_prefix}, $chimera_prefix\n",
		  0;
	}
	if ( $noChimChk == 1 || $noChimChk == 2 ) {
		printL "No deNovo Chimera checking\n", 0;
	}
	if ($mergePreCluster){
		printL "Precluster read merging active :)\n", 0;
	}
	if ( $doBlasting == 1 ) {
		printL "Similarity search with Blast\n", 0;
		if ( !-e $blastBin ) {
			printL "Can't find blast binary at $blastBin\n", 97;
		}
	} elsif ( $doBlasting == 2 ) {
		printL "Similarity search with Lambda\n", 0;
		if ( !-e $lambdaBin ) {
			printL "Can't find LAMBDA binary at $lambdaBin\n", 96;
		}
		if ( !-e $lambdaIdxBin ) {
			printL "Can't find valid labmda indexer executable at $lambdaIdxBin\n",  98;
		}
	}elsif ( $doBlasting == 4 ) {
		printL "Similarity search with VSEARCH\n", 0;
		if ( !-e $VSBinOri ) {
			printL "Can't find VSEARCH binary at $VSBinOri\n", 96;
		}
	}elsif ( $doBlasting == 5 ) {
		printL "Similarity search with USEARCH\n", 0;
		if ( !-e $usBin ) {
			printL "Can't find VSEARCH binary at $usBin\n", 96;
		}
	}
	unless ( $doBlasting < 1 ) {
		printL "ReferenceDatabase=@refDBname\nRefDB location=@TAX_REFDB\n", 0;
		if ( !-e $TAX_REFDB[0] ) {
			printL "RefDB does not exist at loction. Aborting..\n", 103;
		}
	}

	printL "TaxonomicGroup=$organism\n", 0;
	if ( $ClusterPipe == 0 ) {
		printL("PCTID_ERR=$id_OTU_noise\n",0);
	}
	if ($keepUnclassified){
	printL "keeping taxonomic unclassified ${OTU_prefix}'s in matrix\n" ,0;
	} else {
		printL "removing taxonomic unclassified ${OTU_prefix}'s from matrix\n",0 ;
	}

	if ( $custContamCheckDB ne "" ) {
		printL "Custom DB for off-targets: $custContamCheckDB\n", 0;
		if ($keepOfftargets){
			printL "keeping off-target ${OTU_prefix}'s in matrix\n",0;
		} else {
			printL "removing off-target ${OTU_prefix}'s from matrix\n" ,0;
		}
	}
	printL "--------------------------------------------\n", 0;
}


sub getSimBasedTax{
	$doBlasting_pre = lc $doBlasting_pre;
	if ( $doBlasting_pre eq "1" || $doBlasting_pre eq "blast" ) {
		$doBlasting = 1;
	} elsif ( $doBlasting_pre eq "2" || $doBlasting_pre eq "lambda" ) {
		$doBlasting = 2;
	} elsif ( $doBlasting_pre eq "3" || $doBlasting_pre eq "utax" ) {
		$doBlasting = 3;
	} elsif ( $doBlasting_pre eq "4" || $doBlasting_pre eq "vsearch" ) {
		$doBlasting = 4;
	} elsif ( $doBlasting_pre eq "5" || $doBlasting_pre eq "usearch" ) {
		$doBlasting = 5;
	} elsif ( $doBlasting_pre eq "0" ) {
		$doBlasting = 0;
	}
}

sub checkBlastAvaila() {
    if ( !-f $lambdaBin && ( $doBlasting == 2 ) ) {
        printL"Requested LAMBDA based similarity tax annotation; no LAMBDA binary found at $lambdaBin\n",55;
    } elsif ( !-f $blastBin && ( $doBlasting == 1 ) ) {
        printL "Requested blastn based similarity tax annotation; no blastn binary found at $blastBin\n",55;
    } elsif ( !-f $VSBinOri && ( $doBlasting == 4 ) ) {
        printL "Requested VSEARCH based similarity tax annotation; no vsearch binary found at $VSBinOri\n",55;
    } elsif ( !-f $usBin && ( $doBlasting == 5 ) ) {
        printL "Requested USEARCH based similarity tax annotation; no usearch binary found at $usBin\n",55;
    }
}

sub prepLtsOptions{
	
	$lotus_tempDir = $outdir . "/tmpFiles/" if (!defined($lotus_tempDir) || $lotus_tempDir eq "");
	$BlastCores = $uthreads;

	if ($chimCnt){
		$chimCnt = "T";
	} else {
		$chimCnt = "F";
	}

	#still undocumented options: VsearchChimera removePhiX

	#uc/lc some vars
	$platform = lc($platform);
	if ($greengAnno) {
		$pseudoRefOTU   = 1;
		$refDBwanted    = "GG";
		$doBlasting_pre = "2";
		print "Greengenes ID as species ID requested (for integration with software that requires greengenes ID).\nUsing Lambda OTU similarity search and greengenes ref DB\n";
	}    #
	if ( -f $refDBwanted ) {
	}
	elsif ( $refDBwanted =~ m/^GG|SLV|UNITE|HITDB|PR2|BEETAX$/i ) {
		#$refDBwanted = uc($refDBwanted);
	}elsif ( $refDBwanted ne "" ) {
		print "$refDBwanted not recognized. Aborting..\n";
		exit(5);
	}
	if ( $platform eq "hiseq" ) {
		$linesPerFile = 8000000;
	}

	if ( $platform eq "pacbio" ) {
		$dereplicate_minsize_def = 0;
		$ClusterPipe_pre = "CDHIT" if ($ClusterPipe_pre eq "?");
	}
	if ( $dereplicate_minsize !~ m/\D/ && $dereplicate_minsize == -1 ){
		$dereplicate_minsize = $dereplicate_minsize_def ;
	}

	#die $refDBwanted."\n";
	$ampliconType    = uc($ampliconType);
	$organism        = lc($organism);
	$ClusterPipe_pre = uc($ClusterPipe_pre);
	$otuRefDB        = lc $otuRefDB;

	getSimBasedTax();

	if ( $TaxOnly eq "0" && !defined($input) && !defined($outdir) && !defined($mapFile) ) { 
		defined($input) or usage("-i option (input dir/files) is required\n");
		if ( $input =~ m/\*/ ) {
			usage ("\"*\" not supported in input command. Please see documentation on how to set up the mapping file for several input files.");
		}
		defined($outdir) or usage("-o option (output dir) is required\n");
		if ( !defined($mapFile) ) {
			usage("-m missing option (mapping file)\n");
		}
	}
	if ( !defined($sdmOpt) && $TaxOnly eq "0" ) {
		printL("sdm options not set\n","w");
	}
	if ($TaxOnly ne "0" && -f $TaxOnly){
		$TaxOnly =~ m/.*\/(.*)$/;$outdir = $1."/tmp/";;
		printL "Using tmp outdir $outdir\n";
		$rmOutDir=1;
	}
	
	#----- part 2 ----------
		
	my $defDBset = 0;
	if ( $refDBwanted eq "" ) {
		$refDBwanted = "SLV"; #switch to default SIVLA..
		$defDBset    = 1;
	}
	
	#read primary paths to binaries
	readPaths_aligners($lotusCfg);
	#version checks
	#die "$mini2Bin\n";
	my $LCAver = `$LCABin -v`;
	chomp ($LCAver); 
	if ($LCAver !~ m/[\d\.]/ || $LCAver < 0.21){
		printL "LCA is ver $LCAver. Require at least ver 0.21. $LCABin\n",823;
	}
	my $usvstr = `$usBin --version`;
	$usvstr =~ m/usearch v(\d+\.\d+)\.(\d+)/;
	$usearchVer = $1;$usearchsubV = $2;
	
	$usvstr = `$mini2Bin --version`; $usvstr =~ m/(\d\.\d+)/;
	my $miniVer = $1;
	if ($miniVer < 2.17) {die "minimap2 version (found $miniVer at $mini2Bin) too low, expected at least 2.17\n";}

	#if ($usearchVer == 9){printL "Usearch ver 9 currently not supported, please install ver 8.\n",39;}
	if ( $usearchVer > 11 ) {
		printL "Usearch ver $usearchVer is not supported.\n", 55;
	}elsif ( $usearchVer >= 8 && $usearchVer < 9 ) {
		printL"Usearch ver 8 is outdated, it is recommended to install ver 9.\nDownload from http://drive5.com/ and execute \n\"./autoInstall.pl -link_usearch [path to usearch9]\"\n",0;
	}elsif ( $usearchVer < 8 ) {
		printL "Usearch ver 7 is outdated, it is recommended to install ver 9.\nDownload from http://drive5.com/ and execute \n\"./autoInstall.pl -link_usearch [path to usearch9]\"\n",0;
	}
	if ($doLULU && $LULUscript eq "" || !-f $LULUscript){
		printL "Requested LULU matrix corrections, can't find lulu R script at \"$LULUscript\"\n",24;
	}
	if ($Rscript eq "" || !-e $Rscript){
		my $warnLvl = "w";
		$warnLvl=291 if ($doLULU || $ClusterPipe == 7);
        printL "Rscript can't be found, but is required for dada2, LULU and phyloseq\n",$warnLvl;
	}

	#die "$usearchVer $usearchsubV\n";
	if ( $doXtalk == -1 ) {
		if ( $usearchVer >= 11 ) {
			$doXtalk = 0;    #deactivate since useless if not 64-bit uparse
		}  else { $doXtalk = 0;}
	}

	#die "$VSBin\n";
	#which default aligner?
	if ( $doBlasting == -1 ) {
		if ($defDBset) {
			$doBlasting = 0;
		}   else {
			#check which aligner is installed
			my $defaultTaxAlig = 4;
			if ( !-f $VSBin ) {
				if ( -f $lambdaBin ) { $defaultTaxAlig = 2; }
				if ( -f $usBin ) { $defaultTaxAlig = 5; }
				if ( -f $blastBin ) { $defaultTaxAlig = 1; }
				elsif ( !-f $blastBin && $refDBwanted ne "" ) {
					printL "Requested similarity search ($refDBwanted), but no suitable aligner (blast, lambda) binaries found!\n",50;
				}
			}
			if ( substr( $ampliconType, 0, 3 ) eq "ITS" ) {
				$doBlasting  = $defaultTaxAlig;
				$refDBwanted = "UNITE";
			}
			if ( $refDBwanted ne "" ) {
				$doBlasting = $defaultTaxAlig;    #set default to lambda
				printL "RefDB $refDBwanted requested. Setting similarity based search to default Blast option to search $refDBwanted.\n",0;
			}
		}
	}
	elsif ( $doBlasting == 0 && $refDBwanted ne "" ) {
		printL "RefDB $refDBwanted requested, but -taxAligner set to \"0\": therefore RDP classification of reads will be done\n",0;
	}

	readPaths($lotusCfg,$defDBset);

	#die "$refDBwanted\n";
	if ( $doBlasting < 1 ) {
		$doRDPing = 1;
	} else {    #LCA only
		$doRDPing = 0;
	}

	if ( $ClusterPipe_pre eq "CD-HIT" || $ClusterPipe_pre eq "CDHIT" || $ClusterPipe_pre eq "3" ) {
		$ClusterPipe = 3; $clusteringNameStr = "CD-HIT";
		if ( !-e $cdhitBin ) {
			printL "No valid CD-Hit binary found at $cdhitBin\n", 88;
		}
	}elsif ( $ClusterPipe_pre eq "UPARSE" || $ClusterPipe_pre eq "1" ) {
		$ClusterPipe = 1;$clusteringNameStr = "UPARSE";
	}elsif ( $ClusterPipe_pre eq "UNOISE" || $ClusterPipe_pre eq "UNOISE3" || $ClusterPipe_pre eq "6" ) {
		$ClusterPipe = 6; $clusteringNameStr = "UNOISE3";
		${OTU_prefix} = "Zotu";
	}elsif ( $ClusterPipe_pre eq "DADA2" || $ClusterPipe_pre eq "7" ) {
		$ClusterPipe = 7; $clusteringNameStr = "DADA2";
		${OTU_prefix} = "ASV";
	}elsif ( $ClusterPipe_pre eq "SWARM" || $ClusterPipe_pre eq "2" ) {
		$ClusterPipe = 2; $clusteringNameStr = "SWARM";
		if ( !-e $swarmBin ) {
			printL "No valid swarm binary found at $swarmBin\n", 88;
		}
	}elsif ( $ClusterPipe_pre eq "DNACLUST" || $ClusterPipe_pre eq "4" ) {
		$ClusterPipe = 4; $clusteringNameStr = "DNACLUST";
		if ( !-e $dnaclustBin ) {
			printL "No valid DNA clust binary found at $dnaclustBin\n", 88;
		}
	} else {#default pipeline
		$ClusterPipe = 1; $clusteringNameStr = "UPARSE";
		#${OTU_prefix} = "Zotu";
	}
	if ( $platform eq "pacbio" && $ClusterPipe != 3 ) {
		printL("CD-HIT clustering is strongly recommended with PacBio reads (unless you know what you are doing).","w");
	}



	$REFflag = $otuRefDB eq "ref_closed" || $otuRefDB eq "ref_open";
	if ($REFflag) {
		if ( $ClusterPipe != 0  ) {
			printL( "$clusteringNameStr does not support ref DB out clustering\nUse dnaclust instead\n",12);
		}
		if ( $refDBwanted eq "" ) {
			printL "You selected ref based ${OTU_prefix} building, please set -refDB to \"SLV\", \"GG\", \"HITdb\", \"PR2\" or a custom fasta file.\n",22;
		}

	}

	checkBlastAvaila();

	#"LotuS 1.281"
	$selfID =~ m/LotuS (\d\.\d+)/;
	my $sdmVer = checkLtsVer($1);
	if ( $sdmVer < $curSdmV ) {
		printL "Installed sdm version ($sdmVer < $curSdmV) seems to be outdated, please check on \n    lotus2.earlham.ac.uk\nfor the most recent version. Make sure the sdm path in '$lotusCfg' points to the correct sdm binary\n","w";
	}
	$swarmClus_d = int($swarmClus_d);
	if ( $swarmClus_d < 1 ) {
		printL "Please provide as swarm distance an int > 0\n", 29;
	}

	#change automatically lambdaindexer based on mem installed in machine
	if ( 0 && (`cat /proc/meminfo |  grep "MemTotal" | awk '{print \$2}'`) < 16524336 ) {
		$lowMemLambI = 1;
		printL "Less than 16GB Ram detected, switching to mem-friendly workflow\n",0;
	}

	#die();

	if ( !-e $sdmOpt && $TaxOnly eq "0") {
		if ( $sdmOpt eq "" ) {
			$sdmOpt = "configs/sdm_miSeq.txt";
			printL "No sdm Option specified, using standard miSeq sdm options", 0;
		}
		if (!-e $sdmOpt){
			printL "Could not find sdm options file (specified via \"-s $sdmOpt\"). Please make sure this is available.\n Aborting run..\n",33;
		}
	}

	#die $LCABin."\n";
	if ( substr( $ampliconType, 0, 3 ) eq "ITS" ) {
		if ( $organism ne "fungi" && $organism ne "eukaryote" ) {
			$organism = "eukaryote";
			printL	"Setting \"-tax_group\" to \"eukaryote\" as only eukaryote and fungi are supported options for $ampliconType.\n","w";
		}

		if (   $doBlasting == 0
			|| ( !-f $blastBin && !-f $lambdaBin )  || @TAX_REFDB == 0 || !-f $TAX_REFDB[0] )
		{
			my $failedBlastITS = "ITS region was chosen as target; this requires a similarity based taxnomic annotation and excludes RDP tax annotation.\n";
			$failedBlastITS .= "Blast similarity based annotation is not possible due to: ";
			if ( $doBlasting == 0 ) {
				$failedBlastITS .= "Similarity search was not explicitly activated (please use option \"-taxAligner usearch\" or vsearch,lambda,blast).";
			}
			elsif ( !-f $blastBin || !-f $lambdaBin ) {
				$failedBlastITS .=
				  "Neither Lambda nor Blast binary being specified correctly";
			}
			elsif ( @TAX_REFDB == 0 || !-f $TAX_REFDB[0] ) {
				$failedBlastITS .= "Reference DB does not exist ($TAX_REFDB[0]).\n";
			}
			$failedBlastITS .= "\nTherefore LotuS had to abort..\n";
			printL $failedBlastITS, 87;
		}

	#$doBlasting && (-f $blastBin || -f $lambdaBin) && $TAX_REFDB ne "" && -f $TAX_REFDB)
	}
	if ( $noChimChk < 0 || $noChimChk > 3 ) {
		printL "option \"-deactivateChimeraCheck\" has to be between 0 and 3\n", 45;
	}
	if ( $saveDemulti < 0 || $saveDemulti > 3 ) {
		printL "option \"-saveDemultiplex\" has to be between 0 and 3\n", 46;
	}
}



sub help {
    my $opt_name = shift;
    my $exit_code = $opt_name eq 'help' ? 0 : $opt_name;
	my $ver = 0;


# HELP BLOCK LOTUS2

#HELPPARSEBEGIN
my $usage_string = "perl lotus2.pl -i <input fasta|fastq|dir> -o <output_dir> -m/-map <mapping_file>";
my $usage_example = "perl lotus2.pl -i testDir/ -o testOut/ -m testDir/mymap.txt -CL dada2";
my $basic_heading = "Basic Options";
my %basic_options = (
  '-i <file|dir>', 'In case that fastqFile or fnaFile and qualFile were specified in the mapping file, this has to be the directory with input files',
  '-o <dir>', 'Warning: The output directory will be completely removed at the beginning of the LotuS2 run. Please ensure this is a new directory or contains only disposable data!',
  '-m|-map <file>', 'Mapping file'
);

my $further_heading = "Further Options";
my %further_options = (
  '-q <file>', 'input qual file (not defined in case of fastq or input directory)', 
  '-barcode|-MID <file>', 'Filepath to fastq formated file with barcodes (this is a processed mi/hiSeq format). The complementary option in a mapping file would be the column "MIDfqFile"',
  '-s <file>', 'SDM option file, defaults to \"sdm_options.txt\" in current dir',
  '-c <file>', 'LotuS.cfg, config file with program paths',
  '-p <454/miSeq/hiSeq/PacBio>', 'sequencing platform: PacBio, 454, miSeq or hiSeq',
  '-t|-threads <num>', 'number of threads to be used',
  '-tmp|-tmpDir <dir>', 'temporary directory used to save intermediate results'
);

my $workflow_heading = "Workflow Options";
my %workflow_options = (
  '-verbosity <0-3>', 'Level of verbosity from printing all program calls and program output (3) to not even printing errors (0). Default: 1',
  '-saveDemultiplex <0|1|2|3>', '(1) Saves all demultiplexed reads (unfiltered) in the [outputdir]/demultiplexed folder, for easier data upload. (2) Only saves quality filtered demultiplexed reads and continues LotuS2 run subsequently. (3) Saves demultiplexed file into a single fq, saving sample ID in fastq/a header. (0) No demultiplexed reads are saved. (Default: 0)',
  '-taxOnly <file>', 'Skip most of the lotus pipeline and only run a taxonomic classification on a fasta file. E.g. ./lotus2.pl -taxOnly <some16S.fna> -refDB SLV',
  '-redoTaxOnly <0|1>', '(1) Only redo the taxonomic assignments (useful for replacing a DB used on a finished lotus run). (0) Normal lotus run. (Default: 0)',
  '-keepOfftargets <0|1>', '(0)?!?: keep offtarget hits against offtargetDB in output fasta and otu matrix, default 0',
  '-keepTmpFiles <0|1>', '(1) save extra tmp files like chimeric OTUs or the raw blast output in extra dir. (0) do not save these. (Default: 0)',
  '-keepUnclassfied <0|1>', '(1) Includes unclassified OTUs (i.e. no match in RDP/Blast database) in OTU and taxa abundance matrix calculations. (0) does not take these OTUs into account. (Default: 0)',
  '-tolerateCorruptFq <0|1>', '(1) Continue reading fastq files, even if single entries are incomplete (e.g. half of qual values missing). (0) Abort lotus run, if fastq file is corrupt. (Default: 0)'
);

my $taxonomy_heading = "Taxonomy Options";
my %taxonomy_options = (
  '-refDB <SLV|GG|HITdb|PR2|UNITE|beetax>', '(SLV) Silva LSU (23/28S) or SSU (16/18S), (GG greengenes (only SSU available), (HITdb) (SSU, human gut specific), (PR2) LSU spezialized on Ocean environmentas, (UNITE) ITS fungi specific, (beetax) bee gut specific database and tax names. \nDecide which reference DB will be used for a similarity based taxonomy annotation. Databases can be combined, with the first having the highest prioirty. E.g. "PR2,SLV" would first use PR2 to assign OTUs and all unaasigned OTUs would be searched for with SILVA, given that \"-amplicon_type LSU\" was set. Can also be a custom fasta formatted database: in this case provide the path to the fasta file as well as the path to the taxonomy for the sequences using -tax4refDB. See also online help on how to create a custom DB. (Default: GG)',
  '-tax4refDB <file>', 'In conjunction with a custom fasta file provided to argument -refDB, this file contains for each fasta entry in the reference DB a taxonomic annotation string, with the same number of taxonomic levels for each, tab separated.',
  '-amplicon_type <LSU|SSU>', '(LSU) large subunit (23S/28S) or (SSU) small subunit (16S/18S). (Default: SSU)',
  '-tax_group <bacteria|fungi>', '(bacteria) bacterial 16S rDNA annnotation, (fungi) fungal 18S/23S/ITS annotation. (Default: bacteria)',
  '-rdp_thr <0-1>', 'Confidence thresshold for RDP.(Default: 0.8)',
  '-utax_thr <0-1>', 'Confidence thresshold for UTAX. (Default: 0.8)',
  '-taxAligner <0|blast|lambda|utax|vsearch|usearch>', 'Previously doBlast. (0) deavtivated (just use RDP); (1) or (blast) use Blast; (2) or (lambda) use LAMBDA to search against a 16S reference database for taxonomic profiling of OTUs; (3) or (utax): use UTAX with custom databases; (4) or (vsearch) use VSEARCH to align OTUs to custom databases; (5) or (usearch) use USEARCH to align OTUs to custom databases. (Default: 0)',
  '-useBestBlastHitOnly <0|1>', '(1) do not use LCA (lowest common ancestor) to determine most likely taxnomic level (not recommended), instead just use the best blast hit. (0) LCA algorithm. (Default: 0)',
  '-LCA_cover <0-1>', 'Min horizontal coverage of an OTU sequence against ref DB. (Default: 0.5)',
  '-LCA_frac <0-1>', 'Min fraction of reads with identical taxonomy. (Default: 0.9)',
  '-greengenesSpecies <0|1>', '(1) Create greengenes output labels instead of OTU (to be used with greengenes specific programs such as BugBase). (Default: 0)',
  '-ITSx <0|1>', '(1) use ITSx to only retain OTUs fitting to ITS1/ITS2 hmm models; (0) deactivate. (Default: 1)',
  '-itsx_partial <0-100>', 'Parameters for ITSx to extract partial (%) ITS regions as well. (0) deactivate. (Default: 0)',
  '-lulu <0|1>', '(1) use LULU (https://github.com/tobiasgf/lulu) to merge OTUs based on their occurence. (Default: 1)',
  '-buildPhylo <0,1,2,>','(0) do not build OTU phylogeny; (1) use fasttree2; (2) use iqtree2. (Default: 1)',
);

my $clustering_heading = "Clustering Options";
my %clustering_options = (
  '-CL|-clustering <uparse|swarm|cdhit|unoise|dada2>', 'Sequence clustering algorithm: (1) UPARSE, (2) swarm, (3) cd-hit, (6) unoise3, (7) dada2. Short keyword or number can be used to indicate clustering (Default: uparse)',
  '-id <0-1>', 'Clustering threshold for OTUs. (Default: 0.97)',
  '-swarm_distance <0-1> ', 'Clustering threshold for OTUs when using swarm clustering. (Default: 1)',
  '-chim_skew <num>', 'Skew in chimeric fragment abundance (uchime option). (Default: 2)',
  '-count_chimeras', 'Minimum size of dereplicated clustered, one form of noise removal. Can be complex terms like "10:1,3:3" -> meaning at least 10x in 1 sample or 3x in 3 different samples. (Default: 2)',
  '-deactivateChimeraCheck <0|1|2|3>', '(0) do OTU chimera checks. (1) no chimera check at all. (2) Deactivate deNovo chimera check. (3) Deactivate ref based chimera check. (Default: 0)',
  '-readOverlap <num>', 'The maximum number of basepairs that two reads are overlapping. (Default: 300)',
  '-flash_param <string>', 'custom flash parameters, since this contains spaces the command needs to be in parentheses: e.g. -flash_param "-r 150 -s 20". Note that this option completely replaces the default -m and -M flash options (i.e. need to be reinserted, if wanted)]',
  '-endRem <string>', 'DNA sequence, usually reverse primer or reverse adaptor; all sequence beyond this point will be removed from OTUs. This is redundant with the "ReversePrimer" option from the mapping file, but gives more control (e.g. there is a problem with adaptors in the OTU output). (Default: "")',
  '-xtalk <0|1>', '(1) check for crosstalk. Note that this requires in most cases 64bit usearch. (Default: 0)'
);

my $other_heading = "Other Options";
my %other_options = (
  '-v', 'Print LotuS2 version',
  '-check_map <file>', 'Mapping_file: only checks mapping file and exists.',
  '-create_map <file>', 'mapping_file: creates a new mapping file at location, based on already demultiplexed input (-i) dir. E.g. ./lotus2.pl -create_map mymap.txt -i /home/dir_with_demultiplex_fastq',
);
#HELPPARSEEND

my $option_indent = 2;
my $option_width = 24;
my $description_width = 79-$option_width-$option_indent;



sub print_option_pair {
  my $n = scalar($_);
  my $option = $_[0];
  my $description = $_[1];
  my $option_indent = $_[2];
  my $option_width = $_[3];
  my $description_width = $_[4];

  my $indent = ( ' ' x $option_indent );
  my $description_indent = ( ' ' x ($option_width + $option_indent) );
  my $option_spacer = "";

  print "$indent$option";
  if (length $option > $option_width - 2)
  {
    print "\n$description_indent";
  }
  else {
    $option_spacer = ( ' ' x ($option_width - length $option) );
    print "$option_spacer";
  }
  my @description_frags = split / /, $description;
  
  my $line_width = 0;
  foreach my $frag (@description_frags)
  {
    $line_width = $line_width + length $frag;
    if ($line_width > $description_width)
    {
      $line_width = length $frag;
      print "\n$description_indent$frag "
    }
    else
    {
      $line_width = $line_width + 1;
      print "$frag ";
    }
  }
  

  print "\n";
}

######## PRINT HELP HERE ###############

print "Lotus2 usage:\n";
print "$usage_string\n\n";
print "Minimal example:\n$usage_example\n\n";
print "#### OPTIONS ####\n\n";

print "$basic_heading:\n\n";
my $key = "";
foreach $key (keys %basic_options)
{
  print_option_pair($key, $basic_options{$key}, $option_indent, $option_width, $description_width);
}

print "\n\n$further_heading:\n\n";
foreach $key (keys %further_options)
{
  print_option_pair($key, $further_options{$key}, $option_indent, $option_width, $description_width);
}

print "\n\n$workflow_heading:\n\n";
foreach $key (keys %workflow_options)
{
  print_option_pair($key, $workflow_options{$key}, $option_indent, $option_width, $description_width);
}

print "\n\n$taxonomy_heading:\n\n";
foreach $key (keys %taxonomy_options)
{
  print_option_pair($key, $taxonomy_options{$key}, $option_indent, $option_width, $description_width);
}

print "\n\n$clustering_heading:\n\n";
foreach $key (keys %clustering_options)
{
  print_option_pair($key, $clustering_options{$key}, $option_indent, $option_width, $description_width);
}


######## PRINT HELP END #############

#### HELP BLOCK LOTUS2 END #########

#   print "  -pseudoRefOTUcalling [1: create Reference based (open) OTUs, where the chosen reference database (SLV,GG etc) is being used as cluster center. Default: 0]\n";
#print "  -OTUbuild [OTU building strategy: \"ref_closed\", \"ref_open\" or \"denovo\" (default)\n";
	exit($exit_code);
}

sub usage {
    print STDERR @_ if @_;
    help(1);
}

sub parse_duration {
  my $seconds = shift;
  my $hours = int( $seconds / (60*60) );
  my $mins = ( $seconds / 60 ) % 60;
  my $secs = $seconds % 60;
  return sprintf("00:00:%02d", $seconds) if $seconds < 60;
  return sprintf("%02d:%02d:%02d", $hours,$mins,$secs);
}
sub frame {
    my ($txt) = $_[0];
	if (length($txt)==0){return "";}
	my $showT=1; $showT = $_[1] if (@_ > 1);
	my $showFrame=1; $showFrame = $_[2] if (@_ > 2);
	my $dur = " ". parse_duration(time - $start) . "";
    my @txtarrT = split( /\n/, $txt );
	my @txtarr; 
	my $width = `tput cols`;
	#print "$width\n";
	$width = 80 if (!defined $width || $width > 80|| $width < 20);
    my $numOfChar = 10;
	my $txtSpac=  $width - $numOfChar;
	foreach my $s (@txtarrT){
		#last;
		my $cnt=0;
		while (length($s) > $txtSpac && $cnt < 20){
			my $spacIdx = index $s, ' ',$txtSpac-15;
			last if ($spacIdx <= 0);
			push (@txtarr, (substr $s,0,$spacIdx));
			substr ($s,0,$spacIdx+1) = "";
			$cnt++;
		}
		push (@txtarr, $s);
	}

    my $repStr = '-' x $width;#"=========================================================================\n";
    my $ret = "";
	$ret .= $repStr."\n" if ($showFrame==1 || $showFrame==2);
	
    for ( my $i = 0 ; $i < scalar(@txtarr) ; $i++ ) {
		my $pre = "";
		if ($i==0 && $showT){
			$pre = $dur . ' ' x ($numOfChar - length($dur));
		} else {
			$pre = ' ' x $numOfChar;
		}
        $ret .= $pre . $txtarr[$i] . "\n";
    }
    $ret .= $repStr."\n" if ($showFrame==1 || $showFrame==3);
	return $ret;
}

sub firstXseqs($ $ $) {
    my ( $otus, $numS, $out ) = @_;
    open I, "<", $otus;
    open O, ">", $out;
    my $cnt = 0;
    while ( my $line = <I> ) {
        $cnt++ if $line =~ m/^>/;
        last   if $cnt == $numS;
        print O $line;
    }
    close O;
    close I;
    return ($out);
}

sub get16Sstrand($ $) {    #
    my ( $OTUfa, $refDB ) = @_;
    my $ret      = "both";
    my $OTUfa_sh = firstXseqs( $OTUfa, 6, "$lotus_tempDir/otus4blast.tmp" );
    my $cmd = "$blastBin -query $OTUfa_sh -db $refDB -out $lotus_tempDir/blast4dire.blast -outfmt 6 -max_target_seqs 200 -perc_identity 75 -num_threads $BlastCores -strand both;"
      ;                    #-strand plus both
    if ( systemL($cmd) != 0 ) {
        printL "FAILED pre-run blast command:\n$cmd\n", 5;
    }

    #die "$t\n$cmd\n";
    my $plus  = 0;
    my $minus = 0;
    open I, "<", "$lotus_tempDir/blast4dire.blast";
    while ( my $line = <I> ) {
        my @spl = split( "\t", $line );
        if ( $spl[6] < $spl[7] ) {
            if ( $spl[8] > $spl[9] ) {
                $minus++;
            }
            else { $plus++; }
        }
        else {
            if ( $spl[8] > $spl[9] ) {
                $plus++;
            }
            else { $minus++; }
        }
    }
    close I;

    #die("P= $plus M=$minus\n");
    if ( $plus + $minus > 20 ) {
        if ( $plus / ( $plus + $minus ) > 0.9 ) { $ret = "plus"; }
        if ( $plus / ( $plus + $minus ) < 0.1 ) { $ret = "minus"; }
    }

    #die ($plus ." P $minus M \n");
    #printL "Using $ret strand for blast searches.\n",0;
    return ($ret);
}

sub calcHighTax($ $ $ $ $) {
    my ( $hr, $inT, $failsR, $LCAtax, $xtraOut ) = @_;
    printL "\nCalculating higher abundance levels\n", 0;
    my $getHit2DB = 0;
    $getHit2DB = 1 if ( $xtraOut ne "" );
    my ( $hr1, $ar1, $hr2, $tax4RefHR ) =
      readTaxIn( $inT, $LCAtax, 0, $getHit2DB );
    my %Taxo   = %{$hr1};
    my @lvls   = @{$ar1};
    my %hit2db = %{$hr2};
    my %fails  = %{$failsR};
    my %OTUmat = %{$hr};
    my $SEP    = "\t";
    my @avSmps = sort( keys %OTUmat );

    #my @avOTUs = keys %{$OTUmat{$avSmps[0]}};
	my %matOTUs;
	foreach my $sm (keys %OTUmat){
		foreach my $kk (keys %{$OTUmat{$sm}}){
			$matOTUs{$kk} = 0;
		}
	}
	my @matOTUs2 = keys %matOTUs;
	#print "@matOTUs2   ".@matOTUs2."\n"."\n";

    #pseudo OTU ref hits
    if ( $xtraOut ne "" ) {
        open O, ">$xtraOut" or die "Can't open ref ${OTU_prefix} file $xtraOut\n";
        print O "RefOTUs" . $SEP . join( $SEP, @avSmps );
        my @avTax = sort( keys(%hit2db) );
		my %matOTUsTT = %matOTUs; 
        foreach my $ta (@avTax) {
            print O "\n" . $ta;
            my @OTUli = @{ $hit2db{$ta} };
			#first mark which OTUs get hit at all at this level (to count unclassis)
			foreach (@OTUli) {
				$matOTUsTT{$_} = 1;
			}
            foreach my $sm (@avSmps) {
                printL("Sample $sm does not exist\n"), 0 if ( !exists( $OTUmat{$sm} ) );
                my $smplTaxCnt = 0;
                foreach (@OTUli) {
                    if ( !$keepUnclassified && exists( $fails{$_} ) ){
						next;
					}
                    if ( !exists( $OTUmat{$sm}{$_} ) ) {
                        printL( "Sample $sm id $_ does not exist. $ta.\n", 0 );
                    }
                    else {
                        $smplTaxCnt += $OTUmat{$sm}{$_};
                    }
                }
                print O $SEP . $smplTaxCnt;
            }
        }
        close O;
    }

    #standard add up taxonomy on higher levels
    my $cnt = 0;
    my %totCnt;
    my %assCnt;
    my %assTax;
    my %unassTax;
	my @addedUnclOTUs=();
    foreach my $l (@lvls) {
        $cnt++;
        next if ( $cnt == 1 );    #domain doesn't need a file..
        my $lvlSmplCnt       = 0;
        my $assignedCnt      = 0;
        my $taxAssingedCnt   = 0;
        my $taxUnAssingedCnt = 0;
        my $isUnassigned     = 0;
        my @avTax            = sort( keys( %{ $Taxo{$l} } ) );
        my $outF             = $highLvlDir . $l . ".txt";
        open O, ">", $outF;

        #header
        print O $l . $SEP . join( $SEP, @avSmps );
		my %matOTUsTT = %matOTUs; 

        #smplWise counts, summed up to tax
        foreach my $ta (@avTax) {
            print O "\n" . $ta;
            if ( $ta =~ m/\?$/ ) {
                $isUnassigned = 1;
                $taxUnAssingedCnt++;
            }
            else { $taxAssingedCnt++; }
            my @OTUli = @{ $Taxo{$l}{$ta} };
			#first mark which OTUs get hit at all at this level (to count unclassis)
			foreach (@OTUli) {
				$matOTUsTT{$_} = 1;
			}

            #die(join("-",@OTUli)." ".$ta."\n");
			foreach my $sm (@avSmps) {
				printL("Sample $sm does not exist\n"), 0  if ( !exists( $OTUmat{$sm} ) );
				my $smplTaxCnt = 0;
				foreach (@OTUli) {
					next if ( !$keepUnclassified && exists( $fails{$_} ) );
					if ( !exists( $OTUmat{$sm}{$_} ) ) {
						printL( "Sample $sm id $_ does not exist. $l. $ta.\n", 0 );
					}
					else {
						$smplTaxCnt += $OTUmat{$sm}{$_};
					}
				}
				$assignedCnt += $smplTaxCnt if ( !$isUnassigned );
				$lvlSmplCnt  += $smplTaxCnt;
				print O $SEP . $smplTaxCnt;
            }
            $isUnassigned = 0;
			
						#now write out the unclassified ones

        }
		if ($keepUnclassified){
			print O "\n" . "noHit;";
			my @OTUli = ();
			foreach (keys %matOTUsTT){
				push @OTUli, $_ if ($matOTUsTT{$_} == 0);
			}
			push(@addedUnclOTUs , scalar(@OTUli ));
			foreach my $sm (@avSmps) {
				my $smplTaxCnt = 0;
				foreach (@OTUli) {
					if ( !exists( $OTUmat{$sm}{$_} ) ) { printL( "Sample $sm id $_ does not exist.\n", 0 ); 
					} else {$smplTaxCnt += $OTUmat{$sm}{$_};
				}
			}
			$lvlSmplCnt  += $smplTaxCnt;
			print O  $SEP.$smplTaxCnt;
			}
		}
		close O;
        $totCnt{$l}   = $lvlSmplCnt;
        $assCnt{$l}   = $assignedCnt;
        $assTax{$l}   = $taxAssingedCnt;
        $unassTax{$l} = $taxUnAssingedCnt;

        #		printL $l.": ". ($lvlSmplCnt)." ",0;
    }

	if ($keepUnclassified){
		printL "Adding " .join(",",@addedUnclOTUs) . " unclassified OTUs to ".join(",",@lvls) ." levels, respectively\n";
	}

    #print some stats
    $cnt = 0;
    foreach my $l (@lvls) {
        $cnt++;
        next if ( $cnt == 1 );    #domain doesn't need a file..
        if ( $cnt == 2 ) {
            printL
"Total reads in matrix: $totCnt{$l}\nTaxLvl	%AssignedReads	%AssignedTax\n",
              0;
        }
        my $assPerc = 0;
        $assPerc = $assCnt{$l} / $totCnt{$l} if ( $totCnt{$l} > 0 );
        my $assPerc2 = 0;
		my $lsum=( $assTax{$l} + $unassTax{$l} );
        $assPerc2 = $assTax{$l} / $lsum if ( $lsum > 0 );
        printL "$l\t" . int( 100 * $assPerc ) . "\t" . int( 100 * $assPerc2 ) . "\n", 0;
    }
    printL "\n", 0;
    return $tax4RefHR;
}

sub buildTree($ $) {
    my ( $OTUfa, $outdir ) = @_;
	return "" if ($buildPhylo==0);
	my $treePrg = "iqtree2";
	$treePrg = "fasttree" if ($buildPhylo == 1);
	if ( -f $mafftBin && (-f $fasttreeBin || -f $iqTreeBin) ) {
		printL( frame("Building tree ($treePrg) and aligning (mafft) OTUs"),0 );
	} else {
		printL(frame("Skipping tree building and multiple alignment: \nafftBin or fasttreeBin/iqTreeBin are not defined"),0);
	}

	
    my $multAli  = $outdir . "/otuMultAlign.fna";
    my $outTree  = $outdir . "/Tree.tre";
	my $treePrefix = "$lotus_tempDir/IQ";
    my $tthreads = $uthreads;

    #my $cmd = $clustaloBin . " -i $OTUfa -o $multAli --outfmt=fasta --threads=$tthreads --force;";
	my $cmd = "$mafftBin --thread $tthreads --quiet $OTUfa > $multAli";
	#die "$cmd\n";

    if ( $exec == 0 && $onlyTaxRedo == 0 && -f $mafftBin ) {
        if ( !-f $OTUfa ) {
            printL "Could not find ${OTU_prefix} sequence file:\n$OTUfa\n", 5;
        }
		systemL $cmd;
        
        #$citations .= "Clustalo multiple sequence alignments: Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, et al. 2011. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol 7: 539.\n";
		$citations .= "======== Phylogenetic tree building ========\nKatoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30(14):3059-3066. doi:10.1093/nar/gkf436";
    } elsif ($onlyTaxRedo) { printL "Skipping Multiple Alignment step\n", 0; }
    if ( $exec == 0 && $onlyTaxRedo == 0 ) {
		if ( !-f $multAli ) {
			printL "Could not find multiple alignment file:\n$multAli\n", 5;
		}
		if ($buildPhylo==1){
			$cmd = $fasttreeBin . " -nt -gtr -no2nd -spr 4 -log $logDir/fasttree.log -quiet -out $outTree $multAli;";
		} elsif ($buildPhylo==2){
			$cmd = "$iqTreeBin -s $multAli -ntmax $tthreads -pre $treePrefix -seed 678 ";
			$cmd .= "--quiet --fast -m GTR+R10 --alrt 1000\n";#GTR+F+I+G4
			$cmd .= "cp ${treePrefix}.treefile $outTree\n";#cp final good tree over 
		} else {
			printL "Unknown tree option ($buildPhylo)\n",53;
		}
        #die($cmd."\n");
        if ( $exec == 0 ) {
            #printL "Building tree..\n";
            if ( systemL($cmd) != 0 ) {
                printL "FAILED tree building command: " . $cmd . "\n", 5;
            }
            #$citations .="FastTree2 phylogenetic tree construction for OTUs: Price MN, Dehal PS, Arkin AP. 2010. FastTree 2--approximately maximum-likelihood trees for large alignments. ed. A.F.Y. Poon. PLoS One 5: e9490.\n";
			$citations .= "Nguyen, L.-T., Schmidt, H. A., von Haeseler, A. & Minh, B. Q. IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Mol. Biol. Evol. 32, 268–274 (2015)."
        }
    } elsif ($onlyTaxRedo) { 
		printL frame("Skipping Tree building step\n"), 0; 
	}
	return $outTree;
}

sub getGGtaxo($ $) {
    my ( $ggTax, $rDBname ) = @_;
    open TT, "<", $ggTax or printL "Can't open taxonomy file $ggTax\n", 87;

    #my @taxLvls = ("domain","phylum","class","order","family","genus");
    my %ret;
    while ( my $line = <TT> ) {
        chomp $line;
        my @spl = split( "\t", $line );
        my $tmp = $spl[1];
        if ( @spl < 2 ) {
            die( "Taxfile line missing tab separation:\n" . $line . "\n" );
        }
        $tmp =~ s/__;/__\?;/g;
        $tmp =~ s/__unidentified;/__\?;/g;
        $tmp =~ s/s__$/s__\?/g;
        $tmp =~ s/\s*[kpcofgs]__//g;
        $tmp =~ s/\"//g;
        my @sp2 = split( ";", $tmp );
        foreach (@sp2) { s/\]\s*$//; s/^\s*\[//; chomp; }
        my $taxv = join( "\t", @sp2 );

        #die($taxv."\n");
        $ret{ $spl[0] } = \@sp2;
    }
    close TT;
    printL "Read $rDBname taxonomy\n", 0;
    return %ret;
}

sub maxTax($) {
    my ($in)  = @_;
    my @spl22 = @{$in};
    my $cnt   = 0;
    foreach (@spl22) {
        last if ( $_ eq "?" );
        $cnt++;
    }
    return $cnt;
}

sub correctTaxString($ $) {
    my ( $sTax2, $sMaxTaxX ) = @_;
    my @ta  = @{$sTax2};
    my @ta2 = @ta;

    #die "@ta\n".$ta[$sMaxTaxX]." ".$sMaxTaxX."\n";
    for ( my $i = 0 ; $i < $sMaxTaxX ; $i++ ) {
        $ta2[$i] =~ s/^\s+//;
    }
    for ( my $i = $sMaxTaxX ; $i < @ta2 ; $i++ ) {
        $ta2[$i] = "?";
    }
    return \@ta2;
}

sub add2Tree($ $ $) {
    my ( $r1, $r2, $mNum ) = @_;
    my %refT = %{$r1};
    my @cT   = @{$r2};
    my $k    = "";

    #print $cT[3]."\n";
    #my $tmp = join("-",@cT); print $tmp." ".$mNum."\n";
    for ( my $i = 0 ; $i < $mNum ; $i++ ) {
        last if $cT[0] eq "?";
        if ( $i == 0 ) {
            $k = $cT[0];
        }
        else { $k .= ";" . $cT[$i]; }
        if ( exists( $refT{$i}{$k} ) ) {
            $refT{$i}{$k}++;
        }
        else { $refT{$i}{$k} = 1; }
    }
    return \%refT;
}

sub LCA($ $ $) {
    my ( $ar1, $ar2, $maxGGdep ) = @_;
    my @sTax       = @{$ar1};
    my @sMaxTaxNum = @{$ar2};
    if ( scalar(@sTax) == 1 ) {

        #print"early";
        my @tmpX = @{ $sTax[0] };
        my @tmp  = ();
        for ( my $i = 0 ; $i < $sMaxTaxNum[0] ; $i++ ) {
            push( @tmp, $tmpX[$i] );
        }
        for ( my $re = scalar(@tmp) ; $re < $maxGGdep ; $re++ ) {
            push( @tmp, "?" );
        }
        return ( \@tmp );
    }
    my $r1 = {};
    for ( my $i = 0 ; $i < scalar(@sMaxTaxNum) ; $i++ ) {

        #my @temp = split($sTax[$i]);
        #print @{$sTax[$i]}[0]."  $sMaxTaxNum[$i] \n";
        #next if ($sTax[$i] =~ m/uncultured /);
        $r1 = add2Tree( $r1, $sTax[$i], $sMaxTaxNum[$i] );
    }
    my %refT      = %{$r1};
    my $fini      = 0;
    my $latestHit = "";

    #determine which taxa has the highest number of hits
    my $dk;

    #print $sMaxTaxNum[0]."\n";
    my $numHits = int(@sTax) + 1;
    foreach $dk ( sort { $a <=> $b } ( keys %refT ) ) {

        #print $dk." ";
        my @curTaxs = keys %{ $refT{$dk} };
        foreach my $tk (@curTaxs) {

     #if ($dk == 2){print int($LCAfraction*$numHits). " ". $refT{$dk}{$tk}.":";}
     #if ($refT{$dk}{$tk} < $numHits){#need to get active
            if ( $refT{$dk}{$tk} >= int( $LCAfraction * $numHits ) ) {
                $latestHit = $tk;
                $fini      = 0;
                $numHits   = $refT{$dk}{$latestHit};
                last;

                #} #else {#$fini = 1;#last;}
            }
            else {
                $fini = 1;

                #$latestHit = $tk;
            }
        }

        if ($fini) { last; }
    }

    #my $winT = join("\t",@refT);
    #print "LAT ".$latestHit."\n";
    my @ret = split( ";", $latestHit );
    for ( my $re = scalar(@ret) ; $re < $maxGGdep ; $re++ ) {
        push( @ret, "?" );
    }

    #my $maxN = maxTax(\@ret);
    return ( \@ret );    #,$maxN);
}

sub extractTaxoForRefs($ $ $) {
    my ( $OTUrefDBlnk, $BlastTaxR, $GGref ) = @_;
    my %ORDL  = %{$OTUrefDBlnk};
    my %BTaxR = %{$BlastTaxR};     #my %BTaxDepR = %{$BlastTaxDepR};
    my @work  = keys %ORDL;
    if ( !@work ) { return; }

    #open O,">>$taxblastf" or die "Can't open taxOut file $taxblastf\n";
    #TODO: load this only once
    my %GG = %{$GGref};
    foreach my $w (@work) {

        #die $ORDL{$w}." XX ".$GG{$ORDL{$w}}."\n";
        #my @cur = @{$GG{$ORDL{$w}}};
        #$BTaxDepR{$w} = 7;
        $BTaxR{$w} = $GG{ $ORDL{$w} };
    }

    #close O;
    return ( \%BTaxR );    #,\%BTaxDepR);
}

sub rework_tmpLines($ $ $ $ $ $ $ $) {
    my ( $tmpLinesAR, $hr1, $sotu, $splAR, $sID, $sLength, $GGhr, $maxGGdep ) =
      @_;
    my @tmpLines   = @{$tmpLinesAR};
    my $debug_flag = 0;

    #if ($sotu eq "S1__ITS_26"){print "\nYAAYYA\n\n\n"; $debug_flag = 1;}
    #my (%ret,%retD) = (%{$hr1},%{$hr2});

    if ( @tmpLines == 0 || $sID == 0 )
    {    #prob no entry passed inclusion criteria
            #${$hr1}{$sotu} = [];${$hr2}{$sotu} = 0;
            #print $sotu."\n";
        ${$hr1}{$sotu} = [];    #$retD{$sotu} = 0;
        return ($hr1);
    }
    my %ret = %{$hr1};

    #if ($sotu eq "S1__ITS_26"){print "\nYAAYYA\n\n\n";}
    my %GG = %{$GGhr};

    #my @ggk = keys %GG; print @ggk."NN\n";
    my @spl        = @{$splAR};
    my @sTax       = ();
    my @sMaxTaxNum = ();
    my $tolerance  = 1.5;

  #extend tolerance upon lesser hits (to include more spurious, off-target hits)
    if    ( $maxHitOnly == 1 ) { $tolerance = 0; }
    elsif ( $sID == 100 )      { $tolerance = 0.1; }
    elsif ( $sID >= 99.5 )     { $tolerance = 0.25; }
    elsif ( $sID >= 99 )       { $tolerance = 0.5; }
    elsif ( $sID >= 98 )       { $tolerance = 1; }
    elsif ( $sID >= 97 )       { $tolerance = 1.25; }

    #print "XX $sID $tolerance $sLength\n";
    foreach my $lin2 (@tmpLines) {    #just compare if the tax gets any better
            #only hits within 1.5% range of best (first) hit considered
        my @spl2 = @{$lin2};    #split("\t",$lin2);
                                #print "$spl2[2] < ($sID - $tolerance\n";
        if ( $spl2[2] < ( $sID - $tolerance ) )           { next; }
        if ( $spl2[3] < ( $sLength * $lengthTolerance ) ) { next; }
        my $sMax2 = 0;
        foreach (@idThr) {
            if ( $spl2[2] < $_ ) { $sMax2++ }
        }
        $sMax2 = 7 - $sMax2;
        unless ( exists $GG{ $spl2[1] } ) {
            die "Can't find GG entry for $spl2[1]\n";
        }
        my $tTax = $GG{ $spl2[1] };

        #print $tTax." JJ\n";
        my $sMax3 = maxTax($tTax);

        #die "$tTax  $sMax3\n";
        if ( $sMax3 <= $sMax2 ) { $sMax2 = $sMax3; }
        push( @sTax,       $tTax );
        push( @sMaxTaxNum, $sMax2 );
    }

    #print "@sTax\n";
    #print $sID."\n";
    #entry for last OTU with best results etc..
    if ( $sotu ne "" ) {
        die "sTax not defined: LC="
          . @tmpLines
          . "\n@{$tmpLines[0]}\n@{$tmpLines[1]}\n@{$tmpLines[2]}\n@{$tmpLines[3]}\n"
          unless ( @sTax > 0 );
        my ($sTaxX) = LCA( \@sTax, \@sMaxTaxNum, $maxGGdep );

        #print $sotu." ".@{$sTaxX}[0]."\n";
        $ret{$sotu} = $sTaxX;

        #$retD{$sotu} = $sMaxTaxX;
        #if ($debug_flag){print "$sTaxX  $sMaxTaxX $retD{$sotu} $sotu\n";}
    }
    return ( \%ret );    #,\%retD);

}

sub splitBlastTax($ $) {
    my ( $blf, $num ) = @_;
    my $blLines = `wc -l $blf | cut -f1 -d ' ' `;
    my $endL    = int( $blLines / $num );
    if ( $endL < 3000 ) { $endL = 3000; }
    my $subLcnt = $endL;
    my @subf;
    my $fcnt   = 0;
    my $totCnt = 0;
    open I, "<$blf" or die "Can't open to split $blf\n";
    my $lstHit = "";
    my $OO;

    while ( my $l = <I> ) {
        $l =~ m/^(\S+)\s/;

        #my $hit = $1;
        if ( $1 ne $lstHit && $subLcnt >= $endL ) {

            #open new file
            #print "$lstHit  $1 $subLcnt\n";
            $subLcnt = 0;
            $lstHit  = $1;
            close $OO if ( defined $OO );
            open $OO, ">$blf.$fcnt";
            push( @subf, "$blf.$fcnt" );
            $fcnt++;
        }
        else {
            $lstHit = $1;    #has to continue until next time a change occurs..
        }
        print $OO $l;
        $subLcnt++;
        $totCnt++;

    }
    close $OO;
    close I;

    #die $blLines." $blf\n";
    #die "@subf\n$totCnt\n";
    return @subf;
}

sub getTaxForOTUfromRefBlast($ $ $) {
    my ( $blastout, $GGref, $interLMode ) = @_;

    #sp,ge,fa,or,cl,ph
    my %GG       = %{$GGref};
    my @ggk      = keys(%GG);
    my $maxGGdep = scalar( @{ $GG{ $ggk[0] } } );
    @ggk = ();
    open B, "<", $blastout or die "Could not read $blastout\n";
    my $sotu    = "";
    my $sID     = 0;
    my $sLength = 0;

    #my @sTax=(); my @sMaxTaxNum = ();
    my $retRef   = {};    # my $retDRef ={};
    my $cnt      = 0;
    my @tmpLines = ();    #stores Blast lines
    my @spl;              #temp line delim
    my %prevQueries = ();
    my $line2       = "";
    while ( my $line = <B> ) {
        $cnt++;
        chomp $line;
        $line2 = $line;
        my @spl  = split( "\t", $line );
        my $totu = $spl[0];                #line otu
        $totu =~ s/^>//;
        if ( $cnt == 1 ) { $sotu = $totu; }

        #print $line." XX $spl[11] $spl[10]\n"; die () if ($cnt > 50);

#check if this is a 2 read-hit type of match (interleaved mode) & merge subsequently
        if ($interLMode) {
            if ( @tmpLines > 0 && exists( $prevQueries{ $spl[1] } ) ) {
                my @prevHit = @{ $prevQueries{ $spl[1] } };
                die
"something went wrong with the inter matching: $prevHit[1] - $spl[1]\n"
                  unless ( $prevHit[1] eq $spl[1] );
                $prevHit[11] += $spl[11];    #bit score
                $prevHit[3]  += $spl[3];
                $prevHit[5]  += $spl[5];
                $prevHit[4] +=
                  $spl[4];    #alignment length,mistmatches,gap openings
                $prevHit[2] = ( $prevHit[2] + $spl[2] ) / 2;

                #$tmpLines[-1] = \@prevHit;
                @spl = @prevHit;

                #next;
            }
            else { $prevQueries{ $spl[1] } = \@spl; }
        }
        if ( ( $spl[11] < $minBit ) || ( $spl[10] > $minEval ) )
        {                     #just filter out..
            $spl[2] = 0;      #simply deactivate this way...
        }

        if ( $sotu eq $totu ) {
            push( @tmpLines, \@spl );
            if ( $spl[2] > $sID && $spl[3] > $sLength ) {
                $sID     = $spl[2];
                $sLength = $spl[3];
            }
            if ( $spl[3] > ( $sLength * 1.4 ) && $spl[2] > ( $sID * 0.9 ) ) {
                $sID     = $spl[2];
                $sLength = $spl[3];
            }                 #longer alignment is worth it..
                              #print $sID."\n";
        }
        else {
            #print "Maybe\n";
            ($retRef) =
              rework_tmpLines( \@tmpLines, $retRef, $sotu, \@spl, $sID,
                $sLength, \%GG, $maxGGdep )
              if ( $sotu ne "" );
            $sotu = $totu;
            undef @tmpLines;
            undef %prevQueries;
            push( @tmpLines, \@spl );
            $prevQueries{ $spl[1] } = \@spl;
            $sID                    = $spl[2];
            $sLength                = $spl[3];
        }
    }

    #last OTU in extra request
    if ( $sotu ne "" ) {
        my @spl = split( "\t", $line2 );
        ($retRef) =
          rework_tmpLines( \@tmpLines, $retRef, $sotu, \@spl, $sID, $sLength,
            \%GG, $maxGGdep );
    }
    close B;

    #debug
    #my %ret = %{$retRef};	my @tmp = @{$ret{$sotu}};print "\n@tmp  $sotu\n";

    return ($retRef);
}

sub numberOTUs($ $ $) {    #rewrites OTUs with new header (numbered header)
    my ( $inf, $outf, $pref ) = @_;
    open X,  "<", $inf;
    open XX, ">", $outf;
    open L,  ">", $inf . "_convertNames";
    my $wrLine = "";
    my $cnt    = 0;
    my %lnk;
    while ( my $line = <X> ) {
        if ( $line =~ m/^>/ ) {
            print XX ">" . $pref . $cnt . "\n";
            $line =~ m/>(\S+).*/;
            my $tmp = $1;
            $tmp =~ m/(.*);size.*/;
            $lnk{$1} = $pref . $cnt;

            #print $1." CC ".$pref.$cnt."\n";
            print L $1 . "\t" . $pref . $cnt . "\n";
            $cnt++;
        }
        else {
            print XX $line;
        }
    }
    close X;
    close XX;
    close L;

    return \%lnk;
}

sub newUCSizes($ $) {
    my ( $fafil, $ucfil ) = @_;
    my $hr = readFastaHd($fafil);
    my %fa = %{$hr};

    #read UC file
    open UC, "<", $ucfil;
    my $tok;
    my @splu;
    while ( my $line = <UC> ) {
        chomp($line);
        my $cnt = 0;
        if ( $line =~ m/^S/ ) {
            @splu = split( /\t/, $line );
            $tok  = $splu[8];
            $tok =~ s/;size=(\d+);//;

            #die ($line."\n".$1."\n");
            $cnt = $1;
        }
        elsif ( $line =~ m/^H/ ) {
            @splu = split( /\t/, $line );
            $tok  = $splu[9];
            $tok =~ s/;size=\d+;//;
            my $tik = $splu[8];
            $tik =~ m/;size=(\d+);/;
            $cnt = $1;
        }
        if ( !exists( $fa{$tok} ) ) {
            die("Expected fasta seq with head $tok\n");
        }
        $fa{$tok} += $cnt;
    }
    close UC;

    replaceFastaHd( \%fa, $fafil );
}

sub readFastaHd($) {
    my ($fil) = @_;
    open( FAS, "<", "$fil" ) || die("Couldn't open FASTA file $fil.");
    my %Hseq;
    my $line;
    while ( $line = <FAS> ) {
        if ( $line =~ m/^>/ ) {
            chomp($line);
            $line =~ s/;size=\d+;.*//;
            $line = substr( $line, 1 );

            #die $line."\n";
            $Hseq{$line} = 0;
        }
    }
    close FAS;
    return \%Hseq;
}

sub revComplFasta {
    my ($inF) = @_;
    open I, "<$inF";
    open O, ">$inF.tmp";
    while ( my $l = <I> ) {
        chomp $l;
        if ( $l !~ m/^>/ ) { $l = reverse_complement_IUPAC($l); }
        print O $l . "\n";
    }
    close I;
    close O;
    systemL "rm $inF; mv $inF.tmp $inF;";
}

sub reverse_complement_IUPAC {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~
      tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}

sub readFasta($) {
    my ($fil) = @_;
    my %Hseq;
    if ( -z $fil ) { return \%Hseq; }
    open( FAS, "<", "$fil" ) || printL( "Couldn't open FASTA file $fil\n", 88 );
    my $temp="";
    my $line;
    my $trHe = <FAS>;
    chomp($trHe);
    $trHe = substr( $trHe, 1 );
    $trHe =~ s/;size=\d+;.*//;
    $trHe =~ s/\s.*//;

    while ( $line = <FAS> ) {
		chomp $line;
		if ( $line =~ m/^>/ ) {
			#finish old fas`
			$Hseq{$trHe} = $temp;
			#prep new entry
			$trHe = substr( $line, 1 );
			$trHe =~ s/;size=\d+;.*//;
			$trHe =~ s/\s.*//;
			$temp = "";
		} else {
			$temp .= $line;
		}
    }
    $Hseq{$trHe} = $temp;
    close(FAS);
    return \%Hseq;
}

sub writeFasta($ $) {
    my ( $hr, $outFile ) = @_;
    my %fas = %{$hr};
    open O, ">$outFile" or printL "Can't open outfasta $outFile", 99;
    foreach my $k ( keys %fas ) {
        print O ">" . $k . "\n" . $fas{$k} . "\n";
    }
    close O;
}

sub swarm4us_size($) {
    my ($otuf) = @_;
    open I, "<$otuf"      or die "Can't open $otuf\n";
    open O, ">$otuf.ttmp" or die "Can't open $otuf.tmp\n";
    while ( my $line = <I> ) {
        chomp $line;
        if ( $line =~ m/^>/ ) {
            $line =~ s/_(\d+)$/;size=$1;/;
        }
        print O $line . "\n";
    }
    close I;
    close O;

    #	die $otuf."\n";
    systemL "rm $otuf;mv $otuf.ttmp $otuf;";
}

sub replaceFastaHd($ $) {
	my ( $href, $ifil ) = @_;
	my $ofil = $ifil . ".tmp";
	my %fa   = %{$href};
	open O, ">", $ofil;
	open( FAS, "<", "$ifil" ) || die("Couldn't open FASTA file.");
	my %Hseq;
	my $line;
	while ( $line = <FAS> ) {
		if ( $line =~ m/^>/ ) {
			chomp($line);
			$line =~ s/;size=\d+;//;
			$line = substr( $line, 1 );
			if ( !exists( $fa{$line} ) ) { die("can't find key $line\n"); }
			$line = ">" . $line . ";size=" . $fa{$line} . ";\n";
		}
		print O $line;
	}
	close FAS;
	close O;
	rename( $ofil, $ifil ) or die "Unable to overwrite file $ifil\n";
}

sub onelinerSWM($) {
    my ($ifil) = @_;
    my $ofil = $ifil . ".tmp";
    open O, ">", $ofil;
    open( FAS, "<", "$ifil" ) || die("Couldn't open FASTA file.");
    my $seq = "";
    while ( my $line = <FAS> ) {
        chomp($line);
        if ( $line =~ m/^>/ ) {
            $line =~ s/;size=(\d+);.*/_$1/;
            if ( $seq ne "" ) {
                print O $seq . "\n" . $line . "\n";
                $seq = "";
            }
            else {
                print O $line . "\n";
            }
        }
        else {
            $seq .= $line;
        }
    }
    print O $seq;
    close FAS;
    close O;
    rename( $ofil, $ifil ) or die "Unable to overwrite file $ifil\n";
}

sub readRDPtax($) {
    my ($taxf) = @_;    #$avOTUR
    open T, "<", $taxf or die "Can't open $taxf : \n" . $!;
    my @taxLvls = ( "domain", "phylum", "class", "order", "family", "genus" );
    my @Hiera   = ();
    my %Fail;

    #my %AvOTUs = %{$avOTUR};
    my $rankN = -1;
    while ( my $line = <T> ) {
        chomp($line);

        #print $line."\n";
        $line =~ s/"//g;
        my @spl = split( "\t", $line );
        next unless ( $spl[0] =~ m/^${OTU_prefix}/ && @spl > 3 );

        #print $spl[0]."\n";
        if ( $rankN == -1 ) {
            for ( my $i = 3 ; $i < 10 ; $i += 3 ) {
                if ( $spl[$i] =~ m/domain/ ) {
                    $rankN = $i - 1;
                    last;
                }
            }
        }
        next if ( $rankN == -1 );

        #phylum level class
        if ( $spl[ $rankN + 5 ] < $RDPCONF ) {
            $Fail{ $spl[0] } = 1;

            #next;
        }

        #$AvOTUs{$spl[0]} = 1;
        #push(@AvOTUs,$spl[0]);
        my $nhier = "";    #$spl[$rankN]."\t";
        my $tcnt  = 0;     #1,+3
        for ( my $i = $rankN + 0 ; $i <= $rankN + 16 ; $i += 3 ) {
            if ( $taxLvls[$tcnt] eq $spl[ $i + 1 ] ) {
                $tcnt++;
            }
            else { die( "expected " . $taxLvls[$tcnt] . " " . $line . "\n" ); }
            if ( $spl[ $i + 2 ] < $RDPCONF ) {
                $nhier .= "?\t";
            }
            else             { $nhier .= $spl[$i] . "\t"; }
            if ( $i > @spl ) { last; }
        }
        $nhier .= "?\t" . $spl[0];
        push( @Hiera, $nhier );
    }
    close T;

    #write hier
    open H, ">", $RDP_hierFile;
    print H join( "\t", @taxLvls ) . "\tspecies\tOTU\n";
    foreach (@Hiera) { print H $_ . "\n"; }
    close H;

    #print(keys(%Fail)."\n");
    return ( \@Hiera, \%Fail );    #\%AvOTUs,
}

sub blastFmt($ $) {
    my ( $aref, $tD ) = @_;
    $aref = correctTaxString( $aref, $tD );
    my @in = @{$aref};
    return join( "\t", @{$aref} );
}

sub writeUTAXhiera($ $ $) {
    my ( $utout, $avOTUsR, $failsR ) = @_;
    my %fails  = %{$failsR};
    my @avOTUs = @{$avOTUsR};

    #parse utax output
    my %utStr;
    open I, "<$utout" or die "Can;t open $utout\n";
    while (<I>) {
        my @spl   = split /\t/;
        my $k     = $spl[0];
        my @spl2  = split( /,/, $spl[1] );
        my $nline = "";
        my $cnt   = 0;
        foreach my $sus (@spl2) {
            $sus =~ s/^\S://;
            $sus =~ s/\"//g;
            $sus =~ m/(^.*)\(([0-9\.]*)\)/;
            my $sus2 = $1;
            if ( $2 < $utaxConf ) {
                $sus2 = "?";
                if ( $cnt == 0 ) { $fails{$k} = 1; last; }
            }
            $nline .= $sus2 . "\t";
            $cnt++;
        }
        while ( $cnt < 7 ) { $nline .= "?\t"; $cnt++; }

        #die $nline."\n$_";
        $utStr{$k} = $nline;
    }
    close I;

    open H, ">", $SIM_hierFile or die "Can't open $SIM_hierFile\n";
    print H "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOTU\n";

    #print all failed OTUs
    my @faKeys = keys(%fails);
    foreach (@faKeys) {
        if ( exists( $utStr{$_} ) ) {
            delete $fails{$_};
        }
        else {
            print H "?\t?\t?\t?\t?\t?\t?\t" . $_ . "\n";
        }
    }

    #print newly derrived tax
    foreach (@avOTUs) {    #(sort(keys(%{$taxr}))){
        my $failed1 = exists( $fails{$_} );
        if ( !$failed1 ) {    #actually print OTU taxonomy
            die "Can;t find $_ key in utax\n" unless ( exists( $utStr{$_} ) );
            print H $utStr{$_} . "$_\n";
        }
    }

    #die "UTAX: ".@faKeys." unassigned OTUs\n";
    close H;
    return ( \%fails );
}

sub writeBlastHiera($ $ $) {
    my ( $taxr, $avOTUsR, $failsR ) = @_;

    #my %fails = %{$failsR};
    my %fails  = %{$failsR};
    my @avOTUs = @{$avOTUsR};
    my $cnt    = 0;
    open H, ">", $SIM_hierFile or die "Can't open $SIM_hierFile\n";
    print H "${OTU_prefix}\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
    #
    foreach ( keys(%fails) ) {
        if (   exists( ${$taxr}{$_} )
            && ${$taxr}{$_}
            && maxTax( ${$taxr}{$_} ) > 0 )
        {
            delete $fails{$_};
        }
        else {
            print H $_ . "\t?\t?\t?\t?\t?\t?\t?" . "\n";
        }
    }
    foreach (@avOTUs) {    #(sort(keys(%{$taxr}))){
         #if ( exists( ${$taxr}{$_}) && !exists( ${$taxDepr}{$_}) ) {print "entry missing: $_\n";}
         #print maxTax(${$taxr}{$_}). "HH\n";
        my $tdep = 0;
        next if ( exists( $fails{$_} ) );
        $tdep = maxTax( ${$taxr}{$_} )
          if ( exists( ${$taxr}{$_} ) && ${$taxr}{$_} );
        if ( $tdep > 0 ) {
            print H $_ . "\t" . blastFmt( ${$taxr}{$_}, $tdep ) . "\n";
        }
        else {
            print H $_ . "\t?\t?\t?\t?\t?\t?\t?\t" . "\n";
            $fails{$_} = 1;
        }
        $cnt++;
    }
    close H;
    return ( \%fails );
}

sub checkLtsVer($) {
    my ($lver)  = @_;
    my $sdmVstr = `$sdmBin -v`;
    my $sdmVer  = 1;
    if ( $sdmVstr =~ m/sdm (\d\.\d+) \D+/ ) {
        $sdmVer = $1;
    }
    else {
        $sdmVstr = `$sdmBin`;

        #die $sdmVstr;
        if ( $sdmVstr =~ m/This is sdm version (\d\.\d+) \D+/ ) {
            $sdmVer = $1;
        }
    }

    # compare to server version
    unless ( $checkForUpdates == 1 ) { return $sdmVer; }
    #die "LWP:simple package not installed, but required for automatic update checker!\n deactivate \"CheckForUpdates 0\" in $lotusCfg to circumvent the update checks.\n" if ( !$LWPsimple );

    printL "Checking for updates..  ";
    my $url = "http://psbweb05.psb.ugent.be/lotus/lotus/updates/Msg.txt";
    if ( !head($url) ) {
        printL "LotuS server seems to be down!\n";
        return $sdmVer;
    }
    my $updtmpf = get($url);

    my $msg     = "";
    my $hadMsg  = 0;
    my $msgCont = "";
    open( TF, '<', \$updtmpf );
    while (<TF>) { $msg .= $_; }
    close(TF);
    foreach my $lin ( split( /\n/, $msg ) ) {
        my @spl = split /\t/, $lin;
        next if ( @spl == 0 );
        if ( $lver < $spl[0] ) { $msgCont .= $spl[1] . "\n\n" }
        $hadMsg = 1;
    }
    $updtmpf =
      get("http://psbweb05.psb.ugent.be/lotus/lotus/updates/curVer.txt");
    open( TF, '<', \$updtmpf );
    my $lsv = <TF>;
    close(TF);
    chomp $lsv;
    my $msgEnd = "";
    $updtmpf =
      get("http://psbweb05.psb.ugent.be/lotus/lotus/updates/curVerMsg.txt");
    open( TF, '<', \$updtmpf );
    while (<TF>) { $msgEnd .= $_; }
    close(TF);

    $updtmpf =
      get("http://psbweb05.psb.ugent.be/lotus/lotus/updates/UpdateHist.txt");
    my $updates = "";
    open( TF, '<', \$updtmpf );
    $msg = "";
    while (<TF>) { $msg .= $_; }
    close(TF);
    foreach my $lin ( split( /\n/, $msg ) ) {
        my @spl = split /\t/, $lin;
        chomp $lin;
        next if ( @spl < 2 || $spl[0] eq "" );
        $spl[1] =~ m/LotuS (\d?\.\d+)/;
        if ( $lver < $1 ) { $updates .= $spl[0] . "\t" . $spl[1] . "\n" }
    }
    if ( $updates ne "" || $msgCont ne "" ) {
        printL "\n";
        if ( $updates ne "" ) {
            printL
"--------------------------------\nThe following updates are available:\n--------------------------------\n";
            printL $updates;
            printL
              "\n\nCurrent Lotus version is :$lver\nLatest version is: $lsv\n";
        }
        if ( $msgCont ne "" ) {
            printL
"--------------------------------\nThe following messages were on LotuS server:\n--------------------------------\n";
            printL $msgCont;
        }
    }
    else {
        printL "Your LotuS version is up-to-date!\n";
        return $sdmVer;
    }

    if ( $hadMsg || $updates ne "" ) {
        printL
"New LotuS updates discovered (10s wait).\nIf you want to install updates\n  (1) Press \"Ctrl-c\"\n  (2) perl autoInstall.pl\n\n";
        printL
"To deactivate, open $sdmOpt and change \"CheckForUpdates\" to \"0\"\n";
        sleep(10);
    }
    printL "Continuing LotuS run..\n";

    #die;
    return $sdmVer;
}

sub addSeqRun($){
	my ($hr) = @_;
	my %newMap =  %{$hr};
	my %retMap; #second map by SequencingRun, so sorting is easier later
	#SequencingRun
	my @header = @{$newMap{"#SampleID"}};
	my $c=0;
	my $samplTagByFX=0; my $samplTagByDIR = 0; #how to estimate the original miSeq run??
	my %heads = map { $_ => $c++ } @header;
	my %allSRs;
	if (exists($heads{SequencingRun} )){
		my $SRcol = $heads{SequencingRun};
	#	return \%newMap; -> still needs to run through to sort etc
		foreach my $k (sort keys %newMap){
			if ($k eq "#SampleID"){
				$retMap{"#SampleID"}{"#SampleID"} = $newMap{$k};
				next;
			} 
			my $SRid = ${$newMap{$k}}[$SRcol];
			$retMap{$SRid}{$k} = $newMap{$k};
			$allSRs{$SRid} = 1;
		}
		my @SRcats = keys(%allSRs);
		printL "Found \"SequencingRun\" column, with ". scalar(@SRcats) . " categories (@SRcats)\n";
		return (\%newMap,\%retMap);
	}  
	
	
	
	#this route tries to create SequencingRun based on fq info etc
	my $fqFid = -1;
	if (exists($heads{fastqFile})){
		$fqFid = $heads{fastqFile};
	}
	if (exists($heads{fnaFile})){
		die "both \"fnaFile\" and \"fastqFile\" are defined in map. Not allowed .. aborting\n" if ($fqFid != -1);
		$fqFid = $heads{fnaFile};
	}
	if ($fqFid != -1 && exists ($heads{BarcodeSequence}) ){
		$samplTagByFX = 1;
	} else {
		$samplTagByDIR = 1;
	}
	my %allFqs; #stores fastq/dirs to fastqs
	my @SRids = ("a" .. "z", "A" .. "Z", 1 .. 1000);
	my $SRcnt=0; my $maxSRcnt = 0;
	my $smplCnt=0; my $totalSample = scalar(keys %newMap);
	my $splitUp = 100; #how many samples to include in each dada2 batch? (in absence of other information)
	if ($totalSample < 190 && $totalSample > 50){
		$splitUp = $totalSample/2;
	}
	foreach my $k (sort keys %newMap){
		if ($k eq "#SampleID"){
			push(@{$newMap{$k}},"SequencingRun");
			$retMap{"#SampleID"}{"#SampleID"} = $newMap{$k};
			next;
		} 
		if ($samplTagByFX || $samplTagByDIR) {
			my $curFX = @{$newMap{$k}}[$fqFid]; 
			$curFX =~ s/,.*//;
			if ($samplTagByDIR) { #in this case: only use dir but not sample name
				if ($curFX =~ m/^(.*)\/[^\/]+/){
					$curFX = $1;
				} else {$curFX = "";} #important to still delete this, avoid pure files being confused with dirs
			}
			if (!exists($allFqs{$curFX})){
				$allFqs{$curFX} = scalar(keys %allFqs);
			}
			#print "$curFX $SRcnt  $SRids[$SRcnt] \n";
			$SRcnt = $allFqs{$curFX};
		}  else {#no info whatsoever, just ignore
			$smplCnt++;
			if ($smplCnt > $splitUp){ #this inserts a random split to at least ensure dada2 clustering isn't growing too big
				$SRcnt++;$smplCnt=0;
			}
		}
		push(@{$newMap{$k}},$SRids[$SRcnt]);
		$retMap{$SRids[$SRcnt]} {$k} = $newMap{$k};
		$maxSRcnt = $SRcnt if ($SRcnt > $maxSRcnt);
	}
	#print "$maxSRcnt  X\n";
	printL "Can't find \"SequencingRun\" column in map, attempting auto creation with " . ($maxSRcnt+1) ." detected sequencing runs. Please check in $cpMapFile if this is correct\n","w" if ($maxSRcnt > 0);
	return (\%newMap,\%retMap);
}

sub writeMap($$){
	my ($hrM,$outFile) = @_;
	my %mapX = %{$hrM};
	open O,">$outFile" or die "Can't open $outFile\n";
	foreach my $k0 (sort keys %mapX){
		foreach my $k (sort keys (%{$mapX{$k0}}) ){
			print O $k . "\t" . join("\t",@{$mapX{$k0}{$k}}) . "\n";
		}
	}
	close O;
}

sub readMap() {
	my ($createSmplSet) = @_;#tries to automatically derrive sample set
    #find number of samples
    #my @avSMPs = ();
	my $Ltxt = "Reading mapping file\n";
    my %mapH;
    my %combH;
    unless ( open M, "<", $mapFile ) {
        printL( "Couldn't open map file $mapFile: $!\n", 3 );
    }
    local $/ = undef;
    my $mapst = <M>;
    $mapst =~ s/\r\n|\n|\r/\n/g;
    close M;

    #$mapst=~s/\R//g;
    my @mapar = split( /\n/, $mapst );

    #try mac file conversion
    #if (scalar(@mapar) == 1){ @mapar = split(/\r\n/,$mapst);}
    #if (scalar(@mapar) == 1){ @mapar = split(/\r/,$mapst);}
    my $cnt          = 0;
    my $warnTrig     = 0;
    my $fileCol      = -1;
    my $twoFiles     = 0;
    my $CombineCol   = -1;
    my $hasCombiSmpl = 0;
    my $colCnt       = -1;

    foreach my $line (@mapar) {
        $cnt++;
        next if ( $cnt > 1 && ( $line =~ m/^#/ || length($line) < 2 ) );
        chomp($line);
        my @spl = split( "\t", $line );
        if ( @spl == 0 ) {
            $warnTrig = 1;
            printL"Mapping file contains lines that cannot be split by tab seperator (line $cnt):\"$line\"\n","w";
            next;
        }
        if ( $colCnt != -1 && @spl != $colCnt ) {
            for ( my $i = @spl ; $i < $colCnt ; $i++ ) { $spl[$i] = ""; }
        }
        chomp( $spl[0] );
        if ( $spl[0] eq "" ) {
            $warnTrig = 1;
            printL("Empty SampleID in row number $cnt, skipping entry.\n","w");
            next;
        }
        my $smplNms = $spl[0];
        if ( $fileCol != -1 ) {
            if ( $spl[$fileCol] =~ m/,/ ) {
                $Ltxt .= "Switching to paired end read mode\n" if ( $numInput != 2 );
                $numInput = 2;
            }
            elsif ( $numInput == 2 ) {
                printL "Inconsistent file number in mapping. See row with ID $smplNms.\n", 55;
            }
        }
        if ( $CombineCol != -1 && $spl[$CombineCol] ne "" ) {
            $combH{ $spl[$CombineCol] } = $smplNms;
        }
        else {
            $combH{$smplNms} = $smplNms;
        }
        if ( $cnt == 1 ) {
            if ( $line !~ m/^#/ ) {
                printL "First line does not start with \"#\". Please check mapping file for compatibility (http://psbweb05.psb.ugent.be/lotus/documentation.html#MapFile)\n",93;
            }
            my $ccn = 0;
            $colCnt = @spl;
            foreach (@spl) {
                if ( $_ eq "fastqFile" || $_ eq "fnaFile" ) {
                    $Ltxt .= "Sequence files are indicated in mapping file.\n";
                    if ( $fileCol != -1 ) {
                        $Ltxt .= "both fastqFile and fnaFile are given in mapping file, is this intended?\n";
                    }
                    $fileCol = $ccn;
                }
                if ( $_ eq "CombineSamples" ) {
                    $Ltxt .= "Samples will be combined.\n";
                    $CombineCol   = $ccn;
                    $hasCombiSmpl = 1;
                }

                $ccn++;
            }
        }
        if ( $line =~ m/\"/ ) {
            $warnTrig = 0;
            printL("Possible biom incompatibility: Mapping file contains paranthesis characters (\") for sample $smplNms. Lotus is removing this.","w");
        }
        if ( $line =~ m/ / ) {
            $warnTrig = 1;
            printL("Possible biom incompatibility: Mapping file contains spaces (\" \") for sample $smplNms","w");
        }
        if ( $line =~ m/[^\x00-\x7F]/ ) {
            $warnTrig = 1;
            printL("Possible biom incompatibility: Mapping file contains non-ASCII character for sample $smplNms","w");
        }
        $line =~ s/\s+/\t/g;
        if ( $smplNms =~ m/^\s/ || $smplNms =~ m/\s$/ ) {
            printL"SampleID $smplNms contains spaces. Aborting LotuS as this will lead to errors, please fix.\n",5;
        }
        if ( $cnt == 1 ) {    #col names
            if ( $line =~ m/\t\t/ ) {
                printL"Empty column header in mapping file:\n check for double tab chars in line 1:\n$mapFile\n",4;
            }
            if ( $smplNms ne '#SampleID' ) {
                printL"Missing \'#SampleID\' in first line of file $mapFile\n Aborting..\n",65;
            }
			my $hasFwd=0; my $hasRev=0;

            #if ($line =~/\tCombineSamples\t/){$combineSamples=1;} #deprecated
            my $hcnt = 0;
            foreach (@spl) {
                $hcnt++;
                if (m/^\s*$/) {
                    printL"Empty column header in mapping file (column $hcnt)\n$mapFile\n",4;
                }
				if (m/ForwardPrimer/ || m/LinkerPrimerSequence/){$hasFwd=1;}
				if (m/ReversePrimer/){$hasRev=1;}
            }
            if ($hasFwd==0){
				$warnTrig = 1;
				printL("No forward PCR primer for amplicon found in mapping file (column header \"ForwardPrimer\". This might invalidate chimera checks\n","w");
			}
			
			if ( $line =~ m/\t\t/ ) {
                printL"Empty header in mapping file:\n check for double tab chars in line 1:\n$mapFile\n",4;
            }
			
        }
		
        splice( @spl, 0, 1 );

        #print $smplNms." ".$spl[0]."\n";
        for ( my $i = 0 ; $i < @spl ; $i++ ) {
            $spl[$i] =~ s/\"//g;
        }

        #print $smplNms."\n";
        $mapH{$smplNms} = [@spl];

    }

    #print keys %mapH."\n";
    if ( scalar( keys %mapH ) == 0 ) {
        printL("Could not find sample names in mapping file (*nix/win file ending problem?\n",9);
    }

    if ( $warnTrig == 1 ) {
        print "*********\nWarnings for mapping file \n$mapFile \nAbort by pressing Ctrl+c (10 sec wait)\n*********\n";
        sleep(10);
    }
    printL( frame($Ltxt), 0 );
    return ( \%mapH, \%combH, $hasCombiSmpl );
	
}

sub finWarn($) {
    my ($msg) = @_;
    $finalWarnings .= $msg ;
	$finalWarnings .= "\n" unless ($finalWarnings =~ m/\n$/);
    print STDERR "$msg \n";
}

sub readOTUmat($) {
    my ($file) = @_;

    #my @avSmps = sort(keys %OTUmat);
    open I, "<", $file or die "Can't open expected ${OTU_prefix} counts: $file\n";
    my $cnt = -1;
    my @samples;
    my @avOTUs;
    my %OTUm;
    while ( my $line = <I> ) {
        $cnt++;
        chomp($line);
        my @spl = split( "\t", $line );
        if ( $cnt == 0 ) {    #header
            @samples = @spl[ 1 .. $#spl ];

           #die($samples[0]." ".$samples[$#samples]." ".@spl." ".@samples."\n");
            next;
        }
        my $curOT = shift(@spl);
        push( @avOTUs, $curOT );
        for ( my $i = 0 ; $i < @samples ; $i++ ) {
            $OTUm{ $samples[$i] }{$curOT} = $spl[$i];
        }
    }
    close I;
    return ( \%OTUm, \@avOTUs );
}

sub utaxTaxAssign($ $) {
    my ( $query, $taxblastf ) = @_;
    if ( !-f $usBin ) { printL "uearch binary not found: $usBin\n", 81; }
    my $dbfa = $TAX_REFDB[0];
    if ( !-d $dbfa ) { printL "wrong utax input dir: $dbfa\n", 82; }

    #contains path to utax dir
    my $taxLconf = "120";    #250,500,full_length
    if ( $platform eq "454" )    { $taxLconf = 500; }    #454, miSeq, hiSeq
    if ( $platform eq "miseq" )  { $taxLconf = 250; }
    if ( $platform eq "pacbio" ) { $taxLconf = 700; }
    my $utFas  = "$dbfa/fasta/refdb.fa";
    my $utUdb  = "$dbfa/fasta/refdb.$taxLconf.udb";
    my $utConf = "$dbfa/taxconfs/$taxLconf.tc";

    $cmd = "";
    unless ( -e $utUdb ) {
        $cmd .=
          "$usBin -makeudb_utax $utFas -output $utUdb -taxconfsin $utConf;";
    }
    $cmd .= "$usBin -utax $query -db $utUdb -utaxout $taxblastf -strand both;";

    #die $utaCmd."\n";

    if ( $exec == 0 ) {
        printL frame("Assigning taxonomy using UTAX with confidence $taxLconf bp");

        #print $cmd."\n";
        if ( systemL($cmd) ) { printL "UTAX failed:\n$cmd\n", 3; }
    }
    return $taxblastf;
}

sub doDBblasting($ $ $) {
    my ( $query, $DB, $taxblastf ) = @_;
    my $simMethod = "";

#die "$doBlasting\n";
#if ($doBlasting != 3 && $doBlasting && (-f $blastBin || -f $lambdaBin) && $DB ne "" && -f $DB){
    $REFTAX = 1;
    if ( $doBlasting == 1 ) {
        printL "Could not find blast executable.\n$blastBin\n", 33
          unless ( -f $blastBin );
        printL "Could not find blast executable.\n$mkBldbBin\n", 33
          unless ( -f $mkBldbBin );
        $cmd = "$mkBldbBin -in $DB -dbtype 'nucl';";
        unless ( -f $DB . ".nhr" ) {
            if ( systemL($cmd) ) {
                printL "makeBlastDB command failed:\n$cmd\n", 3;
            }
        }
        my $strand = "both";

#if ($exec==0){$strand = get16Sstrand($query,$DB);} #deactivated as speed gain is too moderate for gain
        $cmd = "$blastBin -query $query -db $DB -out $taxblastf -outfmt 6 -max_target_seqs 200 -perc_identity 75 -num_threads $BlastCores -strand $strand ;"
          ;    #-strand plus both minus
        if ( !-s $query ) { $cmd = "touch $taxblastf;"; }
        else {
            $citations .=
"Blast taxonomic similarity search: Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215: 403–10.\n";
        }
        $simMethod = "BLAST";
    
	} elsif ($doBlasting == 4 || $doBlasting == 5){#new default: vsearch 
		my $outCols="query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+ql";
        my $udbDB = $DB . ".vudb";
        $udbDB = $DB . ".udb" if ($doBlasting == 5 );
        if ( !-f $udbDB ) { 
            print "Building UDB database (takes some time the first time)\n";
			if ($doBlasting == 4 ){
				if (systemL("$VSBinOri  --makeudb_usearch $DB -output $udbDB;")){printL "VSEARCH DB building failed\n";}
			} else {
				if (systemL("$usBin  --makeudb_usearch $DB -output $udbDB;")){printL "USEARCH DB building failed\n";}
			}
        }
       if ($doBlasting == 4){ 
		$cmd = "$VSBinOri ";
			$citations .= "VSEARCH taxonomic database search: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584\n" unless ( $citations =~ m/VSEARCH taxonomic database search/ );
	   }
		$simMethod = "VSEARCH";
		if ($doBlasting == 5){
			$cmd = "$usBin "; 
			$simMethod = "USEARCH";
			$citations .= "USEARCH taxonomic database search: \n" unless ( $citations =~ m/USEARCH taxonomic database search/ );
		}
		$cmd .= "--usearch_global $query --db $udbDB  --id 0.75 --query_cov 0.5 -userfields $outCols -userout $taxblastf --maxaccepts 100 --maxrejects 100 -strand both --threads $BlastCores;";
		#--blast6out $taxblastf
		#die "$cmd\n";
    } elsif ( $doBlasting == 2 ) {    #lambda
        printL "Could not find lambda executable.\n$lambdaBin\n", 33   unless ( -f $lambdaBin );
        my $lamVer  = 0.4;
        my $lverTxt = `$lambdaBin --version`;
        if ( $lverTxt =~ m/lambda version: (\d\.\d)[\.0-9]* \(/ ) {$lamVer = $1;}
        ( my $TAX_REFDB_wo = $DB ) =~ s/\.[^.]+$//;
        print $DB. ".fm.txt.concat\n";
        if (   $lamVer > 0.4
            && -f $DB . ".dna5.fm.sa.val"
            && !-f "$DB.dna5.fm.lf.drv.wtc.24" )
        {
            printL "Rewriting taxonomy DB to Lambda > 0.9 version\n", 0;
            systemL "rm -r $DB.*;";
        }
        if ( !-f $DB . ".dna5.fm.sa.val" ) {
            print "Building LAMBDA index anew (may take up to an hour first time)..\n";

            #die($DB.".dna5.fm.sa.val");
            my $xtraLmbdI = "";
            $xtraLmbdI = " --algorithm skew7ext " if ($lowMemLambI);
            my $cmdIdx =  "$lambdaIdxBin -p blastn -t $BlastCores -d $DB $xtraLmbdI;";
            #print $cmdIdx."\n";
            systemL "touch $DB.dna5.fm.lf.drv.wtc.24;";
            if ( systemL($cmdIdx) ) {
                printL( "Lamdba ref DB build failed\n$cmdIdx\n", 3 );
            }
        }

        #OMP_NUM_THREADS = $BlastCores
        print "Starting LAMBDA similarity search..\n";

        #TMPDIR env var.. TODO
        my $tmptaxblastf = "$lotus_tempDir/tax.m8";
		my $outcols = "'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen'";
        $cmd = "$lambdaBin -t $BlastCores -id 75 -nm 200 -p blastn -e 1e-8 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $query -d $DB -o $tmptaxblastf -oc $outcols;";
        $cmd .= "mv $tmptaxblastf $taxblastf;";

        #lambda is not guranteed to return sorted list <- apparently it does
        #$cmd .= "sort $tmptaxblastf > $taxblastf;rm $tmptaxblastf";

        if ( !-s $query ) { $cmd = "rm -f $taxblastf;touch $taxblastf\n"; }
        else { $citations .= "Lambda taxonomic similarity search: Hauswedell H, Singer J, Reinert K. 2014. Lambda: the local aligner for massive biological data. Bioinformatics 30: i349-i355\n" unless ( $citations =~ m/Lambda taxonomic similarity search:/ );
        }
        $simMethod = "LAMBDA";

        #die $cmd."\n";
    } else {
        printL "Unknown similarity comparison program option (-simBasedTaxo): \"$doBlasting\"\n", 98;
    }
    if ($extendedLogs) { $cmd .= "cp $taxblastf $extendedLogD/;"; }
    #die($cmd."\n");
    $duration = time - $start;
	if ( $exec == 0 ) {
		printL frame("Assigning taxonomy against $DB using $simMethod");
		#print $cmd."\n";
		if ( systemL($cmd) ) {
			printL "$simMethod against ref database failed:\n$cmd\n", 3;
		}
	}
	return $taxblastf;
}

sub findUnassigned($ $ $ ) {
    my ( $BTr, $Fr, $outF ) = @_;
    my %BT  = %{$BTr};
    my %Fas = %{$Fr};

    #my @t =keys %Fas; print "$t[0]\n";
    my $cnt    = 0;
    my $dcn    = 0;
    my @kk     = keys %Fas;
    my $totFas = @kk;

    if ( $outF eq "" ) {
        foreach my $k ( keys %BT ) {
            my @curT = @{ $BT{$k} };
            if ( @curT == 0 || $curT[0] eq "?" ) { $dcn++; } #|| $curT[1] eq "?"
            $cnt++;
        }
        if ( @kk == 0 ) {
            printL "Total of "
              . ( $cnt - $dcn )
              . " / $cnt reads have LCA assignments\n", 0;
        }
        else {
            printL "$dcn / $cnt reads failed LCA assignments, checked $totFas reads.\n",
              0;
        }
        return;
    }
    foreach my $k ( keys %BT ) {
        my @curT = @{ $BT{$k} };

        #print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n"
        if ( @curT == 0 || $curT[0] eq "?" ) {    #|| $curT[1] eq "?"){
            delete $BT{$k};
            $dcn++;

            #print ">".$k."\n".$Fas{$k}."\n";
        }    #else {print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" ;}
        else {
            die "Can't find fasta entry for $k\n" unless ( exists $Fas{$k} );
            delete $Fas{$k};
        }
        $cnt++;

        #die if ($cnt ==100);
    }
    @kk = keys %Fas;
    printL "$dcn / $cnt reads failed LCA assignments\nWriting " . @kk ." of previous $totFas reads for next iteration.\n", 0;
    open O, ">$outF" or die "can;t open unassigned fasta file $outF\n";
    foreach my $k (@kk) {
        print O ">" . $k . "\n" . $Fas{$k} . "\n";
    }
    close O;
    return ( $outF, \%BT, $dcn, \%Fas );
}

sub runRDP{
	my $cmd="";

	my $rdpGene = "16srrna";
	$rdpGene = "fungallsu" if ( $organism eq "fungi" || $organism eq "eukaryote" );
	my $msg = "";
	if ( $doRDPing > 0 && $mjar ne "" && -f $mjar ) {
		$msg = "Assigning taxonomy with multiRDP";

		#ampliconType
		$cmd = "java -Xmx1g -jar $mjar --gene=$rdpGene --format=fixrank --hier_outfile=$outdir/hierachy_cnt.tax --conf=0.1 --assign_outfile=$lotus_tempDir/RDPotus.tax $OTUfa;";
		$RDPTAX = 2;
	} elsif ( $doRDPing > 0 ) { #has to do RDP normal
		unless ( $rdpjar ne "" || exists( $ENV{'RDP_JAR_PATH'} ) ) {
			printL"Could not run RDP classifier (required).\nCheck that 'RDP_JAR_PATH' is set as environmental variable and that 'RDPjar' is set in lOTUs.cfg\n",55;
		}
		$msg = "Assigning taxonomy with RDP";
		my $toRDP = $ENV{'RDP_JAR_PATH'};
		if ( $rdpjar ne "" ) { $toRDP = $rdpjar; }
		my $subcmd = "classify";
		$cmd =        "java -Xmx1g -jar " . $toRDP  . " $subcmd -f fixrank -g $rdpGene -h $outdir/hierachy_cnt.tax -q $OTUfa -o $lotus_tempDir/RDPotus.tax -c 0.1;";
		$RDPTAX = 1;
	}
	
	if ( $doRDPing > 0 && $exec == 0 ) {    #ITS can't be classified by RDP
		printL frame($msg), 0;

		#die "XXX  $cmd\n";
		if ( systemL($cmd) ) {
			printL "FAILED RDP classifier execution:\n$cmd\n", 2;
		}
		if ($TaxOnly ne "0" && -e $TaxOnly && -e "$lotus_tempDir/RDPotus.tax"){
			systemL "cp $lotus_tempDir/RDPotus.tax $TaxOnly.rdp";
		}
		$citations .= "RDP ${OTU_prefix} taxonomy: Wang Q, Garrity GM, Tiedje JM, Cole JR. 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Env Microbiol 73: 5261–5267.\n";
	} elsif ( $rdpjar ne "" || exists( $ENV{'RDP_JAR_PATH'} ) ) {
		if ($extendedLogs && -e "$lotus_tempDir/RDPotus.tax" && -d $extendedLogD) { systemL "cp $lotus_tempDir/RDPotus.tax $extendedLogD/;"; }
		if ( $RDPTAX > 0 ) {                #move confusing files
			if ($extendedLogs) {
				systemL "mv $outdir/hierachy_cnt.tax $outdir/cnadjusted_hierachy_cnt.tax $extendedLogD/;";
			} else {
				unlink "$outdir/hierachy_cnt.tax";
				unlink "$outdir/cnadjusted_hierachy_cnt.tax"
				  if ( -e "$outdir/cnadjusted_hierachy_cnt.tax" );
			}
		}
	}
	#RDP finished 
}


sub assignTaxOnly($ $) {
    my ( $taxf, $output ) = @_;    #,$avSmpsARef
                                   #printL $outmat."\n";

    my $hieraR;
    my $failsR = {};

    #RDP taxonomy
    unless ( $doRDPing < 1 || $otuRefDB eq "ref_closed" ) {
        ( $hieraR, $failsR ) = readRDPtax($taxf); #required to rewrite RDP tax..
    }

    my $blastFlag = $doBlasting != 3 && $doBlasting;       #set that blast has to be done

    my $fullBlastTaxR = {};
    my $GGloaded      = 0;
    my %GG;
    if ( $doBlasting == 3 ) {
        my $utaxRaw = "$lotus_tempDir/tax.blast";
        my $utout   = utaxTaxAssign( $taxf, $utaxRaw );
        $failsR = writeUTAXhiera( $utout, $avOTUsR, $failsR );
    }
    elsif ( -f $LCABin && $blastFlag && $maxHitOnly != 1 ) {
        my @blOuts = ();
        for ( my $DBi = 0 ; $DBi < @TAX_REFDB ; $DBi++ ) {
            my $taxblastf = "$lotus_tempDir/tax.$DBi.blast";
            my $DB        = $TAX_REFDB[$DBi];
            my $DBtax     = $TAX_RANKS[$DBi];
            my $blastout  = doDBblasting( $taxf, $DB, $taxblastf );
            push( @blOuts, $blastout );
        }
        my $LCxtr = "";
        if ($pseudoRefOTU) { $LCxtr = "-showHitRead -reportBestHit"; }
		if (-f $TaxOnly) {$SIM_hierFile = $TaxOnly.".hier";}
        my $cmd =  "$LCABin  -i " . join( ",", @blOuts ) . " -r ". join( ",", @TAX_RANKS ) . " -o $SIM_hierFile $LCxtr -LCAfrac $LCAfraction -cover $LCAcover -id ". join( ",", @idThr ) . ";";

        #die $cmd;
        if ( systemL $cmd) { printL "LCA command $cmd failed\n", 44; }
        systemL "cp @blOuts $output\n";
    }
	return $SIM_hierFile;
}

sub makeAbundTable2($ $) {
    my ( $taxf, $retMatR ) = @_;    #,$avSmpsARef
                                    #printL $outmat."\n";

    my $hieraR;
    my $failsR = {};

    #RDP taxonomy
    unless ( $doRDPing < 1 || $otuRefDB eq "ref_closed" ) {
        ( $hieraR, $failsR ) = readRDPtax($taxf); #required to rewrite RDP tax..
    }

    my $blastFlag = $doBlasting != 3 && $doBlasting ;       #set that blast has to be done
    my %fails    = %{$failsR};
    my $curQuery = $OTUfa;

    my $fullBlastTaxR = {};
    my $GGloaded      = 0;
    my %GG;
    if ($REFflag) { #no activated currently
        #$OTUrefSEED
        if ( !$GGloaded ) {
            %GG       = getGGtaxo( $TAX_RANKS[0], $refDBname[0] );
            $GGloaded = 1;
        }

#an extra otu file exists with ref sequences only.. tax needs to be created separately for these & then merged for final output..
#this overwrites all ref OTUs, but de novo OTUs still need to be assigned..
        ($fullBlastTaxR) = extractTaxoForRefs( $OTUrefDBlnk, $fullBlastTaxR, \%GG );
        die "TODO ref assign of denove OTUs\n";
    }
    
	if ( $doBlasting == 3 ) { #utax
        my $utaxRaw = "$lotus_tempDir/tax.blast";
        my $utout   = utaxTaxAssign( $OTUfa, $utaxRaw );
        $failsR = writeUTAXhiera( $utout, $avOTUsR, $failsR );
    
	} elsif ( $blastFlag && $maxHitOnly != 1 ) {
		#TODO: maxHitOnly
        my @blOuts = ();
        for ( my $DBi = 0 ; $DBi < @TAX_REFDB ; $DBi++ ) {
            my $taxblastf = "$lotus_tempDir/tax.$DBi.blast";
            my $DB        = $TAX_REFDB[$DBi];
            my $DBtax     = $TAX_RANKS[$DBi];
            my $blastout  = doDBblasting( $curQuery, $DB, $taxblastf );
            push( @blOuts, $blastout );
        }
        my $LCxtr = "";
        if ($pseudoRefOTU) { $LCxtr = "-showHitRead -reportBestHit"; }
        my $cmd = "$LCABin  -i " . join( ",", @blOuts ) . " -r " . join( ",", @TAX_RANKS ) . " -o $SIM_hierFile $LCxtr -LCAfrac $LCAfraction  -cover $LCAcover -id " . join( ",", @idThr ) . ";";

        #die $cmd;
        if ( systemL $cmd) { printL "LCA command $cmd failed\n", 44; }
        #unlink @blOuts;
    }
    else {    #old Perl implementation
        my $taxblastf = "$lotus_tempDir/tax.blast";
        my $BlastTaxR = {};                  # my $BlastTaxDepR = {};
                                             #read this only once
        my $leftOver  = 111;
        my $fasrA     = readFasta($OTUfa);
        if ($blastFlag) {                    #do blast taxnomy
            my $DBi = 0;
            printL"legacy Perl coded LCA is being used. note that for multi DB tax assignments, this version is not recommended. Please ensure that \'LCA\' program is in your lotus dir or try a clean reinstall if you use multiDB tax assignments.\n","w";
            for ( my $DBi = 0 ; $DBi < @TAX_REFDB ; $DBi++ ) {
                my $DB       = $TAX_REFDB[$DBi];
                my $DBtax    = $TAX_RANKS[$DBi];
                my $blastout = doDBblasting( $curQuery, $DB, $taxblastf );

#my @subf = splitBlastTax($taxblastf,$BlastCores);#paralellize getTaxForOTUfromRefBlast#printL "Running tax assignment in ".@subf." threads..\n";	#my @thrs;
#for (my $i=0;$i<@subf;$i++){$thrs[$i] = threads->create(\&getTaxForOTUfromRefBlast,$subf[$i],\%GG,0);}
#for (my $i=0;$i<@subf;$i++){my $ret = $thrs[$i]->join();$BlastTaxR = {%$BlastTaxR,%$ret};}
#---------  single core way --------------
                %GG        = getGGtaxo( $DBtax, $refDBname[$DBi] );
                $GGloaded  = 1;
                $BlastTaxR = getTaxForOTUfromRefBlast( $blastout, \%GG, 0 );
                if ( $DBi < ( @TAX_REFDB - 1 ) ) {
                    ( $curQuery, $BlastTaxR, $leftOver, $fasrA ) = findUnassigned( $BlastTaxR, $fasrA,$OTUfa . "__U" . $DBi . ".fna" );
                    $taxblastf = "$lotus_tempDir/tax.blast_rem_$DBi";
                }
                else { findUnassigned( $BlastTaxR, $fasrA, "" ); }
                $fullBlastTaxR = { %$fullBlastTaxR, %$BlastTaxR };
                last if ( $leftOver == 0 );
            }
            printL "Assigned @refDBname Taxonomy to ${OTU_prefix}'s\n", 0;
            systemL "rm $OTUfa" . "__U*;" if ( $DBi > 0 );
        }

        #here do the same check on ref seqs
        if ( $REFflag || $blastFlag ) {
            findUnassigned( $fullBlastTaxR, {}, "" );
            ($failsR) = writeBlastHiera( $fullBlastTaxR, $avOTUsR, $failsR );
        }

        #	my @t = @{${$fullBlastTaxR}{"OTU_1"}}; die "@t\n";
        undef $BlastTaxR;
        undef $fullBlastTaxR;

        %fails = %{$failsR};    #my @avOTUs = @{$avOTUsR};
    }
    if ($REFflag) {
        findUnassigned( $fullBlastTaxR, {}, "" );
        ($failsR) = writeBlastHiera( $fullBlastTaxR, $avOTUsR, $failsR );
        die "TODO ref blast\n";
    }

    #get number of samples from aref
    #my @avSMPs = @{$avSmpsARef};
    #sanity check
    #foreach(@avSMPs){die("$_ key not existant") unless (exists($retMat{$_}));}

    return ( \%fails );
}

sub cutUCstring($) {
    my ($aref) = @_;
    my @in = uniq( @{$aref} );

    #my @sa = sort (@in);
    #die("@in\n");
    my @newa = ();
    foreach (@in) {
        my @spl = split($sep_smplID);
        push( @newa, $spl[0] );
    }
    return \@newa;
}

sub mergeUCs($ $) {
    my ( $cref, $dref ) = @_;
    my %cons   = %{$cref};
    my %derep  = %{$dref};
    my $cnt    = 0;
    my $totlen = keys %cons;
    foreach my $k ( keys %cons ) {

        #print $cnt." / ".$totlen."  ".@{$cons{$k}}."\n";
        #if (@{$cons{$k}} == 0){die "$k\n";}
        my @newa = @{ $derep{$k} };  #first add initial seed, should be the same
        die "$k not in derep\n" unless exists( $derep{$k} );
        foreach my $k2 ( @{ $cons{$k} } ) {
            push( @newa, @{ $derep{$k2} } );

            #print "added $k2 to $k\n";
        }
        $cons{$k} = \@newa;
        $cnt++;
    }
    printL "finished merging\n", 0;
    return %cons;
}

sub delineateUCs($ $) {
    my ( $ifi, $mode ) = @_;
    open UC, "<", $ifi;
    my %clus;
    my %expSize = ();
    while ( my $line = <UC> ) {
        chomp($line);
        my @spl     = split( "\t", $line );
        my $hit     = $spl[9];
        my $que     = $spl[8];
        my $hitsize = 1;
        my $quesize = 1;
        if ( $mode >= 1 ) {

            # if $mode==1;
            if ( $mode == 1 ) { $que =~ s/;size=(\d+).*$//; $quesize = $1; }
            $hit =~ s/;size=(\d+).*$//;
            $hitsize = $1;
        }
        if ( $mode == 3 )
        { #switch hit and que, because in UPARSE case, reads are aligned to OTUs (= not backtracing)
                #print($quesize." ".$hitsize."\n");
            my $temp = $hit;
            $hit     = $que;
            $que     = $temp;
            $temp    = $hitsize;
            $hitsize = $quesize;
            $quesize = $temp;
        }

        if ( $spl[0] eq "H" ) {    #in cluster
            if ( exists( $clus{$hit} ) ) {
                push( @{ $clus{$hit} }, $que );
                $expSize{$hit} += $quesize;
            }
            else {
                my @empty = ( $hit, $que );
                $clus{$hit}    = \@empty;
                $expSize{$hit} = $quesize + $hitsize;
            }
        }
        if ( $spl[0] eq "S" ) {    #starts new cluster

            #	print $hit." ".$que."\n"; die();
            if ( !exists( $clus{$que} ) ) {
                my @empty = ($que);
                $clus{$que}    = \@empty;
                $expSize{$que} = $quesize;
            }
        }
        if ( $mode != 2 ) {
            if ( $spl[0] eq "C" )
            { #all "H" should have appeared at this point; cross-check if cluster size is right

                if ( !exists( $clus{$que} ) ) {
                    die("Entry for $que does not exist (but should exist).\n");
                }
                else {
                    my $arraySize = scalar( @{ $clus{$que} } );
                    if ( $spl[2] != $expSize{$que} ) {
                        printL(
                            "$que :: expected "
                              . $spl[2]
                              . " cluster size. found cluster size: "
                              . $expSize{$que}
                              . " from array with S $arraySize.\n",
                            0
                        );
                        foreach ( @{ $clus{$que} } ) {
                            printL $_. "\n", 0;
                        }
                        die();
                    }
                }
            }
        }
    }
    close UC;
    return \%clus;
}

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

sub print_nseq($) {
    my ($tf) = @_;
    my $nseqs = `grep \">\" $tf | wc -l`;
    printL( $nseqs . " reads\n", 0 );
    return ($nseqs);
}

sub printL($ $) {
    my ( $msg, $stat ) = @_;
    if ( defined $stat && lc($stat) eq "w" ) {
		$msg = "WARNING:: $msg";
		#die "$msg";
		print LOG $msg;
        finWarn($msg);
        return;
    }
    print $msg if ($verbosity > 0) ;
    if ( $mainLogFile ne "" ) {
        print LOG $msg;
    }
    if ( defined $stat && $stat ne "0" ) { 
		#print "Error: $msg\n";# if ($verbosity < 1) ; #print in any case, since exiting on this message
		print "\n%@#%@#%@#%@%@#@%#@%#@#%@#%@#%@#@%#@%#@%#@#%@#%@#%@##\nYour LotuS2 run encounterend an error!\nFirst check if the last error occured in a program called by LotuS2 \n\"tail $progOutPut\"\n, if there is an obvious solution (e.g. external program breaking, this we can't fix). To see (and excecute) the last commands by the pipeline, run \n\"tail $cmdLogFile\".\nIn case you decide to contact us on \"https://github.com/hildebra/lotus2/\", please try to include information from these approaches in your message, this will increase our response time. Thank you.\n%@#%@#%@#%@%@#@%#@%#@#%@#%@#%@#@%#@%#@%#@#%@#%@#%@##\n";
		exit($stat); 
	}
}

sub systemL {
    my $cmdF = $_[0];
	my $throwWarn = 1; my $descr="";
	if (@_ > 2){$throwWarn = $_[2];}
	if (@_ > 1){$descr = $_[1];}

    #split up by subcmds..
	my $retStat=0;
	my @cmds = split /;/,$cmdF;
	foreach my $subcmd (@cmds){
		print("[cmd] $subcmd\n") if ($verbosity > 2) ;
		print cmdLOG "[cmd] $subcmd\n";
		my $attach = "";
		$attach = ">> $progOutPut 2>&1" if ($verbosity < 2);
		if ($subcmd =~ m/>/) {$attach = "";}
		my $stat = system("$subcmd $attach");
		if ($stat && $throwWarn){#something went wrong..
			print "$descr CMD failed: $subcmd \nsee $progOutPut for error log\n" ;
			$retStat = $stat;
			exit(9);
		}
	}
    return $retStat;
    # return $ENV{'PIPESTATUS[0]'};
}


sub splitFastas($ $) {
    my ( $inF, $num ) = @_;
    my $protN = cntFastaEntrs($inF);
    my $pPerFile = int( $protN / $num ) + 1;
    my $fCnt     = 0;
    my $curCnt   = 0;
    my @nFiles   = ("$inF.$fCnt");

    open I, "<$inF" or printL "Fatal: Can't open fasta file $inF\n", 59;
    open my $out, ">$inF.$fCnt"
      or printL "Fatal: Can't open output fasta file $inF.$fCnt\n", 59;
    while ( my $l = <I> ) {
        if ( $l =~ m/^>/ ) {
            $curCnt++;
            if ( $curCnt > $pPerFile ) {
                $fCnt++;
                close $out;
                open $out, ">$inF.$fCnt";
                push( @nFiles, "$inF.$fCnt" );
                $curCnt = 0;
            }
        }
        print $out $l;
    }
    close I;
    close $out;
    return \@nFiles;
}

sub extractFastas($ $ $) {
    my ( $in, $hr, $addSize ) = @_;
    my %tars = %{$hr};
    my %fasSubs;
    my $hd = "";
    open I, "<$in" or printL "Fatal: Can't open fasta file $in\n", 59;
    my $use = 0;
    while ( my $l = <I> ) {
        chomp $l;
        if ( $l =~ m/^>/ ) {
            $use = 0;
            $l =~ m/^>(\S+)/;
            $hd = $1;
            if ( exists( $tars{$hd} ) ) {
                if ($addSize) { $hd .= ";size=" . $tars{$hd} . ";"; }
                $use = 1;
                $fasSubs{$hd} = "";
            }
            next;
        }
        $fasSubs{$hd} .= $l if ($use);
    }
    close I;
    return \%fasSubs;
}

sub combine($ $) {
    my ( $fadd, $freal ) = @_;
    if ( $osname eq 'MSWin32' ) {
        systemL("type $fadd  $freal > $freal");
    }
    else {
        systemL("cat $fadd >> $freal");
    }
    unlink("$fadd");
}



sub prefix_find($){
	my ($ar) = @_;
	my @rds = @{$ar}; my @newRds = @rds;
	my $first = $rds[0]; 
	for (my $i=0;$i<@rds; $i++){
		my $second = $rds[$i];
		"$first\0$second" =~ m/^(.*\/).*\0\1/s;
		$first = $1;
	}
	
	print "$first common prefix\n";
	for (my $i=0; $i<@newRds; $i++){
		$newRds[$i] =~ s/^$first//;
		$newRds[$i] .= "/" if ($newRds[$i] !~ m/\/$/ && length($newRds[$i] ) != 0 );
	}
	return (\@newRds,$first);
}


sub autoMap{
	my ($refDirs, $ofile) = @_;
	my $prefix = "SMPL";
	my $smpCnt = 0;
	my $pairedP = "2";#just set default
	my @paired = split /,/,$pairedP;

	if (@paired > 1 ){
		print "Looking for mix of paired and single input files\n";
	}elsif ($paired[0]==2){
		print "Looking for paired input files\n";
	} elsif ($paired[0]==1){
		print "Looking for single read input files\n";
	} else {
		die "Arg 2 has to be supplied and be either \"1\" or \"2\" (single / paired read input in dirs)\n";
	}
	#### regex search string
	my $rawFileSrchStr1 = '(R1_00\d|\.1|_1|\.R1)\.f[^\.]*[aq][\.gz]*$';
	my $rawFileSrchStr2 = '(R2_00\d|\.2|_2|\.R2)\.f[^\.]*[aq][\.gz]*$';
	my $rawSingSrchStr = '.f[^\.]*[aq][\.gz]*$';

	my $map = "#SampleID	fastqFile\n";
	my @preRDs = split(/,/,$refDirs);
	if (@preRDs > 1 && @paired == 1){ #extend to full length
		for (my $i=1;$i<@preRDs;$i++){
			$paired[$i] = $paired[0];
		}
	}
	my ($arRDs,$pathPre) = prefix_find(\@preRDs);
	my @RDs = @{$arRDs};
	my $cnt=0;
	my $paired=2;
	foreach my $RD ( @RDs ){
		#print "DIR: $RD\n";
		opendir(DIR, "$pathPre/$RD") or die $!;
		my @pa2 =(); 
		my @pa1 = ();
		if ($paired == 2){
			@pa2 = sort ( grep { /$rawFileSrchStr2/ && -e "$pathPre/$RD/$_" } readdir(DIR) );	rewinddir DIR;
			@pa1 = sort ( grep { /$rawFileSrchStr1/  && -e "$pathPre/$RD/$_"} readdir(DIR) );	close(DIR);
			#fix reads double assigned, always assume 2nd is correct
			if (@pa1 != @pa2){
				my %h;
				@h{(@pa2)} = undef;
				@pa1 = grep {not exists $h{$_}} @pa1;
			}
			if (@pa1 != @pa2){
				$paired=1;
			}
		}
		
		#in case things are not matching up.. go for single reads
		if ($paired == 1){
			@pa1 = sort ( grep { /$rawSingSrchStr/  && -e "$pathPre/$RD/$_"} readdir(DIR) );	close(DIR);
		}
				
		die "unequal read pairs!: ".@pa1 .",". @pa2."\n" if (@pa1 != @pa2 && $paired == 2);
		print "Found ".@pa1." samples in dir $RD\n";
		print "Assuming unpaired reads, as can not find signature for pair 2 files\n" if (@pa2 == 0 && @pa1 >0 && $paired == 2);
		#add to map
		foreach (my $i=0; $i< @pa1; $i++){
			if (@pa2 > 0){
				$map .= $prefix.$smpCnt."\t".$RD."".$pa1[$i].",".$RD."".$pa2[$i]."\n";
			} else {
				$map .= $prefix.$smpCnt."\t".$RD."".$pa1[$i]."\n";
			}
			$smpCnt ++;
		}
		$cnt++;
		
	}
	open O,">$ofile" or die "Can't open $ofile\n";
	print O $map;
	close O;

	print "Map is $ofile\n";
	print "Please check that all files required are present in map $ofile.\n";
	print "==========================\nTo start analysis:\n./lotus2.pl -m $ofile -i $pathPre/ -o [outdir] [further parameters if desired]\n==========================\n";
	exit(0);

}

sub usearchDerepSort($ ) {
    my ($filtered) = @_;
    my $outputLab = "-output";
    if ( $usearchVer >= 8 && $VSused == 0 ) { $outputLab = "-fastaout"; }

    printL(
"\n =========================================================================\n Dereplicate exact sub-sequences & sort (usearch).\n ",
        0
    );
    $cmd = "$VSBin -derep_fulllength $filtered $outputLab $lotus_tempDir/derep.fa -uc $lotus_tempDir/derep.uc -log $logDir/derep.log -sizeout -threads $uthreads;";

    #die($cmd);
    if ( systemL($cmd) != 0 ) { printL "Failed dereplication\n", 1; }

    my $dcnt          = 1;
    my $bigDerepCnt   = 0;
    my @lotsFiles     = ($filtered);
    my @subDerepFa    = ("$lotus_tempDir/derep.fa");    #my @subDerepUC = ("$lotus_tempDir/derep.uc");
    my @subDerepFaTmp = ();
    while ( -f $filtered . "." . $dcnt ) {
        if ( $dcnt == 1 ) { @subDerepFaTmp = @subDerepFa; @subDerepFa = (); }
        push( @lotsFiles,     $filtered . "." . $dcnt );
        push( @subDerepFaTmp, "$lotus_tempDir/derep.fa.$dcnt" )
          ;    #push(@subDerepUC,"$lotus_tempDir/derep.uc.$dcnt");
        $cmd = "$VSBin -derep_fulllength $filtered.$dcnt $outputLab $subDerepFaTmp[-1] -sizeout -threads $uthreads;"
          ;    #-uc $lotus_tempDir/derep.uc -log $logDir/derep.log
               #die $cmd;
        if ( systemL($cmd) != 0 ) { printL "Failed dereplication\n", 1; }
        if ( -f $filterOutAdd . "." . $dcnt ) {

            #Combining low qual reads..
            combine( $filterOutAdd . "." . $dcnt, $filterOutAdd );
        }
        if ( ( $dcnt + 1 ) % 10 == 0 || !-f $filtered . "." . ( $dcnt + 1 ) )
        {      #intermediary merge.. too big otherwise
            my $firstF = shift(@subDerepFaTmp);

            #print join(" ",@subDerepFaTmp)."\n";
            my $pid =
              open( CMD,
                "| cat " . join( " ", @subDerepFaTmp ) . " >> $firstF" )
              || printL "Merge of derep subfiles failed.\n", 1;
            defined($pid) || printL( "Cannot fork derep of subfiles.\n", 1 );
            close CMD;
            foreach (@subDerepFaTmp) { unlink; }
            @subDerepFaTmp = ();

            my $dereOut = "$lotus_tempDir/derep.inter.$bigDerepCnt.temp";
            $cmd = "$VSBin -derep_fulllength $firstF $outputLab $dereOut -uc $lotus_tempDir/derep.uc -log $logDir/derep.log -sizeout -threads $uthreads;";
            if ( systemL($cmd) != 0 ) {
                printL "Failed intermediate dereplication\n", 1;
            }
            newUCSizes( $dereOut, "$lotus_tempDir/derep.uc" );

            #delete tmp files
            unlink("$lotus_tempDir/derep.uc");
            unlink("$firstF");
            push( @subDerepFa, $dereOut );
            $bigDerepCnt++;

            #consecutive bigger files
            if ( @subDerepFa > 1 ) {
                $firstF = shift(@subDerepFa);
                $pid    = open( CMD,
                    "| cat " . join( " ", @subDerepFa ) . " >> $firstF" );
                defined($pid)
                  || printL( "Cannot fork derep of subfiles.\n", 1 );
                close CMD || printL "Merge of derep subfiles failed.\n", 1;
                $cmd = "$VSBin -derep_fulllength $firstF $outputLab $lotus_tempDir/derep.post.temp -uc $lotus_tempDir/derep.uc -log $logDir/derep.log -sizeout -threads $uthreads;";
                print $cmd. "\n";
                systemL $cmd || printL "Derep Command Failed: \n$cmd\n", 1;
                newUCSizes( "$lotus_tempDir/derep.post.temp", "$lotus_tempDir/derep.uc" );
                foreach (@subDerepFa) { print $_. "\n"; unlink; }
                systemL("mv $lotus_tempDir/derep.post.temp $lotus_tempDir/derep.fa");
                unlink("$lotus_tempDir/derep.pre.temp");
                @subDerepFa = ("$lotus_tempDir/derep.fa");
            }
            elsif ( !-f $filtered . "." . ( $dcnt + 1 ) ) {
                systemL("mv $subDerepFa[0] $lotus_tempDir/derep.fa");
            }
        }
        $dcnt++;
    }
    my $dereplicate_minsizeX = $dereplicate_minsize;
    $dereplicate_minsizeX =~ m/^(\d+)\D?/;
    $dereplicate_minsizeX = $1;
    printL "Reset derep min to $dereplicate_minsizeX, as u/vsearch does not support more complex options\n","w";

    #dereplicate_minsize is 2 by default
    $cmd = "$VSBin -sortbysize  $lotus_tempDir/derep.fa $outputLab $lotus_tempDir/derep.fas -minsize $dereplicate_minsizeX -log $logDir/sortbysize.log;";
    if ( systemL($cmd) != 0 ) { exit(1); }
    unlink("$lotus_tempDir/derep.fa");
    return ("$lotus_tempDir/derep.fas");
}

sub parseSDMlog($) {
    my ($inF) = @_;
    open I, "<", $inF;
    my $totSeq = 0;
    my $SeqLen = 0;
    while (<I>) {
        chomp;
        if (m/^Accepted: (\d+) \(/) { $totSeq = $1; }
        if (m/^\s+\- Seq Length : \d+\/(\d+)\/\d+/) { $SeqLen = $1; }
    }
    close I;

    #die $totSeq." " .$SeqLen."\n";
    return ( $totSeq, $SeqLen );

    #"Accepted: 349515 (0 end-trimmed)","    - Seq Length : 250/250/250"

}

sub removeNonValidCDhits($ $ $) {
    my ( $refDB4otus, $addon, $denNN ) = @_
      ; #idea is to use the addon.clstr to find refSeqs that were used the denNN file (are not in there)

}

sub dnaClust2UC($ $) {
    my ( $in, $out ) = @_;
    print "DNAclust1";
    my %cluster;
    my %clSize;
    my %clus_denovo;
    my %clDNSize;
    open I, "<$in" or printL "Fatal: Can't open dnaclust output $in\n", 38;
    while ( my $line = <I> ) {
        my @spl = split( /\s/, $line );
        next unless ( @spl > 1 );    #/data/falhil/otutmpAN1//clusters.uc
        my $curCl = shift(@spl);

        #print $curCl."\n";
        my $siz = 0;
        foreach (@spl) { m/;size=(\d+);/; $siz += $1; }
        die "DNAclust:: double Cluster detected: $curCl\n"
          if ( exists( $clus_denovo{$curCl} )
            || exists( $clus_denovo{$curCl} ) );
        if ( $curCl =~ m/;size=\d+;/ ) {    #denovo cluster
            $clus_denovo{$curCl} = \@spl;
            $clDNSize{$curCl}    = $siz;
        }
        else {                              #ref cluster
            $cluster{$curCl} = \@spl;
            $clSize{$curCl}  = $siz;
        }
    }
    close I;
    my @refs  = keys %cluster;
    my $CLnUM = 0;

    #print ref clusters in UC format
    #print "uc\n";
    open O, ">$out.ref"
      or printL "Fatal: Can't open dnaclust uc rewrite $in\n", 37;
    foreach my $k (@refs) {
        my @mem = @{ $cluster{$k} };
        my $ref = $k;                  #.";size=".@mem.";";
            #print O "H\t$CLnUM\t170\t100.0\t+\t0\t0\t170M\t$ref\t$ref\n";
        foreach my $mm (@mem) {
            print O "H\t$CLnUM\t170\t100.0\t+\t0\t0\t170M\t$mm\t$ref\n";
        }
        $CLnUM++;
    }
    close O;

    #print denovo clusters in UC format
    open O, ">$out" or printL "Fatal: Can't open dnaclust uc rewrite $in\n", 37;
    foreach my $k ( keys %clus_denovo ) {
        my @mem = @{ $clus_denovo{$k} };
        my $ref = $k;                      #.";size=".@mem.";";
        print O "N\t$CLnUM\t170\t100.0\t+\t0\t0\t170M\t$ref\t$ref\n";
        foreach my $mm (@mem) {
            print O "H\t$CLnUM\t170\t100.0\t+\t0\t0\t170M\t$mm\t$ref\n";
        }
        $CLnUM++;
    }
    close O;

    #die " done\n$out\n";
    return ( \%cluster, \%clSize, \%clus_denovo, \%clDNSize );
}

sub parseLastVUsearch($){
	my ($tarF) = @_;
	my $d2rep = `tail $progOutPut | grep 'Matching unique query sequence' `;
	$d2rep =~ s/Matching unique query sequence://; chomp $d2rep;
	return $d2rep;
}

##################################################
# Core cluster routing on quality filtered files #
##################################################
sub buildOTUs($) {

    my ($outfile) = @_;
    my @UCguide = ( "$lotus_tempDir/finalOTU.uc", "$lotus_tempDir/finalOTU2.uc" );    #,"$lotus_tempDir/otus.uc",1);
    #if ($exec==1){return(\@UCguide);}
    systemL "rm -f $UCguide[0]*;";
    my $refDB4otus = "$TAX_REFDB[0]" if ( @TAX_REFDB > 0 );  #reference database
    #print_nseq("$filterOut");
    my $filtered = "$lotus_tempDir/filtered.fa";

	if ($UPARSEfilter) {
		printL("\n =========================================================================\n Secondary uparse filter\n",0  );
		my $cmd = "$usBin -fastq_filter $filterOut -fastaout $filtered -fastq_trunclen $truncLfwd -threads $uthreads;";

		#-fastq_truncqual $truncQual
		#die($cmd."\n");
		if ( systemL($cmd) != 0 ) { exit(1); }
	}  else {
		$filtered = $filterOut;
	}

	my $derepl = "$lotus_tempDir/derep.fas";    #,$totSeqs,$arL)
	if ($mergePreCluster){
		$derepl = "$lotus_tempDir/derep.merg.fas";
		if ($takeNonMerge && $ClusterPipe != 7 && $ClusterPipe != 1 ){
			system "cat $lotus_tempDir/derep.fas >> $derepl";
		}
	}
	my ( $totSeqs, $SeqLength ) = parseSDMlog("$logDir/demulti.log");
	if ( !$sdmDerepDo ) {
		my ($derepl) = usearchDerepSort($filtered);
	}
	else {
		if ($ClusterPipe != 7 && (!-f $derepl || -z $derepl )) {
			printL "The sdm dereplicated output file was either empty or not existing, aborting lotus.\n$derepl\n",1;
		}

	#my @lotsFiles = ($filtered); my $dcnt=1; while (-f $filtered.".".$dcnt){	push(@lotsFiles,$filtered.".".$dcnt);$dcnt++}
	}
    my $OTUfastaTmp = $outfile;                 #"$lotus_tempDir/uparse.fa";
    my $dnaclustOut = "$lotus_tempDir/clusters_pre.dncl";

    #have to rev seq for dnaclust & cd-hit ref clustering
    if ($REFflag) {
        my $strand = get16Sstrand( $derepl, $refDB4otus );
        print $strand. " strand\n";
        if ( $strand eq "minus" ) {
            printL "reversing 16S for alignment to DB..\n", 0;
            revComplFasta($derepl);
        }
    }
	my $entrMessage=""; my $exitMsg="Finished";
	if ( $ClusterPipe == 1 ) { #UPARSE clustering
		my $maxhot  = 62;
		my $maxdrop = 12;
		#too many files need a more thorough clustering process
		if ( $totSeqs > 12000000 && $totSeqs < 24000000 ) {
			$maxhot  = 72;
			$maxdrop = 15;
		}
		elsif ( $totSeqs >= 24000000 ) {
			$maxhot  = 92;
			$maxdrop = 18;
		}
		my $id_OTUup    = $id_OTU;
		my $outputLab   = "-output";
		my $idLabel     = "-id";
		my $xtraOptions = "";
		if ( $usearchVer >= 10 ) {    #just to control id percentage
			$idLabel  = "";
			$id_OTUup = "";
			if ( $id_OTU != 0.97 ) {
				printL "UPARSE 10 does only support 97% id ${OTU_prefix} clusters\n", "w";
			}
		}
		elsif ( $usearchVer >= 8 ) {
			$idLabel  = "-otu_radius_pct";
			$id_OTUup = 100 - ( $id_OTU * 100 );
			if ( $id_OTUup < 0 || $id_OTUup > 50 ) {
				printL
				  "UPARSE cluster radius $id_OTUup not valid, aborting..\n", 54;
			}
		}
		if ( $usearchVer >= 8 ) {
			$xtraOptions .= " -uparseout $UCguide[0] ";    #-sizeout sizein
			$outputLab = "-fastaout";
		}
		if ( $usearchVer >= 8 && $usearchVer < 9 ) {       #8 specific commands
			$xtraOptions .=" -sizeout -sizein -uparse_maxhot $maxhot -uparse_maxdrop $maxdrop ";
		}
		if ( $noChimChk == 1 || $noChimChk == 2 ) {  #deactivate chimera de novo
			$xtraOptions .= " -uparse_break -999 ";
		}
		$entrMessage = "UPARSE core routine\nCluster at ". 100 * $id_OTU ;
		#"\n =========================================================================\n UPARSE core routine\n Cluster at ". 100 * $id_OTU. "%\n=========================================================================\n",0);

		$cmd = "$usBin -cluster_otus $derepl -otus $OTUfastaTmp $idLabel $id_OTUup -log $logDir/UPARSE.log $xtraOptions;";    # -threads $uthreads"; # -uc ".$UCguide[2]."
		$citations .= "UPARSE ${OTU_prefix} clustering - Edgar RC. 2013. UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nat Methods.\n";
		#die $cmd."\n";
	}
	elsif ( $ClusterPipe == 7 ) { #dada2
		die "incorrect dada2 script defined" unless (-f $dada2Scr);
		$cmd = "$Rscript $defRscriptOpt $dada2Scr $sdmDemultiDir $lotus_tempDir $dada2Seed $uthreads $cpMapFile $derepl;";
		#cleanup of important dada2 files
		$cmd .= "mv -f $lotus_tempDir/*.pdf $logDir;";
		
		$cmd .= "cp $lotus_tempDir/dada2.uc $UCguide[0];";
		$cmd .= "rm -f $OTUfastaTmp;cp $lotus_tempDir/uniqueSeqs.fna $OTUfastaTmp;";
		#$cmd .= "cp $sdmDemultiDir/uniqueSeqs.fna $outdir/primary/;gzip $outdir/primary//uniqueSeqs.fna;";
		$cmd .= "rm -r $sdmDemultiDir;" if ($saveDemulti==0);
		#die "$cmd\n";
		$citations .= "DADA2 ASV clustering - Callahan BJ, McMurdie PJ, Rosen MJ, et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 2016;13:581–3. doi:10.1038/nmeth.3869\n";
		$entrMessage = "DADA2 ASV clustering\ncheck progress at $progOutPut";
	}
	elsif ( $ClusterPipe == 6 ) { #unoise3
		$entrMessage = "UNOISE core routine\n Cluster at ". 100 * $id_OTU . "%";
		#printL("\n =========================================================================\n UNOISE core routine\n Cluster at ". 100 * $id_OTU . "%\n=========================================================================\n",0);

		$cmd = "$usBin -unoise3 $derepl -zotus $OTUfastaTmp -tabbedout $logDir/unoise3_longreport.txt -log $logDir/unoise3.log;";    # -threads $uthreads"; # -uc ".$UCguide[2]."
		$citations .= "UNOISE ASV (zOTU) clustering - R.C. Edgar (2016), UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing, https://doi.org/10.1101/081257 \n";
		#die $cmd."\n";
	} elsif ( $ClusterPipe == 2 ) {

	#prelim de novo OTU filter
		if ( !-e $swarmBin ) {printL "No valid swarm binary found at $swarmBin\n", 88;}
		$entrMessage = "SWARM ${OTU_prefix} clustering\n Cluster with d = ". $swarmClus_d ;
		#printL("\n =========================================================================\n SWARM ${OTU_prefix} clustering\n Cluster with d = ". $swarmClus_d. "\n=========================================================================\n",0);

		#-z: unsearch size output. -u uclust result file
		my $uclustFile = "$lotus_tempDir/clusters.uc";
		my $dofasti    = "-f ";
		if ( $swarmClus_d > 1 ) { $dofasti = ""; }
		$cmd = "$swarmBin -z $dofasti -u $uclustFile -t $uthreads -w $OTUfastaTmp --ceiling 4024 -s $logDir/SWARMstats.log -l $logDir/SWARM.log -o $lotus_tempDir/otus.swarm -d $swarmClus_d < $derepl;";
		$citations .= "swarm v2 ${OTU_prefix} clustering - Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. 2015. Swarm v2: highly-scalable and high-resolution amplicon clustering. PeerJ. DOI: 10.7717/peerj.1420\n";

	} elsif ( $ClusterPipe == 3 ) {
		if ( !-e $cdhitBin ) {printL "No valid CD-Hit binary found at $cdhitBin\n", 88;}
		$entrMessage = "CD-HIT ${OTU_prefix} clustering\n Cluster at ". 100 * $id_OTU ;
		#printL("\n =========================================================================\n CD-HIT ${OTU_prefix} clustering\n Cluster at ". 100 * $id_OTU . "%\n=========================================================================\n", 0 );
		if ($REFflag) {  #$otuRefDB eq "ref_closed" || $otuRefDB eq "ref_open"){
			printL "CD-HIT ref DB clustering not supported!\n", 55;
			#die();
			$cmd = "$cdhitBin-2d -T $uthreads -o $OTUfastaTmp.2 -c $id_OTU -M 0 -i2 $derepl -i $refDB4otus -n 9 -g 1;";#-aL 0.77 -aS 0.98
			if ( $otuRefDB eq "ref_open" ) {
				$cmd .= "$cdhitBin -T $uthreads -i $OTUfastaTmp.2 -c $id_OTU -M 0 -o $OTUfastaTmp.3 -n 9 -g 1 -aS 0.98;";#.3 are the denovo clusters
				die $cmd . "\n";
				$OTUfastaTmp = removeNonValidCDhits( $refDB4otus, "$OTUfastaTmp.2", "$OTUfastaTmp.3" );
			}
		}
		else {           #de novo
			$cmd = "$cdhitBin -T $uthreads -o $OTUfastaTmp -c $id_OTU -G 0 -M 0 -i $derepl -n 9 -g 0 -r 0 -aL 0.0 -aS 0.9;";          #-aL 0.77 -aS 0.98
		}
		$citations .= "CD-HIT ${OTU_prefix} clustering - Fu L, Niu B, Zhu Z, Wu S, Li W. 2012. CD-HIT: Accelerated for clustering the next-generation sequencing data. Bioinformatics 28: 3150–3152.\n";
	}elsif ( $ClusterPipe == 4 ) {    #dnaclust-ref
		my $dnaClOpt = "";
		if ($REFflag) {
			$dnaClOpt .="-l --approximate-filter --predetermined-cluster-centers $refDB4otus ";
			$dnaClOpt .= "--recruit-only " if ( $otuRefDB eq "ref_closed" );
		}
		else {
			printL "DNACLUST de novo clustering not supported in LotuS.\n", 34;
		}

		#ref_closed
		$cmd ="$dnaclustBin -i $derepl -s $id_OTU -t $uthreads --assign-ambiguous $dnaClOpt > $dnaclustOut ;";

		#die $cmd."\n";
		$citations .="DNACLUST - Ghodsi, M., Liu, B., & Pop, M. (2011). DNACLUST: accurate and efficient clustering of phylogenetic marker genes. BMC Bioinformatics, 12, 271. \n";
	}elsif ( $ClusterPipe == 5 ) {    #micca
		die "-CL 5 not supported\n";
		 #"$otuclustBin $derepl -s $id_OTU  --out-clust $otuclust_clust --out-rep $OTUfastaTmp -f fasta -c"; #-d: fast
		 #$citations.= "MICCA";
	}else {
		printL "Unkown \$ClusterPipe $ClusterPipe\n", 7;
	}
	#actual excecution
    if ( $exec == 0 ) {
		printL(frame($entrMessage,1,2),0);
        if ( systemL($cmd) != 0 ) {
            printL( "Failed core ${OTU_prefix} clustering command:\n$cmd\n", 1 );
        }
		
		#post cluster parsing
		my $d2rep;
		if  ( $ClusterPipe == 7 ) {#report within dada2 R script
			$d2rep = `grep 'Found .* ASVs.*(dada2)' $progOutPut`;
			if ($d2rep eq ""){
				$d2rep = `grep 'were chimeric and will be removed (DADA2)' $progOutPut`;
			}
			$exitMsg = "$d2rep";
		}
		printL(frame($exitMsg,1,3),0);
		#die "$progOutPut\n$d2rep\n";
    }
	
	#die $cmd;
    if ( $ClusterPipe == 2 ) {swarm4us_size($OTUfastaTmp); }

    #--------- ref based clustering ------------ 
    my ( $refsCL, $refCLSiz, $denovos, $denoSize );
    if ( $ClusterPipe == 4 ) {    #transcribe dnaclust out to .uc
        ( $refsCL, $refCLSiz, $denovos, $denoSize ) =
          dnaClust2UC( $dnaclustOut, $UCguide[0] );
    }
    if ($REFflag) {
        #1=add correct size tag to each (denovo uchime)
        my $fasRef = extractFastas( $refDB4otus, $refCLSiz, 1 );
        #and get denovo cluster centers separate
        my $fasRefDeno = extractFastas( $refDB4otus, $denoSize, 1 );
        writeFasta( $fasRef,     $OTUfastaTmp . ".ref" );
        writeFasta( $fasRefDeno, $OTUfastaTmp );
    }

#------------------ 2nd part: --------------------
    #backmap reads to OTUs/ASVs/zOTUs
    my $cnt = 0;
    my @allUCs;
    my $userachDffOpt = "--maxhits 1 --query_cov 0.8 --maxrejects 200 -top_hits_only -strand both -id $id_OTU -threads $uthreads";
    my $vsearchSpcfcOpt = "";#" --minqt 0.8 ";
    #$vsearchSpcfcOpt = " --minqt 0.3 " if ( $platform eq "454" );
    $vsearchSpcfcOpt = "" if ( !$VSused );
    $vsearchSpcfcOpt .= " --dbmask none --qmask none";
	my $backMapRep = "";
	
    if ($sdmDerepDo) {
        #only HQ derep need to be mapped, and uparse v8 has this already done
        #my @lotsFiles = ($derepl); #these are dereplicated files
        if (   ( $usearchVer < 8 && $ClusterPipe == 1 )
            || $ClusterPipe == 2 || $ClusterPipe == 6  || $ClusterPipe == 3  )
        {    #usearch8 has .up output instead
            $cmd = "$VSBin --usearch_global ". $derepl. " -db $outfile -uc $UCguide[0] $userachDffOpt $vsearchSpcfcOpt;";    #-threads  $BlastCores";
                   #die $cmd."\n";
            if ( systemL($cmd) != 0 ) {
                printL( "vsearch backmapping command aborted:\n$cmd\n", 1 );
            }
			$backMapRep .= "Backmapping unique reads: ".parseLastVUsearch($progOutPut)."\n"; 
			
        }
    } else {         #map all qual filtered files to OTUs
        if ( $usearchVer >= 8 ) {
            printL "\n\nWarning:: usearch v8 can potentially cause decreased lotus performance in seed extension step, when used with option -highMem 0\n\n", 0;
            sleep(10);
        }
        my @lotsFiles = ($filtered);
        my $dcnt      = 1;
        while ( -f $filtered . "." . $dcnt ) {
            push( @lotsFiles, $filtered . "." . $dcnt );
            $dcnt++;
        }
        foreach my $fil (@lotsFiles) {
            my $tmpUC = "$lotus_tempDir/tmpUCs.$cnt.uct";
            push( @allUCs, $tmpUC );
            $cmd = "$VSBin --usearch_global $fil -db $outfile -uc $tmpUC -strand both -id $id_OTU -threads $uthreads $vsearchSpcfcOpt;";    #-threads  $BlastCores";
            if ( systemL($cmd) != 0 ) { exit(1); }
            if ( $cnt > 0 ) {
                combine( $fil, $filtered );
                unlink($fil);
            }
            $cnt++;
        }
        if (
            systemL( "cat " . join( " ", @allUCs ) . " > " . $UCguide[0] ) != 0 )
        {
            printL "Merge of UC subfiles failed.\n", 1;
            exit(1);
        }
        foreach (@allUCs) { unlink; }
    }
    my @lotsFiles = ($filterOutAdd);
    my $dcnt      = 1;
    while ( -f $filterOutAdd . "." . $dcnt && $dcnt < 100000 ) {
        push( @lotsFiles, $filterOutAdd . "." . $dcnt );
        $dcnt++;
    }

 #if (@lotsFiles > 0){	systemL ("cat ".join(" ",@lotsFiles)." >>$filterOutAdd");}

    #add in mid qual reads
    foreach my $subF (@lotsFiles) {
        if ( -f $subF && !-z $subF ) {
            #make sure there's at least 2 lines
            open T, "<", $subF;
            my $lcnt = 0;
            while (<T>) { $lcnt++; last if ( $lcnt > 5 ); }
            close T;
            if ( $lcnt > 2 ) {    #file contains reads, so map
                #printL frame("Searching with mid qual reads..\n"), 0;
                my $tmpUC = "$lotus_tempDir/add.uc";
                $cmd ="$VSBin -usearch_global $subF -db $OTUfastaTmp -uc ". $tmpUC. " $userachDffOpt $vsearchSpcfcOpt;"; #-threads  $BlastCores";
                if ( -s $OTUfastaTmp ) {
                    if ( systemL($cmd) != 0 ) { printL( "Failed: $cmd\n", 1 ); }
                } else { systemL("touch $tmpUC"); }
				$backMapRep .= "Backmapping  mid qual reads: ".parseLastVUsearch($progOutPut)."\n"; 

                systemL( "cat $tmpUC >> " . $UCguide[0] . ".ADD" );
                unlink $tmpUC;

                #any ref DB to map onto??
                if ( ($REFflag) && -s "$OTUfastaTmp.ref" ) {
                    $cmd = "$VSBin -usearch_global $subF -db $OTUfastaTmp.ref -uc ". $tmpUC . " $userachDffOpt $vsearchSpcfcOpt;";    #-threads  $BlastCores";
                    if ( systemL($cmd) != 0 ) { printL( "Failed: $cmd\n", 1 ); }
                    systemL( "cat $tmpUC >> " . $UCguide[0] . ".ADDREF" );
                    unlink $tmpUC;
                }
            }
        }
    }

    #add in unique abundant reads
    if ( -s $derepl . ".rest" ) {
#push(@lotsFiles,$derepl.".rest"); #these are sdm "uniques" that were too small to map
#my @lotsFiles2 = ($derepl,$derepl.".rest");
        my $restUC = "$lotus_tempDir/rests.uc";
        $cmd =  "$VSBin -usearch_global ". $derepl . ".rest" . " -db $OTUfastaTmp -uc $restUC $userachDffOpt $vsearchSpcfcOpt;";#-threads  $BlastCores";
        if ( -s $OTUfastaTmp ) {
			printL(frame("Backmapping low qual reads to ${OTU_prefix}'s",1,2),0);
			if ( systemL($cmd) != 0 ) { exit(1); }
			printL(frame("Finished backmapping",1,3),0);
        }
        else { systemL("touch $restUC;"); }
        systemL( "cat $restUC >> " . $UCguide[0] . ".REST" );
        unlink $restUC;
        if ( ($REFflag) && -s "$OTUfastaTmp.ref" ) {
            $cmd ="$VSBin -usearch_global " . $derepl . ".rest" . " -db $OTUfastaTmp.ref -uc " . $restUC. " $userachDffOpt $vsearchSpcfcOpt;";    #-threads  $BlastCores";
            if ( systemL($cmd) != 0 ) { printL( "Failed: $cmd\n", 1 ); }
            systemL( "cat $restUC >> " . $UCguide[0] . ".RESTREF;" );
            unlink $restUC;
			$backMapRep .= "Backmapping  low count unique reads: ".parseLastVUsearch($progOutPut)."\n"; 

        }
    }
	
	printL frame($backMapRep),0;

#if (systemL("cat ".join(" ",@allUCs)." >> ".$UCguide[0]) != 0){printL "Merge of UC subfiles failed.\n",1;};
#foreach (@allUCs){unlink;}
    return ( \@UCguide );
}

#eg. $CONT_REFDB_PHIX
#removes all seqs from otusFA, that match at 95% to refDB
sub derepBash($) {

    #SWARM derep way
    my ($filtered) = @_;
    die("DERPRECATED bash derep\n");
    my $derepCmd = 'grep -v "^>" ' . $filtered . ' | \
	grep -v [^ACGTacgt] | sort -d | uniq -c | \
	while read abundance sequence ; do
		hash=$(printf "${sequence}" | sha1sum)
		hash=${hash:0:40}
		printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
	done | sort -t "_" -k2,2nr -k1.2,1d | \
	sed -e \'s/\_/\n/2\' > ' . "$lotus_tempDir/derep.fa";    #amplicons_linearized_dereplicated.fasta
    die $derepCmd . "\n";
    open( CMD, $derepCmd ) || printL( "Can't derep in swarm:\n$derepCmd\n", 1 );
    close CMD;
    return ("$lotus_tempDir/derep.fa");
}

sub swarmClust($) {
    my ($outfile) = @_;
    die("DEPRECATED swarm function\n");

    my $swarmThreads = $uthreads;
    my $swPath       = "/home/falhil/bin/swarm-dev/";
    my $swarmBin     = "$swPath/swarm";

    #my $swarmBreakBin = "python $swPath/scripts/swarm_breaker.py";
    my @UCguide = ( "$lotus_tempDir/finalOTU.uc", 2 );    #,"$lotus_tempDir/otus.uc",1);
    if ( $exec == 1 ) { return ( \@UCguide ); }

    my $filtered = $filterOut;
    my $derepl   = "$lotus_tempDir/derep.fas";            #,$totSeqs,$arL)
    if ( !-f $derepl || -z $derepl ) {
        printL "The sdm dereplicated output file was either empty or not existant, aborting lotus.\n$derepl\n", 1;
    }

    my ( $totSeqs, $SeqLength ) = parseSDMlog("$logDir/demulti.log");
    my @lotsFiles = ($filtered);
    my $dcnt      = 1;
    while ( -f $filtered . "." . $dcnt ) {
        push( @lotsFiles, $filtered . "." . $dcnt );
        $dcnt++;
    }
    if ( !$sdmDerepDo ) {
        my ($dereplFi) = usearchDerepSort($filtered);
        onelinerSWM($derepl);
        $derepl = $dereplFi;
    }

    return ( \@UCguide );
}

sub annotateFaProTax{
	my ($hir,$FaPr) = @_;
	return unless (defined($FaPr));
	
	 $citations .=
	"FaProTax (functional abundances based on OTUs) - Louca, S., Parfrey, L.W., Doebeli, M. (2016) - Decoupling function and taxonomy in the global ocean microbiome. Science 353:1272-1277\n";
	
}


sub systemW {
    my ($cc) = @_;
    my $ret = systemL $cc;
    if ($ret) { printL "Failed command $cc\n", 99; }
    printL("[sys]: $cc",0);
    return $ret;
}

