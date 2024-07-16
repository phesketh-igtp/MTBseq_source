#!/usr/bin/env perl

package TBresi;

use strict;
use warnings;
use File::Copy;
use TBtools;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

##needed for resi --> add to requirements!
use Statistics::R;
use Getopt::Long;
use File::Basename;
use File::Slurp;

$VERSION    =  1.0.0;
@ISA        =  qw(Exporter);
@EXPORT     =  qw(tbresi, tbresisummary, tbcombinedresi);

sub tbresi {
   # get parameter and input from front-end.
   my $logprint                  =  shift;
   my $VAR_dir                   =  shift;
   my $POS_OUT                   =  shift;
   my $STRAIN_OUT                =  shift;
   my $ref                       =  shift;
   my $refg                      =  shift;
   my $micovf                    =  shift;
   my $micovr                    =  shift;
   my $mifreq                    =  shift;
   my $miphred20                 =  shift;
   my $date_string               =  shift;
   #my $all_vars                  =  shift;
   my $all_vars                  =  1; # set to 1 for fetching all genome positions.
   my $snp_vars                  =  0; #deactivated the option
   my $lowfreq_vars              =  shift;
   my @position_tables           =  @_;
   my $genes                     =  {};
   my $annotation                =  {};
   my $phylo_positions_homolka   =  {};
   my $phylo_positions_coll      =  {};
   my $phylo_positions_beijing   =  {};
   my $phylo_positions           =  {};
   my $positions_data            =  {};
   my $output_file               =  "Strain_Classification.tab";
   my %check_up;
   # save already detected strains.
   if(-f "$STRAIN_OUT/$output_file") {
      open(IN,"$STRAIN_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $output_file: TBstrains.pm line: ", __LINE__ , " \n";
      <IN>;
      while(<IN>) {
         my $line       =  $_;
         $line          =~ s/\015?\012?$//;
         my @fields     =  split(/\t/);
         $check_up{$fields[1]."_".$fields[2]}   =  $fields[0];
      }
   }
   # parse the genomic sequence for determining substitutions.
   print $logprint "<INFO>\t",timer(),"\tStart parsing $ref...\n";
   my $genome  =  parse_fasta($logprint,$VAR_dir,$ref);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $ref!\n";
   print $logprint "<INFO>\t",timer(),"\tStart parsing $refg...\n";
   # parse annotation.
   parse_annotation($logprint,$VAR_dir,$genes,$annotation,$refg);
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $refg!\n";
   print $logprint "<INFO>\t",timer(),"\t","Start loading classification schemes...\n";
   # parse lineage classification.
   parse_classification($phylo_positions_homolka,$phylo_positions_coll,$phylo_positions_beijing);
   print $logprint "<INFO>\t",timer(),"\t","Finsihed loading classification schemes!\n";
   print $logprint "<INFO>\t",timer(),"\t","Start combining classification schemes...\n";
   # this hash will contain the data we need.
   foreach my $pos (keys %$phylo_positions_homolka) {
      $phylo_positions->{$pos}->{0}    =     $phylo_positions_homolka->{$pos};
   }
   foreach my $pos (keys %$phylo_positions_coll) {
      $phylo_positions->{$pos}->{0}    =     $phylo_positions_coll->{$pos};
   }
   foreach my $pos (keys %$phylo_positions_beijing) {
      $phylo_positions->{$pos}->{0}    =     $phylo_positions_beijing->{$pos};
   }
   print $logprint "<INFO>\t",timer(),"\t","Finished combining classification schemes!\n";
   unless(-f "$STRAIN_OUT/$output_file") {
      print $logprint "<INFO>\t",timer(),"\t","Start writing $output_file...\n";
      open(OUT,">$STRAIN_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't create $output_file: TBstrains.pm line: ", __LINE__ , " \n";
      my $header     =     "Date\tSampleID\tLibraryID\tFullID";
      $header        .=    "\tHomolka species\tHomolka lineage\tHomolka group\tQuality";
      $header        .=    "\tColl lineage (branch)\tColl lineage_name (branch)\tColl quality (branch)";
      $header        .=    "\tColl lineage (easy)\tColl lineage_name (easy)\tColl quality (easy)";
      $header        .=    "\tBeijing lineage (easy)\tBeijing quality (easy)";
      $header        .=    "\n";
      print OUT $header;
      close(OUT);
      print $logprint "<INFO>\t",timer(),"\t","Finished writing $output_file!\n";
   }
   open(OUT,">>$STRAIN_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't create $output_file: TBstrains line: ", __LINE__ , " \n";
   # start logic...
   foreach my $file (sort { $a cmp $b } @position_tables) {
      print $logprint "<INFO>\t",timer(),"\t","Start parsing $file...\n";
      $file       =~ /(\S+)\.gatk_position_table\.tab$/;
      my $id      =  $1;
      my @sample  =  split(/_/,$id);
      # check if strains classification already exists.
      if(exists $check_up{"\'$sample[0]"."_"."\'$sample[1]"}) {
         print $logprint "<INFO>\t",timer(),"\t","Skipping $file. Lineage classification already existing!\n";
         next;
      }
      # parse position table file to hash.
      my $position_table            =  {};
      my $phylo_position_table      =  {};
      parse_position_table($logprint,$POS_OUT,$file,$micovf,$micovr,$miphred20,$mifreq,$position_table);
      # skip every position not included in phylo positions.
      foreach my $pos (keys %$phylo_positions) {
         foreach my $index (keys %{$phylo_positions->{$pos}}) {
            $phylo_position_table->{$pos}->{$index}   =  $position_table->{$pos}->{$index};
         }
      }
      $position_table               =  {};
      print $logprint "<INFO>\t",timer(),"\t","Finished parsing $file!\n";
      print $logprint "<INFO>\t",timer(),"\t","Start fetching variants for $id...\n";
      my $variants                  =  {};
      my $statistics                =  {};
      call_variants($logprint,$phylo_position_table,$variants,$statistics,$micovf,$micovr,$miphred20,$mifreq,$annotation,$genes,$genome,$all_vars,$snp_vars,$lowfreq_vars);
      $statistics                   =  {};
      $phylo_position_table         =  {};
      $positions_data->{$id}        =  $variants;
      print $logprint "<INFO>\t",timer(),"\t","Finished fetching variants for $id!\n";
      print $logprint "<INFO>\t",timer(),"\t","Start writing lineage result for $id...\n";
      foreach my $id (sort { $a cmp $b } keys %$positions_data) {
         my $lineage                =  "unknown";
         my $species                =  "unknown";
         my $lineage_name           =  "unknown";
         my $quality_homolka        =  "good";
         my $quality_coll           =  "good";
         my $quality_beijing        =  "good";
         my $IDpositions            =  {};
         foreach my $pos (keys %$phylo_positions_homolka) {
            my $allel1              =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
            my $freq1               =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
            my $count1              =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
            $quality_homolka        =  "ugly"   unless($allel1 =~ /[AGCTagct]/);
            $quality_homolka        =  "bad"    unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
            $IDpositions->{$pos}    =  $allel1;
         }
         foreach my $pos (keys %$phylo_positions_coll) {
            my $allel1              =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
            my $freq1               =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
            my $count1              =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
            $quality_coll           =  "ugly"   unless($allel1 =~ /[AGCTagct]/);
            $quality_coll           =  "bad"    unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
            $IDpositions->{$pos}    =  $allel1;
         }
         foreach my $pos (keys %$phylo_positions_beijing) {
            my $allel1              =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[2];
            my $freq1               =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[6];
            my $count1              =  (split(/\t/, $positions_data->{$id}->{$pos}->{0}))[11];
            $quality_beijing        =  "ugly"   unless($allel1 =~ /[AGCTagct]/);
            $quality_beijing        =  "bad"    unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded frequency and coverage.
            $IDpositions->{$pos}    =  $allel1;
         }
         my $quality                =  0;
         my $outline                =  "\'$date_string\t\'$sample[0]\t\'$sample[1]\t'$id";
         ($species,$lineage_name)   =  specificator_homolka($IDpositions);
         $quality                   =  $quality_homolka;
         $lineage                   =  translate_homolka2coll($lineage_name);
         $species                   =  "\'" . $species;
         $lineage                   =  "\'" . $lineage;
         $lineage_name              =  "\'" . $lineage_name;
         $quality                   =  "\'" . $quality;
         $outline                   .= "\t$species\t$lineage\t$lineage_name\t$quality";
         ($species,$lineage)        =  specificator_coll_branch($IDpositions);
         $quality                   =  $quality_coll;
         $lineage_name              =  translate_coll2homolka($lineage);
         $species                   =  "\'" . $species;
         $lineage                   =  "\'" . $lineage;
         $lineage_name              =  "\'" . $lineage_name;
         $quality                   =  "\'" . $quality;
         $outline                   .= "\t$lineage\t$lineage_name\t$quality";
         ($species,$lineage)        =  specificator_coll_easy($IDpositions);
         $quality                   =  $quality_coll;
         $lineage_name              =  translate_coll2homolka($lineage);
         $species                   =  "\'" . $species;
         $lineage                   =  "\'" . $lineage;
         $lineage_name              =  "\'" . $lineage_name;
         $quality                   =  "\'" . $quality;
         $outline                   .= "\t$lineage\t$lineage_name\t$quality";
         ($species,$lineage)        =  specificator_beijing_easy($IDpositions);
         $quality                   =  $quality_beijing;
         $lineage_name              =  $lineage;
         $lineage                   =  "\'" . $lineage;
         $quality                   =  "\'" . $quality;
         $outline                   .= "\t$lineage\t$quality";
         $outline                   .= "\n";
         print OUT $outline;
      }
      $positions_data   =  {};
      $variants         =  {};
      print $logprint "<INFO>\t",timer(),"\t","Finished writing lineage result for $id!\n";
   }
   close(OUT);
   undef(%check_up);
}


sub tbresisummary {

	my $cutoff                   =      shift;
	my $qcutoff                  =      shift;
	my $fcutoff                  =      shift;
	my $pcutoff                  =      shift;
	my $mincutoff                =      shift;
	my @resi_files               =      shift;
	my $resi_file                =       {};
	my $file                     =       {};
	my $empty					 =		 {};	
	my @line_number				 =		 @_;
	my $line                     =       {};
	my @ID                       =       @_;
	my $outputfile               =       {};
	my %check_up;



  foreach my $resi_file (@resi_files) { #check whether the file is empty
	  $empty = "full";
	  my @line_number = read_file("$resi_file");
	  my $lines = @line_number;
	  if ($lines == 1){
		  $empty = "empty";
	  }
	  next if ($empty eq "empty");
	  
	$resi_file=~/^(.+).tab/ or die "strange file format: $resi_file\n";
	my $file=$1;
	open(Fout,">${file}_summary.tab") or die "\n\ncannot write output file\n\n\n";    
	print Fout "SampleID\tLibID\tINH\tFreq_INH\tRMP\tFreq_RMP\tSM\tFreq_SM\tEMB\tFreq_EMB\tPZA\tFreq_PZA\tMFX\tFreq_MFX\tLFX\tFreq_LFX\tCFZ\tFreq_CFZ\tKAN\tFreq_KAN\tAMK\tFreq_AMK\tCPR\tFreq_CPR\tETH/PTH\tFreq_ETH/PTH\tLZD\tFreq_LZD\tBDQ\tFreq_BDQ\tCS\tFreq_CS\tPAS\tFreq_PAS\tDLM\tFreq_DLM\tPrediction\n";
    @ID=split("_",$file);
    my %mutations;
    open (Fin,"<$resi_file") or die "\n\ncannot open $file\n\n\n";
	my $R = Statistics::R->new();
       $R->startR;
       $R->send('library("readr")');#R-Library for reading in tab-seperated files
        
       $R->run(qq'table<-read_delim(file="${resi_file}", "\t", escape_double =F, trim_ws=T);'); #Read in Table with resistance-calls 
       @names=("INH-R (INH)","RIF-R (RMP)","SM-R (SM)","EMB-R (EMB)","PZA-R (PZA)","MFX-R (MFX)","LFX-R (LFX)","CFZ-R (CFZ)","KAN-R (KAN)","AMI-R (AMK)","CAP-R (CPR)","ETH-R (ETH)","LZD-R (LZD)","BDQ-R (BDQ)","CS-R (CS)","PAS-R (PAS)","DEL-R (DEL)");
            foreach $a (@names){#loop over all antibiotics
    
			$R->run(qq'name<-"$a"
                cutoff<-as.numeric("$cutoff")
                qcutoff<-as.numeric("$qcutoff")
                fcutoff<-as.numeric("$fcutoff")
				pcutoff<-as.numeric("$pcutoff")
				mincutoff<-as.numeric("$mincutoff")
                short<-as.character(unlist(strsplit("$a"," ",fixed=TRUE))[1])
                anti<-as.character(unlist(strsplit("$a"," ",fixed=TRUE))[2])
                antibiotic<-gsub("[()]","", anti)');
			$name=$R->get('short');
			$antibiotic=$R->get('antibiotic');
        
			print "Search for mutations mediating resistance to $antibiotic\n";   
			$R->run('table1<-table[which(table$CovFor >=cutoff & table$CovRev >=cutoff & table$Qual20 >=qcutoff & table$Freq >=fcutoff & table$Freq <=25 & table$Cov >= mincutoff & (table$Qual20/(table$CovRev + table$CovFor)*100) >= pcutoff),]
				table2<-table[which(table$CovFor >=cutoff & table$CovRev >=cutoff & table$Qual20 >=qcutoff & table$Freq >=fcutoff & table$Freq >25 & table$Cov >= mincutoff),]
				table3<-rbind(table1, table2)
                idx<-table3[grep(short,table3$better_resi),]
                count<-dim(idx)[1]'); #subtable with all lines, that have an entry in the better_resi column for the specific drug with "-R" and fullfilling the thresholds
			$dim=$R->get('count');
        

            if($dim >= 1){#check if there are any resistance calls
            $R->run(q'mut<-1
                     freq<-1
                if(dim(idx)[1]>1){
                    for(i in 1:dim(idx)[1]){
                        
                            mut[i]<-idx$Mutation_Annotation[i];
                            freq[i]<-round(idx$Freq[i],digits=2);
                    }
                    if(length(mut)==2){mut<-paste(mut[1],mut[2],sep="; ")
                                       freq<-paste(freq[1],freq[2],sep="; ")
                    }else if (length(mut)==3){mut<-paste(mut[1],mut[2],mut[3],sep="; ")
                                              freq<-paste(freq[1],freq[2], freq[3],sep="; ")
                    }else if (length(mut)==4){mut<-paste(mut[1],mut[2],mut[3], mut[4],sep="; ")
                                            freq<-paste(freq[1],freq[2], freq[3], freq[4],sep="; ")
                    }else if (length(mut)==5){mut<-paste(mut[1],mut[2],mut[3],mut[4], mut[5],sep="; ")
                                              freq<-paste(freq[1],freq[2], freq[3], freq[4],freq[5],sep="; ")
                    }else if (length(mut)==6){mut<-paste(mut[1],mut[2],mut[3],mut[4], mut[5], mut[6],sep="; ")
                                            freq<-paste(freq[1],freq[2], freq[3], freq[4],freq[5],freq[6],sep="; ")
                    }else {mut<-paste(mut[1],mut[2],mut[3],mut[4], mut[5], mut[6],"and others",sep="; ")
                    freq<-paste(freq[1],freq[2],freq[3],freq[4], freq[5], freq[6],"and others",sep="; ")}
                }else{
                        mut<-idx$Mutation_Annotation
                        freq<-round(idx$Freq, digits=2)
                    }');
            
            $mutations{$antibiotic}=$R->get('mut');
            $mutations{$antibiotic.'_Freq'}=$R->get('freq');
           
            }

            else{$mutations{$antibiotic}="-";
                $mutations{$antibiotic.'_Freq'}="-";}

        
        }
    if ($mutations{INH} ne "-" && $mutations{RMP} ne "-" && ($mutations{MFX} ne "-" || $mutations{LFX} ne "-") && ($mutations{BDQ} ne "-" || $mutations{LZD} ne "-" )){$prediction = "XDR";
    }elsif ($mutations{INH} ne "-" && $mutations{RMP} ne "-" && (($mutations{MFX} ne "-" || $mutations{LFX} ne "-"))){$prediction = "preXDR";
    }elsif ($mutations{INH} ne "-" && $mutations{RMP} ne "-"){$prediction = "MDR";
    }elsif ($mutations{RMP} ne "-") {$prediction = "RR";
    }elsif ($mutations{INH} eq "-" && $mutations{RMP} eq "-" && $mutations{SM} eq "-" && $mutations{EMB} eq "-" && $mutations{PZA} eq "-" && $mutations{MFX} eq "-" && $mutations{LFX} eq "-" && $mutations{CFZ} eq "-" && $mutations{KAN} eq "-" && $mutations{AMK} eq "-" && $mutations{CPR} eq "-" && $mutations{ETH} eq "-" && $mutations{LZD} eq "-" && $mutations{BDQ} eq "-" && $mutations{CS} eq "-" && $mutations{PAS} eq "-" && $mutations{DEL} eq "-"){$prediction = "S";
    }else{$prediction = "nonMDR";}
    print Fout "$ID[0]\t$ID[1]\t$mutations{INH}\t$mutations{INH_Freq}\t$mutations{RMP}\t$mutations{RMP_Freq}\t$mutations{SM}\t$mutations{SM_Freq}\t$mutations{EMB}\t$mutations{EMB_Freq}\t$mutations{PZA}\t$mutations{PZA_Freq}\t$mutations{MFX}\t$mutations{MFX_Freq}\t$mutations{LFX}\t$mutations{LFX_Freq}\t$mutations{CFZ}\t$mutations{CFZ_Freq}\t$mutations{KAN}\t$mutations{KAN_Freq}\t$mutations{AMK}\t$mutations{AMK_Freq}\t$mutations{CPR}\t$mutations{CPR_Freq}\t$mutations{ETH}\t$mutations{ETH_Freq}\t$mutations{LZD}\t$mutations{LZD_Freq}\t$mutations{BDQ}\t$mutations{BDQ_Freq}\t$mutations{CS}\t$mutations{CS_Freq}\t$mutations{PAS}\t$mutations{PAS_Freq}\t$mutations{DEL}\t$mutations{DEL_Freq}\t$prediction\n";
    $R->stopR();
    close Fin;
    close Fout;
 }
}
sub tbcombinedresi{
	my $output_file                =  "Strain_Resistance.tab";
	my $RESI_OUT                   =  shift;
	my @resisum_files              =  shift;
	
	if(-f "$RESI_OUT/$output_file") {
      open(IN,"$RESI_OUT/$output_file") || die print $logprint "<ERROR>\t",timer(),"\tCan't open $output_file: TBresi.pm line: ", __LINE__ , " \n";
      <IN>;
      while(<IN>) {
         my $line       =  $_;
         $line          =~ s/\015?\012?$//;
         my @fields     =  split(/\t/);
         $check_up{$fields[1]."_".$fields[2]}   =  $fields[0];
      }
   }
}
1;
