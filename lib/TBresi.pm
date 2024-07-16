#!/usr/bin/env perl

package TBresi;

use strict;
use warnings;
use File::Copy;
use File::Basename;
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
@EXPORT     =  qw(tbresi);

sub tbresi {
   # get parameter and input from front-end.
   my $logprint                  =  shift;
   my $VAR_dir                   =  shift;
   my $CALL_OUT                  =  shift;
   my $RESI_OUT                  =  shift;
   my $date_string               =  shift;
   my $resi_list_master          =  shift;
   my @truecodon_files           =  @_;
   my $minion                    =     "";
   my $res_hash = {};
   my $change_hash ={};
   

print "$resi_list_master\n";
   print $logprint ("<INFO>\t",timer(),"\tNo resistance file $resi_list_master. Will skip resistance annotation.\n") unless(-f $resi_list_master);
   #if($resi_list_master            eq   ''||"NONE"   ) { die "\n[ERROR]\t",timer(),"\tNeed to provide a resistance file. Use --help for usage information\n";}
   my $res_table_date = "custom";
   $res_table_date = $1 if ($resi_list_master) =~/^.*(\d\d\d\d.\d\d.\d\d)/;
   #print "$res_table_date\n";
   
   ###aminoacid changes###
   my %codon = ('Ser'=>'S',
			'Phe'=>'F',
			'Leu'=>'L',
			'Tyr'=>'Y',
			'_'=>'_',
			'Cys'=>'C',
			'Trp'=>'W',
			'Leu'=>'L',
			'Pro'=>'P',
			'His'=>'H',
			'Gln'=>'Q',
			'Arg'=>'R',
			'Ile'=>'I',
			'Met'=>'M',
			'Thr'=>'T',
			'Asn'=>'N',
			'Lys'=>'K',
			'Val'=>'V',
			'Ala'=>'A',
			'Asp'=>'D',
			'Glu'=>'E',
			'Gly'=>'G'
			);
   print $logprint "<INFO>\t",timer(),"\tStart parsing $resi_list_master...\n";
   if ($resi_list_master ne ""||"NONE"){
   open(Fin, "<", $resi_list_master) or die "\n[ERROR]\t",timer(),"\tUnable to open $resi_list_master\n";
   <Fin>;
   while (my $line = <Fin>){
   $line          =~  s/\015?\012?$//;
   my @line = split("\t",$line);
   my $pos = $line[0];
   ($pos = $pos) =~s/[i|\+]\d+$//;
   my $ref_allel = $line[1];
   my $mut_allel = $line[2];
   #if ($mut_allel =~ /^INS_/){print "$mut_allel\n";}
   $mut_allel =~s/^INS_//;
   #if ($mut_allel =~ /A|G|T|C/){print "$mut_allel\n";}
   my $type = $line[3];
   my $affects = $line[4];
   my $gene = $line[5];
   my $orientation = $line[6];
   my $gene_name = $line[7];
   if ($gene_name eq "gidB"){$gene_name="gid"};
   my $change = $line[8];
   my $predict_phyl = $line[9];
   my $predict = $line[10];
   my $comment = $line[13];
   my $mic = $line[14];
   my $whotier = $line[19];

	if ($pos ne "NA"){
		$change_hash->{$pos}->{$mut_allel} =   {change       =>   $change,
												predict      =>   $predict,
												predict_phyl =>   $predict_phyl,
												affects      =>   $affects,
												mut_allel    =>   $mut_allel,
												gene_name    =>   $gene_name,
												gene         =>   $gene,
												comment      =>   $comment};
	}
	$res_hash->{$gene}->{$change} = {predict => $predict,
									predict_phyl => $predict_phyl,
									gene_name => $gene_name,
									affects => $affects,
									comment => $comment};
	}
   print $logprint "<INFO>\t",timer(),"\tFinished parsing $resi_list_master!\n";
   close (Fin);
   }

   foreach my $variant_file (@truecodon_files) {
      print $logprint ("<INFO>\t",timer(),"\tNo variant file $variant_file. Will skip resistance annotation.\n") unless(-f "$CALL_OUT/$variant_file");
      print $logprint ("<INFO>\t",timer(),"\tStart calling resistance for $variant_file...\n");
		if ($variant_file =~ /.*gatk_position_true-codon-variants.*/){
			(my $out_file = basename($variant_file)) =~ s/\..*$//g;
			my $output_mode = $1 if ($variant_file) =~/^.*outmode(\d\d\d)/;
			open(OUT,">","$RESI_OUT/${out_file}.gatk_position_true-codon-variants_outmode${output_mode}_${res_table_date}_res.tsv") or die "\n<ERROR>\t",timer(),"\tUnable to create $${out_file}.gatk_position_true-codon-variants_outmode${output_mode}_${res_table_date}_res.tsv\n";

			open (IN, "<", "$CALL_OUT/$variant_file") or die "\n[ERROR]\t",timer(),"\tUnable to open $variant_file\n";
			my $header = <IN>;
			$header =~ s/\015?\012?$//;
			$header .= "\tbetter_resi\tBenign\tMutation_Annotation\tComment";
			print OUT "$header\n";

			while (my $line = <IN>) {
			chomp($line);
			next unless $line;
			next if $line =~ /^[A-Za-z]/;
			$line =~ s/\015?\012?$//;
			my @line = split("\t", $line);
			my $pos			= $line[0];		$pos 			= 0 unless($pos);
			my $pos_insertion = $pos;
			$pos_insertion =~s/[i|\+]\d+$//;
			my $ref			= $line[1];		$ref 			= 0 unless($ref);
			my $type		= $line[2];		$type 			= 0 unless($type);
			next if ($minion and $type eq "Del");
			my $allel		= $line[3];		$allel 			= 0 unless($allel);
			my $cov_forward = $line[4];		$cov_forward 	= 0 unless($cov_forward);
			my $cov_reverse = $line[5];		$cov_reverse 	= 0 unless($cov_reverse);
			my $qual_20		= $line[6];		$qual_20 		= 0 unless($qual_20);
			my $freq1		= $line[7];		$freq1 			= 0 unless($freq1);
			my $coverage	= $line[8];		$coverage 		= 0 unless($coverage);
			my $subs		= $line[9];		$subs 			= 0 unless($subs);
			my $subst		= $subs;
			$subst=~s/\s.*$//;
			#print "$subst\n";
			if (($subst ne "-") and $subst !~ /^\s*$/){
				my $old_aa;
				my $new_aa;
					if (substr($subst, 0, 3) !~ /\d+/) {
						$old_aa = substr( $subst, 0, 3 );
						# print "$old_aa\n";
						$old_aa = $codon{$old_aa};
						#print "$old_aa\n";
						}
					else{$old_aa = substr($subst, 0, 1)}

					if (substr($subst, -3) !~ /\d+/){
						$new_aa = substr( $subst, -3);
						$new_aa = $codon{$new_aa};
					}
					else{$new_aa = substr($subst, -1)}

				$subst =~ s/[^0-9]//g;
				$subst = $old_aa.$subst.$new_aa;
			}
		my $gene		= $line[10];	$gene 			= 0 unless($gene);
		my $gene_name	= $line[11];	$gene_name 		= 0 unless($gene_name);
		if ($gene_name eq "-"){$gene_name = $gene};
		if ($gene_name eq "gidB"){$gene_name = "gid"};
		my $annotation	= $line[12];	$annotation 	= 0 unless($annotation);
		#my $resistance	= $line[11];	$resistance 	= 0 unless($resistance);
		#my $phylo		= $line[12];	$phylo 			= 0 unless($phylo);
		my $region		= $line[13];	$region 		= 0 unless($region);
		my $warning		= $line[14];	$warning 		= "-" unless($warning);

		my $gene_position;
		my $wt_allel = $ref;
		my $subst_allel = $allel;
		#print "$subst_allel\n";
		my $better_res = "-";
		my $better_res_change = "-";
		my $benigninfo = "-";
		my $res_comment = "-";

#if ($gene_name eq "rrl"){print Dumper(\$res_hash->{$gene_name});}
			if (exists($res_hash->{$gene}->{$subst})){
			#if ($gene_name eq "rrl"){print Dumper(\$res_hash->{$gene_name});}
				if ($res_hash->{$gene}->{$subst}->{predict} ne "benign" and $res_hash->{$gene}->{$subst}->{predict} ne "likely benign"){
					$better_res = $res_hash->{$gene}->{$subst}->{predict};
					$gene_name = $res_hash->{$gene}->{$subst}->{gene_name};
					$better_res_change = $gene_name." ".$subst;
					$res_comment = $res_hash->{$gene}->{$subst}->{comment};
				}
				elsif ($res_hash->{$gene}->{$subst}->{predict} eq "benign" or $res_hash->{$gene}->{$subst}->{predict} eq "likely benign"){
					$benigninfo = $res_hash->{$gene}->{$subst}->{predict};
					$gene_name = $res_hash->{$gene}->{$subst}->{gene_name};
					$better_res_change = $gene_name." ".$subst;
					$res_comment = $res_hash->{$gene}->{$subst}->{comment};
				}
			}

			elsif (exists($change_hash->{$pos_insertion}->{$subst_allel})){
				if ($change_hash->{$pos_insertion}->{$subst_allel}->{predict} ne "benign" and $change_hash->{$pos_insertion}->{$subst_allel}->{predict} ne "likely benign" and ($subst_allel=~ /^((?!NA).)*[A|T|C|G]$/gm or $subst_allel eq "GAP")){
					$better_res = $change_hash->{$pos_insertion}->{$subst_allel}->{predict};
					$gene_name = $change_hash->{$pos_insertion}->{$subst_allel}->{gene_name};
					$better_res_change = $gene_name." ".$change_hash->{$pos_insertion}->{$subst_allel}->{change};
					$res_comment = $change_hash->{$pos_insertion}->{$subst_allel}->{comment};
				}
				elsif ($change_hash->{$pos_insertion}->{$subst_allel}->{predict} eq "benign" or $change_hash->{$pos_insertion}->{$subst_allel}->{predict} eq "likely benign" and ($subst_allel=~ /^((?!NA).)*[A|T|C|G]$/gm or $subst_allel eq "GAP")){
					#if ($pos_insertion == 1474959){print Dumper(\$change_hash->{$pos_insertion});print "$subst_allel\n";print "$pos_insertion\n";}
					$benigninfo = $change_hash->{$pos_insertion}->{$subst_allel}->{predict};
					$gene_name = $change_hash->{$pos_insertion}->{$subst_allel}->{gene_name};
					$better_res_change = $gene_name." ".$change_hash->{$pos_insertion}->{$subst_allel}->{change};
					$res_comment = $change_hash->{$pos_insertion}->{$subst_allel}->{comment};
				}
			}
			###general rules for dealing with Del, Ins and stop codons in certain genes
			
			###pncA reverse
			elsif ($gene_name eq "pncA" and $type eq "Del"){
			$gene_position = (2289241-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res = "PZA-R";
			}
			elsif ($gene_name eq "pncA" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (2289241-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res = "PZA-R";
				}
			}
			elsif ($gene_name eq "pncA" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="PZA-R";
				$better_res_change = $subst;
				}
				elsif($subst=~/(^\w).*(\w$)/ and $1 ne $2){
				$better_res ="PZA-R#";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			###katG reverse
			elsif ($gene_name eq "katG" and $type eq "Del"){
			#$res = 1;
			$gene_position = (2156111-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="INH-R";
			}
			elsif ($gene_name eq "katG" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (2156111-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="INH-R";
				}
			}
			elsif ($gene_name eq "katG" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="INH-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			###ethA reverse
			elsif ($gene_name eq "ethA" and $type eq "Del"){
			$gene_position = (4327473-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="ETH-R";
			}
			elsif ($gene_name eq "ethA" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (4327473-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="ETH-R";
				}
			}
			elsif ($gene_name eq "ethA" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="ETH-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			##gidB reverse - gene name = gid in variant file
			elsif (($gene_name eq "gid" or $gene_name eq "gidB") and $type eq "Del"){
			$gene_position = (4408202-$pos)+1;
			$wt_allel =~ tr/AGTC/tcag/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="SM-R";
			}
			elsif (($gene_name eq "gid" or $gene_name eq "gidB") and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = (4408202-$pos_insertion)+1;
				$subst_allel =~ tr/AGTC/tcag/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="SM-R";
				}
			}
			elsif (($gene_name eq "gid" or $gene_name eq "gidB") and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="SM-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## rpoB forward
			elsif ($gene_name eq "rpoB" and $type eq "Del"){
			$gene_position = ($pos-759807)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="RIF-R";
			}
			elsif ($gene_name eq "rpoB" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-759807)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="RIF-R";
				}
			}
			elsif ($gene_name eq "rpoB" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="RIF-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## ald forward
			elsif ($gene_name eq "ald" and $type eq "Del"){
			$gene_position = ($pos-3086820)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="CS-R";
			}
			elsif ($gene_name eq "ald" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-3086820)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="CS-R";
				}
			}
			elsif ($gene_name eq "ald" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="CS-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## Rv0678 forward
			elsif ($gene eq "Rv0678" and $type eq "Del"){
			$gene_position = ($pos-778990)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="BDQ-R,CFZ-R";
			}
			elsif ($gene eq "Rv0678" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-778990)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="BDQ-R,CFZ-R";
				}
			}
			elsif ($gene eq "Rv0678" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="BDQ-R,CFZ-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			
			## ddn  Rv3547 forward
			elsif ($gene_name eq "ddn" and $type eq "Del"){
			$gene_position = ($pos-3986844)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="DEL-R";
			}
			elsif ($gene_name eq "ddn" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-3986844)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="DEL-R";
				}
			}
			elsif ($gene_name eq "ddn" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="DEL-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}
			## tlyA Rv1694 forward
			elsif ($gene_name eq "tlyA" and $type eq "Del"){
			$gene_position = ($pos-1917940)+1;
			$wt_allel =~ tr/AGTC/agtc/;
			$better_res_change = $gene_name." ".$gene_position."_del_".$wt_allel;
			$better_res ="CAP-R";
			}
			elsif ($gene_name eq "tlyA" and $type eq "Ins"){
			my $subst_allel = $allel;
				if ($subst_allel=~/A|G|T|C/){
				$gene_position = ($pos_insertion-1917940)+1;
				$subst_allel =~ tr/AGTC/agtc/;
				$better_res_change = $gene_name." ".$gene_position."_ins_".$subst_allel;
				$better_res ="CAP-R";
				}
			}
			elsif ($gene_name eq "tlyA" and $type eq "SNP"){
				if ($subst=~/.*_$/){
				$better_res ="CAP-R";
				$better_res_change = $gene_name." ".$subst;
				}
			}

   print OUT "$pos\t$ref\t$type\t$allel\t$cov_forward\t$cov_reverse\t$qual_20\t$freq1\t$coverage\t$subs\t$gene\t$gene_name\t$annotation\t$region\t$warning\t$better_res\t$benigninfo\t$better_res_change\t$res_comment\n";

			}
   close (IN);
   close (OUT);
   @truecodon_files    =  ();
   $res_hash           =  {};
   $change_hash        =  {};
		}
	}
}
1;