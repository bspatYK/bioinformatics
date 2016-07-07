#!/usr/bin/evn perl -w
use strict;
use warnings;

my $binsts_file = $ARGV[0];
my $sample  = $ARGV[1];
my $bin_size = $ARGV[2];
my $step = $ARGV[3];
my $sr_GC_file = $sample."_".$bin_size."K_".$step."_sr_GC.statastics";
my $gn_GC_file = $sample."_".$bin_size."K_".$step."_gGC.statastics";
my $GC_corrected_sts  = $sample."_".$bin_size."K_corrected1.statastics";
my $GC_corrected_GC_sts  = $sample."_".$bin_size."K_".$step."_GC_corrected_GC.statastics";

my $out_dir = (@ARGV == 5) ? "$ARGV[4]" : "$ENV{'PWD'}";
open SRGCSTS, ">$out_dir/$sr_GC_file";
open GNGCSTS, ">$out_dir/$gn_GC_file";
open (GCcorrectedsts, ">$out_dir/$GC_corrected_sts") or die "Cannot open $out_dir/$GC_corrected_sts, $!";
open GCcorrectedGCsts, ">$out_dir/$GC_corrected_GC_sts";
print SRGCSTS "GC_content\tbin_count\tSum\tMean\tMin\tQ1\tMedian\tQ3\tMax\treads_count_every_bin\n";
print GNGCSTS "GC_content\tbin_count\tSum\tMean\tMin\tQ1\tMedian\tQ3\tMax\treads_count_every_bin\n";
print GCcorrectedsts "binID\tchr\tchr_region\tCorrected_Relative_Reads_Number1\n";
print GCcorrectedGCsts "GC_content\tbin_count\tSum\tMean\tMin\tQ1\tMedian\tQ3\tMax\treads_count_every_bin\n";

my %reads_count_of_sr_GC = ();
my %reads_count_of_gGC = ();
my $raw_available_bin_count = 0;
my $total_mapped_reads =0;
my %w_of_bins = ();
my @raw_mapped_reads_number = ();

open (BINSTS, "$binsts_file") or die "Cannot open $binsts_file, $!";

while (<BINSTS>){
	chomp;
#binID	chr	chr_region	all_region	reads_count	reads_GC	reads_bases	reads_GC_content	gGC_content
#1	chr1	1-1000000	1-1000000	1515	75855	148691	0.510151925805866	0.468094680839666
#2	chr1	1000001-2000000	1000001-2000000	3948	215109	384850	0.558942445108484	0.568914
	next if (/^binID/);
	my @line = split /\t/;
	#ignore no reads bins
	next if ($line[4] == 0);
	#ignore bins with 0 GC bases
	next if ($line[5] == 0);
	#ignore bins with 100%  GC bases
	next if ($line[7] == 1);
	#for autosome
	next if ($line[1] !~ /^chr\d+$/);
	$raw_available_bin_count ++;
	$total_mapped_reads += $line[4];
	push (@raw_mapped_reads_number, $line[4]);
	my $reads_step = (int ($line[7] / $step) + 1) * $step;
	my $genome_step = (int ($line[8] / $step) + 1) * $step;
	$reads_count_of_sr_GC{$reads_step} .= "$line[4]"."|";
	$reads_count_of_gGC{$genome_step} .= "$line[4]"."|";
	}
close BINSTS or die "Cannot close $binsts_file, $!";

my $rmmrn1 = $total_mapped_reads / $raw_available_bin_count;
my @raw_mapped_reads_number_statistics = statistics(@raw_mapped_reads_number);
my $rmmrn = $raw_mapped_reads_number_statistics[7];


for (my $i = 0; $i <= 1 / $step; $i += 1){
	my $j = $i * $step;
	if (! $reads_count_of_sr_GC{$j}){
		$reads_count_of_sr_GC{$j} = 0 ;
	}else{
		chop $reads_count_of_sr_GC{$j};
		}
	if (! $reads_count_of_gGC{$j}){
		$reads_count_of_gGC{$j} = 0 ;
	}else{
		chop $reads_count_of_gGC{$j};
		}
		my @reads_count_of_sr_GC = split /\|/, $reads_count_of_sr_GC{$j};
		my @sr_statistics = statistics(@reads_count_of_sr_GC);
		my $sr_count   = (($sr_statistics[0] == 1) && ($reads_count_of_sr_GC[0] == 0)) ? 0 : $sr_statistics[0] ;
		my $sr_sum     = $sr_statistics[1];
		my $sr_mean    = $sr_statistics[2];
		my $sr_min     = $sr_statistics[5];
		my $sr_Q1      = $sr_statistics[6];
		my $sr_median  = $sr_statistics[7];
		my $sr_Q3      = $sr_statistics[8];
		my $sr_max     = $sr_statistics[9];
		print SRGCSTS "$j\t$sr_count\t$sr_sum\t$sr_mean\t$sr_min\t$sr_Q1\t$sr_median\t$sr_Q3\t$sr_max\t$reads_count_of_sr_GC{$j}\n";
		my @reads_count_of_gGC = split /\|/, $reads_count_of_gGC{$j};
		my @gs = statistics(@reads_count_of_gGC);
		my $genome_count  = (($gs[0] == 1) && ($reads_count_of_gGC[0] == 0)) ? 0 : $gs[0] ;
		my $genome_sum    = $gs[1];
		my $genome_mean   = $gs[2];
		my $genome_min    = $gs[5];
		my $genome_Q1     = $gs[6];
		my $genome_median = $gs[7];
		my $genome_Q3     = $gs[8];
		my $genome_max    = $gs[9];
		print GNGCSTS "$j\t$genome_count\t$genome_sum\t$genome_mean\t$genome_min\t$genome_Q1\t$genome_median\t$genome_Q3\t$genome_max\t$reads_count_of_gGC{$j}\n";
		if ($sr_mean + $genome_mean == 0){
			$w_of_bins{$j} = 0;
		}elsif ($sr_mean * $genome_mean == 0){
			$w_of_bins{$j} = $rmmrn1 / ($sr_mean + $genome_mean) ;
		}else{
			$w_of_bins{$j} = $rmmrn1 / (($sr_mean + $genome_mean) / 2);
			}
		my @corrected_reads_count_of_sr_GC  = ();
		for my $rc(@reads_count_of_sr_GC){
			my $rc_corrected = $rc * $w_of_bins{$j};
			push @corrected_reads_count_of_sr_GC, $rc_corrected;
			}
		my $corrected_reads_count_of_sr_GC = join '|', @corrected_reads_count_of_sr_GC;
		print $corrected_reads_count_of_sr_GC."\n";
		my @cs = statistics(@corrected_reads_count_of_sr_GC);
		my $corrected_count  = (($cs[0] == 1) && ($corrected_reads_count_of_sr_GC[0] == 0)) ? 0 : $cs[0];
		my $corrected_sum    = $cs[1];
		my $corrected_mean   = $cs[2];
		my $corrected_min    = $cs[5];
		my $corrected_Q1     = $cs[6];
		my $corrected_median = $cs[7];
		my $corrected_Q3     = $cs[8];
		my $corrected_max    = $cs[9];
		print GCcorrectedGCsts "$j\t$corrected_count\t$corrected_sum\t$corrected_mean\t$corrected_min\t$corrected_Q1\t$corrected_median\t$corrected_Q3\t$corrected_max\t$corrected_reads_count_of_sr_GC\n";		
	}
close SRGCSTS;
close GNGCSTS;
close GCcorrectedGCsts;

my $corrected_total_mapped_reads = 0;
my $corrected_available_bin_count = 0;
my @corrected_mapped_reads_number = ();

open (BINSTS, "$binsts_file") or die "Cannot open $binsts_file, $!";
while (<BINSTS>){
	chomp;
#binID	chr	chr_region	all_region	reads_count	reads_GC	reads_bases	reads_GC_content	gGC_content
#1	chr1	1-1000000	1-1000000	1515	75855	148691	0.510151925805866	0.468094680839666
#2	chr1	1000001-2000000	1000001-2000000	3948	215109	384850	0.558942445108484	0.568914
	next if (/^binID/);
	my @line = split /\t/;
	#ignore no reads bins
	next if ($line[4] == 0);
	#ignore bins with 0 GC bases
	next if ($line[5] == 0);
	#ignore bins with 100%  GC bases
	next if ($line[7] == 1);
	#for autosome
	next if ($line[1] !~ /^chr\d+$/);
	my $reads_step = (int ($line[7] / $step) + 1) * $step;
	my $w = exists $w_of_bins{$reads_step} ? $w_of_bins{$reads_step} : 0 ;
	$corrected_total_mapped_reads += $w * $line[4];
	$corrected_available_bin_count ++ if ($w * $line[4] > 0);
	push (@corrected_mapped_reads_number, $w * $line[4]) if ($w * $line[4] > 0);
	}
close BINSTS or die "Cannot close $binsts_file, $!";

my $corrected_mean_mapped_reads_number = $corrected_total_mapped_reads / $corrected_available_bin_count;
my @corrected_mapped_reads_number_statistics = statistics(@corrected_mapped_reads_number);
my $cmmrn = $corrected_mapped_reads_number_statistics[7];

open (BINSTS, "$binsts_file") or die "Cannot open $binsts_file, $!";
while (<BINSTS>){
	chomp;
#binID	chr	chr_region	all_region	reads_count	reads_GC	reads_bases	reads_GC_content	gGC_content
#1	chr1	1-1000000	1-1000000	1515	75855	148691	0.510151925805866	0.468094680839666
#2	chr1	1000001-2000000	1000001-2000000	3948	215109	384850	0.558942445108484	0.568914
	next if (/^binID/);
	my @line = split /\t/;
	my $reads_step = (int ($line[7] / $step) + 1) * $step;
	my $w = exists $w_of_bins{$reads_step} ? $w_of_bins{$reads_step} : 0 ;
	my $corrected_reads_count = $w * $line[4];
	my $RRN = $line[4] / $rmmrn1;
	my $corrected_RRN = $corrected_reads_count / $cmmrn ;
	print GCcorrectedsts "$line[0]\t$line[1]\t$line[2]\t$corrected_RRN\n";
	}
close BINSTS or die "Cannot close $binsts_file, $!";
close GCcorrectedsts or die "Cannot close $out_dir/$GC_corrected_sts, $!";


sub statistics{
	my @array = sort {$a <=> $b} @_;
	my $min = $array[0];
	my $max = $array[-1];
	my $count = @array;
	my $sum = 0;
	foreach $_(@array){
		$sum += $_;
		}
		my $mean = $sum / $count;
	my $a = 0;
	foreach $_(@array){
		$a += ($_ - $mean) ** 2;
		}
	my $S2 = $a / $count;
	my $SD = $S2 ** 0.5;
	my $median = 0;
	my $Q1 = 0;
	my $Q3 = 0;
	if ($count % 2 == 0){
		my $pos1 = $count/2;
		my $pos2 = $count/2 + 1;
		$median = ($array[$pos1-1] + $array[$pos2-1])/2;
		if (($count / 2) % 2 == 0){
			my $pos1 = $count/2/2;
			my $pos2 = $count/2/2 + 1;
			my $pos3 = $count/2/2 * 3;
			my $pos4 = $count/2/2 * 3 + 1;
			$Q1 = ($array[$pos1-1] + $array[$pos2-1])/2;
			$Q3 = ($array[$pos3-1] + $array[$pos4-1])/2;
		}else{
			my $pos1 = ($count/2 + 1)/2;
			my $pos2 = $count/2 + ($count/2 + 1)/2;
			$Q1 = $array[$pos1-1];
			$Q3 = $array[$pos2-1];
			}
	}else{
		my $pos1 = ($count + 1)/2;
		$median = $array[$pos1-1];
		if ((($count - 1) / 2) % 2 == 0){
			my $pos1 = ($count - 1)/2/2;
			my $pos2 = ($count - 1)/2/2 + 1;
			my $pos3 = ($count - 1)/2/2 * 3 + 1;
			my $pos4 = ($count - 1)/2/2 * 3 + 1 + 1;
			$Q1 = ($array[$pos1-1] + $array[$pos2-1])/2;
			if ($count == 1){
				$Q3 = $median;
			}else{
				$Q3 = ($array[$pos3-1] + $array[$pos4-1])/2;
				}
		}else{
			my $pos1 = ($count + 1)/2/2;
			my $pos2 = ($count + 1)/2/2*3;
			$Q1 = $array[$pos1-1];
			$Q3 = $array[$pos2-1];
			}
		}
		my @sts = ();
		push @sts, $count, $sum, $mean, $S2, $SD, $min, $Q1, $median, $Q3, $max;
		return @sts;
	} 