#!/usr/bin/env perl 
# =============================================================================
#
#         FILE: seeker.pl
#
#        USAGE: ./seeker.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 08/06/2018 03:55:51 PM
#     REVISION: ---
# =============================================================================

use strict;
use warnings;
use utf8;

#use lib 'lib';
#use Info;

use Data::Dumper;
use Storable;
use Storable 'dclone';
use Parallel::ForkManager;
use File::Basename;



# ==============================================================================

# Verify arguments number.
die "Input arguents missing" unless @ARGV >= 2;

my $outputs;

# Testing on inputs and outputs.
die "Input file doesn't exist" unless ( -e -f -r "$ARGV[0]");
die "You don't have permission to write on directory " . dirname($ARGV[1]) . "." if ( ! -d -r dirname($ARGV[1]));
if ("$ARGV[1]") {$outputs = $ARGV[1];} else {$outputs = "out_neopipe.tsv";}


# Input an output variables.
my $bam = $ARGV[0];
my $dbin = "/home/scratch60/lbuzzo/RTC/tmp/.ref.perldb";
my $dbout = $outputs;

# Set number of jobs.
my $jobs = $ARGV[2] ? $ARGV[2] : 1;

# Set the maximun number of process.
my $pm = Parallel::ForkManager->new($jobs);


# Some runtime variables.
my $read_size = 101;
my $fragment_average_size = 300;
my $fragment_std = 50;
my $range = $fragment_average_size + 3*$fragment_std - $read_size;
my $Kmer_size = 1000;

my $pos_before_gene_start = 0;
my $pos_before_gene_end = 0;
my $pos_beyond_gene_start = 0;
my $pos_beyond_gene_end = 0;

# Main table.
my %table;
my %genome;
my $criteria;


# Run in parent right after finishing child process.
$pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $table_ref) = @_;
		while (my ($key, $value) = each %$table_ref) {
			$table{$key} = $value;
		}
		
		# Debug code.
		#print Dumper $table_ref;
	}
);


# The criteria to search on each chromossome.
sub runner {
	
	# Inner variables.
	my ($file, $table, @options) = @_;
	
	# Debug code.
	#print "$file $criteria $table\n";
	#return;
	
	# Options must e a hash.
	my $chr = $options[0];
	
	
	# Open reading pipe on $bam with 'samtools view ...'.
	open my $ph, "-|" => "samtools view $file $chr 2> /dev/null"
		or die "Can't open pipe for file $file: $!.";
	
	while (<$ph>) {
		chomp;
		my @read = split;
		
		my $Kmer = int($read[3]/1000);
		
		# Debug code.
		
		
		# Skip abnormal chromosome names.
		next if (length($read[2]) > 5 || length($read[6]) > 5);
		
		# Debug code.
		#print "@read\n";
		
		# Search for abonormal read. 
		if ( ($read[6] eq "=" && abs($read[8]) > 10000) || ($read[6] ne "=" && $read[6] ne "*") ) {
			
			# Map the second read.
			my $chr_rtc;
			
			if ($read[6] eq "=") {
			$chr_rtc = $read[2];
			}
			else {
			$chr_rtc = $read[6];
			}
			
			
			# Insert read positions in table.
			if ($read[5] =~ qr/^${read_size}M$/) {
				push @{ $table->{$chr}{$chr_rtc}{$Kmer}{'regular'} } => \@read;
				print "$chr\t$Kmer\t$read[3]\texpected\tregular\t$read[5]\t$chr_rtc\t$read[7]\t$read[0]\t$read[9]\n";
			}
			else {
				push @{ $table->{$chr}{$chr_rtc}{$Kmer}{'border'} } => \@read;
				print "$chr\t$Kmer\t$read[3]\texpected\tborder\t$read[5]\t$chr_rtc\t$read[7]\t$read[0]\t$read[9]\n";
			}
		}
		
		#Emergency stop
		#last;
	} # while.
	
	# Close reading pipe.
	close $ph;
	
	return;
}



# ==============================================================================

# Open a pipe on $bam to 'samtools view -H ...' and extract 'chr' pattern.
open my $ph, "-|" => "samtools view -H $bam 2> /dev/null | tail -n3 | head -n1 | cut -f2"
	or die "Can't open pipe for file $bam: $!.";

my $pattern = <$ph>;
if ($pattern) {
	$pattern = substr $pattern, 3, 3;
}
else {
	die "Couldn't read $bam on pipe.";
}

close $ph;


# Header line to print.
my $header = "chr1\tkmer\tposition1\texpected\ttype\tmatch\tchr2\tposition2\tid\tseq";
print "$header\n";


my $counter = 0;

# Iterate over all chromosomes.
while ( $counter < 24 ) {
	
	# Adjust iterations through chromosomes.
	$jobs = 24 - $counter if ($counter + $jobs >= 24);
	
	
	# Dispatch parallel jobs.
	DATA_LOOP:
	foreach my $jid (1..$jobs) {
		# Forks and returns the pid for the child:
		my $pid = $pm->start and next DATA_LOOP;
		
		
		# The main hash to store data.
		my %table;
		
		# Some pseudo-global vars.
		my $sfx;
		# Take correct suffix for 'chr' pattern.
		if ( $counter + $jid > 22 ) {
			$sfx = "X";
			$sfx = "Y" if ($counter + $jid == 24);
		}
		else {
			$sfx = $counter + $jid;
		}
		
		# Create chromossome pattern.
		my $chr = $pattern . $sfx;
		
		
		# The main search criteria.
		runner($bam, \%table, $chr);
		
		# Coverage calculator.
		#coverage();
		
		
		# Terminate the child process.
		$pm->finish(0, \%table); # Terminates the child process
	} # foreach of the children.
	
	# Wait for child processes.
	$pm->wait_all_children;
	
	$counter += $jobs;
}





# ==============================================================================

# Save table on output.
store \%table, $dbout
	or die "Can't store $dbout in file: $!.";


# Print with Dumper.
#print Dumper \%table;


# Simple print.
print STDERR "Table hash: " . (keys %table) . ".\n";
