package Info;

#===============================================================================
#
#         FILE: Info.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 08/06/2018 07:52:42 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
 

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


1;
