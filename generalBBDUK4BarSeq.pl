############################################################################################################################################################
# generalBBDUK4BarSeq.pl - A wrapper script over the BBDuk tool for NGS Bar-seq analysis.
# Copyright (C) 2015-present Antonio Claudio Bello Ribeiro and contributors
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#	Contact details: antonio.ribeiro@abdn.ac.uk or at Centre for Genome-Enabled Biology and Medicine (CGEBM) 23 St Machar Drive - Aberdeen AB24 3RY - United Kingdom
############################################################################################################################################################
use strict;
use warnings;
#use DateTime;
#use Cwd;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

# Quit unless we have the correct number of command-line args (reference: http://alvinalexander.com/perl/perl-command-line-arguments-read-args)
my $num_args = $#ARGV + 1;

if ($num_args != 4)
{
	print "\nUsage: generalBBDUK4BarSeq.pl <absolute path to the gzipped file of reads (typically the R1 dataset) to be processed by the code> <mismatch rate allowed when looking for the barcode sequence (0 to 3)> <absolute path to the list of barcodes in FASTA format> <output tabulated file>\n";
	# Reports will be merged later with UNIX xargs
	exit;
}

# File(s) input
print"\nStarting the processing...\n\n";
my $file = $ARGV[0]; # Path to the sample's reads file
my $mismatchRate = $ARGV[1]; # The aimed mismatch rate being used in the experiment (for usage in the bbduk third-party script later)
my $list_file = $ARGV[2]; # The FASTA-like file with the list of barcode sequences being searched in the experiment (for usage in the bbduk third-party script later)
my $output_filename = $ARGV[3]; # The filename of the aimed output file

# Printing the \@ARGV array for debugging purposes if needed (reference: http://www.perlmonks.org/?node_id=714601)
print "Info: ARGV is: @ARGV\n";

unless (-e $list_file) { # Testing whether the file exists or not
	die "$0: File not found: $list_file";
}

if(($mismatchRate < 0) || ($mismatchRate > 3)){
	print "Sorry, I'm not prepared to deal with this mismatch rate input value.... Please, restart the program and provide any number between 0 and 3.\n";
	die;
}

print "\nParsing the barcode list file and populating the data structures...\n";
open (IN, $list_file) or die "$0: ERROR: $list_file: $!"; # Testing whether the file can be opened

# References: https://www.biostars.org/p/90936/ and http://seqanswers.com/forums/showthread.php?t=24777
my @barcodes = ();

my %hash_FASTA_Seqs = ();
my $hash_FASTA_Seqs = '';
local $/ = '>';
my $junk = <IN>;
my %hash_barcode_Lengths = ();
my $hash_barcode_Lengths = '';

while(my $record = <IN>) {
	chomp $record;
	my ($defLine, @seqLines) = split /\n/, $record;
	print "\n\$defLine: $defLine\n";
	my $sequence = join('',@seqLines);
	print "Seq Length: ", length($sequence), "\n";
	$hash_barcode_Lengths{$defLine} = length($sequence);
	$hash_FASTA_Seqs{$defLine} = $sequence;
	push @barcodes, $defLine;
}
close IN;

# Testing the contents of the hash (reference: http://www.tutorialspoint.com/perl/perl_hashes.htm)
my @keys_hash_FASTA_Seqs = keys %hash_FASTA_Seqs;
my $size_hash_FASTA_Seqs = @keys_hash_FASTA_Seqs;
print "\nHash size for hash_FASTA_Seqs: $size_hash_FASTA_Seqs\n";
my @values_hash_FASTA_Seqs = values %hash_FASTA_Seqs;
my $size_hash_FASTA_Seqs_values = @values_hash_FASTA_Seqs;
print "\nHash values size for hash_FASTA_Seqs: $size_hash_FASTA_Seqs_values\n";

print "\n[DEBUG]################################################################################################\n\n";
print "@barcodes\n";
print "\n[DEBUG]################################################################################################\n";

#############################################
# String manipulation helpers and extensions:
#############################################
my $underscore = "_";
my $extBBDuk = ".BBDukED";
my $extBBDukED_NOTmatched = ".NOTmatched";
my $extBBDukED_Matched = ".Matched";
my $extFQ = ".fq";
my $gString;

###########################
# Global arrays and hashes:
###########################
my %gHoHoTallyMatched; # Hash to try to store, per FASTQ file/sample, the $tallyMatched value for each barcode and which will be used later
my %gHoPropsANYBarcode; # Hash to try to store, per FASTQ file/sample, the "first key metric" (i.e. the proportions of reads with internal barcode (any barcode) relative to total reads (all reads irrespective of whether they have internal barcode or not) in each sample
my %gHoHoBarcodesProportions; # Hash to try to store the proportion of each barcode per FASTQ file/sample
my @gAoHeaders = ("Barcode related information (below) | Samples (right)");
my @gAoRealBasenames = ();
my @gAo1stRow = ("Proportion of reads with ANY barcodes relative to total reads");

#########################
# Other global variables:
#########################
my $real_basename;
my $barcode;
my $sampleHeader;
my $without_extension;
my $testString1;
my $testString2;
#########################

# A basic approach to try to prevent inadvertent usage of the other read file if present in the directory
chomp(my $compressedFileBasename = `basename $file`);
print "\n";
print "\$compressedFileBasename: $compressedFileBasename\n";
if ($compressedFileBasename =~ /_L001_R1_001.fastq.gz/) {
		print "It matches R1\n\n";
		$testString1="_L001_R1_001.fastq.gz";
		$testString2="_L001_R1_001.fastq";
	} elsif($compressedFileBasename =~ /_L001_R2_001.fastq.gz/) {
		print "It matches R2\n\n";
		$testString1="_L001_R2_001.fastq.gz";
		$testString2="_L001_R2_001.fastq";
	} else {
		print "Can't recognise this read file... Aborting...";
		die;
	}

#if ($file =~ /\.fastq.gz$/i) { # Previous simpler approach considering that only R1 reads would be present
if ($file =~ /$testString1$/i) {
	# This might be improved in a later version to prevent wasting time with a gunzip operation especially if dealing with a massive number of files... For now, since typically BarSeq-related samples are not so big in size and for my own debugging/development time purposes, just including a double-check step to allow the comparison of the obtained number of reads prior and after the gunzipping stage
	chomp(my $priorRawLineCount = `zcat $file | wc -l`);
	print "\$file $file (raw lines):\t$priorRawLineCount\n";
	chomp(my $priorSeqsLineCount = $priorRawLineCount/4);
	print "\$file $file (seqs lines):\t$priorSeqsLineCount\n\n";
	print "I'm trying to gunzip the file $file ...\n\n";
	($without_extension = $file) =~ s/\.[^.]+$//; # Trying to get the "basename" of the file (reference: https://stackoverflow.com/questions/3667859/remove-file-extension-and-path-from-a-string-in-perl )
	print "The uncompressed version of the file $file will be named as: $without_extension\n\n";
	system("gunzip $file");
	opendir(DIR, ".");
	my @internalFiles_2 = readdir(DIR);
	close DIR;
	foreach my $internalFile (@internalFiles_2) {
		#if ($internalFile =~ /\.fastq$/i) { # Previous simpler approach considering that only R1 reads would be present
		if ($internalFile =~ /$testString2$/i) {
			print "I have, now, the following FASTQ file(s) under the directory: $internalFile ... Trying to get the respective tallying information...\n\n";
			chomp(my $rawLineCount = `cat $internalFile | wc -l`);
			print "\$internalFile $internalFile (raw lines):\t$rawLineCount\n";
			chomp(my $seqsLineCount = $rawLineCount/4);
			print "\$internalFile $internalFile (seqs lines):\t$seqsLineCount\n";
			($real_basename = $internalFile) =~ s/\.[^.]+$//;
			print "\nThe real basename of the file is: $real_basename\n\n";
			# This might be improved in a later version to prevent wasting time with a gunzip operation especially if dealing with a massive number of files... For now, since typically BarSeq-related samples are not so big in size and for my own debugging/development time purposes, just including a double-check step to allow the comparison of the obtained number of reads prior and after the gunzipping stage
			print "Trying to compare, now, the read numbers before and after the gunzip operation...\n\n";
			if ($priorSeqsLineCount == $seqsLineCount) {
				print "Read numbers obtained prior and after gunzipping are EQUAL for dataset $real_basename ... We're good to proceed!\n";
			} else {
				die "Read numbers obtained prior and after gunzipping are NOT EQUAL for dataset $real_basename ... Shutting down...\n";
			}
			# Trying to populate the array of headers
			# Test of modification to make the header of the sample follow the "basename" of the output report file aiming for more succinctness
			($sampleHeader = `basename $output_filename _report.txt`) =~ s/\R//g;
			#push @gAoHeaders, $real_basename;
			push @gAoHeaders, $sampleHeader;
			print "\n[DEBUG]#### \$sampleHeader: $sampleHeader\n";
			push @gAoRealBasenames, $real_basename;
			# Zeroing the $SUM_TALLY_MATCHED variable for the first iteration within the FASTQ file/sample
			my $SUM_TALLY_MATCHED = 0;
			print "\n[DEBUG]#### My \$SUM_TALLY_MATCHED value before the first barcode of file $internalFile is $SUM_TALLY_MATCHED\n\n";
			foreach $barcode (@barcodes) {
				`sh bbduk.sh -Xmx4g "in=$internalFile" "out=$real_basename$underscore$barcode$extBBDuk$extBBDukED_NOTmatched$extFQ" "outm=$real_basename$underscore$barcode$extBBDuk$extBBDukED_Matched$extFQ" "k=$hash_barcode_Lengths{$barcode}" "literal=$hash_FASTA_Seqs{$barcode}" "maskmiddle=f" "hammingdistance=$mismatchRate" `;
				print "\n";
				my $outm = $real_basename.$underscore.$barcode.$extBBDuk.$extBBDukED_Matched.$extFQ;
				my $out = $real_basename.$underscore.$barcode.$extBBDuk.$extBBDukED_NOTmatched.$extFQ;
				# Trying to get and report the raw number of lines for the "matched" file for the given barcode
				chomp(my $rawLines4TallyingMatched = `cat $outm | wc -l`);
				print "\t$outm (raw lines):\t$rawLines4TallyingMatched\n";
				# Trying to get and report the "payload" number of lines for the "matched" file for the given barcode
				chomp(my $tallyMatched = $rawLines4TallyingMatched/4);
				print "\t$outm (seqs lines):\t$tallyMatched\n";
				# Trying to store the $tallyMatched information in a hash so it can be used later in the calculation of the "second metric"
				$gHoHoTallyMatched{$real_basename}{$barcode} = $tallyMatched;
				# Trying to accumulate the total number of "matched" lines in each "matched" file created by the iterations of the "bbduk" script
				$SUM_TALLY_MATCHED = add($SUM_TALLY_MATCHED,$tallyMatched);
				# Trying to report the current value of the SUM_TALLY_MATCHED accumulator for eventual debugging purposes
				print "\n[DEBUG]################################################################################################\n";
				print "\tCurrent SUM_TALLY_MATCHED accumulator value after barcode $barcode is:\t$SUM_TALLY_MATCHED\n";
				print "[DEBUG]################################################################################################\n\n";
				# Trying to calculate and report a basic percentage of the presence of the barcode within the total lines of the file (just for debugging purposes)
				chomp(my $pctTallyMatched = $tallyMatched/$seqsLineCount);
				print "\t$outm (%):\t$pctTallyMatched\n";
				print "\t$outm:\t$tallyMatched\t$pctTallyMatched\n";
				print "\n";
				# Trying to get and report the information for the "NOT matched" figures (just for debugging purposes)
				chomp(my $rawLines4TallyingNOTMatched = `cat $out | wc -l`);
				print "\t$out (raw lines):\t$rawLines4TallyingNOTMatched\n";
				chomp(my $tallyNOTMatched = $rawLines4TallyingNOTMatched/4);
				print "\t$out (seqs lines):\t$tallyNOTMatched\n";
				chomp(my $pctTallyNOTMatched = $tallyNOTMatched/$seqsLineCount);
				print "\t$out (%):\t$pctTallyNOTMatched\n";
				print "\t$out:\t$tallyNOTMatched\t$pctTallyNOTMatched\n";
				print "\n";
				print "[DEBUG]############### JUST CHECKED A GIVEN BARCODE AT THE FIRST MAIN ITERATION OF THE CODE ###############\n";
				print "\n";
				print "[DEBUG] Now, trying to get rid of the intermediate files generated by BBDUK in order to save storage space (on the fly)...\n";
				print "[DEBUG] Removing $outm file...\n";
				`rm $outm`;
				print "[DEBUG] Removing $out file...\n";
				`rm $out`;
				print "[DEBUG] ###################################################################################################\n";
			}
			print "\n";
			print "[DEBUG]################ GETTING TO THE FIRST 'MAIN METRIC' COMPUTATION STAGE ################\n";
			print "\n";
			# "First key metric" to look for: 1) proportion of reads with internal barcode (any barcode) relative to total reads (all reads irrespective of whether they have internal barcode or not) in each sample... I can get it now because I've iterated of each barcode per FASTQ file and has now the final calculation for the given $SUM_TALLY_MATCHED value
			print "[DEBUG]###################################################################\n";
			print "\tFinal \$SUM_TALLY_MATCHED for file $internalFile is $SUM_TALLY_MATCHED ...\n";
			print "[DEBUG]###################################################################\n\n";
			chomp(my $proportionOfReadsWithInternalBarcode_ANY_relative2TotalReads = $SUM_TALLY_MATCHED/$seqsLineCount);
			print "$internalFile \(proportion of reads with internal barcode relative to total reads present in the file\):\t$proportionOfReadsWithInternalBarcode_ANY_relative2TotalReads\n";
			$gHoPropsANYBarcode{$real_basename} = $proportionOfReadsWithInternalBarcode_ANY_relative2TotalReads; # Trying to populate the hash with the proportion of ANY barcoded reads per sample
			push @gAo1stRow, $proportionOfReadsWithInternalBarcode_ANY_relative2TotalReads; # Trying to populate the "printing" array for the first row with the information gotten for the proportion of ANY barcoded reads per sample
			# "Second key metric" to look for: 2) within a sample, the proportion of each different internal barcode relative to the total number of reads with an internal barcode (reads without a barcode are excluded from this total)
			# Taking care of the cases where, eventually, $SUM_TALLY_MATCHED value == 0... Reference: http://forums.devshed.com/perl-programming-6/perl-illegal-division-zero-376242.html... $z is $x divided by $y, unless $y is zero, in which case $z should be zero too
			print "\n";
			print "[DEBUG]################ GETTING TO THE SECOND 'MAIN METRIC' COMPUTATION STAGE ################\n";
			print "\n";
			foreach $barcode (@barcodes) {
				my $tallyMatched42ndTime = $gHoHoTallyMatched{$real_basename}{$barcode};
				print "[DEBUG]\$real_basename:$real_basename\t\$barcode:$barcode\t\$tallyMatched42ndTime:$tallyMatched42ndTime\n";
				my $pctTallyMatched4The2ndTime = 0;
				if ($SUM_TALLY_MATCHED == 0) {
					$pctTallyMatched4The2ndTime = 0;
				} else {
					chomp($pctTallyMatched4The2ndTime = $tallyMatched42ndTime/$SUM_TALLY_MATCHED);
				}
				print "$real_basename$underscore$barcode$extBBDuk$extBBDukED_Matched$extFQ (%):\t$pctTallyMatched4The2ndTime\n";
				# Trying to store the $pctTallyMatched4The2ndTime information in a hash so it can be used later in the printing report
				$gHoHoBarcodesProportions{$real_basename}{$barcode} = $pctTallyMatched4The2ndTime;
				print "$real_basename$underscore$barcode$extBBDuk$extBBDukED_Matched$extFQ:\t$tallyMatched42ndTime\t$pctTallyMatched4The2ndTime\n";
				print "\n[DEBUG]##########################################################################################################\n";
			}
		}
	}
}
#system("gzip *fastq");
system("gzip $without_extension");
print "\n[DEBUG]########################################################################################################################\n";
print "[DEBUG] I've finished what I had to do in this directory... You should now have the barcode countings done for the original compressed FASTQ file $file properly reported... Hopping back to the upper level folder...\n";
print "[DEBUG]##########################################################################################################################\n\n";

# Trying out the final report PRINTING... Reference: http://www.perlmonks.org/?node_id=1070897
print "\n#############################\n";
print "Printing the output report...\n";
print "#############################\n\n";
open(RAWFILE, '>', $output_filename) or die "Cannot open $output_filename: $!";
print RAWFILE join("\t", @gAoHeaders), "\n"; # Reference: http://www.perlmonks.org/?node_id=668523
print RAWFILE join("\t", @gAo1stRow), "\n";

# Trying to prepare the printing rows for each barcode
foreach $barcode (@barcodes) {
	print RAWFILE "$barcode";
	foreach my $realBasenameElement (@gAoRealBasenames) {
		$gString = "\t$gHoHoBarcodesProportions{$realBasenameElement}{$barcode}";
		print RAWFILE $gString;
	}
	print RAWFILE "\n";
}

###############
# Sub-routines:
###############

sub add { # Reference: https://stackoverflow.com/questions/20971475/adding-two-number-in-perl
	(my $x, my $y) = @_;
	my $res = $x + $y ;
	return $res ;
}

print "\nJOB DONE!!!";
print " \n";
