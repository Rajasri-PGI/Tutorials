#!/usr/bin/perl

use strict;
use warnings;

print "\n=========================================\n";
print "      PROTEOME SEQUENCE READER\n";
print "=========================================\n";

# Ask user for FASTA file
print "\nEnter the FASTA proteome file path: ";
chomp(my $file = <STDIN>);

# Check if file exists
unless (-e $file) {
    die "\nERROR: File '$file' does not exist.\n";
}

# Ask user for Protein ID or keyword
print "Enter Protein ID or Protein Name: ";
chomp(my $search = <STDIN>);

open(my $FH, "<", $file) or die "Cannot open file: $!";

my $header = "";
my $sequence = "";
my $found = 0;

while (my $line = <$FH>) {

    chomp($line);

    if ($line =~ /^>/) {

        # Print previous protein if matched
        if ($header =~ /\Q$search\E/i) {

            $found = 1;

            print "\nProtein Found\n";
            print "-----------------------------------------\n";
            print "$header\n";
            print "Sequence Length : " . length($sequence) . " amino acids\n\n";
            print "$sequence\n";

            last;
        }

        $header = $line;
        $sequence = "";
    }

    else {

        $sequence .= $line;

    }

}

# Check last protein in file
if (!$found && $header =~ /\Q$search\E/i) {

    print "\nProtein Found\n";
    print "-----------------------------------------\n";
    print "$header\n";
    print "Sequence Length : " . length($sequence) . " amino acids\n\n";
    print "$sequence\n";

    $found = 1;
}

close($FH);

if (!$found) {

    print "\nProtein not found in the proteome.\n";

}

print "\nProgram Finished.\n";
