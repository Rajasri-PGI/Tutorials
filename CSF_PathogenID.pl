#!/usr/bin/perl

###############################################################################
# CSF-PathogenID v1.0
#
# Standalone Perl Proteomics Search Engine
#
#
#
###############################################################################

use strict;
use warnings;

use Getopt::Long;
use POSIX qw(strftime);
use File::Path qw(make_path);
use MIME::Base64 qw(decode_base64);
use Compress::Zlib;
use XML::LibXML::Reader;

use constant {
    XML_READER_TYPE_ELEMENT     => 1,
    XML_READER_TYPE_TEXT        => 3,
    XML_READER_TYPE_END_ELEMENT => 15,
};
###############################################################################
# CONFIGURATION
###############################################################################

my %CONFIG = (

    precursor_tolerance_ppm => 20,

    fragment_tolerance_da   => 0.02,

    enzyme                  => 'Trypsin',

    missed_cleavages        => 2,

    min_peptide_length      => 7,

    max_peptide_length      => 45,

    output                  => "Results",

);

###############################################################################
# COMMAND LINE
###############################################################################

my $mzml  = "";
my $fasta = "";
my $output = "";

GetOptions(

    "mzml=s"   => \$mzml,
    "fasta=s"  => \$fasta,
    "output=s" => \$output

);

die "\nUsage:\nperl CSF_PathogenID.pl --mzml sample.mzML --fasta database.fasta --output Results\n\n"

unless($mzml && $fasta);

$CONFIG{output}=$output if($output ne "");

make_path($CONFIG{output}) unless(-d $CONFIG{output});

###############################################################################
# LOGGING
###############################################################################

open(my $LOG, ">", "$CONFIG{output}/Search.log") or die $!;

sub logmsg{

    my($msg)=@_;

    my $time=strftime("%Y-%m-%d %H:%M:%S",localtime);

    print $LOG "[$time] $msg\n";

    print "[$time] $msg\n";

}

###############################################################################
# GLOBAL VARIABLES
###############################################################################

my %AA;

my %PROTEINS;

my %PEPTIDES;

my %RESULTS;

my %STATISTICS;
my %PEPTIDE_INDEX;
my %PROTEIN_RESULTS;
my $accession;
###############################################################################
# MAIN
###############################################################################

main();

exit;

###############################################################################
# MAIN PROGRAM
###############################################################################

sub main{

    logmsg("==========================================");

    logmsg("CSF-PathogenID Started");

    logmsg("==========================================");

    initialise_amino_acids();

    read_fasta();

	digest_proteins();

    read_mzml();

	infer_proteins();

	classify_pathogens();

	generate_reports();

    logmsg("Finished");

}

###############################################################################
# PLACEHOLDER FUNCTIONS
###############################################################################


sub initialise_amino_acids{

%AA=(

A=>71.03711,

R=>156.10111,

N=>114.04293,

D=>115.02694,

C=>103.00919,

E=>129.04259,

Q=>128.05858,

G=>57.02146,

H=>137.05891,

I=>113.08406,

L=>113.08406,

K=>128.09496,

M=>131.04049,

F=>147.06841,

P=>97.05276,

S=>87.03203,

T=>101.04768,

W=>186.07931,

Y=>163.06333,

V=>99.06841

);

logmsg("Loaded ".scalar(keys %AA)." amino acid masses.");



}

###############################################################################
# READ FASTA DATABASE
###############################################################################

sub read_fasta {

    logmsg("Reading FASTA database...");

    open(my $FASTA, "<", $fasta)
        or die "Cannot open FASTA file: $fasta\n";

    my $accession  = "";
    my $description = "";
    my $species     = "Unknown";
    my $sequence    = "";

    my $protein_count = 0;

    while (my $line = <$FASTA>) {

        chomp($line);

        next if $line =~ /^\s*$/;

        ##############################################################
        # New FASTA header
        ##############################################################

        if ($line =~ /^>/) {

            # Save previous protein
            if ($accession ne "") {

                $sequence =~ s/\s+//g;
                $sequence = uc($sequence);

                $PROTEINS{$accession}{description} = $description;
                $PROTEINS{$accession}{species}     = $species;
                $PROTEINS{$accession}{sequence}    = $sequence;
                $PROTEINS{$accession}{length}      = length($sequence);
                $PROTEINS{$accession}{mw}          = 0;   # calculated later

                $protein_count++;
            }

            ##########################################################
            # Reset variables
            ##########################################################

            $sequence = "";
            $species  = "Unknown";

            ##########################################################
            # Parse UniProt header
            ##########################################################

            # Example:
            # >sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606

            if ($line =~ /^>\w+\|([^|]+)\|(.+)$/) {

                $accession  = $1;
                $description = $2;

            }
            else {

                $line =~ s/^>//;

                ($accession) = split(/\s+/, $line);

                $description = $line;

            }

            ##########################################################
            # Species extraction
            ##########################################################

            if ($description =~ /OS=([^=]+?)\sOX=/) {

                $species = $1;

                $species =~ s/\s+$//;

            }

        }

        ##############################################################
        # Sequence lines
        ##############################################################

        else {

            $line =~ s/\s+//g;

            $line = uc($line);

            # Remove invalid characters
            $line =~ s/[^A-Z]//g;

            $sequence .= $line;

        }

    }

    ##############################################################
    # Save last protein
    ##############################################################

    if ($accession ne "") {

        $PROTEINS{$accession}{description} = $description;
        $PROTEINS{$accession}{species}     = $species;
        $PROTEINS{$accession}{sequence}    = $sequence;
        $PROTEINS{$accession}{length}      = length($sequence);
        $PROTEINS{$accession}{mw}          = 0;

        $protein_count++;

    }

    close($FASTA);

    ##############################################################
    # Statistics
    ##############################################################

    $STATISTICS{proteins_loaded} = $protein_count;

    logmsg("Proteins loaded : $protein_count");

    ##############################################################
    # Display first few proteins
    ##############################################################

    my $shown = 0;

    foreach my $id (sort keys %PROTEINS) {

        logmsg(
            sprintf(
                "%-12s Length=%5d Species=%s",
                $id,
                $PROTEINS{$id}{length},
                $PROTEINS{$id}{species}
            )
        );

        last if (++$shown >= 5);

    }

    logmsg("FASTA parsing completed.");


}
print scalar(keys %PROTEINS), "\n";
###############################################################################
# STREAMING MZML PARSER
###############################################################################

sub read_mzml {

    logmsg("Reading mzML file...");

    my $reader = XML::LibXML::Reader->new(location => $mzml)
        or die "Cannot open mzML file: $mzml\n";

    my %spectrum;
    my %binary;

    my $spectrum_count = 0;

    while ($reader->read()) {

        my $type = $reader->nodeType();
        my $name = $reader->name();

        ###############################################################
        # START OF A SPECTRUM
        ###############################################################

        if ($type == XML_READER_TYPE_ELEMENT && $name eq "spectrum") {

            %spectrum = ();
            %binary   = ();

            $spectrum{id} =
                $reader->getAttribute("id") // "";

            $spectrum{index} =
                $reader->getAttribute("index") // "";

            next;
        }

        ###############################################################
        # CV PARAMETERS
        ###############################################################

        if ($type == XML_READER_TYPE_ELEMENT && $name eq "cvParam") {

            my $acc =
                $reader->getAttribute("accession") // "";

            my $value =
                $reader->getAttribute("value") // "";

            if ($acc eq "MS:1000511") {

                $spectrum{ms_level} = $value;

            }

            elsif ($acc eq "MS:1000744") {

                $spectrum{precursor_mz} = $value;

            }

            elsif ($acc eq "MS:1000514") {

                $binary{array_type} = "mz";

            }

            elsif ($acc eq "MS:1000515") {

                $binary{array_type} = "intensity";

            }

            elsif ($acc eq "MS:1000523") {

                $binary{precision} = 64;

            }

            elsif ($acc eq "MS:1000521") {

                $binary{precision} = 32;

            }

            elsif ($acc eq "MS:1000574") {

                $binary{compression} = "zlib";

            }

            elsif ($acc eq "MS:1000576") {

                $binary{compression} = "none";

            }

        }

        ###############################################################
        # BINARY BLOCK
        ###############################################################

        elsif ($type == XML_READER_TYPE_ELEMENT &&
               $name eq "binary") {

            unless ($reader->isEmptyElement()) {

                $reader->read();

                if ($reader->nodeType() ==
                    XML_READER_TYPE_TEXT) {

                    $binary{base64} = $reader->value();

                }

            }

        }

        ###############################################################
        # END OF ONE binaryDataArray
        ###############################################################

        elsif ($type == XML_READER_TYPE_END_ELEMENT &&
               $name eq "binaryDataArray") {

            if (defined $binary{array_type}) {

                if ($binary{array_type} eq "mz") {

                    $spectrum{mz_binary} = { %binary };

                }
                else {

                    $spectrum{intensity_binary} = { %binary };

                }

            }

            %binary = ();

        }

        ###############################################################
        # END OF ONE SPECTRUM
        ###############################################################

        elsif ($type == XML_READER_TYPE_END_ELEMENT &&
               $name eq "spectrum") {

            process_spectrum(\%spectrum);

            $spectrum_count++;

            if ($spectrum_count % 500 == 0) {

                logmsg("$spectrum_count spectra processed...");

            }

            %spectrum = ();

        }

    }

    $STATISTICS{spectra} = $spectrum_count;

    logmsg("Total spectra : $spectrum_count");

}

###############################################################################
# PROCESS ONE SPECTRUM
###############################################################################

sub process_spectrum {

    my ($spec) = @_;

    ###############################################################
    # Basic validation
    ###############################################################

    return unless defined $spec;
    return unless ref($spec) eq "HASH";

    return unless exists $spec->{mz_binary};
    return unless exists $spec->{intensity_binary};

    ###############################################################
    # Only search MS/MS spectra
    ###############################################################

    if (defined $spec->{ms_level}) {

        return unless $spec->{ms_level} == 2;

    }

    ###############################################################
    # Decode m/z array
    ###############################################################

    my $mz_array = decode_binary(
        $spec->{mz_binary}
    );

    ###############################################################
    # Decode intensity array
    ###############################################################

    my $intensity_array = decode_binary(
        $spec->{intensity_binary}
    );

    return unless defined $mz_array;
    return unless defined $intensity_array;

    ###############################################################
    # Convert arrays into peak list
    ###############################################################

    my $peaks = extract_peaks(
        $mz_array,
        $intensity_array
    );

    return unless @$peaks;

    ###############################################################
    # Store decoded peaks
    ###############################################################

    $spec->{peaks} = $peaks;

    ###############################################################
    # Spectrum statistics
    ###############################################################

    $spec->{peak_count} = scalar(@$peaks);

    ###############################################################
    # Perform spectrum search
    ###############################################################

    my $result = match_spectrum($spec);

    ###############################################################
    # Save identification
    ###############################################################

    if (defined $result) {

        $RESULTS{$spec->{id}} = $result;

        $STATISTICS{identified_spectra}++;

    }

    ###############################################################
    # Count processed spectra
    ###############################################################

    $STATISTICS{processed_spectra}++;

}

###############################################################################
# DECODE BINARY DATA
###############################################################################

sub decode_binary {

    my ($binary) = @_;

    ###############################################################
    # Validate input
    ###############################################################

    return undef unless defined $binary;
    return undef unless ref($binary) eq "HASH";
    return undef unless exists $binary->{base64};

    ###############################################################
    # Decode Base64
    ###############################################################

    my $buffer = decode_base64($binary->{base64});

    ###############################################################
    # Decompress if zlib-compressed
    ###############################################################

    if (defined $binary->{compression} &&
        $binary->{compression} eq "zlib") {

        my $inflated = Compress::Zlib::uncompress($buffer);

        unless (defined $inflated) {

            logmsg("ERROR: Failed to decompress binary array.");

            return undef;

        }

        $buffer = $inflated;
    }

    ###############################################################
    # Determine precision
    ###############################################################

    my $precision = 64;

    if (defined $binary->{precision}) {

        $precision = $binary->{precision};

    }

    ###############################################################
    # Unpack floating-point values
    ###############################################################

    my @values;

    if ($precision == 32) {

        my $count = length($buffer) / 4;

        @values = unpack("f<*", $buffer);

        if (scalar(@values) != $count) {

            logmsg("WARNING: 32-bit unpack count mismatch.");

        }

    }
    elsif ($precision == 64) {

        my $count = length($buffer) / 8;

        @values = unpack("d<*", $buffer);

        if (scalar(@values) != $count) {

            logmsg("WARNING: 64-bit unpack count mismatch.");

        }

    }
    else {

        logmsg("Unsupported precision: $precision");

        return undef;

    }

    ###############################################################
    # Remove invalid numbers
    ###############################################################

    @values = grep {

        defined($_)
        && $_ eq $_          # not NaN
        && abs($_) != 9**9**9

    } @values;

    ###############################################################
    # Return decoded array
    ###############################################################

    return \@values;

}


###############################################################################
# EXTRACT PEAKS
###############################################################################

sub extract_peaks {

    my ($mz_array, $intensity_array) = @_;

    ###############################################################
    # Validate input
    ###############################################################

    return [] unless defined $mz_array;
    return [] unless defined $intensity_array;

    return [] unless ref($mz_array) eq "ARRAY";
    return [] unless ref($intensity_array) eq "ARRAY";

    my $n = @$mz_array;

    return [] unless $n;
    return [] unless $n == @$intensity_array;

    ###############################################################
    # Determine base peak intensity
    ###############################################################

    my $base_peak = 0;

    foreach my $i (@$intensity_array) {

        next unless defined $i;

        $base_peak = $i if $i > $base_peak;

    }

    $base_peak = 1 if $base_peak == 0;

    ###############################################################
    # Build peak list
    ###############################################################

    my @peaks;

    for (my $i = 0; $i < $n; $i++) {

        my $mz  = $mz_array->[$i];
        my $int = $intensity_array->[$i];

        next unless defined $mz;
        next unless defined $int;

        next unless $mz > 0;
        next unless $int > 0;

        my $norm = ($int / $base_peak) * 100.0;

        push @peaks, {

            mz => $mz,

            intensity => $int,

            normalized => $norm

        };

    }

    ###############################################################
    # Sort by m/z
    ###############################################################

    @peaks = sort {

        $a->{mz} <=> $b->{mz}

    } @peaks;

    ###############################################################
    # Remove duplicate m/z values
    ###############################################################

    my @filtered;

    my $previous = -1;

    foreach my $p (@peaks) {

        if (abs($p->{mz} - $previous) < 1e-6) {

            next;

        }

        push @filtered, $p;

        $previous = $p->{mz};

    }

    ###############################################################
    # Statistics
    ###############################################################

    $STATISTICS{total_peaks} += scalar(@filtered);

    ###############################################################
    # Return peak list
    ###############################################################

    return \@filtered;

}

###############################################################################
# IN-SILICO TRYPSIN DIGESTION
###############################################################################

###############################################################################
# DIGEST PROTEINS INTO TRYPTIC PEPTIDES
###############################################################################

sub digest_proteins {

    logmsg("Digesting proteins...");

    ###############################################################
    # Reset peptide database
    ###############################################################

    %PEPTIDES      = ();
    %PEPTIDE_INDEX = ();

    $STATISTICS{proteins_digested} = 0;
    $STATISTICS{unique_peptides}   = 0;

    ###############################################################
    # Loop through every protein loaded from FASTA
    ###############################################################

    foreach my $accession (sort keys %PROTEINS) {

        ###########################################################
        # Skip invalid entries
        ###########################################################

        next unless exists $PROTEINS{$accession};

        my $sequence = $PROTEINS{$accession}{sequence};

        next unless defined $sequence;
        next if $sequence eq "";

        ###########################################################
        # Remove whitespace/newlines
        ###########################################################

        $sequence =~ s/\s+//g;

        ###########################################################
        # Debug output
        ###########################################################

        print "Digesting $accession (Length = ",
              length($sequence), " aa)\n";

        ###########################################################
        # Count proteins digested
        ###########################################################

        $STATISTICS{proteins_digested}++;

        ###########################################################
        # Perform trypsin digestion
        ###########################################################

        my @fragments = split /(?<=[KR])(?!P)/, $sequence;

        my $position = 1;

        ###########################################################
        # Process each peptide
        ###########################################################

        foreach my $fragment (@fragments) {

            next unless defined $fragment;

            $fragment =~ s/\s+//g;

            #######################################################
            # Ignore very short peptides
            #######################################################

            next if length($fragment) < $CONFIG{min_peptide_length};

            #######################################################
            # Ignore very long peptides
            #######################################################

            next if length($fragment) > $CONFIG{max_peptide_length};

            #######################################################
            # Calculate peptide mass
            #######################################################

            my $mass = calculate_peptide_mass($fragment);

            #######################################################
            # Store peptide
            #######################################################

            $PEPTIDES{$fragment}{mass} = $mass;

            push @{ $PEPTIDES{$fragment}{proteins} }, $accession;

            push @{ $PEPTIDES{$fragment}{locations} }, $position;

            #######################################################
            # Build precursor-mass index
            #######################################################

            my $bin = int($mass);

            push @{ $PEPTIDE_INDEX{$bin} }, $fragment;

            #######################################################
            # Update peptide position
            #######################################################

            $position += length($fragment);

        }

    }

    ###############################################################
    # Remove duplicate protein assignments
    ###############################################################

    foreach my $pep (keys %PEPTIDES) {

        my %seen;

        @{$PEPTIDES{$pep}{proteins}} =
            grep { !$seen{$_}++ }
            @{$PEPTIDES{$pep}{proteins}};

    }

    ###############################################################
    # Statistics
    ###############################################################

    $STATISTICS{unique_peptides} = scalar(keys %PEPTIDES);

    logmsg("Proteins digested : $STATISTICS{proteins_digested}");
    logmsg("Unique peptides   : $STATISTICS{unique_peptides}");
    logmsg("Mass bins         : " . scalar(keys %PEPTIDE_INDEX));

}
print "Proteins digested = $STATISTICS{proteins_digested}\n";

print "Peptides = ", scalar(keys %PEPTIDES), "\n";

print "Mass bins = ", scalar(keys %PEPTIDE_INDEX), "\n";
###############################################################################
# CALCULATE PEPTIDE MONOISOTOPIC MASS
###############################################################################

sub calculate_peptide_mass {

    my ($peptide) = @_;

    ###############################################################
    # Validate peptide
    ###############################################################

    return undef unless defined $peptide;

    $peptide = uc($peptide);

    return undef unless length($peptide);

    ###############################################################
    # Mass of water
    ###############################################################

    my $mass = 18.01056;

    ###############################################################
    # Add amino acid masses
    ###############################################################

    foreach my $aa (split //, $peptide) {

        unless (exists $AA{$aa}) {

            logmsg("Unknown amino acid '$aa' in peptide $peptide");

            return undef;

        }

        $mass += $AA{$aa};

    }

    ###############################################################
    # Round to six decimal places
    ###############################################################

    $mass = sprintf("%.6f", $mass);

    return $mass;


}

###############################################################################
# GENERATE THEORETICAL FRAGMENTS (b and y ions)
###############################################################################

sub generate_fragments {

    my ($peptide) = @_;

    ###############################################################
    # Validate input
    ###############################################################

    return undef unless defined $peptide;

    $peptide = uc($peptide);

    my $length = length($peptide);

    return undef if $length < 2;

    ###############################################################
    # Constants (Monoisotopic Masses)
    ###############################################################

    my $PROTON = 1.007276466812;
    my $H2O    = 18.010564684;
    my $H      = 1.00782503223;
    my $OH     = 17.00273965177;

    ###############################################################
    # Prefix and suffix cumulative masses
    ###############################################################

    my @prefix_mass;
    my @suffix_mass;

    my $mass = 0;

    # Prefix masses
    for (my $i = 0; $i < $length; $i++) {

        my $aa = substr($peptide, $i, 1);

        unless (exists $AA{$aa}) {

            logmsg("Unknown amino acid '$aa'");

            return undef;

        }

        $mass += $AA{$aa};

        $prefix_mass[$i] = $mass;

    }

    # Suffix masses
    $mass = 0;

    for (my $i = $length - 1; $i >= 0; $i--) {

        my $aa = substr($peptide, $i, 1);

        $mass += $AA{$aa};

        $suffix_mass[$i] = $mass;

    }

    ###############################################################
    # Generate ions
    ###############################################################

    my @fragments;

    for (my $i = 0; $i < $length - 1; $i++) {

        ##############################################
        # b-ion
        ##############################################

        my $b_mass =
            $prefix_mass[$i] + $PROTON;

        push @fragments, {

            ion      => "b".($i+1),

            number   => $i+1,

            type     => "b",

            charge   => 1,

            mz       => sprintf("%.6f",$b_mass)

        };

        ##############################################
        # y-ion
        ##############################################

        my $y_number = $length - $i - 1;

        my $y_mass =
            $suffix_mass[$i+1] +
            $H2O +
            $PROTON;

        push @fragments, {

            ion      => "y".$y_number,

            number   => $y_number,

            type     => "y",

            charge   => 1,

            mz       => sprintf("%.6f",$y_mass)

        };

    }

    ###############################################################
    # Sort by m/z
    ###############################################################

    @fragments = sort {

        $a->{mz} <=> $b->{mz}

    } @fragments;

    ###############################################################
    # Return fragment list
    ###############################################################

    return \@fragments;

}

###############################################################################
# MATCH EXPERIMENTAL SPECTRUM AGAINST PEPTIDE DATABASE
###############################################################################

sub match_spectrum {

    my ($spec) = @_;

    ###############################################################
    # Validate input
    ###############################################################

    return undef unless defined $spec;
    return undef unless ref($spec) eq "HASH";

    return undef unless exists $spec->{peaks};
    return undef unless exists $spec->{precursor_mz};

    my $precursor = $spec->{precursor_mz};

    ###############################################################
    # Approximate precursor mass
    # (Assumes charge = +1 for Version 1)
    ###############################################################

    my $PROTON = 1.007276466812;

    my $precursor_mass = $precursor - $PROTON;

    ###############################################################
    # Candidate peptide retrieval
    ###############################################################

    my $bin = int($precursor_mass);

    my @candidate_peptides;

    for my $b ($bin-1 .. $bin+1) {

        next unless exists $PEPTIDE_INDEX{$b};

        push @candidate_peptides,
            @{ $PEPTIDE_INDEX{$b} };

    }

    return undef unless @candidate_peptides;

    ###############################################################
    # Search candidates
    ###############################################################

    my $best_score = -1;

    my $best_match;

    foreach my $peptide (@candidate_peptides) {

        next unless exists $PEPTIDES{$peptide};

        my $pep_mass = $PEPTIDES{$peptide}{mass};

        ###########################################################
        # Precursor tolerance (ppm)
        ###########################################################

        my $ppm =
            abs($pep_mass - $precursor_mass)
            / $pep_mass
            * 1000000;

        next if $ppm > $CONFIG{precursor_tolerance_ppm};

        ###########################################################
        # Generate theoretical fragments
        ###########################################################

        my $fragments =
            generate_fragments($peptide);

        next unless defined $fragments;

        ###########################################################
        # Score spectrum
        ###########################################################

        my $score =
            score_psm(
                $spec->{peaks},
                $fragments
            );

        next unless defined $score;

        ###########################################################
        # Keep best hit
        ###########################################################

        if ($score > $best_score) {

            $best_score = $score;

            $best_match = {

                peptide => $peptide,

                protein =>
                    $PEPTIDES{$peptide}{proteins}[0],

                proteins =>
                    $PEPTIDES{$peptide}{proteins},

                precursor_mass =>
                    $pep_mass,

                observed_mass =>
                    $precursor_mass,

                ppm_error =>
                    sprintf("%.3f",$ppm),

                score =>
                    sprintf("%.2f",$score),

                fragment_count =>
                    scalar(@$fragments)

            };

        }

    }

    ###############################################################
    # Return best identification
    ###############################################################

    return $best_match;


}

###############################################################################
# PEPTIDE-SPECTRUM MATCH SCORING
###############################################################################

sub score_psm {

    my ($experimental, $theoretical) = @_;

    ###############################################################
    # Validate input
    ###############################################################

    return undef unless defined $experimental;
    return undef unless defined $theoretical;

    return undef unless ref($experimental) eq "ARRAY";
    return undef unless ref($theoretical) eq "ARRAY";

    my $tolerance = $CONFIG{fragment_tolerance_da};

    my $matched_fragments = 0;
    my $matched_intensity = 0;

    ###############################################################
    # Compare theoretical ions with experimental peaks
    ###############################################################

    foreach my $ion (@$theoretical) {

        my $target = $ion->{mz};

        foreach my $peak (@$experimental) {

            my $diff = abs($peak->{mz} - $target);

            if ($diff <= $tolerance) {

                $matched_fragments++;

                $matched_intensity +=
                    $peak->{normalized};

                last;

            }

        }

    }

    ###############################################################
    # No matches
    ###############################################################

    return 0 if $matched_fragments == 0;

    ###############################################################
    # Calculate score
    ###############################################################

    my $coverage =
        $matched_fragments / scalar(@$theoretical);

    my $intensity_score =
        $matched_intensity / 100.0;

    my $score =
          ($matched_fragments * 10)
        + ($coverage * 100)
        + $intensity_score;

    ###############################################################
    # Round score
    ###############################################################

    return sprintf("%.2f", $score);

}

###############################################################################
# PROTEIN INFERENCE
###############################################################################

sub infer_proteins {

    logmsg("Performing protein inference...");

    %PROTEIN_RESULTS = ();

    ###############################################################
    # Process all identified spectra
    ###############################################################

    foreach my $spectrum (keys %RESULTS) {

        my $psm = $RESULTS{$spectrum};

        next unless defined $psm;
        next unless exists $psm->{proteins};

        my $score = $psm->{score};

        foreach my $protein (@{$psm->{proteins}}) {

            #######################################################
            # Protein score
            #######################################################

            $PROTEIN_RESULTS{$protein}{score} += $score;

            #######################################################
            # Number of spectra
            #######################################################

            $PROTEIN_RESULTS{$protein}{spectra}++;

            #######################################################
            # Store peptide
            #######################################################

            $PROTEIN_RESULTS{$protein}{peptides}
                {$psm->{peptide}}++;

            #######################################################
            # Best peptide
            #######################################################

            if (
                !exists $PROTEIN_RESULTS{$protein}{best_score}
                ||
                $score >
                $PROTEIN_RESULTS{$protein}{best_score}
            ) {

                $PROTEIN_RESULTS{$protein}{best_score}
                    = $score;

                $PROTEIN_RESULTS{$protein}{best_peptide}
                    = $psm->{peptide};

            }

        }

    }

    ###############################################################
    # Calculate statistics
    ###############################################################

    my $protein_count = 0;

    foreach my $protein (keys %PROTEIN_RESULTS) {

        my @unique =
            keys %{ $PROTEIN_RESULTS{$protein}{peptides} };

        $PROTEIN_RESULTS{$protein}{unique_peptides}
            = scalar(@unique);

        ###########################################################
        # Protein annotation
        ###########################################################

        if (exists $PROTEINS{$protein}) {

            $PROTEIN_RESULTS{$protein}{description}
                = $PROTEINS{$protein}{description};

            $PROTEIN_RESULTS{$protein}{species}
                = $PROTEINS{$protein}{species};

            $PROTEIN_RESULTS{$protein}{length}
                = $PROTEINS{$protein}{length};

        }

        $protein_count++;

    }

    ###############################################################
    # Statistics
    ###############################################################

    $STATISTICS{identified_proteins}
        = $protein_count;

    logmsg("Proteins identified : $protein_count");

}

###############################################################################
# CLASSIFY PROTEINS BY SPECIES
###############################################################################

sub classify_pathogens {

    logmsg("Classifying identified proteins...");

    my %counts = (

        Human      => 0,
        Bacterial  => 0,
        Viral      => 0,
        Fungal     => 0,
        Parasite   => 0,
        Archaeal   => 0,
        Other      => 0,
        Unknown    => 0

    );

    foreach my $protein (keys %PROTEIN_RESULTS) {

        my $species =
            $PROTEIN_RESULTS{$protein}{species} // "Unknown";

        my $class = "Unknown";

        ###############################################################
        # Human
        ###############################################################

        if ($species =~ /Homo\s+sapiens/i) {

            $class = "Human";

        }

        ###############################################################
        # Viruses
        ###############################################################

        elsif (
            $species =~ /virus/i ||
            $species =~ /coronavirus/i ||
            $species =~ /influenza/i ||
            $species =~ /adenovirus/i ||
            $species =~ /herpes/i ||
            $species =~ /papilloma/i ||
            $species =~ /hepatitis/i ||
            $species =~ /HIV/i
        ) {

            $class = "Viral";

        }

        ###############################################################
        # Bacteria
        ###############################################################

        elsif (

            $species =~ /Escherichia/i ||
            $species =~ /Staphylococcus/i ||
            $species =~ /Streptococcus/i ||
            $species =~ /Pseudomonas/i ||
            $species =~ /Klebsiella/i ||
            $species =~ /Acinetobacter/i ||
            $species =~ /Enterococcus/i ||
            $species =~ /Salmonella/i ||
            $species =~ /Mycobacterium/i ||
            $species =~ /Bacillus/i ||
            $species =~ /Clostridium/i ||
            $species =~ /Helicobacter/i ||
            $species =~ /Vibrio/i ||
            $species =~ /Campylobacter/i
        ) {

            $class = "Bacterial";

        }

        ###############################################################
        # Fungi
        ###############################################################

        elsif (

            $species =~ /Candida/i ||
            $species =~ /Aspergillus/i ||
            $species =~ /Cryptococcus/i ||
            $species =~ /Saccharomyces/i ||
            $species =~ /Penicillium/i ||
            $species =~ /Fusarium/i
        ) {

            $class = "Fungal";

        }

        ###############################################################
        # Parasites
        ###############################################################

        elsif (

            $species =~ /Plasmodium/i ||
            $species =~ /Leishmania/i ||
            $species =~ /Trypanosoma/i ||
            $species =~ /Toxoplasma/i ||
            $species =~ /Giardia/i ||
            $species =~ /Entamoeba/i
        ) {

            $class = "Parasite";

        }

        ###############################################################
        # Archaea
        ###############################################################

        elsif (

            $species =~ /Methan/i ||
            $species =~ /Archae/i
        ) {

            $class = "Archaeal";

        }

        ###############################################################
        # Other
        ###############################################################

        elsif ($species ne "Unknown") {

            $class = "Other";

        }

        ###############################################################
        # Store classification
        ###############################################################

        $PROTEIN_RESULTS{$protein}{classification} = $class;

        $counts{$class}++;

    }

    ###############################################################
    # Save statistics
    ###############################################################

    $STATISTICS{human_proteins}      = $counts{Human};
    $STATISTICS{bacterial_proteins}  = $counts{Bacterial};
    $STATISTICS{viral_proteins}      = $counts{Viral};
    $STATISTICS{fungal_proteins}     = $counts{Fungal};
    $STATISTICS{parasite_proteins}   = $counts{Parasite};
    $STATISTICS{archaeal_proteins}   = $counts{Archaeal};
    $STATISTICS{other_proteins}      = $counts{Other};
    $STATISTICS{unknown_proteins}    = $counts{Unknown};

    ###############################################################
    # Logging
    ###############################################################

    logmsg("Human proteins      : $counts{Human}");
    logmsg("Bacterial proteins  : $counts{Bacterial}");
    logmsg("Viral proteins      : $counts{Viral}");
    logmsg("Fungal proteins     : $counts{Fungal}");
    logmsg("Parasite proteins   : $counts{Parasite}");
    logmsg("Archaeal proteins   : $counts{Archaeal}");
    logmsg("Other proteins      : $counts{Other}");
    logmsg("Unknown proteins    : $counts{Unknown}");
}

###############################################################################
# GENERATE REPORTS
###############################################################################

sub generate_reports {

    logmsg("Generating reports...");

    ###############################################################
    # Protein Report
    ###############################################################

    open(my $PROT, ">", "$CONFIG{output}/Protein_Report.csv")
        or die "Cannot create Protein_Report.csv\n";

    print $PROT join(",",
        "Protein",
        "Description",
        "Species",
        "Classification",
        "Score",
        "Spectra",
        "Unique_Peptides",
        "Best_Peptide"
    ), "\n";

    foreach my $protein (
        sort {
            ($PROTEIN_RESULTS{$b}{score} // 0)
                <=>
            ($PROTEIN_RESULTS{$a}{score} // 0)
        } keys %PROTEIN_RESULTS
    ) {

        print $PROT join(",",

            $protein,

            '"' .
            ($PROTEIN_RESULTS{$protein}{description} // "") .
            '"',

            '"' .
            ($PROTEIN_RESULTS{$protein}{species} // "") .
            '"',

            $PROTEIN_RESULTS{$protein}{classification} // "",

            sprintf("%.2f",
                $PROTEIN_RESULTS{$protein}{score} // 0),

            $PROTEIN_RESULTS{$protein}{spectra} // 0,

            $PROTEIN_RESULTS{$protein}{unique_peptides} // 0,

            $PROTEIN_RESULTS{$protein}{best_peptide} // ""

        ), "\n";

    }

    close($PROT);

    ###############################################################
    # PSM Report
    ###############################################################

    open(my $PSM, ">", "$CONFIG{output}/PSM_Report.csv")
        or die "Cannot create PSM_Report.csv\n";

    print $PSM join(",",
        "Spectrum",
        "Peptide",
        "Protein",
        "Observed_Mass",
        "Theoretical_Mass",
        "PPM_Error",
        "Score"
    ), "\n";

    foreach my $spec (sort keys %RESULTS) {

        my $r = $RESULTS{$spec};

        next unless defined $r;

        print $PSM join(",",

            $spec,

            $r->{peptide} // "",

            ref($r->{proteins}) eq "ARRAY"
                ? join(";", @{$r->{proteins}})
                : "",

            $r->{observed_mass} // "",

            $r->{precursor_mass} // "",

            $r->{ppm_error} // "",

            $r->{score} // ""

        ), "\n";

    }

    close($PSM);

    ###############################################################
    # Summary Report
    ###############################################################

    open(my $SUM, ">", "$CONFIG{output}/Summary.txt")
        or die "Cannot create Summary.txt\n";

    print $SUM "=========================================\n";
    print $SUM "CSF-PathogenID Summary Report\n";
    print $SUM "=========================================\n\n";

    print $SUM "Proteins Loaded        : ",
        ($STATISTICS{proteins_loaded} // 0), "\n";

    print $SUM "Proteins Digested      : ",
        ($STATISTICS{proteins_digested} // 0), "\n";

    print $SUM "Unique Peptides        : ",
        ($STATISTICS{unique_peptides} // 0), "\n";

    print $SUM "Spectra Processed      : ",
        ($STATISTICS{processed_spectra} // 0), "\n";

    print $SUM "Identified Spectra     : ",
        ($STATISTICS{identified_spectra} // 0), "\n";

    print $SUM "Identified Proteins    : ",
        ($STATISTICS{identified_proteins} // 0), "\n";

    print $SUM "\n";

    print $SUM "Human Proteins         : ",
        ($STATISTICS{human_proteins} // 0), "\n";

    print $SUM "Bacterial Proteins     : ",
        ($STATISTICS{bacterial_proteins} // 0), "\n";

    print $SUM "Viral Proteins         : ",
        ($STATISTICS{viral_proteins} // 0), "\n";

    print $SUM "Fungal Proteins        : ",
        ($STATISTICS{fungal_proteins} // 0), "\n";

    print $SUM "Parasite Proteins      : ",
        ($STATISTICS{parasite_proteins} // 0), "\n";

    print $SUM "Archaeal Proteins      : ",
        ($STATISTICS{archaeal_proteins} // 0), "\n";

    print $SUM "Other Proteins         : ",
        ($STATISTICS{other_proteins} // 0), "\n";

    print $SUM "Unknown Proteins       : ",
        ($STATISTICS{unknown_proteins} // 0), "\n";

    close($SUM);

    ###############################################################
    # Log completion
    ###############################################################

    logmsg("Protein report  : $CONFIG{output}/Protein_Report.csv");
    logmsg("PSM report      : $CONFIG{output}/PSM_Report.csv");
    logmsg("Summary report  : $CONFIG{output}/Summary.txt");

    logmsg("All reports successfully generated.");


}

