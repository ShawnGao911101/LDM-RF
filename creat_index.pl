#!/usr/bin/env perl

use strict;
use warnings;

use IO::Handle;
use Getopt::Long;
use Carp;
use File::Find;
use File::Basename;
use Storable;
use List::Util qw(sum max);
use Math::CDF qw(pbinom);

my $options = check_params();

print "Reading hash...\n";
my $unique_data_ref   = retrieve( $options->{'hash'} );
my $unique_keys_total = keys %{$unique_data_ref};
print "unique data keys total [$unique_keys_total]\n";

generate_ids($unique_data_ref);

#----------------------------------------------------------------
sub generate_ids {
	my ($unique_data_ref) = @_;

	my $seqs_bino_ids_file      = 'seqs_bino_ids.csv';
	my $unique_records_ids_file = 'unique_records_ids.csv';

	open my $SEQS_BINO_fh, '>', $seqs_bino_ids_file
	  or croak("Can't open output file [$seqs_bino_ids_file]: $!.\n");

	open my $REC_IDS_fh, '>', $unique_records_ids_file
	  or croak("Can't open output file [$unique_records_ids_file]: $!.\n");

	foreach my $unique_record ( keys %$unique_data_ref ) {
		my $id = $unique_data_ref->{$unique_record}{'id'};
		print $REC_IDS_fh "$id,$unique_record\n";

		foreach my $seq ( keys %{ $unique_data_ref->{$unique_record}{'seqs'} } )
		{
			print $SEQS_BINO_fh "$seq,$id" . '_bino' . "\n";
		}
	}

	close $SEQS_BINO_fh
	  or croak("Failed to close [$seqs_bino_ids_file]: $!\n");

	close $REC_IDS_fh
	  or croak("Failed to close [$unique_records_ids_file]: $!\n");
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "hash:s" );

	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0")   if $options{'man'};
	exec("pod2usage $0")
	  if !( $options{'hash'} );

	return \%options;
}

__DATA__

=head1 NAME

   creat_index.pl

=head1 COPYRIGHT

   Copyright (c) 2019 Jiri Stiller and Shang Gao, CSIRO Primary Industry, Australia.
   All rights reserved.
   Intended for CSIRO internal use only.
   
=head1 DESCRIPTION
   
   This scirpt generates the one-to-many links between unique record and GBS tags from a hash.

=head1 SYNOPSIS

   creat_index.pl -hash <file> [-help] [-man]
   
      -help     Displays basic usage information
      -man      Displays more detailed information 
      -hash     Input hash
   
   Example:   
      creat_index.pl -hash all_hash
      
   Results:
      seqs_bino_ids.csv
      unique_records_ids.csv
   
=cut
