#!/usr/bin/env perl

use strict;
use warnings;

use IO::Handle;
use Getopt::Long;
use Carp;
use File::Find;
use Storable;

my $options = check_params();

print "Generating data\n";
my @data_files;
my %data;
my %unique_records;
retrieve_data_files( \@data_files );
merge_data( \@data_files, \%data );
make_unique_data( \%data, \%unique_records );
print "Done\n";

#----------------------------------------------------------------

sub make_unique_data {
	my ( $data_ref, $unique_data_ref ) = @_;

	my $unique_records_id = 0;

  SEQ:
	foreach my $seq ( keys %{$data_ref} ) {
		my (@records) = sort keys %{ $data_ref->{$seq} };

		if ( $options->{'limit'} ) {
			my $total_records = @records;
			next SEQ if ( $total_records < $options->{'limit'} );
		}

		my $records = join( '_', @records );

		if ( !$unique_data_ref->{$records}{'id'} ) {
			$unique_records_id++;
			$unique_data_ref->{$records}{'id'} = $unique_records_id;
		}

		$unique_data_ref->{$records}{'seqs'}{$seq}++;
	}

	my $data_keys_total = keys %{$data_ref};
	print "\t... data keys total [$data_keys_total]\n";

	my $unique_data_keys_total = keys %{$unique_data_ref};
	print "\t... unique data keys total [$unique_data_keys_total]\n";

	print "\n\tSaving unique hash... \n";
	store( $unique_data_ref, 'all_hash' );
	print "\t\t... Done \n";

	if ( $options->{'split'} ) {
		print "\n\tSaving hashes split from unique hash... \n";

		my $split_hashes_ref =
		  split_hash( $options->{'split'}, $unique_data_ref );

		my $hash_count = 0;

		foreach my $hash_ref (@$split_hashes_ref) {
			$hash_count++;

			if ( $hash_count == 1 ) {
				my $split_keys_total =
				  keys %{$hash_ref};
				print "\t... each split hash keys total [$split_keys_total]\n";
			}

			my $hash_file = $hash_count . '_hash';
			store( $hash_ref, $hash_file );
		}

		print "\t\t... Done \n";
	}
}

sub split_hash {
	my ( $x, $hash ) = @_;

	my @keys = keys %$hash;
	my @hashes;

	while ( my @subset = splice( @keys, 0, $x ) ) {
		push @hashes, { map { $_ => $hash->{$_} } @subset };
	}

	return \@hashes;
}

sub retrieve_ids {
	my ($ids_ref) = @_;

	open my $IDS_fh, '<', $options->{'ids'}
	  or croak("Can't open input file [$options->{'ids'}]: $!.\n");

	while ( my $line = <$IDS_fh> ) {
		chomp $line;
		my ( $id_dump, $id_snp ) = split( /\t/, $line );
		$ids_ref->{$id_dump} = $id_snp;
	}

	close $IDS_fh
	  or croak("Failed to close [$options->{'ids'}]: $!\n");
}

sub merge_data {
	my ( $data_files_ref, $data_ref ) = @_;

	my $count = 0;
	my %ids   = ();

	retrieve_ids( \%ids );

	foreach my $file_name (@$data_files_ref) {
		$count++;

		if ( !( $count % 100 ) ) {
			print "   *** Processed [$count] dump files\n";
		}

		open my $DATA_fh, '<', $file_name
		  or croak("Can't open input file [$file_name]: $!.\n");

		$file_name =~ s/.*\///xms;
		$file_name =~ s/_dump\z//xms;
		$file_name = $ids{$file_name};

	  DATA_LINE:
		while ( my $line = <$DATA_fh> ) {

			my @values = split( /\t/, $line );
			my ( $seq, $val ) = @values[ 0, 1 ];
			$data_ref->{$seq}{$file_name} = $val;
		}

		close $DATA_fh or croak("Failed to close [$file_name]: $!\n");
	}
}

sub retrieve_data_files {
	my ($data_files_ref) = @_;

	find(
		{
			wanted => sub {
				push @$data_files_ref, $File::Find::name
				  if ( (m/_dump\z/) );
			},
		},
		$options->{'dir'}
	);
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "dir:s", "ids:s", "split:i",
		"limit:i", "hash+" );

	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0")   if $options{'man'};
	exec("pod2usage $0")
	  if !( $options{'dir'} && $options{'ids'} );

	return \%options;
}

__DATA__

=head1 NAME

   create_unique_records.pl

=head1 COPYRIGHT

   Copyright (c) 2019 Jiri Stiller and Shang Gao, All rights reserved.
   
=head1 DESCRIPTION
   
    This script creates an integrated unique record (UR) hash and split UR hashes for parallel computing.

=head1 SYNOPSIS

   create_unique_records.pl -dir <dir> -ids <file> [-limit <number>] [-split <number>] [-help] [-man]
   
      -help		Displays basic usage information
      -man		Displays more detailed information 
      -dir      Input dir containg result files of k-mer analysis with se
      -split    Split hash into subhashes according to keys number specified
      -limit    Minimum number of genotypes within unique record key
      -ids      Replace genotype ID to match SNP matrix if needed
				Change_ID.sed:
				ERR2331144      ERX2382092
				ERR2331145      ERX2382093
				ERR2331146      ERX2382094
					..              ..
					..              ..

   
   Example:   
      create_unique_records.pl -dir /input_dir -ids /Change_ID.sed -split 1000 -limit 10
      
   Results:
      all_hash
      1_hash
      ..
      
=cut
