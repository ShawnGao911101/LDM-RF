#!/usr/bin/env perl

use strict;
use warnings;

use IO::Handle;
use Getopt::Long;
use Carp;
use File::Find;
use List::Util qw(sum max min);

my $options = check_params();

my @bino_files;
retrieve_bino_files( \@bino_files );
process_bino_files( \@bino_files );

#----------------------------------------------------------------
sub process_bino_files {
	my ($bino_files_ref) = @_;

	my $machine_file = $bino_files[0];
	$machine_file =~ s/.*\///xms;
	$machine_file =~ s/bino/machine/xms;
	
#	my $Pvalue_file    = $machine_file . '_Pvalue';
#	my $power_one_file = $machine_file . '_power_one';
	my $sum_file       = $machine_file . '_sum';

#	open my $Pvalue_fh, '>', $Pvalue_file
#	  or croak("Can't open output file [$Pvalue_file]: $!.\n");

#	open my $power_one_fh, '>', $power_one_file
#	  or croak("Can't open output file [$power_one_file]: $!.\n");

	open my $sum_fh, '>', $sum_file
	  or croak("Can't open output file [$sum_file]: $!.\n");
	  
#	print_header($Pvalue_fh);
#	print_header($power_one_fh);
	print_header($sum_fh);

	foreach my $bino_file (@$bino_files_ref) {
#		process_type( $bino_file, $machine_file, $Pvalue_fh,    'Pvalue' );
#		process_type( $bino_file, $machine_file, $power_one_fh, 'power_one' );
		process_type( $bino_file, $machine_file, $sum_fh,       'sum' );
	}

#	close $Pvalue_fh    or croak("Failed to close [$Pvalue_file]: $!\n");
#	close $power_one_fh or croak("Failed to close [$power_one_file]: $!\n");
	close $sum_fh       or croak("Failed to close [$sum_file]: $!\n");
}

sub process_type {
	my ( $bino_file, $machine_file, $fh, $type ) = @_;

	my %chrs        = ();
	my %chrs_count  = ();
	my @chrs_values = ();
	my $best_snp_id;
	my %best_chr_values  = ();
	my %chrs_total_lines = ();

	open my $BINO_fh, '<', $bino_file
	  or croak("Can't open input file [$bino_file]: $!.\n");

	my $count      = 0;
	my $type_count = 0;

	my $type_column =
	    ( $type eq 'Pvalue' )    ? 9
	  : ( $type eq 'power_one' ) ? 10
	  :                            4;

  BINO_LINE:
	while ( my $line = <$BINO_fh> ) {
		$count++;

		if ( $count == 1 ) {
			next BINO_LINE;
		}

		chomp $line;
		my @values = split( /\t/, $line );
		my ( $calc_type, $snp_id, $type_value ) = @values[ 0, 1, $type_column ];

		next BINO_LINE if !( $calc_type eq $type );

		$type_count++;
		my ( $chr, $pos ) = split( '_', $snp_id );
		$chr =~ s/_ .+ \z//xms;
		$best_snp_id = $snp_id if !$best_snp_id;

		if ( !$chrs_count{$chr} && ( @chrs_values < 5 ) ) {
			$chrs_count{$chr}++;
			push @chrs_values, $type_value;
		}

		if ( $best_snp_id =~ m/$chr/xms ) {
			$best_chr_values{$chr}{'count'}++;
			push @{ $best_chr_values{$chr}{'pos'} }, $pos;
		}

		my $chr_number = $chr;
		$chr_number =~ s/S//xms;
		$chrs_total_lines{$chr_number}++;
	}

	my $total_lines          = $type_count;
	my ($key)                = keys %best_chr_values;
	my $total_lines_best_chr = $best_chr_values{$key}{'count'};
	my @total_lines          = ();

	foreach my $key ( sort keys %chrs_total_lines ) {
		push @total_lines, $chrs_total_lines{$key};
	}
		
	my $chrs_total_lines = join( "\t", @total_lines );


        my $best_val         = $chrs_values[0];
        my $sec_val          = $chrs_values[1];


        my $odds_best        = $chrs_values[0] / (2 - $chrs_values[0] + 0.0000000001 );
        my $odds_sec         = $chrs_values[1] / (2 - $chrs_values[1] + 0.0000000001 );
        my $middle_vlalue    = ( $chrs_values[3] + $chrs_values[4] ) / 2;
        my $odds_middle      = $middle_vlalue / ( 2 - $middle_vlalue + 0.0000000001 );
        my $best_second      = $odds_best / $odds_sec;
        my $best_middle      = $odds_best / $odds_middle;
        my $best_second_org      = $chrs_values[0] / $chrs_values[1];
        my $best_middle_org      =
          $chrs_values[0] / ( ( $chrs_values[3] + $chrs_values[4] ) / 2 );


        my $best_chr         = $best_snp_id;
        $best_chr =~ s/_ .+ \z//xms;
        my $max_dist = max( @{ $best_chr_values{$best_chr}{'pos'} } );
        my $min_dist = min( @{ $best_chr_values{$best_chr}{'pos'} } );
        my $distance = $max_dist - $min_dist;

        my $bino_value = $bino_file;
	$bino_value =~ s/.*\///xms;
	my $best_value = $chrs_values[0];

	print $fh
"$type\t$bino_value\t$total_lines\t$best_snp_id\t$best_value\t$chrs_total_lines\t$best_second\t$best_middle\t$distance\t$best_val\t$sec_val\t$middle_vlalue\t$best_second_org\t$best_middle_org\n";

	close $BINO_fh or croak("Failed to close [$bino_file]: $!\n");
}

sub retrieve_bino_files {
	my ($data_files_ref) = @_;

	find(
		{
			wanted => sub {
				push @$data_files_ref, $File::Find::name
				  if ( (m/_bino\z/) );
			},
		},
		$options->{'dir'}
	);
}

sub print_header {
	my ($fh) = @_;

	print $fh
"type\tuni_id\ttot_lines\tbest_snp\tbest_value\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tbest_second\tbest_middle\tdist\tbest_val\tsec_val\tmid_val\tbest_second_org\tbest_middle_org\n";
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "dir:s" );

	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0")   if $options{'man'};
	exec("pod2usage $0")
	  if !( $options{'dir'} );

	return \%options;
}

__DATA__

=head1 NAME

   creat_ML_features.pl

=head1 COPYRIGHT

   Copyright (c) 2019 Jiri Stiller and Shang Gao, All rights reserved.
   
=head1 DESCRIPTION
   
   This script generates features used for training machine learning models.

=head1 SYNOPSIS

   creat_ML_features.pl -dir <dir> [-help] [-man]
   
      -help		Displays basic usage information
      -man		Displays more detailed information 
      -dir      Directory containing SUM files
   
   Example:   
      creat_ML_features.pl -dir /input_dir
      
   Results:
      2_machine_Pvalue
      2_machine_power_one
      2_machine_sum
      ..
      ..
      ..
   
=cut
