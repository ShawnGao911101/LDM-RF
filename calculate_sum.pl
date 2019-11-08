#!/usr/bin/env perl

use strict;
use warnings;

use IO::Handle;
use Getopt::Long;
use Carp;
use File::Find;
use Storable;
use List::Util qw(sum max);
use Math::GSL::Randist qw(gsl_ran_binomial_pdf);
#use PDL::Stats::Distr;
#use Math::CDF qw(pbinom);


my $options = check_params();

print "Reading hash...\n";
my $unique_data_ref   = retrieve( $options->{'hash'} );
my $unique_keys_total = keys %{$unique_data_ref};
print "unique data keys total [$unique_keys_total]\n";
my %bino_files;
process_unique_records( $unique_data_ref, \%bino_files );

#----------------------------------------------------------------
sub process_unique_records {
	my ( $unique_data_ref, $bino_files_ref ) = @_;

	my @records_order = ();
	my %SNP           = ();
	retrieve_records( \@records_order );
	retrieve_snp_values( \%SNP );

  UNIQUE_RECORD_LINE:
	foreach my $unique_record ( keys %$unique_data_ref ) {
		my $unique_record_id = $unique_data_ref->{$unique_record}{'id'};
		my $bino_file        = $unique_record_id . '_bino';

		open my $BINO_fh, '>', $bino_file
		  or croak("Can't open output file [$bino_file]: $!.\n");

		print_header($BINO_fh);

		my @PHE       = ();
		my %bino_data = ();
		generate_phe_values( \$unique_record, \@records_order, \@PHE );
		my $phe = join( "", @PHE );
		merge_snp_phe_values( \%bino_data, \%SNP, \$phe );
		sort_snp_column( $BINO_fh, \%bino_data, 'Pvalue' );
		sort_snp_column( $BINO_fh, \%bino_data, 'power_one' );
		sort_snp_column( $BINO_fh, \%bino_data, 'sum' );

		close $BINO_fh
		  or croak("Failed to close [$bino_file]: $!\n");
	}
}

sub retrieve_records {
	my ($records_order_ref) = @_;

	open my $DATA_fh, '<', $options->{'rec'}
	  or croak("Can't open input file [$options->{'rec'}]: $!.\n");

	my @records_order;

	while ( my $line = <$DATA_fh> ) {
		chomp $line;
		@$records_order_ref = split( /\t/, $line );
	}

	close $DATA_fh or croak("Failed to close [$options->{'rec'}]: $!\n");
}

sub retrieve_snp_values {
	my ($SNP_ref) = @_;

	open my $SNP_fh, '<', $options->{'snp'}
	  or croak("Can't open input file [$options->{'snp'}]: $!.\n");

	my $snp_count = 0;

  SNP:
	while ( my $snp_line = <$SNP_fh> ) {
		$snp_count++;
		next SNP if ( $snp_count == 1 );

		chomp $snp_line;
		my ( $snp_id, $snp ) = ( $snp_line =~ /\A (\S+) \t (.+) \z/xms );
		$snp =~ tr/\t//d;
		$SNP_ref->{$snp_id} = $snp;
	}

	close $SNP_fh
	  or croak("Failed to close [$options->{'snp'}]: $!\n");
}

sub print_header {
	my ($fh) = @_;

	print $fh "sort\t snp_id\t",
	  join( "\t",
		qw(one_calc zero_calc sum_one_zero_calc one_match zero_match all_match p_one_match bino_Pvalue power_one)
	  ),
	  "\n";
}

sub generate_phe_values {
	my ( $unique_record_ref, $records_order_ref, $PHE_ref ) = @_;

	my @records = split( '_', $$unique_record_ref );
	my %values;
	foreach my $record (@records) {
		$values{$record}++;
	}

	foreach my $ordered_record (@$records_order_ref) {

		if ( $values{$ordered_record} ) {

			# leave SNP letters unchanged when or'ed below
			push @$PHE_ref, "@";
		}
		else {
			# change SNP letters to lowercase
			push @$PHE_ref, " ";
		}
	}
}

sub merge_snp_phe_values {
	my ( $bino_data_ref, $SNP_ref, $phe_ref ) = @_;

  SNP_LINE:
	while ( my ( $snp_id, $snp ) = each %$SNP_ref ) {
		my %data = ();

		# combine snp and phe using bitwise or
		$snp |= $$phe_ref;

		# tr/x// returns the number of occurrences of x.
		$data{A1} = $snp =~ tr/A//;
		$data{T1} = $snp =~ tr/T//;
		$data{G1} = $snp =~ tr/G//;
		$data{C1} = $snp =~ tr/C//;
		$data{N1} = $snp =~ tr/N//;
		$data{A0} = $snp =~ tr/a//;
		$data{T0} = $snp =~ tr/t//;
		$data{G0} = $snp =~ tr/g//;
		$data{C0} = $snp =~ tr/c//;
		$data{N0} = $snp =~ tr/n//;

		my (
			$flag,              $one_calc,    $zero_calc,
			$sum_one_zero_calc, $one_match,   $zero_match,
			$all_match,         $p_one_match, $binomial_Pvalue,
			$power_one
		) = calculate_PK_values( \%data );

		next SNP_LINE if ( $flag == 1 );

		my $line = join( "\t",
			$snp_id,            $one_calc,    $zero_calc,
			$sum_one_zero_calc, $one_match,   $zero_match,
			$all_match,         $p_one_match, $binomial_Pvalue,
			$power_one )
		  . "\n";
		$bino_data_ref->{$snp_id}{'line'}      = $line;
		$bino_data_ref->{$snp_id}{'Pvalue'}    = "$binomial_Pvalue";
		$bino_data_ref->{$snp_id}{'power_one'} = $power_one;
		$bino_data_ref->{$snp_id}{'sum'}       = $sum_one_zero_calc;
	}
}

sub calculate_PK_values {
	my ($data_ref) = @_;

	my @one_values = ();
	push @one_values, $data_ref->{'A1'};
	push @one_values, $data_ref->{'T1'};
	push @one_values, $data_ref->{'G1'};
	push @one_values, $data_ref->{'C1'};
	push @one_values, $data_ref->{'N1'};

	my $one_max = max(@one_values);
	my $one_sum = sum(@one_values);

	return (1) if ( $one_sum == 0 );

	my @zero_values = ();
	push @zero_values, $data_ref->{'A0'};
	push @zero_values, $data_ref->{'T0'};
	push @zero_values, $data_ref->{'G0'};
	push @zero_values, $data_ref->{'C0'};
	push @zero_values, $data_ref->{'N0'};

	my $zero_sum = sum(@zero_values);

	return (1) if ( $zero_sum == 0 );

	my $zero_correspondig_to_one_max_key            = '';
	my @zero_values_corresponding_to_one_max_values = ();

  KEY:
	foreach my $key ( keys %$data_ref ) {

		if (   ( $data_ref->{$key} == $one_max )
			&& ( $key =~ m/1\z/xms )
			&& ( $key =~ m/[ATGCN]/xms ) )
		{
			$zero_correspondig_to_one_max_key = $key;
			$zero_correspondig_to_one_max_key =~ s/1/0/xms;
			push @zero_values_corresponding_to_one_max_values,
			  $data_ref->{$zero_correspondig_to_one_max_key};
		}
	}

	my $zero_max_values_count = @zero_values_corresponding_to_one_max_values;
	my $zero_value_correspondig_to_one_max_value =
	  max(@zero_values_corresponding_to_one_max_values);
	my $one_calc = $one_max / $one_sum;
	my $zero_calc =
	  ( $zero_sum - $zero_value_correspondig_to_one_max_value ) / $zero_sum;
	my $one_match         = $one_max;
	my $zero_match = $zero_sum - $zero_value_correspondig_to_one_max_value;
	my $all_match  = $one_match + $zero_match;
	my $sum_one_zero_calc = ($one_match + $zero_match) / ($one_sum + $zero_sum);
	my $p_one_match =
	  ( $one_max + $zero_value_correspondig_to_one_max_value ) /
	  $options->{'psize'};

#	my $binomial_Pvalue = pbinom( $one_match, $all_match, $p_one_match );
#	my $binomial_Pvalue = pmf_binomial($one_match, $all_match, $p_one_match);
	my $binomial_Pvalue = gsl_ran_binomial_pdf($one_match,$p_one_match,$all_match);
	my $power_one = $p_one_match**$one_match;

	return (
		0,                  $one_calc,    $zero_calc,
		$sum_one_zero_calc, $one_match,   $zero_match,
		$all_match,         $p_one_match, $binomial_Pvalue,
		$power_one
	);
}

sub sort_snp_column {
	my ( $fh, $bino_data_ref, $column ) = @_;

	my %chrs       = ();
	my @snp_sorted = ();

	if ( $column eq 'sum' ) {
		@snp_sorted = sort {
			$bino_data_ref->{$b}{$column} <=> $bino_data_ref->{$a}{$column}
		} keys %{$bino_data_ref};
	}
	else {
		@snp_sorted = sort {
			$bino_data_ref->{$a}{$column} <=> $bino_data_ref->{$b}{$column}
		} keys %{$bino_data_ref};
	}

  SNP:
	foreach my $snp (@snp_sorted) {
		my $chr = $snp;
		$chr =~ s/_ .+ \z//xms;
		$chrs{$chr}++;
		my $total_chrs = keys %chrs;
		print $fh $column . "\t" . $bino_data_ref->{$snp}{'line'};
		last SNP if ( $total_chrs == 8 );
	}
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "rec:s", "snp:s", "hash:s",
		"psize:i" );

	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0")   if $options{'man'};
	exec("pod2usage $0")
	  if !($options{'rec'}
		&& $options{'snp'}
		&& $options{'hash'}
		&& $options{'psize'} );

	return \%options;
}

__DATA__

=head1 NAME

   calculate_sum.pl

=head1 COPYRIGHT

   Copyright (c) 2019 Jiri Stiller and Shang Gao, All rights reserved.
   
=head1 DESCRIPTION
   
   This scirpt computes and descending sorts the statistic ‘SUM’ for each UR in split hash against whole SNP matrix.

=head1 SYNOPSIS

   calculate_sum.pl -rec <file> -hash <file> -snp <file> -psize <number> [-help] [-man]
   
      -help		Displays basic usage information
      -man		Displays more detailed information 
	  -rec		A list of all genotype IDs (must be same with the IDs in SNP matrix)
	  -hash		The uniq_record hash table created by "create_unique_records.pl"
	  -snp		The SNP matrix. (Hapmap format without annotation head lines and columns)
				eample:
				SNP_ID    Genotype1  Genotype2  Genotype3
				S1_42269      T          T          T
				S1_42284      G          G          G
				S1_51494      G          G          G 

	  -psize	The size of population used 
   
   Example:   
      calculate_sum.pl -rec records.csv -hash 1_hash -snp snp_matrix -psize 1140
      
=cut
