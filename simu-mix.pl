#!/usr/bin/env perl

use strict;
use warnings;

my $script_dir = $0;
$script_dir =~ s#[^/]*$##;

my $mix_program = "${script_dir}mixer";

my ($dat_file, $snp_file, $out_prefix);
my @pop_names;
my @pop_files;
my %hash_pop_names;


&parse_args();
&check_pops_match_dat(); # ensure populations in $dat_file match @pop_names

# generate a breakpoints filename
my $bp_file = "$out_prefix.bp";

# call simulator
print "Calling simulator:\n";
my $ret = system("$mix_program $dat_file $snp_file $bp_file");

if ($ret != 0) {
  # call failed: clean up
  &clean_up();
  die "Failed to run simulator program $mix_program\n";
}

print "\nDone simulating; producing results using breakpoints file\n";

my @person_bp;
my @pop_indices;
&parse_break_points(\@person_bp, \@pop_indices);
&print_haplotypes(\@person_bp, \@pop_indices);

&clean_up();

print "done.\n";

################################################################################
# subroutines

sub parse_args() {
  if (@ARGV < 7) {
    print "Too few arguments.\n\n";
    &print_usage();
  }

  $dat_file = $ARGV[0];
  $snp_file = $ARGV[1];
  $out_prefix = $ARGV[2];
  for(my $i = 3; $i < @ARGV; $i+=2) {
    if ($ARGV[$i] !~ /^-/) {
      print "Bad argument: $ARGV[$i]\n\n";
      &print_usage(); # should have a hyphen before the population name
    }
    my $pop_name = $ARGV[$i];
    $pop_name =~ s/^-//;  # remove the hyphen
    push @pop_names, $pop_name;
    push @pop_files, $ARGV[$i+1];
    $hash_pop_names{ $pop_name } = $#pop_names;
  }
}

sub check_pops_match_dat() {
  open DAT, "$dat_file" or die "Couldn't open $dat_file: $!\n";
  my $header = <DAT>;
  close DAT;

  $header =~ s/^\s+//; # remove any leading whitespace
  
  # get fields from header:
  my @fields = split /\s+/, $header;

  # remove first two fields which are number of samples and admixed population
  # label
  shift @fields;
  shift @fields;

  my $num_pops = @pop_names;
  my $num_pop_fields = @fields;
  if ($num_pops != $num_pop_fields) {
    print "Error: user supplied arguments for $num_pops populations, but dat ";
    print "file contains\n";
    print "$num_pop_fields populations\n";
    exit;
  }

  my %pop_names_seen;
  while (@fields) {
    my $cur = shift @fields;
    if (not defined $hash_pop_names{$cur}) {
      print "Error: dat file specifies population $cur, which wasn't supplied ";
      print "as an argument\n";
      exit;
    }
    elsif (defined $pop_names_seen{$cur}) {
      print "Error: dat file contains label for population $cur more than once";
      print "\n";
      exit;
    }
    $pop_names_seen{$cur} = 1; # saw this population
  }
}


sub print_usage() {
  print "Usage:\n";
  print "  $0 [in.dat] [in.snp] [out.phgeno] [POPULATION FILES]\n\n";
  print "POPULATION FILES is a list of arguments, one for each population ";
  print "specified in\n";
  print "in.dat.  The format is -POP_LABEL [in.phgeno].\n\n";
  print "Example using files in 'examples' subdirectory:\n";
  print "  $0 AA.dat YRI.snp AA.phgeno -CEU ceu.phgeno -YRI yri.phgeno\n\n";
  print "The SNP and phased genotype files are expected to be in ";
  print "Eigenstrat format.\n\n";
  exit;
}

# Useful if we need a version of this that doesn't read all the genotype data
# into memory:
sub parse_break_points {
  my $person_bp_ref = shift;
  my $pop_indices_ref = shift;

  open BP, "$bp_file" or
      die "Couldn't open simulator output file $bp_file: $!\n";

  my $header = <BP>;
  chomp $header;
  # indicates which population label corresponds to each index number
  @{$pop_indices_ref} = split /\s+/, $header;

  # read each line from the file, which describes the breakpoints for each
  # individual
  my $total_bps = 0;
  my $num_indivs = 0;
  while (my $line = <BP>) {
    my @indiv_bps = split /\s+/, $line;
    push @{$person_bp_ref}, \@indiv_bps;
    
    # for averaging the number of breakpoints:
    $total_bps += @indiv_bps;
    $num_indivs++;
  }

  my $ave = $total_bps / $num_indivs;
  print "Average number of break points per individual: $ave\n";

  close BP;
}


sub print_haplotypes {
  my $person_bp_ref = shift;
  my $pop_indices_ref = shift;

  # read in the haplotypes corresponding to each of the populations:
  my @pop_haps_lists;
  my $nMarkers = -1;
  foreach my $pop (@{$pop_indices_ref}) {
    my $pop_num = $hash_pop_names{ $pop };
    my $pop_file = $pop_files[ $pop_num ];

    my @haps;
    my $ret = open POPIN, "$pop_file";
    if (!$ret) {
      &clean_up();
      die "Couldn't open $pop_file for $pop: $!\n";
    }
    while (my $hap = <POPIN>) {
      chomp $hap;
      push @haps, $hap;
    }
    close POPIN;
    push @pop_haps_lists, \@haps;
    if ($nMarkers == -1) {
      $nMarkers = @haps;
    }
    elsif ($nMarkers != @haps) {
      &clean_up();
      die "Error: different number of markers between population files!\n";
    }
  }

  # open output file for admixed samples
  my $ret = open OUT, ">$out_prefix.phgeno";
  if (!$ret) {
    &clean_up();
    die "Couldn't open output file $out_prefix.phgeno: $!\n";
  }

  # indicates whether a sample is being copied from:
  my @copying_sample;
  for(my $p = 0; $p < @{$pop_indices_ref}; $p++) {
    my $nSamples = length $pop_haps_lists[$p]->[0];
    for(my $i = 0; $i < $nSamples; $i++) {
      $copying_sample[$p][$i] = 0;
    }
  }

  # indicates which sample the admixed individual is copying from:
  my @copying_ind;
  for(my $i = 0; $i < @{$person_bp_ref}; $i++) {
    $copying_ind[$i] = -1; # initially
  }

  for (my $marker = 0; $marker < $nMarkers; $marker++) {

    # process breakpoints for each sample, identifying a suitable source
    # sample to copy from (i.e., one that isn't being copied from already)
    for (my $ind = 0; $ind < @{$person_bp_ref}; $ind++) { # for each sample...
      # All the unprocessed breakpoints for $indiv
      my $indiv_bps_ref = $person_bp_ref->[$ind];

      # get the next breakpoint downstream of $marker:
      my ($pop_num, $end_pos) = split /:/, $indiv_bps_ref->[0];

      my $dont_copy_id = -1;
      my $new_copy_id_needed = 0;
      if ($end_pos < $marker) {
	my $old_pop_num = $pop_num;
	my $old_copy_id = $copying_ind[ $ind ];

	# Recombination break point!
	do {
	  shift @{$indiv_bps_ref}; # remove the stale bp
	  if (@{$indiv_bps_ref} == 0) {
	    die "Empty bp list!! Have marker = $marker\n";
	  }
	  ($pop_num, $end_pos) = split /:/, $indiv_bps_ref->[0]; # get next bp
	  # Note: can sometimes get multiple breakpoints at the same marker
	  # this is arguably a bug in the mixing program, but it's easy enough
	  # to just skip any such breakpoint: they define an empty region
	  # between the previous breakpoint and the next breakpoint that occurs
	  # at the same marker.
	} while ($end_pos < $marker);

	# ensure that we don't "recombine" to the sample source sample:
	if ($old_pop_num == $pop_num) {
	  $dont_copy_id = $copying_ind[ $ind ];
	}
	# no longer copying from the old sample
	$copying_sample[ $old_pop_num ][ $old_copy_id ] = 0;

	$new_copy_id_needed = 1;
      }

      if ($new_copy_id_needed || $copying_ind[ $ind ] == -1) {

	my %searched_ids;

	# find new sample to copy from:
	my $new_copy_id;
	do {
	  $new_copy_id = int( rand( scalar @{$copying_sample[ $pop_num ]} ) );
	  $searched_ids{ $new_copy_id } = 1;
	} while (($copying_sample[ $pop_num ][ $new_copy_id ] ||
			  $new_copy_id == $dont_copy_id) &&
			  (scalar keys %searched_ids !=
				      scalar @{$copying_sample[ $pop_num ]} ));

        if ($copying_sample[$pop_num][ $new_copy_id ]) {
	  die "Error: items out of sync" unless (scalar keys %searched_ids ==
					scalar @{$copying_sample[ $pop_num ]} );
	  my $num_haplotypes = @{$copying_sample[ $pop_num ]};
	  my $pop_name = $pop_names[$pop_num];
	  die "Error: used all $num_haplotypes haplotypes in pop $pop_name\n";
	}

	# found source that's not being copied from:
	$copying_sample[$pop_num][$new_copy_id] = 1;
	$copying_ind[ $ind ] = $new_copy_id;
      }
      die "Error End pos = $end_pos, marker = $marker\n" if($end_pos < $marker);

      # haplotype alleles for population <$pop_num> at <$marker>; will copy a
      # single allele corresponding to the sample that <$ind> is being copied
      # from
      my $allele_line = $pop_haps_lists[$pop_num]->[$marker];

      my $genotype = substr $allele_line, $copying_ind[ $ind ], 1;
      print OUT "$genotype";
    }

    print OUT "\n";
  }

  close OUT;

  # Note: all the genotype files we opened will close on exit
}


sub clean_up {
  # don't delete this file anymore
#  unlink $break_points_file;
}
