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
  print "Eigenstrat format.\n";
  print "To simulate N haplotypes, 2N source haplotypes are required from ";
  print "each population\n\n";
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

  # open source population genotype files; store them in @pop_geno_fhs in the
  # same order as they are in @{$pop_indices_ref}
  my @pop_geno_fhs;
  my @test;
  foreach my $pop (@{$pop_indices_ref}) {
    my $pop_num = $hash_pop_names{ $pop };
    my $pop_file = $pop_files[ $pop_num ];
    my $ret = open my $fh, "$pop_file";
    if (!$ret) {
      &clean_up();
      die "Couldn't open $pop_file for $pop: $!\n";
    }
    push @pop_geno_fhs, $fh;
  }

  # open output file for admixed samples
  my $ret = open OUT, ">$out_prefix.phgeno";
  if (!$ret) {
    &clean_up();
    die "Couldn't open output file $out_prefix.phgeno: $!\n";
  }

  # We have two chromosomes for each population so that within-population
  # recombination breakpoints can be represented.  To introduce these
  # recombinations, we will copy from alternate copy numbers (chromosomes);
  # thus, these values will be 0 or 1. We have one copy number value per sample:
  my @copy_num_to_use;
  foreach (@{$person_bp_ref}) {
    push @copy_num_to_use, 0;
  }

  my $fh0 = $pop_geno_fhs[0];
  my $locus = 0;
  while (my $line = readline($fh0)) { # while there are still data in $fh0...
    my @locus_genos;
    push @locus_genos, $line;

    # Read in genotypes for all the source populations; have data for $fh0
    # (i.e., $pop_geno_fhs[0]), now read for the other populations (starting
    # from index 1, not 0):
    for (my $i = 1; $i < @pop_geno_fhs; $i++) {
      my $fh = $pop_geno_fhs[$i];
      $line = readline($fh);
      push @locus_genos, $line;
    }


    # Get the genotype for each sample and print it
    for (my $ind = 0; $ind < @{$person_bp_ref}; $ind++) { # for each sample...
      # All the unprocessed breakpoints for $indiv
      my $indiv_bps_ref = $person_bp_ref->[$ind];

      # get the next breakpoint downstream of $locus:
      my ($pop_num, $end_pos) = split /:/, $indiv_bps_ref->[0];
      if ($end_pos < $locus) {
	# Recombination break point!
	shift @{$indiv_bps_ref}; # remove the stale bp
	($pop_num, $end_pos) = split /:/, $indiv_bps_ref->[0]; # get next bp
	# swap the chromsome copy number the sample will use
	$copy_num_to_use[$ind] = ($copy_num_to_use[$ind] + 1) % 2;
      }
      die if ($end_pos < $locus);

      # all genotypes at this locus for the populaton to be copied from; will
      # copy a single genotype from this string (i.e., list)
      my $geno_line = $locus_genos[$pop_num];

      # sample $ind uses chromosome number 2*$ind or 2*$ind+1
      my $indiv_to_copy = 2 * $ind + $copy_num_to_use[$ind];
      if (length $geno_line <= $indiv_to_copy) {
	print "Have only " . (length($geno_line)+1) . " haplotypes in ";
	print "population " . $pop_names[$pop_num] . "\n";
	print "Need " . 2 * @{$person_bp_ref} . ". Dying.\n";
	exit 1;
      }
      my $genotype = substr $geno_line, $indiv_to_copy, 1;
      print OUT "$genotype";
    }

    print OUT "\n";

    $locus++; # going to next locus
  }

  close OUT;

  # Note: all the genotype files we opened will close on exit
}


# Old code from when I thought the genotype files were transposed with one
# sample per line (rather than one locus per line -- what it really is)
#sub parse_break_points_print_results {
#  open BP, "$break_points_file" or
#      die "Couldn't open simulator output file $break_points_file: $!\n";
#
#  my $header = <BP>;
#  chomp $header;
#  # indicates which population label corresponds to each index number
#  my @pop_order = split /\s+/, $header;
#
#  # open source population genotype files
#  my @pop_geno_fhs;
#  foreach my $pop (@pop_order) {
#    my $pop_num = $hash_pop_names{ $pop };
#    my $pop_file = $pop_files[ $pop_num ];
#    my $ret = open FILE, "$pop_file";
#    if (!$ret) {
#      &clean_up();
#      die "Couldn't open genotype file $pop_file: $!\n";
#    }
#    push @pop_geno_fhs, *FILE;
#  }
#
#  # open output file for admixed samples
#  my $ret = open OUT, ">$out_prefix.phgeno";
#  if (!$ret) {
#    &clean_up();
#    die "Couldn't open output file $out_prefix.phgeno: $!\n";
#  }
#
#  while (my $indiv_bps = <BP>) {
#    # remove trailing whitespace:
#    $indiv_bps =~ s/\s+$//;
#
#    # read two lines from each of the source populations for use with this
#    # sample; @genos will be array of arrays, with two elements for each
#    # population
#    my @genos;
#    foreach my $fh (@pop_geno_fhs) {
#      my $line1 = <$fh>;
#      my $line2 = <$fh>;
#      push @genos, [$line1, $line2];
#    }
#
#    my @bps = split /\s+/, $indiv_bps;
#    my $start_pos = 0;
#    # We have two genotypes for each population; to introduce breakpoints,
#    # we alternate which genotype we copy from for each break point.
#    my $copy_num_to_use = 0;
#    foreach my $bp (@bps) {
#      my ($pop_num, $end_pos) = split /:/, $bp;
#
#      my $copy_geno = $genos[$pop_num][$copy_num_to_use];
#      my $segment = substr $copy_geno, $start_pos, ($end_pos - $start_pos + 1);
#      print OUT $segment;
#
#      $start_pos = $end_pos + 1;
#      $copy_num_to_use = ($copy_num_to_use + 1) % 2;
#    }
#    print "\n";
#  }
#}

sub clean_up {
  # don't delete this file anymore
#  unlink $break_points_file;
}
