#!/usr/bin/env perl

use strict;
use warnings;

my $saw_first_line = 0; # read first line of bp file?
my $num_markers = -1;

while (<>) {
  if (!$saw_first_line) {
    $saw_first_line = 1;
    next; # first line gives pop labels
  }

  my @bps = split /\s/;

  my $cur_marker_num = 0;
  foreach my $bp (@bps) {
    my ($pop, $end_marker) = split /:/, $bp;
    for( ; $cur_marker_num <= $end_marker; $cur_marker_num++) {
      print "$pop";
    }
  }

  if ($num_markers < 0) {
    $num_markers = $cur_marker_num;
  }
  elsif ($num_markers != $cur_marker_num) {
    die "Saw different number of markers  on different lines?!\n";
  }

  print "\n";
}
