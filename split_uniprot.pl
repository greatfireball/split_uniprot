#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my $file = shift;
my $seqio_object = Bio::SeqIO->new(-file => $file);
my $seq_object   = $seqio_object->next_seq;

print join(";", $seq_object->get_keywords()), "\n";

my @classification = $seq_object->species()->classification();
