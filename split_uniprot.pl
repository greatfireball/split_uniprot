#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my $file = shift;
my $seqio_object = Bio::SeqIO->new( -file => $file, -format => 'swiss' );

my @types =
  qw(bacteria archaea viruses eukaryota_not_metazoa eukaryota_and_metazoa);

my %filehandles = ();

foreach my $type (@types) {
    my $fasta_file_name = $type . ".fasta";
    $filehandles{$fasta_file_name} =
      Bio::SeqIO->new( -file => $fasta_file_name, -format => "fasta" );
    my $sp_file_name = $type . ".sp";
    $filehandles{$sp_file_name} =
      Bio::SeqIO->new( -file => $sp_file_name, -format => "swiss" );
}

# go through all sequences
while ( my $seq_object = $seqio_object->next_seq ) {
    print join( ";", $seq_object->get_keywords() ), "\n";

    my @classification = $seq_object->species()->classification();
}
