#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my $file = shift;
my $seqio_object = Bio::SeqIO->new( -file => "<" . $file, -format => 'swiss' );

my @types =
  qw(bacteria archaea viruses eukaryota_not_metazoa eukaryota_and_metazoa);

my %filehandles = ();

foreach my $type (@types) {
    $filehandles{$type} = {
        fasta => Bio::SeqIO->new(
            -file   => ">" . $type . ".fasta",
            -format => "fasta"
        ),
        swiss =>
          Bio::SeqIO->new( -file => ">" . $type . ".sp", -format => "swiss" )
    };
}

# go through all sequences
while ( my $seq_object = $seqio_object->next_seq ) {

    # first check if the keyword "Complete proteome" is present
    next
      unless ( grep { $_ =~ /Complete proteome/i }
        $seq_object->get_keywords() );

    my @classification = $seq_object->species()->classification();

    foreach my $type (@types) {
        if ( test_by_string( $type, @classification ) ) {
            print "Valid for type: '$type'\n";
        }
        else {
            print "NOT Valid for type: '$type'\n";
        }
    }
}

sub test_by_string {
    my ( $type, @classification ) = @_;

    # test if type contains an not_.+
    my @not_allowed = $type =~ s/(_?not_[^_]+)//g;

    # all and can be deleted
    $type =~ s/and_//g;

    # the rest of the items are required
    my @required = split( /_+/, $type );

    # prepare hashs for faster access
    my %req;
    my %forbidden;
    @req{@required}          = map { 0 } (@required);
    @forbidden{@not_allowed} = map { 0 } (@not_allowed);

    # got though the classification and increase the hash values
    foreach my $item (@classification) {
        $item = lc($item);
        $forbidden{$item}++ if ( exists $forbidden{$item} );
        $req{$item}++       if ( exists $req{$item} );
    }

    # check that no forbidden item was marked as present
    return undef if ( grep { $forbidden{$_} > 0 } ( keys %forbidden ) );

    # and return true if all required items are marked as present
    return 1 if ( ( grep { $req{$_} > 0 } ( keys %req ) ) == ( keys %req ) );
}
