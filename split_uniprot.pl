#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;
use IO::Uncompress::Gunzip qw($GunzipError);

use Term::ProgressBar 2.00;

my $file = shift;
my $z    = new IO::Uncompress::Gunzip $file
  or die "gunzip failed: $GunzipError\n";

my $filesize = 169981337736;

my $progress = Term::ProgressBar->new(
    {
        name  => 'Reading',
        count => $filesize,
        ETA   => 'linear',
    }
);

$progress->max_update_rate(1);
my $next_update = 0;

my $num_datasets = 0;

my $seqio_object = Bio::SeqIO->new( -fh => $z, -format => 'swiss' );

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

    $num_datasets++;

    $next_update = $progress->update( $z->tell() )
      if $z->tell() > $next_update;

    if ( $num_datasets % 100000 == 0 ) {
        $progress->message( sprintf "Number of data sets %d", $num_datasets );
    }

    # first check if the keyword "Complete proteome" is present
    unless ( grep { $_ =~ /Complete proteome/i } $seq_object->get_keywords() ) {

        #	print "Skipping entry ".$z->tell()."\n";
        next;
    }

    my @classification = $seq_object->species()->classification();

    foreach my $type (@types) {
        if ( test_by_string( $type, @classification ) ) {

            #           print "Valid for type: '$type'\n";
            # store the sequence in the subset
            foreach my $file ( keys %{ $filehandles{$type} } ) {
                $filehandles{$type}{$file}->write_seq($seq_object);
            }
        }
        else {
            #            print "NOT Valid for type: '$type'\n";
        }
    }
}

$progress->update($filesize)
  if $filesize >= $next_update;

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
