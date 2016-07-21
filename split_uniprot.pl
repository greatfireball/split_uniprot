#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my $verbose = 0;
eval {require Term::ProgressBar; };

unless ($@) { $verbose = 1; }


my $file = shift;
my $chunknumber = shift;
my $chunksize = shift;

#my $z    = new IO::Uncompress::Gunzip $file
#  or die "gunzip failed: $GunzipError\n";

open( my $z, "<", $file ) || die "Unable to open file: $!";

my $filesize = $chunksize;
my $startoffset = ($chunknumber-1)*$chunksize;
my $endoffset = $startoffset+$chunksize;

seek($z, $startoffset, 0) || die "$!";

# find the next end-block
# go through the input data and extract complete swissprot records
while ( !eof($z) && tell($z) < $endoffset) {
    my $line = <$z>;
    next unless ($line =~ /^\/{2}/);
}


my $progress;

if ($verbose)
{
    $progress = Term::ProgressBar->new(
	{
	    name  => 'Reading',
	    count => $filesize,
	    ETA   => 'linear',
	}
	);

    $progress->max_update_rate(1);
}

my $next_update = 0;

my $num_datasets = 0;

my @types =
  qw(bacteria archaea viruses eukaryota_not_metazoa eukaryota_and_metazoa);

my %filehandles = ();

foreach my $type (@types) {
    $filehandles{$type} = {
        fasta => Bio::SeqIO->new(
            -file   => ">" . sprintf("chunk%05d_%s.fasta", $chunknumber, $type),
            -format => "fasta"
        ),
        swiss =>
          Bio::SeqIO->new( -file => ">" . sprintf("chunk%05d_%s.sp", $chunknumber, $type), -format => "swiss" )
    };
}

my $inputblock = "";

# go through the input data and extract complete swissprot records
while ( !eof($z) && tell($z) < $endoffset) {

    # collect a complete swissprot block
    $inputblock = "";
    my $wanted_set =
      0;    # does the block contains Complete proteome or Reference proteome?

    while (<$z>) {
        $inputblock .= $_;
        $wanted_set = 1
          if ( $_ =~ /^KW.+Complete proteome|^KW.+Reference proteome/ );

        last if ( $_ =~ /^\/{2}/ );
    }

    if ($verbose)
    {
	$next_update = $progress->update( tell($z)-$startoffset )
	    if tell($z)-$startoffset > $next_update;
    }

    $num_datasets++;

    if ( $num_datasets % 100000 == 0 ) {
	if ($verbose)
	{
	    $progress->message( sprintf "Number of data sets %d" . '@' . "%d",
				$num_datasets, tell($z) );
	}
    }

    # check if the seq contains the keyword ($wanted_set is set)
    next unless ($wanted_set);

    open( my $tmp, "<", \$inputblock ) || die "$!";
    my $seqio_object = Bio::SeqIO->new( -fh => $tmp, -format => 'swiss' );

    # parse that sequence
    while ( my $seq_object = $seqio_object->next_seq ) {

        # first check if the keyword "Complete proteome" is present
        unless ( grep { $_ =~ /Complete proteome/i }
            $seq_object->get_keywords() )
        {

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
}

if ($verbose)
{
    $progress->update($filesize)
	if $filesize >= $next_update;
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

__END__
# Bacteria (2)
wget -O bacteria.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:2)&fil=&format=fasta&force=yes'
wget -O bacteria_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:2) keyword:1185&fil=&format=fasta&force=yes'

# Metazoa (33208)
wget -O metazoa.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:33208)&fil=&format=fasta&force=yes'
wget -O metazoa_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:33208) keyword:1185&fil=&format=fasta&force=yes'

# Viridiplantae (33090)
wget -O viridiplantae.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:33090)&fil=&format=fasta&force=yes'
wget -O viridiplantae_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:33090) keyword:1185&fil=&format=fasta&force=yes'

# Fungi (4751)
wget -O fungi.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:4751)&fil=&format=fasta&force=yes'
wget -O fungi_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:4751) keyword:1185&fil=&format=fasta&force=yes'

# Archaea (2157)
wget -O archaea.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:2157)&fil=&format=fasta&force=yes'
wget -O archaea_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:2157) keyword:1185&fil=&format=fasta&force=yes'

# Viruses (10239)
wget -O viruses.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:10239)&fil=&format=fasta&force=yes'
wget -O viruses_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:10239) keyword:1185&fil=&format=fasta&force=yes'

# Eukaryota (2759)
wget -O other_eukaryota.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:2759 -taxonomy:4751 -taxonomy:33090 -taxonomy:33208)&fil=&format=fasta&force=yes'
wget -O other_eukaryota_ref.fasta.gz 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=proteomes:(taxonomy:2759 -taxonomy:4751 -taxonomy:33090 -taxonomy:33208) keyword:1185&fil=&format=fasta&force=yes'

