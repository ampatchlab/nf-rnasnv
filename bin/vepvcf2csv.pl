#!/usr/bin/env perl
################################################################################
#
#    vepvcf2tsv.pl
#
#    Parses an Ensembl VEP VCF into comma-separated values
#
#    Copyright (C) 2019 Stephen Kazakoff
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

use strict;
use warnings;

use Getopt::Long;
use IO::Zlib;
use Text::CSV;
use Text::ParseWords;
use URI::Escape;

my $input_file;
my $output_file;

GetOptions(
  'input_file|i=s' => \$input_file,
  'output_file|o=s' => \$output_file,
);

die "ERROR: No input file specified\n" unless $input_file;
die "ERROR: No output file specified\n" unless $output_file;

die "ERROR: Could not find input file: $input_file\n" unless -e $input_file;
die "ERROR: Output file already exists: $output_file\n" if -e $output_file;

my $csq_field_name = 'CSQ';

my $in_fh = IO::Zlib->new($input_file, 'rb');
my $out_fh = IO::Zlib->new($output_file, 'wb9');

my $meta = parse_headers($in_fh);

my $headers = $meta->{headers};
my $csv = Text::CSV->new({binary => 1, sep_char => ",", always_quote => 1, eol => "\n" });

my @info_headers = @{ $headers->{INFO} };
my @format_headers = @{ $headers->{FORMAT} };
my @csq_headers = @{ $headers->{$csq_field_name} };

my @out_col_names = qw(CHROM POS ID REF ALT QUAL FILTER);
push(@out_col_names, map { sprintf("INFO:%s", $_) } @info_headers);
for my $sample (@{ $meta->{samples} }) {
  push(@out_col_names, map { sprintf("%s/%s", $sample, $_) } @format_headers);
}
push(@out_col_names, @csq_headers);

$csv->column_names(@out_col_names);
$csv->print($out_fh, \@out_col_names);

while (my $line = <$in_fh>) {
  chomp $line;

  my %data;
  @data{@{ $meta->{column_names} }} = split(/\t/, $line, -1);

  my %out = map { $_ => $data{$_} } qw(CHROM POS ID REF ALT QUAL FILTER);

  # parse INFO fields
  my @csq_records;
  for my $key_value (split(/;/, $data{INFO})) {
    my ($k, $v) = split(/=/, $key_value, 2);
    if ($k eq $csq_field_name) {
      @csq_records = split(/,/, $v);
    }
    else {
      my $info_col_name = sprintf("INFO:%s", $k);
      $out{$info_col_name} = defined $v ? $v : 'True';
    }
  }

  # parse FORMAT fields
  my @format_fields = split(/:/, $data{FORMAT});
  for my $sample (@{ $meta->{samples} }) {
    my @sample_format_fields = map { sprintf("%s/%s", $sample, $_) } @format_fields;
    @out{@sample_format_fields} = map { $_ eq '' ? undef : $_ } split(/:/, $data{$sample}, -1);
  }

  # parse CSQ fields
  for my $rec (@csq_records) {
    @out{@csq_headers} = map { s/&/,/g; $_ eq '' ? undef : uri_unescape($_) } split(/\|/, $rec, -1);

    # write out the results
    $csv->print_hr($out_fh, \%out);
  }
}

sub parse_headers {
  my $fh = shift;

  my %metadata;

  while (my $line = <$fh>) {
    chomp $line;

    $line =~ s/^(#*)//;
    my $hash_count = length($1);

    if ($hash_count == 2) { # metadata-info lines
      if ($line =~ /^(INFO|FORMAT)=<(.*)>/) {
        my %fields = map { split(/=/, $_, 2) } grep { /=/ } quotewords(',', 0, $2);

        # skip any superfluous VEP 'custom' INFO fields
        next if exists($metadata{headers}{$csq_field_name}) && $1 eq 'INFO';

        if ($1 eq 'INFO' && $fields{ID} eq $csq_field_name) {
          $fields{Description} =~ m/Format: (.*)/;
          $metadata{headers}{$csq_field_name} = [split(/\|/, $1)];
        }
        else {
          push(@{ $metadata{headers}{$1} }, $fields{ID});
        }
      }

      next;
    }

    if ($hash_count == 1) { # column names
      my @cols = split(/\t/, $line);
      $metadata{column_names} = \@cols;
      my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = @cols;
      $metadata{samples} = \@samples;
    }

    last;
  }

  return \%metadata;
}
