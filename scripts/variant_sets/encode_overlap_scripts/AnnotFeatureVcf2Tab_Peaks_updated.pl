#!/bin/perl
##
## Original version by TEL 2016
## Updated for ENCODE data EDF 2020-01
##
##	build a table from annotation overlaps
##
##	arguments: all_sites.vcf files_celltypes_assays_list >out

use strict;
use IO::Uncompress::Gunzip;



open my $VCF, '-|', '/usr/bin/gunzip', '-c', "$ARGV[0]" 
	or die "can't open input 1\n";
open(FILES, "$ARGV[1]") 
	or die "can't open input 2\n";

my $intfiles = {};
while (<FILES>) {
	chomp;
	my @ts = split(/\t/);
	$intfiles -> {$ts[1]} -> {$ts[2]} = $ts[0];
}
close FILES;

my $data = {};
my @columns;
foreach my $cell (keys(%{$intfiles})) { 
	#print $cell;
	foreach my $assay (keys(%{$intfiles->{$cell}})) {
		my $id = join("_", $assay, $cell);
		push(@columns, $id);
		my $intersect_file = join("/","peak_intersect_vcfs",$intfiles->{$cell}->{$assay});
		#my $orig_file = $intfiles->{$cell}->{$assay};
		#my @spl_1 = split(/\//, $orig_file);
		#my @spl = split(/\./, @spl_1[1] );
		#my $intersect_file = join(".", "peak_intersect/homo_sapiens.GRCh38", $cell, $assay, @spl[0], "peaks.20200116.bed");
		open my $INT, '-|', '/usr/bin/gunzip', '-c',  $intersect_file
			or die "can't open $intersect_file\n";
		print STDERR "opened ", $intersect_file, "\n";
		
		while ( <$INT> ) {
			chomp;
			my @ts = split(/\t/);
			$data -> {$ts[2]} -> {$id} = ();
		}
		close $INT;
	}
}
print join("\t", "#CHR", "POS", "VARIANT", @columns), "\n";

print STDERR "writing output file", "\n";
while ( <$VCF> ) {
	chomp;
	if ($_ =~ /^#/) {next}
	my @ts = split(/\t/);
	print $ts[0], "\t", $ts[1], "\t", $ts[2];
	foreach my $col (@columns) {	
		if (exists($data -> {$ts[2]} -> {$col})) {
			print "\t", "1";
		} else {
			print "\t", "0";
		}
	}
	print "\n";

}
close $VCF;


