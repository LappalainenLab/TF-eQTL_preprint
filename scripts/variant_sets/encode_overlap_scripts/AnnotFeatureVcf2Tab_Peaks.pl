#build a table from annotation overlaps
#arguments: all_sites.vcf files_celltypes_assays_list >out

use strict;



open(VCF, "$ARGV[0]") or die "can't open input 1\n";
open(FILES, "$ARGV[1]") or die "can't open input 2\n";

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
#	print $cell;
	foreach my $assay (keys(%{$intfiles->{$cell}})) {
		my $id = join("-", $cell, $assay);
		push(@columns, $id);
		my $orig_file = $intfiles->{$cell}->{$assay};
		my @spl = split(/\./, $orig_file);
		my $intersect_file = join(".", "Peak_intersect/homo_sapiens.GRCh38", $cell, $assay, @spl[scalar(@spl)-5], "peaks.20161111.bed");
		open(INT, $intersect_file) or die "can't open $intersect_file\n";
		print STDERR "opened ", $intersect_file, "\n";
		
		while (<INT>) {
			chomp;
			my @ts = split(/\t/);
			$data -> {$ts[2]} -> {$id} = ();
		}
		close INT;
	}
}
print join("\t", "VARIANT", @columns), "\n";

while (<VCF>) {
	chomp;
	if ($_ =~ /^#/) {next}
	my @ts = split(/\t/);
	print $ts[2];
	foreach my $col (@columns) {	
		if (exists($data -> {$ts[2]} -> {$col})) {
			print "\t", "1";
		} else {
			print "\t", "0";
		}
	}
	print "\n";

}
close VCF;


