#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);

my %opts;
GetOptions(\%opts,"i=s","o=s","h");

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h})) {
    print "Usage: perl $Script -i fragments.tsv.gz -o /path/to/prefix\n";
    print "Example: perl $Script -i input.gz -o /results/M8.filter4\n";
    exit;
}

my $umifile = $opts{i};
my $base_prefix = $opts{o};
my $outfile = "${base_prefix}_result.txt";

# ---------- Step 0: 预扫描所有 fragments ----------
print "[", scalar(localtime), "] Pre-scanning fragments...\n";
open (PRE, "gzip -dc $umifile |") || die "Cannot open $umifile: $!";
my %bc_all = ();
while(<PRE>){
    next if (m/^\#/);
    my ($chr, $start, $end, $bc, $num) = split;
    next unless ($chr =~ /^chr(?:[1-9]|1[0-9]|X)$/);  # 只保留核基因组
    $bc_all{$bc} += $num;
}
close PRE;

my $total_barcodes = scalar(keys %bc_all);
print "Total barcodes: $total_barcodes\n";

# ---------- Step 1: 加载全部有效 reads 用于抽样 ----------
print "[", scalar(localtime), "] Loading valid fragments...\n";
open (IN, "gzip -dc $umifile |") || die $!;
my ($total_reads, $index, @simu_arr, @idx_to_bc) = (0, 0, (), ());
while(<IN>){
    next if (m/^\#/);
    my ($chr, $start, $end, $bc, $num) = split;
    next unless ($chr =~ /^chr(?:[1-9]|1[0-9]|X)$/);
    $index++;
    $total_reads += $num;
    push @idx_to_bc, $bc;
    push @simu_arr, ($index) x $num;
}
close IN;

# ---------- Step 2: 随机打乱 ----------
print "[", scalar(localtime), "] Shuffling...\n";
shuffle(\@simu_arr);

# ---------- Step 3: 逐步抽样，计算统计量 ----------
print "[", scalar(localtime), "] Sampling...\n";
my ($percent, $current_uniq, @is_seen, %hstat, %bc_counts) = (12.5, 0, (), (), ());
for (my $i=0; $i<@simu_arr; $i++){
    my $f_idx = $simu_arr[$i];
    if (!defined $is_seen[$f_idx]){
        $is_seen[$f_idx] = 1;
        $current_uniq++;
        $bc_counts{$idx_to_bc[$f_idx-1]}++;
    }
    if ($i+1 >= $total_reads * $percent / 100 || $i+1 == @simu_arr){
        my @counts = values %bc_counts;
        # 计算 median nFrags，基于当前有 reads 的 barcode，不再补零到 total_barcodes
        my $median = calculate_median_from_array(\@counts);
        $hstat{$percent} = [$current_uniq, $i+1, (1-$current_uniq/($i+1))*100, $median];
        $percent += 12.5;
        last if $percent > 100.1;
    }
}

# ---------- Step 4: 输出结果 ----------
open (OUT, ">$outfile");
print OUT "Percent\tUnique_Fragments\tTotal_Fragments\tPercent_Duplicates\tMedian_nFrags\n";
foreach my $p (sort {$a<=>$b} keys %hstat){
    printf OUT "%.1f\t%d\t%d\t%.2f\t%.2f\n", $p, @{$hstat{$p}};
}
close OUT;

print "[", scalar(localtime), "] Done. Result: $outfile\n";

# ========== Subroutines ==========
sub calculate_median_from_array {
    my $arr = shift;
    return 0 unless @$arr;
    my @sorted = sort {$a<=>$b} @$arr;
    my $len = scalar @sorted;
    if ($len % 2 == 0) {
        return ($sorted[$len/2 - 1] + $sorted[$len/2]) / 2;
    } else {
        return $sorted[int($len/2)];
    }
}

sub shuffle {
    my $a = shift;
    for (my $i = @$a - 1; $i > 0; $i--){
        my $j = int(rand($i + 1));
        @$a[$i,$j] = @$a[$j,$i];
    }
}