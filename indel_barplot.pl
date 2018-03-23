#!/usr/bin/perl

#***************************************************************************************************
# FileName:
# Creator: Chen Y.L. <shenyulan@genomics.cn>
# Create Time: Wed Feb  4 16:19:36 CST 2015

# Description:
# CopyRight:
# vision: 0.1
# ModifyList:
#   Revision:
#   Modifier:
#   ModifyTime:
#   ModifyReason:
#***************************************************************************************************
#use strict;
#use warnings;

my $usage=<<usage;
    Usage: perl $0 <OT|outDir> <OT|insert_size.table> <IN|tag:indel-vcf file 1> <IN|tag:indel-vcf file 2> ...
    Example: perl $0 indel_plot insert_size.table HiSeq:$Bin/141101_I649_FCC5C90ACXX_L6_HUMntgEBAAAAAA-25.brecal.HC.indel.vcf CG_1Ad:GS82110-FS3-ON4610.sort.rmdup.gatkHC.mbq5.indel.vcf
usage

die($usage) unless @ARGV >1;

my ($outDir, $table, @vcfs) = @ARGV;
`mkdir $outDir` unless -e $outDir;

my @insert_sizes = ('<-10', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '>10' );
my ($color_code, $i) = ('',0);
open TB, "> $outDir/$table" or die $!;
open TBR, "> $outDir/$table.R" or die $!;
foreach(@vcfs){
    chomp;
    $i ++;
    $color_code .= "\'$i\',";
    my ($tag,$vcf) = (split /:/);
#    my @insertSizes = (split /\n/,`grep -v '^#' $vcf | awk 'length(\$4) > 1 || length(\$5) >1 {print length(\$5)-length(\$4)}'| sort -nr`);
    map{
        if ($_ < -10){
            $insert_sizes{$tag}{'<-10'}++;
        }elsif($_ > 10){
            $insert_sizes{$tag}{'>10'}++;
        }else{
            $insert_sizes{$tag}{$_}++;
        }
    }(split /\n/,`less $vcf | grep -v '^#' | awk 'length(\$4) > 1 || length(\$5) >1 {print length(\$5)-length(\$4)}'`);
}

map{print TB "\t$_";}(keys %insert_sizes);
print TB "\n";
foreach my $size (@insert_sizes){
    my $tag_size;
    foreach my $tag (keys %insert_sizes){
        $tag_size .= "\t$insert_sizes{$tag}{$size}";
    }
    print TB  "$size$tag_size\n";
}
$color_code =~ s/,$//;
print TBR "pdf(file=\"$outDir/$table.pdf\",w=10,h=8)\nm<-read.table(\"$outDir/$table\")\nbarplot(t(m),beside=T,col=c($color_code),legend=T,xlab='indel size',ylab='indel number')\naxis(1,seq(1,110,5),labels=F)\n";

system("/opt/blc/genome/biosoft/R/bin/R CMD BATCH $outDir/$table.R");
my $convert=(-f "/usr/bin/convert")?"/usr/bin/convert":"convert";
system("$convert -density 64 $outDir/$table.pdf $outDir/$table.png");
