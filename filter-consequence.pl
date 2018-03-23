#! /usr/bin/perl
use File::Basename;
use threads;
use DBI;
use FindBin qw($Bin);
use Getopt::Std;
use lib ("$FindBin::RealBin");

do "$Bin/resource.pm";

my ($vcf, $output, $config) = @ARGV;

my ($dbs, $db_count);

my $result_dir = dirname $output;
my $result_file = basename $output;

my $md5vcf = `md5sum $vcf | /bin/awk '{print \$1}'`;
chomp $md5vcf;

my (%pd_threshold) = &parseconfig($config);

foreach (keys %pd_threshold){
	$dbs .= "$_,";
	$db_count ++;
}

$dbs =~ s/,$//g;
#print "dbs: $dbs\n";

#for my $ssig_chr (@chrs)
#for my $ssig_chr (1..22, 'X', 'Y')
for my $ssig_chr (21)
{
    my $pid=fork();
    die "Cannot fork: $!" if (! defined $pid);
    if($pid ==0)
    {
#		screen: ($vcf, $sig_chr, $out_dir, $out_file) = @_;
        &screen($vcf, $ssig_chr, $result_dir, "$ssig_chr.$md5vcf");
        exit;
    }
    else
    {
        push (@child_pids, $pid);
    }
} 

for my $pid (@child_pids)
{
   waitpid($pid, 0);
}

#for my $ssig_chr (@chrs)
#for my $ssig_chr (1..22, 'X', 'Y')
for my $ssig_chr (21)
{
#    `sort -t "  " -k2 -n $result_dir/$ssig_chr.$md5vcf.tsv >> ${out_result}.tsv;rm $result_dir/$ssig_chr.$md5vcf.tsv` if -e "$result_dir/$ssig_chr.$md5vcf.tsv";
#    `sort -t "	" -k2 -n $result_dir/$ssig_chr.$md5vcf.screen >> $output.screen && rm $result_dir/$ssig_chr.$md5vcf.screen` if -e "$result_dir/$ssig_chr.$md5vcf.screen";
#    `sort -t "	" -k2 -n $result_dir/$ssig_chr.$md5vcf.highAF >> $output.highAF && rm $result_dir/$ssig_chr.$md5vcf.highAF` if -e "$result_dir/$ssig_chr.$md5vcf.highAF";
    `sort -t "	" -k2 -n $result_dir/$ssig_chr.$md5vcf.CC >> $output.CC && rm $result_dir/$ssig_chr.$md5vcf.CC` if -e "$result_dir/$ssig_chr.$md5vcf.CC";
#    `sort -t "	" -k2 -n $result_dir/$ssig_chr.$md5vcf.dbsnp >> $output.dbsnp && rm $result_dir/$ssig_chr.$md5vcf.dbsnp` if -e "$result_dir/$ssig_chr.$md5vcf.dbsnp";
}

#open VCF, "$vcf" or die $!;

sub parseconfig {
	my ($config) = @_;
	my ($flag, $flag_pd, $config_app, $config_threshold);
	&mkopen('CFG', $config, "<");
	while(<CFG>){
		chomp;
		next if /^#/;
		next if /^$/;
		
#print "$_\n";
		$flag_pd = 1 if /^FilterScreenStack:Consequence:begin/;
		$flag = 1 if /^\[Consequence:default\]/;
		next if /^\[Consequence:default\]/;
		$flag = 2 if /^\[Consequence:user defined\]/;
		next if /^\[Consequence:user defined\]/;
		$flag_pd = 0 if /^FilterScreenStack:Consequence:end/;

#print "flag: $flag\n";

		if($flag == 1){
			$config_default .= "$_,";
		}elsif($flag == 2){
			$config_user_defined .= "$_,";
		}else{
#			die "config erro: Population data, please check!\n";
		}
	}

	$config_app = $config_user_defined if $flag_pd == 0 and $flag == 2;
	$config_app = $config_default if $flag_pd == 0 and $flag == 1;

	map{
		s/\s//g;
		my ($db, $threshold) = ($1, $2) if /(\S+)\(>(\S+)\)/;
		$config_threshold{$db} = $threshold if $db =~ /\S/;
#		print "$db => $threshold\n" if $db =~ /\S/;
	}(split /,/, $config_app);
	
	return (%config_threshold);
}

sub screen {
	my ($vcf, $sig_chr, $out_dir, $out_file) = @_;
	my ($v_num, $v_snv_num, $v_indel_num, $v_highaf_num, $v_snv_highaf_num, $v_indel_highaf_num) = (0,0,0,0);
	&mkopen('VCF', $vcf, "<");
#	&mkopen('OT', "$out_file.screen", ">");
	&mkopen('OT', "$out_file.AF", ">");
#	&mkopen('OT2', "$out_file.dbsnp", ">");
	while(<VCF>){
		chomp;
#		print OT "$_\n" if /^#/;
		next if /^##/;
		next if /^#/;
		my ($chr, $pos, $rs, $ref, $alt, $qua, $filter, $info, $format, @sample_infos) = (split /\t/);
		$chr =~ s/chr//;
		next unless $chr eq $sig_chr;

		my $sample_info = (join "\t", @sample_infos);
		foreach my $allele (split /,/, $alt){
		$v_num  ++;

		my ($temp_v_snv_num, $temp_v_indel_num, $ref_dbpp_query_results) = &filter_concequence($chr, $pos, $ref, $allele, $dbs);
		$v_snv_num += $temp_v_snv_num;
		$v_indel_num += $temp_v_indel_num;
		my @dbpp_query_results = @$ref_dbpp_query_results;

		my ($is_pass, $i, $af_info, $db_num) = ('PASS', 0, '', 0);
		foreach (split /,/, $dbs){
			if ($dbpp_query_results[$i] > $pd_threshold{$_}){
				$af_info .= "$_=$dbpp_query_results[$i];";
				$is_pass = 'HighAF';
				$db_num ++;
			}
			$i ++;
		}

		$af_info = '.' unless $af_info =~ /\w+/;
		$af_info =~ s/;$//;

		if ($is_pass eq 'HighAF'){
			$v_highaf_num ++ ;
			if(length($ref) == 1 && length($allele) == 1){
				$v_snv_highaf_num ++;
			}else{
				$v_indel_highaf_num ++;
			}
		}
		print OT "chr$chr\t$pos\t$rs\t$ref\t$allele\t$is_pass\t$db_num\t$af_info\n";
		}
	}
	print "chr$sig_chr, total allele: $v_num (snv: $v_snv_num | indel: $v_indel_num ), high AF allele:$v_highaf_num (snv: $v_snv_highaf_num | indel: $v_indel_highaf_num )\n";
	close VCF;
}

sub filter_concequence {
	my ($chr, $pos, $ref, $allele, $dbs) = @_;
	my ($db_file, $dbh, $db_ary_ref, @query_results) = ('', 0, 0);

        $db_file="$_wfdata_{dbsnp}/dbsnp.chr$chr.db";

	my ($key, $chr_r, $pos_r, $end, $ref_r, $allele_r, $vartype, $rs_r, $gene,$c_feature, $c_biotype, $c_symbol, $c_consequence, $c_impact, $c_hgvsc, $c_hgvsp, $c_abb_aac, $c_ei, $c_domains, $c_swissprot, $c_trembl, $c_uniparc, $c_sift, $c_polyphen, $m_feature, $m_biotype, $m_symbol, $m_consequence, $m_impact, $m_hgvsc, $m_hgvsp, $m_abb_aac, $m_ei, $m_domains, $m_swissprot, $m_trembl, $m_uniparc, $m_sift, $m_polyphen, $d_feature, $d_biotype, $d_symbol, $d_consequence, $d_impact, $d_hgvsc, $d_hgvsp, $d_abb_aac, $d_ei, $d_domains, $d_swissprot, $d_trembl, $d_uniparc, $d_sift, $d_polyphen, $all_transcript_info, $forward_sequence_context, $reverse_sequence_context, $forward_up50, $forward_down50, $reverse_up50, $reverse_down50);

        if(-f $db_file){
                $dbh = DBI->connect("dbi:SQLite:dbname=$db_file","","");
                $db_ary_ref=$dbh->selectall_arrayref("SELECT $dbs FROM dbpp WHERE key37 = '$chr:$pos:$ref:$allele'");
		 map{
                    foreach my $i (0..60){
                        $query_results[$i] = @$_[$i];
                    }
                } @$ary_ref;
                $dbh->disconnect();
        }else{
                $rs_r = '.' unless -f $db_file;
        }

	($key, $chr_r, $pos_r, $end, $ref_r, $allele_r, $vartype, $rs_r, $gene,$c_feature, $c_biotype, $c_symbol, $c_consequence, $c_impact, $c_hgvsc, $c_hgvsp, $c_abb_aac, $c_ei, $c_domains, $c_swissprot, $c_trembl, $c_uniparc, $c_sift, $c_polyphen, $m_feature, $m_biotype, $m_symbol, $m_consequence, $m_impact, $m_hgvsc, $m_hgvsp, $m_abb_aac, $m_ei, $m_domains, $m_swissprot, $m_trembl, $m_uniparc, $m_sift, $m_polyphen, $d_feature, $d_biotype, $d_symbol, $d_consequence, $d_impact, $d_hgvsc, $d_hgvsp, $d_abb_aac, $d_ei, $d_domains, $d_swissprot, $d_trembl, $d_uniparc, $d_sift, $d_polyphen, $all_transcript_info, $forward_sequence_context, $reverse_sequence_context, $forward_up50, $forward_down50, $reverse_up50, $reverse_down50) = @query_results unless $rs_r eq '.';

	return (@query_results);
}

sub mkopen {
    my ($fh,$file,$mode)=@_;
    if( $file =~ /\.gz$/ ){
        open $fh ,"gzip -dc $file |" or die "Cannot open input file $file" unless $mode =~ />/;
        open $fh ,"| gzip > $file" or die "Cannot open output file $file" if $mode eq '>';
        open $fh ,"| gzip >> $file" or die "Cannot open output file $file" if $mode eq '>>';
    }elsif( $file =~ /\.bz2$/ ){
        open $fh ,"bunzip2 -dc $file |" or die "Cannot open input file $file" unless $mode =~ />/;
        open $fh ,"| bzip2 > $file" or die "Cannot open output file $file" if $mode eq '>';
        open $fh ,"| bzip2 >> $file" or die "Cannot open output file $file" if $mode eq '>>';
    }else{
       open $fh, $file or die "Cannot open input file $file" unless $mode =~ />/;
        open $fh, "> $file" or die "Cannot open output file $file" if $mode eq '>';
        open $fh, ">> $file" or die "Cannot open output file $file" if $mode eq '>>';
	}
}
