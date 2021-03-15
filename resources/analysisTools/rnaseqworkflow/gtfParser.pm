package gtfParser;

use strict;

#######################################################################################################################
# AUTHOR: Nagarajan Paramasivam
# DATE: 26th Feb 2021
# EMAIL: n.paramasivam@dkfz.de

# Made fore HIPO2 RNAseq workfow
# Perl package containing functions to parse GTF files

# An example gene line in the GTF file (tab separated)
# chr1    HAVANA  gene    11869   14409   .       +       .       gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
#######################################################################################################################

sub parseGTFFile($$){

  my ($file_fh, $return_names) = @_;

  my %data          = ();
  my %data_genename = ();
  my $count         = 0;

  while (! eof($file_fh)) {

   my $line = readline($file_fh) 
      || die "Error reading from '$file_fh': $!";
  
    chomp $line;
    my @line_split = split ("\t", $line);

    if ($line_split[2] eq "gene") {
      
      my @comment = split (";",$line_split[8]);
      my ($gene_id, $gene_name);

      for(my $i = 0; $i <= $#comment; $i++){

        if($comment[$i] =~ /gene_id/){
          ($gene_id) = $comment[$i] =~ /\"(.*)\"/;
        }
        if($comment[$i] =~ /gene_name/){
          ($gene_name) = $comment[$i] =~ /\"(.*)\"/;
        }
      }

      if(defined $gene_id && $gene_id !~ /^$/){
        $data{$gene_id} = join("\t", 
                                @line_split[0,3,4], 
                                $gene_id, 
                                ".",
                                $line_split[6],
                                $gene_name
                              );

        $data_genename{$gene_id} = $gene_name;
        $count++;
      }
    }
  }

  if($return_names == 1){
    return(\%data, \$count, \%data_genename);
  } 
  else{
    return(\%data, \$count)
  }
}

1;