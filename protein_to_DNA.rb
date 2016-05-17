require "open-uri"

module Transc_Transl
  # This module contains useful hashes
  # used in the Protein_sequence class

  # Notes:
  # The codons chosen are based on frequency in the human genome relative 
  # to synonymous codons.
  # GAC for aspartate is the default for B (aspartate or asparagine)
  # UGA (stop) is used for selenocystine.

  @@protein_DNA_codon = {
    :A => 'GCC', # Alanine, Ala
    :B => 'GAC', # aspartate or asparagine, Asx ASPARTATE DEFAULT
    :C => 'TGC', # Cysteine, Cys
    :D => 'GAC', # Aspartic Acid, Asp
    :E => 'GAG', # Glutamic Acid, Glu
    :F => 'TTC', # Phenylalanine, Phe
    :G => 'GGC', # Glycine, Gly
    :H => 'CAC', # Histidine, His
    :I => 'ATC', # Isoleucine, Ile
    :K => 'AAG', # Lysine, Lys
    :L => 'CTG', # Leucine, Leu
    :M => 'ATG', # Methionine, Met
    :N => 'AAC', # Asparagine, Asn
    :P => 'CCC', # Proline, Pro     
    :Q => 'CAG', # Glutamine, Gln
    :R => 'AGA', # Arginine, Arg
    :S => 'AGC', # Serine, Ser
    :T => 'ACC', # Threonine, Thr
    :U => 'TGA', #  selenocystine
    :V => 'GTG', # Valine, Val
    :W => 'TGG', # Tryptophan, Trp
    :Y => 'TAC', # Tyrosine, Tyr
    :- => '-',   # gap of indeterminate length
    :* => 'TGA', # Stop
    :X => 'unk'  # 'X' may be any amino acid
  }

  @@fasta_legals = @@protein_DNA_codon.map{|key, value| key.to_s }

end  

class Protein_sequence
  include Transc_Transl # includes module defined above

  def initialize
    # takes ARGV[1] as the outfile
    # sets the outfile to standard out if none given at command line
    ARGV[1] ? @outFile = open(ARGV[1], 'w') : @outFile = $stdout
    @protein_seq = get_FASTA
  end

  def get_FASTA
    # opens ARGV[0] as the input file.  can be a text webpage, or a local file
    # calling a webpage requires "open uri"
    protein_seq = ''
    open(ARGV[0], "r") do |f|
      @outFile.write f.gets #this is the doc line that should lead a FASTA file
      f.each_line { |line| protein_seq << line.chomp }
    end
    protein_seq
  end

  def is_FASTA?(seq)
    # takes a string and checks if it is FASTA legal
    chi = seq.chars.detect {|c| !@@fasta_legals.include?(c)}
    @outFile.write "Error! Unknown amino acid #{chi}\n" if chi
    !chi
  end

  def protein_to_codon
    # takes the @protein_seq from init and maps it to 
    # DNA codons from Transc_Transl
    # result is split to 80 character lines and written to the output
    if is_FASTA?(@protein_seq)
      a = @protein_seq.chars
      codons = a.map{|aminoacid| @@protein_DNA_codon[aminoacid.to_sym].chars}      
      @outFile.write codons.join.scan(/.{1,80}/).join("\n").downcase + "\n"
    else
      puts "The input must be a string, and a protein sequence in FASTA format."
      puts "For further feedback see #{ARGV[1]}" if ARGV[1] 
    end  
    @outFile.close
  end

end

puts   
Protein_sequence.new.protein_to_codon # calls the protein_to_codon method.

__END__

This program takes a protein sequence and returns a DNA sequence that 
would code for the protein represented.

tested on ruby 2.0.0p645 (2015-04-13 revision 50299) [universal.x86_64-darwin15]

usage: ruby protein_to_DNA.rb fileLoc outputfile
fileLoc may be a FASTA format protein sequence text file or a web pointer to 
a FASTA format protein sequence text file (such as found on the UniProt
protein sequence repository.)

A protein sequence in FASTA format begins with a single-line description, 
followed by lines of sequence data.
Example:
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
KEFTPPVQAAYQKVVAGVANALAHKYH

sample command lines:

$ruby protein_to_DNA.rb http://www.uniprot.org/uniprot/P11310.fasta outfile.txt
$cat outfile.txt
>sp|P11310|ACADM_HUMAN Medium-chain specific acyl-CoA dehydrogenase, mitochondrial OS=Homo sapiens GN=ACADM PE=1 SV=1
atggccgccggcttcggcagatgctgcagagtgctgagaagcatcagcagattccactggagaagccagcacaccaaggc
caacagacagagagagcccggcctgggcttcagcttcgagttcaccgagcagcagaaggagttccaggccaccgccagaa
agttcgccagagaggagatcatccccgtggccgccgagtacgacaagaccggcgagtaccccgtgcccctgatcagaaga
gcctgggagctgggcctgatgaacacccacatccccgagaactgcggcggcctgggcctgggcaccttcgacgcctgcct
gatcagcgaggagctggcctacggctgcaccggcgtgcagaccgccatcgagggcaacagcctgggccagatgcccatca
tcatcgccggcaacgaccagcagaagaagaagtacctgggcagaatgaccgaggagcccctgatgtgcgcctactgcgtg
accgagcccggcgccggcagcgacgtggccggcatcaagaccaaggccgagaagaagggcgacgagtacatcatcaacgg
ccagaagatgtggatcaccaacggcggcaaggccaactggtacttcctgctggccagaagcgaccccgaccccaaggccc
ccgccaacaaggccttcaccggcttcatcgtggaggccgacacccccggcatccagatcggcagaaaggagctgaacatg
ggccagagatgcagcgacaccagaggcatcgtgttcgaggacgtgaaggtgcccaaggagaacgtgctgatcggcgacgg
cgccggcttcaaggtggccatgggcgccttcgacaagaccagacccgtggtggccgccggcgccgtgggcctggcccaga
gagccctggacgaggccaccaagtacgccctggagagaaagaccttcggcaagctgctggtggagcaccaggccatcagc
ttcatgctggccgagatggccatgaaggtggagctggccagaatgagctaccagagagccgcctgggaggtggacagcgg
cagaagaaacacctactacgccagcatcgccaaggccttcgccggcgacatcgccaaccagctggccaccgacgccgtgc
agatcctgggcggcaacggcttcaacaccgagtaccccgtggagaagctgatgagagacgccaagatctaccagatctac
gagggcaccagccagatccagagactgatcgtggccagagagcacatcgacaagtacaagaac

$ruby protein_to_DNA.rb http://www.uniprot.org/uniprot/P68871.fasta
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2
atggtgcacctgacccccgaggagaagagcgccgtgaccgccctgtggggcaaggtgaacgtggacgaggtgggcggcga
ggccctgggcagactgctggtggtgtacccctggacccagagattcttcgagagcttcggcgacctgagcacccccgacg
ccgtgatgggcaaccccaaggtgaaggcccacggcaagaaggtgctgggcgccttcagcgacggcctggcccacctggac
aacctgaagggcaccttcgccaccctgagcgagctgcactgcgacaagctgcacgtggaccccgagaacttcagactgct
gggcaacgtgctggtgtgcgtgctggcccaccacttcggcaaggagttcaccccccccgtgcaggccgcctaccagaagg
tggtggccggcgtggccaacgccctggcccacaagtaccac

(after saving P68871.fasta to a local file:)
$ruby protein_to_DNA.rb P68871.fasta outfile.txt
$cat outfile.txt
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2
atggtgcacctgacccccgaggagaagagcgccgtgaccgccctgtggggcaaggtgaacgtggacgaggtgggcggcga
ggccctgggcagactgctggtggtgtacccctggacccagagattcttcgagagcttcggcgacctgagcacccccgacg
ccgtgatgggcaaccccaaggtgaaggcccacggcaagaaggtgctgggcgccttcagcgacggcctggcccacctggac
aacctgaagggcaccttcgccaccctgagcgagctgcactgcgacaagctgcacgtggaccccgagaacttcagactgct
gggcaacgtgctggtgtgcgtgctggcccaccacttcggcaaggagttcaccccccccgtgcaggccgcctaccagaagg
tggtggccggcgtggccaacgccctggcccacaavgtaccac

given a test fail file fail_example.fasta:
>sp|P0
MVHLJTPEEKS

$ruby protein_to_DNA.rb fail_example.fasta outfile.txt

The input must be a string, and a protein sequence in FASTA format.
For further feedback see outfile.txt
$cat outfile.txt
>sp|P0
Error! Unknown amino acid J

  





