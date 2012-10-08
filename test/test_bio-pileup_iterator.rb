require 'helper'

class TestBioPileupIterator < Test::Unit::TestCase
  #enable debug logging for pileup_iterator
  def setup
    log_name = 'bio-pileup_iterator'
    Bio::Log::CLI.logger('stderr')
    #Bio::Log::CLI.configure(log_name) # when commented out no debug is printed out
  end
  
  def test_pileup_parsing
    line = "contig00001\t199\tA\t4\t.$...$\t>a^>"
    #contig00001\t200\tT\t2\t..\taR"
    piles = Bio::DB::PileupIterator.new(line).to_a
    pileup = piles[0]
    reads = piles[0].reads
  
    assert_equal 'A', reads[0].sequence
    assert_equal 4, reads.length
    assert_kind_of Bio::DB::Pileup, pileup
  end
  
  def test_2_pileup_columns
    line = "contig00001\t199\tA\t4\t.$...$\t>a^>\ncontig00001\t200\tT\t2\t..\taR"
    piles = Bio::DB::PileupIterator.new(line).to_a
  
    pileup = piles[0]
    reads = piles[0].reads
    reads2 = piles[1].reads
  
    assert_equal 'A', piles[0].ref_base
    assert_equal 'T', piles[1].ref_base
    assert_equal 4, reads.length
    assert_equal 2, reads2.length
    assert_equal 'AT', reads2[0].sequence
  end
  
  def test_fwd_rev
    line = "contig00001\t199\tA\t4\t.$,..$\t>a^>\ncontig00001\t200\tT\t2\t,.\taR"
    piles = Bio::DB::PileupIterator.new(line).to_a
  
    pileup = piles[0]
    reads = piles[0].reads
    reads2 = piles[1].reads
  
    assert_equal 4, reads.length
    assert_equal 2, reads2.length
    assert_equal 'AT', reads2[0].sequence
    assert_equal '-', reads2[0].direction
    assert_equal '+', reads2[1].direction
  end
  
  def test_deletion
    line = "contig00001\t199\tA\t4\t.-1T...$\t>a^>\ncontig00001\t200\tT\t2\t*..\taR"
    piles = Bio::DB::PileupIterator.new(line).to_a
  
    pileup = piles[0]
    reads = piles[0].reads
    reads2 = piles[1].reads
  
    assert_equal 'A*', reads[0].sequence
    assert_equal Hash.new, reads[0].insertions
  end
  
  def test_substitution
    line = "contig00001\t199\tA\t4\t.G..$\t>a^>"
    piles = Bio::DB::PileupIterator.new(line).to_a
  
    pileup = piles[0]
    reads = piles[0].reads
  
    assert_equal 'A', reads[0].sequence
    assert_equal 'G', reads[1].sequence
    assert_equal 'A', reads[0].sequence
  end
  
  def test_substitution_with_insertion
    line = "contig00001\t199\tA\t4\tG-1T..$.\t>a^>\ncontig00001\t200\tT\t2\t*..\taR"
    piles = Bio::DB::PileupIterator.new(line).to_a
  
    pileup = piles[0]
    reads = piles[0].reads
    reads2 = piles[1].reads
  
    assert_equal 2, piles.length
    assert_equal 4, reads.length
    assert_equal 3, reads2.length
    assert_equal 'G*', reads[0].sequence
    assert_equal 'AT', reads[1].sequence
    assert_equal 'A', reads[2].sequence
    assert_equal 'AT', reads[3].sequence
  end
  
  def test_start_read_warning_of_deletion_next
    line = "contig00001\t8\tG\t4\t..,^],-1g\ta!U!\n"+
    "contig00001\t9\tg\t4\t..,*\ta!aU"
    piles = Bio::DB::PileupIterator.new(line).to_a
  
    pileup = piles[0]
    reads = piles[0].reads
    reads2 = piles[1].reads
  end
  
  def test_star_then_insert
    line = "contig00001\t23\tC\t40\t.*+1G..\t~~~~\n"
    piles = Bio::DB::PileupIterator.new(line).to_a
    
    pileup = piles[0]
    reads = piles[0].reads
    
    assert_equal 4, reads.length
    assert_equal({}, reads[0].insertions)
    assert_equal '*', reads[1].sequence
    assert_equal({23 => 'G'}, reads[1].insertions)
  end
  
  def test_star_finishing_a_read
    line = "contig00001\t717\tC\t47\t,$.$,$*$,$,$,$,$,$,$,$*$*$*$*$,$.$*$*$*$*$.$.$.$,$,$*$,$,$,$,$.$.$*$*$.$,$,$,$,$.$,$.$,$*$,$,$\t0..~2-.-.,#~~~~+,~~~~+**,!~!!!!!!~~(((((((((~!!\n"
    piles = Bio::DB::PileupIterator.new(line).to_a
    assert_equal '*', piles[0].reads[3].sequence
  end
  
  def test_start_finishing_a_read
    line = "contig00002\t1\tC\t47\t^],$\t~\n"
    piles = Bio::DB::PileupIterator.new(line).to_a
    assert_equal 'C', piles[0].reads[0].sequence
  end
  
  def test_start_with_a_gap
    line = "contig00075\t503\tT\t24\t,^]*\tU\n"
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal 'T', piles[0].reads[0].sequence
    assert_equal '*', piles[0].reads[1].sequence
  end
  
  def test_start_then_insert_then_end
    line = "contig00075\t503\tG\t24\t^].+1T$^].\t~~\n"
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal 'G', piles[0].reads[0].sequence
    assert_equal({503 => 'T'}, piles[0].reads[0].insertions)
    assert_equal 'G', piles[0].reads[1].sequence
  end
  
  def test_star_then_insert2
    line = "contig00075\t503\tG\t24\t,*+1g.\t~~\n"
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal 'G', piles[0].reads[0].sequence
    assert_equal '*', piles[0].reads[1].sequence
    assert_equal 'G', piles[0].reads[2].sequence
  end
  
  def test_start_with_gap_then_insertion
    line = "contig00075\t503\tG\t24\t,,.^]*+1g\tE~\n"+
    "contig00075\t504\tA\t24\t,,.,\tE~\n"
    
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal 'GA', piles[0].reads[0].sequence
    assert_equal 'GA', piles[0].reads[1].sequence
    assert_equal 'GA', piles[0].reads[2].sequence
    assert_equal '*A', piles[0].reads[3].sequence
    assert_equal({503 => 'g'}, piles[0].reads[3].insertions)
  end
  
  def test_double_insertion
    line = "contig00075\t503\tG\t24\t*+2gg\tE\n"
    
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal({503 => 'gg'}, piles[0].reads[0].insertions)
  end
  
  def test_non_perfect_starting_read
    line = "contig00075\t503\tG\t24\t^].*+2gg\tE\n"
    
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal '+', piles[0].reads[0].direction
    assert_equal 'G', piles[0].reads[0].sequence
    assert_equal '*', piles[0].reads[1].sequence
  end
  
  def test_non_matching_finish
    line = "contig00002\t6317\tC\t2\ta$.\t!B\n"+
     "contig00002\t6318\tT\t1\t.\tA\n"
  
    
    piles = Bio::DB::PileupIterator.new(line)
    piles = piles.to_a
    assert_equal 2, piles[0].reads.length
    assert_equal 'a', piles[0].reads[0].sequence
    assert_equal 'CT', piles[0].reads[1].sequence
  end
  
  def test_insertion_then_mismatch
    line = "contig00044\t867\tC\t6\t,,,,,.\t!:!!:=\n"+
     "contig00044\t868\tG\t6\tt,+1ttt,.\t!A!!C9\n"
  
    piles = Bio::DB::PileupIterator.new(line)
    
    piles = piles.to_a
    assert_equal 6, piles[0].reads.length
    assert_equal 'Ct', piles[0].reads[0].sequence
    assert_equal 'CG', piles[0].reads[1].sequence
    hash = {868=>'t'}
    assert_equal hash, piles[0].reads[1].insertions
    assert_equal 'Ct', piles[0].reads[2].sequence
  end
  
  def test_some_beta_testing_bug
    #<Bio::DB::Pileup:0x00000009d5e8d8 @ref_name="contig03007", @pos=658, @ref_base="a", @coverage=14.0, @read_bases="gg+1cG-1Ag+1c*+1aG-1A*+1a****g**+1A", @read_quals="!!!!~!~!!!!!!~", @ref_count=nil, @non_ref_count_hash=nil, @non_ref_count=nil> (Exception)
    line = "contig03007\t658\ta\t14\tgg+1cG-1Ag+1c*+1aG-1A*+1a****g**+1A\t!!!!~!~!!!!!!~\n"
    piles = Bio::DB::PileupIterator.new(line).to_a #parse, it should fail otherwise
    assert_equal 14, piles[0].coverage
  end
  
  def test_unexpected_equals
    line = "contig00032\t264\tg\t10\t*+1=$.*+1c,.*+1c...,\t~^~^^~^^^W\n"
    piles = Bio::DB::PileupIterator.new(line).to_a #parse, it should fail otherwise
    assert_equal 10, piles[0].coverage
    assert_equal 10, piles[0].reads.length
  end
  
  def test_optional_mapping_quality
    line = "gi|308171891|ref|NC_014551.1|\t2\tA\t2\t^:,^~,\t!!\n"
    piles = Bio::DB::PileupIterator.new(line).to_a #parse, it should fail otherwise
    assert_equal 2, piles[0].coverage
    assert_equal 2, piles[0].reads.length
    assert_equal 'A', piles[0].reads[0].sequence
    assert_equal '-', piles[0].reads[0].direction
    
    line_without_mapping_quality = "gi|308171891|ref|NC_014551.1|\t2\tA\t2\t^,^,\t!!\n"
    piles = Bio::DB::PileupIterator.new(line).to_a #parse, it should fail otherwise
    assert_equal 2, piles[0].coverage
    assert_equal 2, piles[0].reads.length
    assert_equal 'A', piles[0].reads[0].sequence
    assert_equal '-', piles[0].reads[0].direction
  end
  
  def test_when_read_mapping_quality_is_dot
    lines = "gi|308171891|ref|NC_014551.1|\t2\tA\t2\t^:,^~,\t!!\n"+
      "gi|308171891|ref|NC_014551.1|\t3\tT\t4\t,,^!.^.,\t!!!!\n"+ # This is the line that is really being tested
      "gi|308171891|ref|NC_014551.1|\t4\tT\t4\t,,.,\t!!!!\n"
    piles = Bio::DB::PileupIterator.new(lines).to_a #parse, it should fail otherwise
    assert_equal 2, piles[0].coverage
    assert_equal 4, piles[2].coverage
    assert_equal 2, piles[0].reads.length
    assert_equal 'ATT', piles[0].reads[0].sequence
    assert_equal '-', piles[0].reads[0].direction
  end
  
  def test_n
    lines = "gi|308171891|ref|NC_014551.1|\t111565\tN\t7\taaAaAAA\t~~~~~~~\n"
    piles = Bio::DB::PileupIterator.new(lines).to_a #parse, it should fail otherwise
    assert_equal 1, piles.length
    assert_equal 7, piles[0].coverage
    assert_equal 'a', piles[0].reads[1].sequence
    assert_equal 'N', piles[0].ref_base
  end
end
