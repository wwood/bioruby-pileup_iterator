require 'helper'

class TestBioPileupIterator < Test::Unit::TestCase
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
    # piles.log.level = Bio::Log::DEBUG
    piles = piles.to_a
    assert_equal 'T', piles[0].reads[0].sequence
    assert_equal '*', piles[0].reads[1].sequence
  end

  def test_start_then_insert_then_end
    line = "contig00075\t503\tG\t24\t^].+1T$^].\t~~\n"
    piles = Bio::DB::PileupIterator.new(line)
    # piles.log.level = Bio::Log::DEBUG
    piles = piles.to_a
    assert_equal 'G', piles[0].reads[0].sequence
    assert_equal({503 => 'T'}, piles[0].reads[0].insertions)
    assert_equal 'G', piles[0].reads[1].sequence
  end

  def test_star_then_insert2
    line = "contig00075\t503\tG\t24\t,*+1g.\t~~\n"
    piles = Bio::DB::PileupIterator.new(line)
    # piles.log.level = Bio::Log::DEBUG
    piles = piles.to_a
    assert_equal 'G', piles[0].reads[0].sequence
    assert_equal '*', piles[0].reads[1].sequence
    assert_equal 'G', piles[0].reads[2].sequence
  end
  
  def test_start_with_gap_then_insertion
    line = "contig00075\t503\tG\t24\t,,.^]*+1g\tE~\n"+
    "contig00075\t504\tA\t24\t,,.,\tE~\n"
    
    piles = Bio::DB::PileupIterator.new(line)
    # piles.log.level = Bio::Log::DEBUG
    piles = piles.to_a
    assert_equal 'GA', piles[0].reads[0].sequence
    assert_equal 'GA', piles[0].reads[1].sequence
    assert_equal 'GA', piles[0].reads[2].sequence
    assert_equal '*A', piles[0].reads[3].sequence
    assert_equal({503 => 'g'}, piles[0].reads[3].insertions)
  end
  
  def test_double_insertion
    line = "contig00075\t503\tG\t24\t*+1gg\tE\n"
    
    piles = Bio::DB::PileupIterator.new(line)
    # piles.log.level = Bio::Log::DEBUG
    piles = piles.to_a
    assert_equal({503 => 'gg'}, piles[0].reads[0].insertions)
  end
end
