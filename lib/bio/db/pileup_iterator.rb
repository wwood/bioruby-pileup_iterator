

class Bio::DB::Pileup
  # Bio::DB::PileupIterator::PileupRead objects that occur at this position
  attr_accessor :reads
  
  def num_deletions
    return read_bases.gsub(/[^\*]/,'').length
  end
end

class Bio::DB::Pileup
#contig00091 10  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB
attr_accessor :ref_name, :pos, :ref_base, :quality, :read_bases, :qualities
  def initialize(line)
    (@ref_name, @pos, @ref_base, @quality, @read_bases, @qualities) = line.split(/\s+/)
    @pos = @pos.to_i
  end
  
  def coverage
    quality.length
  end
end


class Bio::DB::PileupIterator
  include Enumerable
  
  def initialize(io)
    @io = io
  end

  # Iterates through the positions of the a pileup, returning an instance of Bio::DB::Pileup complete with an instance variable @reads, an Array of Bio::DB::PileupRead objects.
  #
  # Known problems:
  # * Doesn't record start or ends of each read
  # * Doesn't lookahead to determine the sequence of each read (though it does give the preceding bases)
  # * Gives no information with mismatches
  def each
    current_ordered_reads = []
    log = Bio::Log::LoggerPlus['bio-pileup_iterator']
    
    @io.each_line do |line|
      log.debug "new current_line: #{line.inspect}"
      pileup = Bio::DB::Pileup.new(line.strip)
      current_read_index = 0
      reads_ending = []

      bases = pileup.read_bases
      log.debug "new column's read_bases: #{bases.inspect}"
      log.debug "pileup entry parsed: #{pileup.inspect}"
      while bases.length > 0
        log.debug "bases remaining: #{bases}    ------------------------"
        
        # Firstly, what is the current read we are working with
        current_read = current_ordered_reads[current_read_index]
        # if adding a new read
        if current_read.nil?
          log.debug 'adding a new read'
          current_read = PileupRead.new
          current_ordered_reads.push current_read
        else
          log.debug 'reusing a read'
        end
        matches = nil

        # Now, parse what the current read is
        if matches = bases.match(/^([ACGTNacgtn\.\,])([\+\-])([0-9]+)([ACGTNacgtn]+)(\${0,1})/)
          log.debug "matched #{matches.to_s} as insertion/deletion"
          # insertion / deletion
          if matches[1] == '.'
            raise if !current_read.direction.nil? and current_read.direction != PileupRead::FORWARD_DIRECTION
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          elsif matches[1] == ','
            raise if !current_read.direction.nil? and current_read.direction != PileupRead::REVERSE_DIRECTION
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          else
            # Could sanity check the direction here by detecting case, but eh
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          end
          
          # record the insertion
          if matches[2] == '+'
            current_read.add_insertion pileup.pos, matches[3], matches[4]
          end
          
          if matches[5].length > 0
            log.debug "Ending this read"
            # end this read
            reads_ending.push current_read_index
          end
        # currently I don't care about indels, except for the direction, so I'll leave it at that for now

        # end of the read
        elsif matches = bases.match(/^([\.\,])\$/)
          log.debug "matched #{matches.to_s} as end of read"
          # regular match in some direction, end of read
          if matches[1]=='.' # if forwards
            raise if current_read.direction and current_read.direction != PileupRead::FORWARD_DIRECTION
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          else # else must be backwards, since it can only be , or .
            raise if current_read.direction and current_read.direction != PileupRead::REVERSE_DIRECTION
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          end
          log.debug "current read after deletion: #{current_read.inspect}"
          reads_ending.push current_read_index
          
        # regular match continuuing onwards
        elsif matches = bases.match(/^\./)
          log.debug "matched #{matches.to_s} as forward regular match"
          # regular match in the forward direction
          raise if !current_read.direction.nil? and current_read.direction != PileupRead::FORWARD_DIRECTION
          current_read.direction = PileupRead::FORWARD_DIRECTION
          log.debug "before adding this base, current sequence is '#{current_read.sequence}'"
          current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          log.debug "after adding this base, current sequence is '#{current_read.sequence}', ref_base: #{pileup.ref_base}"
        elsif matches = bases.match(/^\,/)
          log.debug "matched #{matches.to_s} as reverse regular match"
          # regular match in the reverse direction
          if !current_read.direction.nil? and current_read.direction != PileupRead::REVERSE_DIRECTION
            error_msg = "Unexpectedly found read a #{current_read.direction} direction read when expecting a positive direction one. This suggests there is a problem with either the pileup file or this pileup parser. Current pileup column #{pileup.inspect}, read #{current_read.inspect}, chomped until #{bases}"
            log.error error_msg
            raise Exception, error_msg
          end
          current_read.direction = PileupRead::REVERSE_DIRECTION
          current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
        
        # starting a new read (possibly with a gap), with an accompanying insertion/deletion
        elsif matches = bases.match(/^\^\]([ACGTNacgtn\.\,\*])([\+\-])([0-9]+)([ACGTNacgtn]+)(\${0,1})/)
          if matches[1] == '.'
            log.debug 'forward match starting a read'
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          elsif matches[1] == ','
            log.debug 'reverse match starting a read'
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          elsif matches[1] == '*'
            log.debug 'starting a read with a gap'
            # leave direction unknown at this point
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          elsif matches[1] == matches[1].upcase
            log.debug 'forward match starting a read, warning of insertion next'
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          else
            log.debug 'forward match starting a read, warning of insertion next'
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          end
          
          # record the insertion
          if matches[2] == '+'
            current_read.add_insertion pileup.pos, matches[3], matches[4]
          end
          
          if matches[5].length > 0
            log.debug "Ending this read"
            # end this read
            reads_ending.push current_read_index
          end

          
        # regular match, starting a new read
        elsif matches = bases.match(/^\^\]([ACGTNacgtn\.\,\*])(\${0,1})/)
          if matches[1] == '.'
            log.debug 'forward match starting a read'
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          elsif matches[1] == ','
            log.debug 'reverse match starting a read'
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          elsif matches[1] == '*'
            log.debug 'gap starting a read'
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          elsif matches[1] == matches[1].upcase
            log.debug 'forward match starting a read, warning of insertion next'
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          else
            log.debug 'forward match starting a read, warning of insertion next'
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          end
          if matches[2].length > 0
            log.debug "Ending this read, even though it started here too.. it happens.."
            # end this read
            reads_ending.push current_read_index
          end

          
        elsif matches = bases.match(/^\*([\+\-])([0-9]+)([ACGTNacgtn=]+)(\${0,1})/)
          log.debug 'gap then insert/delete found'
          # gap - should already be known from the last position
          current_read.sequence = "#{current_read.sequence}*"
          if matches[4].length > 0
            log.debug "Ending this read"
            # end this read
            reads_ending.push current_read_index
          end
                    
          # record the insertion
          if matches[1] == '+'
            current_read.add_insertion pileup.pos, matches[2], matches[3]
          end

          
        elsif matches = bases.match(/^\*(\${0,1})/)
          log.debug 'gap found'
          # gap - should already be known from the last position
          current_read.sequence = "#{current_read.sequence}*"
          if matches[1].length > 0
            log.debug "Ending this read"
            # end this read
            reads_ending.push current_read_index
          end
          
        elsif matches = bases.match(/(^[ACGTNacgtn])/)
          log.debug 'regular mismatch found'
          # simple mismatch
          current_read.sequence = "#{current_read.sequence}#{matches[1]}"
        end
        log.debug "current read's sequence: #{current_read.sequence}"
        
        #raise Exception, "implement mismatch parsing here!!!"
        raise Exception, "Unexpected Pileup format bases, starting here: #{bases}, from #{pileup.inspect}" if matches.nil?
        
        #remove the matched part from the base string for next time
        bases = bases[matches.to_s.length..bases.length-1]

        current_read_index += 1
      end

      # Create a new copy of the array and yield that, otherwise when things get deleted they get removed from the yielded array as well (which is unwanted)
      yielded_array = Array.new(current_ordered_reads)
      pileup.reads = yielded_array
      log.debug "Number of reads yielded: #{pileup.reads.length}"
      yield pileup
      
      # Remove reads that ended. In reverse order since removing the last ones first doesn't mess with the indices beforehand in the array
      reads_ending.reverse.each do |i|
        log.debug "Deleting read of index #{i} (total reads #{current_ordered_reads.length}): #{current_ordered_reads[i].inspect}"
        current_ordered_reads.delete_at i
      end
      log.debug "Ended up with #{current_ordered_reads.length} reads that should be present next time"
    end
  end

  class PileupRead
    # Directions relative to reference
    FORWARD_DIRECTION = '+'
    REVERSE_DIRECTION = '-'

    # sequence is always in the direction of the start of the reference to the end - only @direction gives direction information
    attr_accessor :direction, :sequence
    
    # A hash of recorded insertions. Key of hash is the position in the consensus that is has been added to in the alignment, value is the bases that have been inserted
    attr_reader :insertions
    
    def initialize
      @sequence = ''
      @insertions = {}
    end
    
    def read
      @sequence[@sequence.length-2..@sequence.length-1]
    end
    
    def add_insertion(position, insertion_length, insertion_bases)
      insertions[position] = insertion_bases
    end
  end
end