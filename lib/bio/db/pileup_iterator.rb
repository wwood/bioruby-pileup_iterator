require 'pp'

class Bio::DB::Pileup
  # Bio::DB::PileupIterator::PileupRead objects that occur at this position
  attr_accessor :reads
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
  # * Doesn't record the mapping quality of each read
  def each
    current_ordered_reads = []
    log = Bio::Log::LoggerPlus['bio-pileup_iterator']
    logging = true
    
    @io.each_line do |line|
      pileup = Bio::DB::Pileup.new(line.strip)
      current_read_index = 0
      reads_ending = []

      bases = pileup.read_bases
      log.debug "new column's read_bases: #{bases.inspect}" if log.debug?
      log.debug "pileup entry parsed: #{pileup.inspect}" if log.debug?
      while bases.length > 0
        log.debug "==== new read within a single pileup being parsed. Starting with #{bases}" if log.debug?
        
        # Firstly, what is the current read we are working with
        current_read = current_ordered_reads[current_read_index]
        # if adding a new read
        if current_read.nil?
          log.debug 'adding a new read: '+bases if log.debug?
          current_read = PileupRead.new
          current_ordered_reads.push current_read
        end
        matches = nil
        
        # if starting, remove it
        log.debug "before read start removal, pileup is #{bases}, read is #{current_read}" if log.debug?
        matched_string = ''
        if bases[0]=='^'
          # Match the ^ and the mapping quality
          matched_string += bases[0..1]
          bases = bases[2...bases.length]
        end
        log.debug "after read start removal, pileup is #{bases}" if log.debug?
        
        # next expect the actual base bit
        if matches = bases.match(/^([ACGTNacgtn\.\,\*])/)
          matched_string += bases[0]
          if matches[1] == '.'
            if !current_read.direction.nil? and current_read.direction != PileupRead::FORWARD_DIRECTION
              pp current_read
              raise "Unexpectedly found direction #{current_read.direction}, expected #{PileupRead::FORWARD_DIRECTION}, in starting at '#{bases}'(EndOfLine) in '#{line}', in the read above"
            end
            current_read.direction = PileupRead::FORWARD_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          elsif matches[1] == ','
            if !current_read.direction.nil? and current_read.direction != PileupRead::REVERSE_DIRECTION
              pp current_read
              raise "Unexpectedly found direction #{current_read.direction}, expected #{PileupRead::REVERSE_DIRECTION}, in starting at '#{bases}'(EndOfLine) in '#{line}', in the read above"
            end
            current_read.direction = PileupRead::REVERSE_DIRECTION
            current_read.sequence = "#{current_read.sequence}#{pileup.ref_base}"
          else
            # Could sanity check the direction here by detecting case, but eh
            current_read.sequence = "#{current_read.sequence}#{matches[1]}"
          end
          # remove the matched base
          bases = bases[1...bases.length]
        else
          raise Exception, "Expected a character corresponding to a base, one of '[ACGTNacgtn.,]'. Starting here: #{bases}, from #{pileup.inspect}"
        end
        log.debug "after regular position removal, pileup is #{bases}" if log.debug?
        
        # then read insertion or deletion in the coming position(s)
        if matches = bases.match(/^([\+\-])([0-9]+)/)
          matched_length = matches[1].length+matches[2].length
          bases = bases[matched_length...bases.length]
          matched_string += matches[1]+matches[2]
          log.debug "after removal of bases leading up to an insertion/deletion, pileup is #{bases}" if log.debug?
          
          regex = /^([ACGTNacgtn=]{#{matches[2].to_i}})/
          log.debug "insertion/deletion secondary regex: #{regex.inspect}" if log.debug?
          last_matched = bases.match(regex)
          
          if last_matched.nil?
            raise Exception, "Failed to parse insertion. Starting here: #{bases}, from #{pileup.inspect}"
          else
            bases = bases[last_matched[1].length...bases.length]
            if matches[1]=='+'
              # record the insertion
              current_read.add_insertion pileup.pos, matches[2], last_matched[1]
            elsif matches[1]=='-'
              #currently deletions are not recorded, slipped to future
            end
          end
        end
        log.debug "after indel removal, pileup is now #{bases}" if log.debug?
        
        # Then read an ending read
        if bases[0]=='$'
          reads_ending.push current_read_index
          matched_string += '$'
          bases = bases[1...bases.length]
        end
        
        log.debug "Matched '#{matched_string}', now the bases are '#{bases}'" if log.debug?
        
        current_read_index += 1
      end

      # Create a new copy of the array and yield that, otherwise when things get deleted they get removed from the yielded array as well (which is unwanted)
      yielded_array = Array.new(current_ordered_reads)
      pileup.reads = yielded_array
      #log.debug "Number of reads yielded: #{pileup.reads.length}"
      yield pileup
      
      # Remove reads that ended. In reverse order since removing the last ones first doesn't mess with the indices beforehand in the array
      reads_ending.reverse.each do |i|
        #log.debug "Deleting read of index #{i} (total reads #{current_ordered_reads.length}): #{current_ordered_reads[i].inspect}"
        current_ordered_reads.delete_at i
      end
      #log.debug "Ended up with #{current_ordered_reads.length} reads that should be present next time"
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