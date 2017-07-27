#!/usr/bin/awk -f
# Parses the MD tag from SAM files and adds the `TC:z:<pos>` tag for reads
# that contain T>C mutations at the given position.
# Can be used in a pipe from samtools view, but can not be piped back in to
# samtools easily, as it's missing the header!

### We need these helper functions to convert Phred from ASCII to a number

function _ord_init(    low, high, i, t)
{
    low = sprintf("%c", 7) # BEL is ascii 7
    if (low == "\a") {    # regular ascii
        low = 0
        high = 127
    } else if (sprintf("%c", 128 + 7) == "\a") {
        # ascii, mark parity
        low = 128
        high = 255
    } else {        # ebcdic(!)
        low = 0
        high = 255
    }

    for (i = low; i <= high; i++) {
        t = sprintf("%c", i)
        _ord_[t] = i
    }
}

function ord(str,    c)
{
    # only first character is of interest
    c = substr(str, 1, 1)
    return _ord_[c]
}

function chr(c)
{
    # force c to be numeric by adding 0
    return sprintf("%c", c + 0)
}

## END ASCII-to-int helper functions

function parseMdTag(md_str, read_seq, pos, flag, mut_array, mutKeeper) {
    # parseMdTag does as the name says, it takes the string from md_str,
    # iterates over it, and saves the position of each mismatch to the array
    # mutKeeper, where each mismatch is stored and adressable with the index
    # [P<position>]
    # Meanwhile, mut_array counts all mismatches in this read.
    regex="[A,C,G,T]"
    complement["A"]="T"
    complement["C"]="G"
    complement["G"]="C"
    complement["T"]="A"
    complement["N"]="N"
    seqLen=length(read_seq)
    if (match(md_str, regex)) {
        # m_len contains the number of matching bases
        m_len=substr(md_str, 1, RSTART-1)
        # ref_base contains the ref_seq at position of mismatch (m_len+1)
        ref_base=substr(md_str, RSTART, RLENGTH)
        # remain contains the remainder of the MD-flag, is passed on recursively.
        remain=substr(md_str, RSTART+RLENGTH)
        read_base=substr(read_seq, pos+m_len, 1)
        if (flag == 16) {
            ref_base=complement[ref_base]
            mut_base=complement[read_base]
        } else if (flag == 0) {
                mut_base=read_base
        }
        curMutPos=pos + m_len

        mutKeeper["P" curMutPos] = ref_base"-"mut_base
        # printf "\tR:%d-%s:%s\tmK:%s", curMutPos, ref_base"-"mut_base, remain, mutKeeper["P" + curMutPos]
        # !!! If read maps to neg. strand, read_seq position in SAM will be reverse complement !!!
        # !!! Here, ref_base refers to the reference base on the pos or neg strand.
        # !!! mut_base is base that actual read has instead at this position.
        # printf "  F:%s Match:%s Mut:%s-%s -- Seq_substr:%s\n", flag, m_len, ref_base, mut_base, substr(read_seq, pos, m_len)
        # printf "  F:%s Match:%s Mut:%s-%s", flag, m_len, ref_base, mut_base
        if (!(mut_base == "") && (mut_base != ref_base) && (mut_base != "N")) {
            mut_array[ref_base"-"mut_base]+=1
        }
        if (match(remain, regex)) {
            # Recursion galore, required for several mutations
            parseMdTag(remain, read_seq, pos+m_len+1, flag, mut_array, mutKeeper)
        }
    }
    return
}

{
    cigar=0
    md=0
    seq=0
    map_flag=0
    mut_array[""]=0
    split("", tcTags) # Initialise empty tcTags array
}
{
  _ord_init()
  if ($1 ~ /^@/) {
    print
  }
  else {
    map_flag=$2
    cigar=$6
    seq=$10
    regex="[A,C,G,T]"
    seqLen=length(seq)
    # Initate the array mutKeeper for every sequence.
    for (i=1;i<=seqLen;i++) {
        mutKeeper["P" i]=0
    }
    # Get MD Tag by iterating over accessory info columns
    for (i=11;i<=NF;i++) {
        # We only need to start at column 11, as this is where the accessory info
        # starts. If the beginning of column i matches MD:Z
        if ($i ~ /^MD:Z/) {
            md=substr($i,6,length($i))
        }
    }
    printf $0
    if (match(md,regex)) {
      # If there is an MD-Tag parse it and store in mut_array
      # This allows counting current mutations
      # mutKeeper is used to get positional information of where
      # the mutaions are
      parseMdTag(md, seq, 1, map_flag, mut_array, mutKeeper)
      for (i=1;i<=seqLen;i++) {
        # Do we have T>C mutations?
        seqPos = "P" i
        if (i <= 18) {
          if (mutKeeper[seqPos] == "T-C") {
            # Get phred score for this base
            phBaseQ = substr($11,i,1)
            baseQ = ord(phBaseQ) - 33
            # Only display T>C conversion if it meets the cutoff of 27
            # as discussed with Veronika
            if(baseQ > 27) {
              # printf "\tTC:i:%d", i
              tcTags["good"] = tcTags["good"] "," i
            } else {
              # printf "\tTN:i:%d", i
              tcTags["bad"] = tcTags["bad"] "," i
            }
          }
        }
      }
      if (length(tcTags) > 0) {
        if (length(tcTags["good"]) > 0) {
          printf "\tTC:Z:%s", substr(tcTags["good"],2,length(tcTags["good"]))
        }
        if (length(tcTags["bad"]) > 0) {
          printf "\tTN:Z:%s", substr(tcTags["bad"],2,length(tcTags["bad"]))
        }
      }
    }
    printf "\n"
  }
}
