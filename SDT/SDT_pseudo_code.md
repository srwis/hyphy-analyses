#pseudo code
1) define single substitution Q matrix
  a)
2) add each together



check how many changes are within one codon at position i
  if single nucleotide change:
    check if in the last position:
      if it is in the last positon:
        check if the first nucleotide of the next codon is a change:
          if the first nucleotide of codon i+1 is not a change:
            treat it like a typical single nucleotide changes
          else if the first nucleotide of codon i+1 is a change:
            check if the second nucleotide of codon i+1 is a change:
              if it is not a change:
                treat like double hit across codons
              if it is a change:
                treat like triple hit across codons
      else is not in last position:
        treat like single nucleotide change

  if a double nucleotide change in codon at position i
    check if in last two positions:
      if in last two positions:
        check if the first nucleotide of the next codon is a change:
          if the first nucleotide of codon i+1 is not a change:
            treat it like a typical double nucleotide changes
          else if the first nucleotide of codon i+1 is a change:
            treat like a triple hit across codons
      else if not in last two positions:
        treat like typical double hits

  if triple nucleotide change in codon at position i:
    treat like typical triple hit
