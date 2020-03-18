#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
#    start = rna_sequence.find('AUG')
#    protein = ""
#    i = start
#    while i < len(rna_sequence)-2:
#        codon = rna_sequence[i:i+3]
#        a = genetic_code
#        if a == '*':
#            break
#        i += 3
#        protein.append(a)
#    return ''.join(protein)

    rna_sequence=rna_sequence.upper()
    protein = ""
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i: i+3]
        amino_acid = genetic_code.get(codon, '*')
        if amino_acid != '*':
            protein += amino_acid
        else:
            break
    return protein

#if genetic_code == "UAA" "UAG" "UGA" :

#    start = rna_sequence.find('AUG')
#    peptide = []
#    i = start
#    while i < len(rna_sequence)-2:
#        codon = rna_sequence[i:i+3]
#        a = genetic_code[codon]
#        if a == '*' :
#            break
#        i += 3
#        peptide.append(a)
#    return ''.join(peptide)

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
##the find method will give you the first index position (int) of the codon you're referring to 
#    start_codon = rna_sequence.find('AUG')
#
#    ##You create a string of the starting nucleotide (A in this case) and rest of the sequence 
#    seq_start = rna_sequence[int(start_codon):]
#    ##This selects the first nucleotide in codon/starting with with zero, selects for multiples of 3; n are multiples of 3 starting with 0
#    protein_seq = ''
#    for n in (0, len(seq_start), 3):
#        if seq_start[n:n+3] in genetic_code:
#            protein_seq += genetic_code[seq_start[n:n+3]]
#    return protein_seq

    rna_sequence=rna_sequence.upper()
    start = rna_sequence.find('AUG')
    protein = []
    i = start
    while i < len(rna_sequence)-2:
        codon = rna_sequence[i:i+3]
        amino_acid = genetic_code[codon]
        if amino_acid == '*':
            break
        i += 3
        protein.append(amino_acid)
    return ''.join(protein)

#    rna_sequence=rna_sequence.upper()
#    protein = ""
#    for i in range(0, len(rna_sequence), 3):
#        codon = rna_sequence[i: i+3]
#        amino_acid = genetic_code.get(codon, '*')
#        if amino_acid != '*':
#            protein += amino_acid
#        else:
#            break
#    return protein

#    for i in range(0, len(rna_sequence), 3):
#        codon = rna_sequence[i: i+3]
#        amino_acid = genetic_code.get(codon, '*')
#        if amino_acid != '*':
#            protein += amino_acid
#        else:
#            break
#    return protein

#    peptide = []
#    i = start
#    while i < len(rna_sequence)-2:
#        codon = rna_sequence[i:i+3]
#        a = genetic_code[codon]
#        if a == '*' :
#            break
#        i += 3
#        peptide.append(a)
#    return ''.join(peptide)

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    >>> get_reverse('ATGC')
    'CGTA'
    """

    sequence = sequence.upper()
    return sequence[::-1]

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'UACG'
    >>> get_reverse('ATGC')
    'TACG'
    """

    sequence = sequence.upper()
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    return ''.join([complement[base] for base in sequence])

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    >>> reverse_and_complement('ATGC')
    'GCAT'
    """
    sequence = sequence.upper()
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    pass


if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
