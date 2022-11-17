alphabet = {"A" : 0, "C" : 1, "G" : 2, "T" : 3}
complement = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}

def kmer2str(val, k):
    """ Transform a kmer integer into its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'G', 'T']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def stream_kmers(text, k):
    alphabet = {"A" : 0, "C" : 1, "G" : 2, "T" : 3}
    complement = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}
    list_kmer = []
    kmer = 0
    rkmer = 0
    alph = ['C','G','T']
    
    for i in range(k-1):
        #print(i)
        nucl = text[i]
        if nucl not in alph:
            nucl = 'A'
        kmer = kmer << 2
        kmer += alphabet[nucl]
        
        rkmer = rkmer >> 2
        rkmer += (alphabet[complement[nucl]] << (k - 1) * 2)
    
    mask = (1 << (k - 1) * 2) - 1
    for nt in text[k-1:]:
        if nt not in alph:
            nt = 'A'
        kmer = kmer & mask
        kmer = kmer << 2
        kmer = kmer + alphabet[nt]
        
        rkmer = rkmer >> 2
        rkmer = rkmer + (alphabet[complement[nt]] << (k - 1) * 2)
        
        list_kmer.append(min(kmer,rkmer))
    
    return list_kmer
