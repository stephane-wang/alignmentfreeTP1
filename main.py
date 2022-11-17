from loading import load_directory
from kmers import stream_kmers, kmer2str
import collections


def similarity(A, inter, B):
    
    return inter / (A + inter) , inter / (B + inter)


def jaccard(A, inter, B):
    
    return inter / (A + inter + B)


def listtostr(listseq):
    newstr = ""
    for i in range(len(listseq)):
        newstr+=listseq[i]
    
    return newstr


def intersect(Kmer1,Kmer2):
    A = []
    B = []
    K1 = collections.Counter(Kmer1)
    K2 = collections.Counter(Kmer2)
    inter = K1 & K2
    A = K1 - inter
    B = K2 - inter
    
    return A, inter, B


if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    kmer_dict = {}
    
    for i in range(len(files)):
        text = listtostr(files[filenames[i]])
        #print(text,stream_kmers(text,k))
        kmer_dict[filenames[i]] = stream_kmers(text,k)
        kmer_dict[filenames[i]].sort()
    
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            #print(kmer_dict[filenames[i]],"\n",kmer_dict[filenames[j]])
            A, inter, B = intersect(kmer_dict[filenames[i]], kmer_dict[filenames[j]])
            #print(A,"\n",inter,"\n",B)
            lA = sum(A.values())
            lB = sum(B.values())
            linter = sum(inter.values())
            print(filenames[i], filenames[j], jaccard(lA, linter, lB), similarity(lA, linter, lB))
