import sys
import math
import json

# https://jgeisler0303.github.io/deBruijnDecode/#decoderTest

from logger import *
from deBruijnJsTools import *

class DeBruijn():
    def __init__(self, vocab_size, kmer_size):
        self.vocab_size   = vocab_size
        self.kmer_size    = kmer_size
        self.vocab_matrix = [0 for j in range(vocab_size*kmer_size)]
        self.dbsequence   = []

        print_log(f"DeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} vocab_matrix: {join_list(self.vocab_matrix)} dbsequence: {join_list(self.dbsequence)}")
        
        self.generate(1, 1)

        print_log(f"DeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} vocab_matrix: {join_list(self.vocab_matrix)} dbsequence: {join_list(self.dbsequence)}")

    def generate(self, t, p):
        print_debug(f"DeBruijn :: generate :: t: {t} p: {p} c: {self.vocab_size} kmer_size: {self.kmer_size}")
        
        if t > self.kmer_size:
            if (self.kmer_size % p) == 0:
                for frame in range(p):
                    self.dbsequence.append(self.vocab_matrix[j+1])
                    print_debug(f"DeBruijn :: generate :: t: {t} p: {p} vocab_size: {self.vocab_size} kmer_size: {self.kmer_size} frame: 0 vocab_matrix: {join_list(self.vocab_matrix)} dbsequence: {join_list(self.dbsequence)}")

        else:
            self.vocab_matrix[t] = self.vocab_matrix[t-p]
            self.generate(t+1, p)

            for frame in range(self.vocab_matrix[t-p]+1, self.vocab_size):
                print_debug(f"DeBruijn :: generate :: t: {t} p: {p} vocab_size: {self.vocab_size} kmer_size: {self.kmer_size} frame: {frame} vocab_matrix: {join_list(self.vocab_matrix)} dbsequence: {join_list(self.dbsequence)}")
                self.vocab_matrix[t] = frame
                self.generate(t+1, t)

def genDecodableDeBruijn(T, K, L, vocab_size, kmer_size):
    print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size}\n\tT: {join_list(T)}\n\tK: {join_list(K)}\n\tL: {join_list(L)}")

    if kmer_size <= 2:
        db          = DeBruijn(vocab_size, 2)
        dbsequence  = db.dbsequence
        dbsequence_ = dbsequence + dbsequence[0:2]

        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbsequence: {join_list(dbsequence)} dbsequence_: {join_list(dbsequence_)}")
        
        for i in range(vocab_size ** 2):
            word         = dbsequence_[i:i+2]
            decoded_char = wordIndex(word, vocab_size)
            print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbsequence: {join_list(dbsequence)} dbsequence_: {join_list(dbsequence_)} i: {i} word: {word} decoded_char: {decoded_char}")
            T[decoded_char] = i

    else:
        kmer_size_ = kmer_size-1
        _,_,_,dbsequence_    = genDecodableDeBruijn(T, K, L, vocab_size, kmer_size_)
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbsequence_: {join_list(dbsequence_)}")
        
        ones_pos     = findOnes(dbsequence_, kmer_size_)
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} ones_pos: {ones_pos}")
        
        dbsequence__ = dbsequence_[0:ones_pos] + dbsequence_[ones_pos+1:]
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} ones_pos: {ones_pos} dbsequence__: {join_list(dbsequence__)} [{len(dbsequence__)}]")

        dbsequence__hat = operator_D_inv(dbsequence__, 0, vocab_size)
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} ones_pos: {ones_pos} dbsequence__hat: {join_list(dbsequence__hat)}  [{len(dbsequence__hat)}]")

        dbs_pos = (vocab_size-1) * ((vocab_size ** kmer_size_) - 1) + ones_pos
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbs_pos: {dbs_pos}")
        
        dbs_val = dbsequence__hat[dbs_pos]
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbs_val: {dbs_val}")
        
        dbsequence = rho(vocab_size, dbs_pos, dbs_val, dbsequence__hat)
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbsequence: {join_list(dbsequence)} [{len(dbsequence)}]")

        lidx = kmer_size - 3 
        tidx = ((vocab_size ** kmer_size_)-1)
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} lidx : {lidx}")
        print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} tidx : {tidx}")
        
        L[lidx] = dbsequence[0:tidx]

    K[kmer_size-2] = findOnes(dbsequence, kmer_size)
    
    print_debug(f"genDecodableDeBruijn ::\n\tdbsequence: {dbsequence} ({len(dbsequence)})\n\tT: {T} ({len(T)})\n\tK: {K} ({len(K)})\n\tL: {L} ({len(L)})")

    return T, K, L, dbsequence

def vocabToDeBruijn(vocab, kmer_size):
    print_debug(f"vocabToDeBruijn :: vocab: {vocab} kmer_size: {kmer_size}")

    vocab_size = len(vocab)
    T          = [None] * (vocab_size ** 2)
    K          = [None] * (kmer_size-1)
    L          = [None] * (kmer_size-2)

    print_debug(f"vocabToDeBruijn :: vocab: {vocab} vocab_size: {vocab_size} kmer_size: {kmer_size} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)}")

    T, K, L, dbsequence = genDecodableDeBruijn(T, K, L, vocab_size, kmer_size)

    print_debug(f"vocabToDeBruijn :: vocab: {vocab} vocab_size: {vocab_size} kmer_size: {kmer_size} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)} dbsequence: {join_list(dbsequence)}")

    return T, K, L, dbsequence, vocab_size


def decodeDeBruijnWord(T, K, L, vocab_size, word):
    print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size}")

    word_min_size = 2
    kmer_size     = len(word)

    print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size}")

    if (kmer_size == word_min_size):
        decoded_char  = wordIndex(word, vocab_size)
        decoded_char_ = T[decoded_char]
        print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} decoded_char_: {decoded_char_} decoded_char: {decoded_char}")
        return decoded_char_
  
    word_D  = operator_D(word, vocab_size)
    i       = kmer_size-word_min_size-1
    k       = K[i]
    p       = (vocab_size-1) * ((vocab_size ** (kmer_size-1)) - 1) + k
    allOnes = True

    for i in range(kmer_size-1):
        if word_D[i] != 1:
            allOnes = False
    
    print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} word_D: {word_D} i: {i} k: {k} p: {p} allOnes: {allOnes}")

    if allOnes:
        lv           = L[kmer_size-word_min_size-1][k]
        e            = lv + ((vocab_size-1) ** 2)
        decoded_char = p  + mod(word[0]-e, vocab_size)
        print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} word_D: {word_D} lv: {lv} e: {e} decoded_char: {decoded_char} allOnes: true")

    else:
        decoded_char_ = decodeDeBruijnWord(T, K, L, vocab_size, word_D)

        if (decoded_char_>k):
            decoded_char_ -= 1

        lv           = L[kmer_size-word_min_size-1][decoded_char_]
        e            = lv
        decoded_char = decoded_char_ + (((vocab_size ** (kmer_size-1)) - 1) * mod(e-word[0], vocab_size))
        
        print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} word_D: {word_D} lv: {lv} e: {e} decoded_char: {decoded_char} allOnes: false decoded_char_: {decoded_char_}")

        if (decoded_char<0) or (decoded_char>p-1):
            decoded_char = decoded_char + vocab_size

    print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} decoded_char: {decoded_char}")

    return decoded_char

def genDeBruijnDecodeMatrix(T, K, L, dbsequence, vocab_size, kmer_size):
    matrix_size   = (vocab_size ** kmer_size)
    decode_matrix = [None] * matrix_size

    print_log(f"genDeBruijnDecodeMatrix ::\n\tT: {join_list(T)} ({len(T)})\n\tK: {join_list(K)} ({len(K)})\n\tL: {L}\n\tJ: {join_list(decode_matrix)} ({len(decode_matrix)})\n\tdbsequence: {join_list(dbsequence)} ({len(dbsequence)})\n\tvocab_size: {vocab_size}\n\tkmer_size: {kmer_size}")

    for matrix_pos in range(matrix_size):
        word        = [None] * kmer_size
        matrix_pos_ = matrix_pos
        print_debug(f"genDeBruijnDecodeMatrix :: matrix_pos: {matrix_pos} matrix_pos_ {matrix_pos_} word {word} ")
        
        for kmer_pos in range(kmer_size-1, -1, -1):
            pos_frame      = matrix_pos_ % vocab_size
            word[kmer_pos] = pos_frame
            matrix_pos_    = math.floor(matrix_pos_ / vocab_size)
            print_debug(f"genDeBruijnDecodeMatrix :: matrix_pos: {matrix_pos} matrix_pos_ {matrix_pos_} word {word} kmer_pos: {kmer_pos} pos_frame: {pos_frame}")

        decoded_char_ = findWord(dbsequence, word)
        decoded_char  = decodeDeBruijnWord(T, K, L, vocab_size, word)
        print_debug(f"genDeBruijnDecodeMatrix :: matrix_pos: {matrix_pos} matrix_pos_ {matrix_pos_} word {word} decoded_char: {decoded_char}")
        assert decoded_char_==decoded_char, f"decoded_char_ == decoded_char. decoded_char_: {decoded_char_} decoded_char: {decoded_char}"
        decode_matrix[matrix_pos] = decoded_char

    print_debug(f"genDeBruijnDecodeMatrix :: decode_matrix: {decode_matrix}")
    
    return decode_matrix

def main(vocab, kmer_size):
    T, K, L, dbsequence, vocab_size = vocabToDeBruijn(vocab, kmer_size)
    decode_matrix                   = genDeBruijnDecodeMatrix(T, K, L, dbsequence, vocab_size, kmer_size)
    dbsequence_str                  = "".join([vocab[l] for l in dbsequence])

    print_info(f" vocabulary    : {vocab} ({len(vocab)})")
    print_info(f" T             : {T} ({len(T)})")
    print_info(f" K             : {K} ({len(K)})")
    print_info(f" L             : {L} ({len(L)})")
    for l in range(len(L)):
        print_info(f"   L[{l}]        : {L[l]} ({len(L[l])})")
    print_info(f" dbsequence    : {dbsequence} ({len(dbsequence)})")
    print_info(f" dbsequence_str: {dbsequence_str}")
    print_info(f" vocab_size    : {vocab_size}")
    print_info(f" kmer_size     : {kmer_size}")
    print_info(f" decode_matrix : {decode_matrix} ({len(decode_matrix)})")

def test():
    vocab       = "ACGT"
    kmer_size   = 3
    vocab_size_ = len(vocab)
    L_          = [[0,0,0,1,1,3,3,2,3,1,2,1,3,1,0]]
    K_          = [7,45]
    T_          = [0,2,4,15,1,7,9,6,3,8,12,11,5,10,13,14]
    dbsequence_ = [0,0,0,1,1,3,3,2,3,1,2,1,3,1,0,3,3,3,0,0,2,2,1,2,0,1,0,2,0,3,2,2,2,3,3,1,1,0,1,3,0,3,1,3,2,1,1,1,2,2,0,0,3,0,1,2,3,0,2,3,2,0,2,1]

    T, K, L, dbsequence, vocab_size = vocabToDeBruijn(vocab, kmer_size)

    assert vocab_size_ == vocab_size, f"{vocab_size_} != {vocab_size}"
    assert T_ == T, f"{T_} != {T}"
    assert K_ == K, f"{K_} != {K}"
    assert L_ == L, f"{L_} != {L}"
    assert dbsequence_ == dbsequence, f"{L_} != {L}"

    dbsequence_str = "".join([vocab[l] for l in dbsequence])

    print_info(f"test :: T             : {T} ({len(T)})")
    print_info(f"test :: K             : {K} ({len(K)})")
    print_info(f"test :: L             : {L} ({len(L)})")
    for l in range(len(L)):
        print_info(f"test ::  L[{l}]         : {L[l]} ({len(L[l])})")
    print_info(f"test :: dbsequence    : {dbsequence} ({len(dbsequence)})")
    print_info(f"test :: dbsequence_str: {dbsequence_str}")
    print_info(f"test :: vocab_size    : {vocab_size}")
    print_info(f"test :: kmer_size     : {kmer_size}")

    decode_matrix = genDeBruijnDecodeMatrix(T, K, L, dbsequence, vocab_size, kmer_size)

    print_info(f"test :: decode_matrix : {decode_matrix} ({len(decode_matrix)})")

    dbsequence_str_circ = dbsequence_str + dbsequence_str[:kmer_size-1] # circular
    keys   = {}
    for kpos in range(len(dbsequence_str)):
        kmer = dbsequence_str_circ[kpos:kpos+kmer_size]
        karr = [vocab.index(p) for p in kmer]
        d    = decodeDeBruijnWord(T, K, L, vocab_size, karr)
        print(kpos, kmer, karr, d)
        keys[kmer] = d

    #TODO: generalize
    for st in range(len(vocab)):
        for nd in range(len(vocab)):
            for rd in range(len(vocab)):
                kmer = vocab[st] + vocab[nd] + vocab[rd]
                kid  = 4**2*st + 4**1*nd + 4**0*rd
                karr = [st, nd, rd]
                d    = decodeDeBruijnWord(T, K, L, vocab_size, karr)
                assert keys[kmer] == d
                print(kmer, kid, d)

    print_info("all tests passed")


if __name__ == '__main__':
    test()

    vocab     =     sys.argv[1]
    kmer_size = int(sys.argv[2])
    main(vocab, kmer_size)
