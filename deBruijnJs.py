import sys
import os
import math
import json

# https://jgeisler0303.github.io/deBruijnDecode/#decoderTest

"""
4_7
(16384)
real    12m18.553s
user    11m49.953s
sys     00m21.078s

real    4m55.889s
user    4m53.109s
sys     0m02.125s

real    0m0.164s
user    0m0.078s
sys     0m0.063s

"""

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
        # print_debug(f"DeBruijn :: generate :: t: {t} p: {p} c: {self.vocab_size} kmer_size: {self.kmer_size}")
        
        if t > self.kmer_size:
            if (self.kmer_size % p) == 0:
                for frame in range(p):
                    self.dbsequence.append(self.vocab_matrix[frame+1])
                    # print_debug(f"DeBruijn :: generate :: t: {t} p: {p} vocab_size: {self.vocab_size} kmer_size: {self.kmer_size} frame: 0 vocab_matrix: {join_list(self.vocab_matrix)} dbsequence: {join_list(self.dbsequence)}")

        else:
            self.vocab_matrix[t] = self.vocab_matrix[t-p]
            self.generate(t+1, p)

            for frame in range(self.vocab_matrix[t-p]+1, self.vocab_size):
                # print_debug(f"DeBruijn :: generate :: t: {t} p: {p} vocab_size: {self.vocab_size} kmer_size: {self.kmer_size} frame: {frame} vocab_matrix: {join_list(self.vocab_matrix)} dbsequence: {join_list(self.dbsequence)}")
                self.vocab_matrix[t] = frame
                self.generate(t+1, t)

def genDecodableDeBruijn(T, K, L, vocab_size, kmer_size):
    # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size}\n\tT: {join_list(T)}\n\tK: {join_list(K)}\n\tL: {join_list(L)}")

    if kmer_size <= 2:
        db          = DeBruijn(vocab_size, 2)
        dbsequence  = db.dbsequence
        dbsequence_ = dbsequence + dbsequence[0:2]

        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} db.dbsequence: {join_list(db.dbsequence)} dbsequence_: {join_list(dbsequence_)}")
        
        for i in range(vocab_size ** 2):
            word         = dbsequence_[i:i+2]
            decoded_char = wordIndex(word, vocab_size)
            # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} db.dbsequence: {join_list(db.dbsequence)} dbsequence_: {join_list(dbsequence_)} i: {i} word: {word} decoded_char: {decoded_char}")
            T[decoded_char] = i

    else:
        kmer_size_        = kmer_size-1
        _,_,_,dbsequence_ = genDecodableDeBruijn(T, K, L, vocab_size, kmer_size_)
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbsequence_: {join_list(dbsequence_)}")
        
        ones_pos     = findOnes(dbsequence_, kmer_size_)
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} ones_pos: {ones_pos}")
        
        dbsequence__ = dbsequence_[0:ones_pos] + dbsequence_[ones_pos+1:]
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} ones_pos: {ones_pos} dbsequence__: {join_list(dbsequence__)} [{len(dbsequence__)}]")

        dbsequence__hat = operator_D_inv(dbsequence__, 0, vocab_size)
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} ones_pos: {ones_pos} dbsequence__hat: {join_list(dbsequence__hat)}  [{len(dbsequence__hat)}]")

        dbs_pos = (vocab_size-1) * ((vocab_size ** kmer_size_) - 1) + ones_pos
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbs_pos: {dbs_pos}")
        
        dbs_val = dbsequence__hat[dbs_pos]
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbs_val: {dbs_val}")
        
        dbsequence = rho(vocab_size, dbs_pos, dbs_val, dbsequence__hat)
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} dbsequence: {join_list(dbsequence)} [{len(dbsequence)}]")

        lidx = kmer_size - 3 
        tidx = ((vocab_size ** kmer_size_)-1)
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} lidx : {lidx}")
        # print_debug(f"genDecodableDeBruijn :: vocab_size: {vocab_size} kmer_size: {kmer_size} tidx : {tidx}")
        
        L[lidx] = dbsequence[0:tidx]

    K[kmer_size-2] = findOnes(dbsequence, kmer_size)
    
    # print_debug(f"genDecodableDeBruijn ::\n\tdbsequence: {dbsequence} ({len(dbsequence)})\n\tT: {T} ({len(T)})\n\tK: {K} ({len(K)})\n\tL: {L} ({len(L)})")

    return T, K, L, dbsequence

def vocabToDeBruijn(vocab, kmer_size):
    # print_debug(f"vocabToDeBruijn :: vocab: {vocab} kmer_size: {kmer_size}")

    vocab_size = len(vocab)

    db_file    = f"{vocab_size}_{kmer_size}.vars.json"

    if os.path.exists(db_file):
        print(f"reading vars")
        with open(db_file, 'r') as fhd:
            print(f" reading vars from {db_file}")
            T, K, L, len_t, len_k, len_l, len_dbsequence, vocab_size_, kmer_size_, dbsequence = json.load(fhd)
            assert vocab_size == vocab_size_
            assert kmer_size  == kmer_size_
            assert len_t == len(T)
            assert len_k == len(K)
            assert len_l == len(L)
            assert len_dbsequence == len(dbsequence)
            return T, K, L, vocab_size, dbsequence

    else:
        print(f"generating vars")

        T = [None] * (vocab_size ** 2)
        K = [None] * (kmer_size-1)
        L = [None] * (kmer_size-2)

        # print_debug(f"vocabToDeBruijn :: vocab: {vocab} vocab_size: {vocab_size} kmer_size: {kmer_size} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)}")

        T, K, L, dbsequence = genDecodableDeBruijn(T, K, L, vocab_size, kmer_size)

        # print_debug(f"vocabToDeBruijn :: vocab: {vocab} vocab_size: {vocab_size} kmer_size: {kmer_size} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)} dbsequence: {join_list(dbsequence)}")

        with open(db_file, 'w') as fhd:
            print(f" saving vars to {db_file}")
            json.dump([T, K, L, len(T), len(K), len(L), len(dbsequence), vocab_size, kmer_size, dbsequence], fhd, indent=1)

        return T, K, L, vocab_size, dbsequence





def decodeDeBruijnWord(T, K, L, vocab_size, word):
    # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size}")

    word_min_size = 2
    kmer_size     = len(word)

    # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size}")

    if (kmer_size == word_min_size):
        decoded_char  = wordIndex(word, vocab_size)
        decoded_char_ = T[decoded_char]
        # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} decoded_char_: {decoded_char_} decoded_char: {decoded_char}")
        return decoded_char_
  
    word_D   = operator_D(word, vocab_size)
    word_pos = kmer_size-word_min_size-1
    kv       = K[word_pos]
    pv       = (vocab_size-1) * ((vocab_size ** (kmer_size-1)) - 1) + kv
    allOnes  = True

    for kp in range(kmer_size-1):
        if word_D[kp] != 1:
            allOnes = False
    
    # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} word_D: {word_D} word_pos: {word_pos} kv: {kv} pv: {pv} allOnes: {allOnes}")

    if allOnes:
        lv           = L[kmer_size-word_min_size-1][kv]
        e            = lv + ((vocab_size-1) ** 2)
        decoded_char = pv + mod(word[0]-e, vocab_size)
        # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} word_D: {word_D} lv: {lv} e: {e} decoded_char: {decoded_char} allOnes: true")

    else:
        decoded_char_ = decodeDeBruijnWord(T, K, L, vocab_size, word_D)

        if decoded_char_ > kv:
            decoded_char_ -= 1

        lv           = L[kmer_size-word_min_size-1][decoded_char_]
        e            = lv
        decoded_char = decoded_char_ + (((vocab_size ** (kmer_size-1)) - 1) * mod(e-word[0], vocab_size))
        
        # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} word_min_size: {word_min_size} kmer_size: {kmer_size} word_D: {word_D} lv: {lv} e: {e} decoded_char: {decoded_char} allOnes: false decoded_char_: {decoded_char_}")

        if (decoded_char<0) or (decoded_char>pv-1):
            decoded_char = decoded_char + vocab_size

    # print_debug(f"decodeDeBruijnWord :: word: {join_list(word)} vocab_size: {vocab_size} decoded_char: {decoded_char}")

    return decoded_char

def genDeBruijnDecodeMatrix(dbsequence, vocab_size, kmer_size):
    db_file       = f"{vocab_size}_{kmer_size}.matrix.json"
    matrix_size   = (vocab_size ** kmer_size)

    if os.path.exists(db_file):
        print(f"reading matrix {db_file}")

        with open(db_file, 'r') as fhd:
            print(f" reading matrix from {db_file}")
            decode_matrix, matrix_size_, vocab_size_, kmer_size_ = json.load(fhd)
            assert vocab_size  == vocab_size_
            assert kmer_size   == kmer_size_
            assert matrix_size == matrix_size_
            assert matrix_size == len(decode_matrix)
            return decode_matrix

    else:
        print(f"generating matrix")

        decode_matrix = [None] * matrix_size

        print_log(f"genDeBruijnDecodeMatrix ::\n\tJ: {join_list(decode_matrix)} ({len(decode_matrix)})\n\tdbsequence: {join_list(dbsequence)} ({len(dbsequence)})\n\tvocab_size: {vocab_size}\n\tkmer_size: {kmer_size}")

        for matrix_pos in range(matrix_size):
            word        = [None] * kmer_size
            matrix_pos_ = matrix_pos
            # print_debug(f"genDeBruijnDecodeMatrix :: matrix_pos: {matrix_pos} matrix_pos_ {matrix_pos_} word {word} ")
            
            for kmer_pos in range(kmer_size-1, -1, -1):
                pos_frame      = matrix_pos_ % vocab_size
                word[kmer_pos] = pos_frame
                matrix_pos_    = math.floor(matrix_pos_ / vocab_size)
                # print_debug(f"genDeBruijnDecodeMatrix :: matrix_pos: {matrix_pos} matrix_pos_ {matrix_pos_} word {word} kmer_pos: {kmer_pos} pos_frame: {pos_frame}")

            decoded_char_ = findWord(dbsequence, word)
            # decoded_char  = decodeDeBruijnWord(T, K, L, vocab_size, word)
            # # print_debug(f"genDeBruijnDecodeMatrix :: matrix_pos: {matrix_pos} matrix_pos_ {matrix_pos_} word {word} decoded_char: {decoded_char}")
            # assert decoded_char_ == decoded_char, f"decoded_char_ == decoded_char. decoded_char_: {decoded_char_} decoded_char: {decoded_char}"
            decode_matrix[matrix_pos] = decoded_char_

        with open(db_file, 'w') as fhd:
            print(f" saving matrix to {db_file}")
            json.dump([decode_matrix, matrix_size, vocab_size, kmer_size], fhd, indent=0)

    # print_debug(f"genDeBruijnDecodeMatrix :: decode_matrix: {decode_matrix}")
    
    return decode_matrix


def main(vocab, kmer_size):
    T, K, L, vocab_size, dbsequence = vocabToDeBruijn(vocab, kmer_size)
    dbsequence_str                  = "".join([vocab[l] for l in dbsequence])

    seqfile = f"{vocab}_{kmer_size}.seq"
    with open(seqfile, 'w') as fhd:
        print(f" saving sequence to {seqfile}")
        fhd.write(dbsequence_str)

    T_ = str(T) if len(T) < 45 else str(T[:20]).strip("]") + " ... " + str(T[-20:]).strip("[")
    K_ = str(K) if len(K) < 45 else str(K[:20]).strip("]") + " ... " + str(K[-20:]).strip("[")
    L_ = [
        (str(l) if len(l) < 45 else str(l[:20]).strip("]") + " ... " + str(l[-20:]).strip("["))
        for l in L
        ]
    dbsequence_     = str(dbsequence)     if len(dbsequence)     < 45 else str(dbsequence[:20])    .strip("]") + " ... " + str(dbsequence[-20:])    .strip("[")
    dbsequence_str_ = str(dbsequence_str) if len(dbsequence_str) < 45 else str(dbsequence_str[:20]).strip("]") + " ... " + str(dbsequence_str[-20:]).strip("[")

    print_info(f" vocabulary    : {vocab} ({len(vocab)})")
    print_info(f" T             : {T_} ({len(T)})")
    print_info(f" K             : {K_} ({len(K)})")
    print_info(f" L             : ({len(L)})")
    for l in range(len(L)):
        print_info(f"   L[{l}]        : {L_[l]} ({len(L[l])})")
    print_info(f" dbsequence    : {dbsequence_} ({len(dbsequence)})")
    print_info(f" dbsequence_str: {dbsequence_str_}")
    print_info(f" vocab_size    : {vocab_size}")
    print_info(f" kmer_size     : {kmer_size}")

    decode_matrix  = genDeBruijnDecodeMatrix(dbsequence, vocab_size, kmer_size)
    decode_matrix_ = str(decode_matrix) if len(decode_matrix) < 25 else str(decode_matrix[:10]).strip("]") + " ... " + str(decode_matrix[-10:]).strip("[")

    print_info(f" decode_matrix : {decode_matrix_} ({len(decode_matrix)})")

def test():
    vocab          = "ACGT"
    kmer_size      = 3
    vocab_size_    = len(vocab)
    T_             = [0,2,4,15,1,7,9,6,3,8,12,11,5,10,13,14]
    K_             = [7,45]
    L_             = [[0,0,0,1,1,3,3,2,3,1,2,1,3,1,0]]
    dbsequence_    = [0,0,0,1,1,3,3,2,3,1,2,1,3,1,0,3,3,3,0,0,2,2,1,2,0,1,0,2,0,3,2,2,2,3,3,1,1,0,1,3,0,3,1,3,2,1,1,1,2,2,0,0,3,0,1,2,3,0,2,3,2,0,2,1]
    decode_matrix_ = [0,1,18,50,24,2,53,37,26,61,19,57,51,40,28,14,63,36,25,13,35,45,46,3,22,9,47,54,38,11,42,4,49,23,60,27,62,44,21,10,48,20,30,31,55,7,58,32,17,52,56,39,12,34,8,41,59,43,29,6,16,33,5,15]

    T, K, L, vocab_size, dbsequence = vocabToDeBruijn(vocab, kmer_size)

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

    decode_matrix = genDeBruijnDecodeMatrix(dbsequence, vocab_size, kmer_size)
    assert decode_matrix_ == decode_matrix

    print_info(f"test :: decode_matrix : {decode_matrix} ({len(decode_matrix)})")

    dbsequence_str_circ = dbsequence_str + dbsequence_str[:kmer_size-1] # circular
    keys   = {}
    for kpos in range(len(dbsequence_str)):
        kmer = dbsequence_str_circ[kpos:kpos+kmer_size]
        karr = [vocab.index(p) for p in kmer]
        dwor = decodeDeBruijnWord(T, K, L, vocab_size, karr)
        print(kpos, kmer, karr, dwor)
        keys[kmer] = dwor

    #TODO: generalize
    for st in range(len(vocab)):
        for nd in range(len(vocab)):
            for rd in range(len(vocab)):
                kmer = vocab[st] + vocab[nd] + vocab[rd]
                kid  = 4**2*st + 4**1*nd + 4**0*rd
                karr = [st, nd, rd]
                dwor = decodeDeBruijnWord(T, K, L, vocab_size, karr)
                assert keys[kmer] == dwor
                print(kmer, kid, dwor)

    print_info("all tests passed")


if __name__ == '__main__':
    test()

    vocab     =     sys.argv[1]
    kmer_size = int(sys.argv[2])
    main(vocab, kmer_size)
