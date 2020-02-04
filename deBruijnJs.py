import sys
import math
import json

# https://jgeisler0303.github.io/deBruijnDecode/#decoderTest

from logger import *
from deBruijnJsTools import *

class DeBruijn():
    def __init__(self, c, n):
        self.c = c
        self.n = n
        self.a = [0 for j in range(c*n)]
        self.s = []

        print_log(f"DeBruijn :: c: {c} n: {n} a: {join_list(self.a)} s: {join_list(self.s)}")
        
        self.generate(1, 1)

        print_log(f"DeBruijn :: c: {c} n: {n} a: {join_list(self.a)} s: {join_list(self.s)}")

    def generate(self, t, p):
        print_debug(f"DeBruijn :: generate :: t: {t} p: {p} c: {self.c} n: {self.n}")
        
        if t > self.n:
            if (self.n%p) == 0:
                for j in range(p):
                    self.s.append(self.a[j+1])
                    print_debug(f"DeBruijn :: generate :: t: {t} p: {p} c: {self.c} n: {self.n} j: 0 a: {join_list(self.a)} s: {join_list(self.s)}")

        else:
            self.a[t] = self.a[t-p]
            self.generate(t+1, p)

            for j in range(self.a[t-p]+1, self.c):
                print_debug(f"DeBruijn :: generate :: t: {t} p: {p} c: {self.c} n: {self.n} j: {j} a: {join_list(self.a)} s: {join_list(self.s)}")
                self.a[t] = j
                self.generate(t+1, t)

def genDecodableDeBruijn(T, K, L, c, n):
    print_debug(f"genDecodableDeBruijn :: c: {c} n: {n}\n\tT: {join_list(T)}\n\tK: {join_list(K)}\n\tL: {join_list(L)}")

    if n <= 2:
        db = DeBruijn(c, 2)
        t  = db.s
        # t_ = "".join([str(x) for x in t[0:2]])
        t_ = t + t[0:2]

        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} t: {join_list(t)} t_: {join_list(t_)}")
        
        for i in range(c ** 2):
            k   = t_[i:i+2]
            idx = wordIndex(k, c)
            print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} t: {join_list(t)} t_: {join_list(t_)} i: {i} k: {k} idx: {idx}")
            T[idx] = i

    else:
        n_      = n-1
        _,_,_,s = genDecodableDeBruijn(T, K, L, c, n_)
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} s    : {join_list(s)}")
        
        k = findOnes(s, n_)
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} k    : {k}")
        
        s_    = s[0:k] + s[k+1:]
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} k    : {k} s_   : {join_list(s_)} [{len(s_)}]")

        s_hat = operator_D_inv(s_, 0, c)
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} k    : {k} s_hat: {join_list(s_hat)}  [{len(s_hat)}]")

        p = (c-1) * ((c ** n_) - 1) + k
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} p    : {p}")
        
        e = s_hat[p]
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} e    : {e}")
        
        t = rho(c, p, e, s_hat)
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} t    : {join_list(t)} [{len(t)}]")

        lidx = n-3 
        tidx = ((c ** n_)-1)
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} lidx : {lidx}")
        print_debug(f"genDecodableDeBruijn :: c: {c} n: {n} tidx : {tidx}")
        L[lidx] = t[0:tidx]

    K[n-2] = findOnes(t, n)
    
    print_debug(f"genDecodableDeBruijn ::\n\tt: {t} ({len(t)})\n\tT: {T} ({len(T)})\n\tK: {K} ({len(K)})\n\tL: {L} ({len(L)})")

    return T, K, L, t

def createBruijnDecoder(T, K, L, w, c):
    print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c}")

    r = 2
    n = len(w)

    print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n}")

    if (n==r):
        idx = wordIndex(w, c)
        v   = T[idx]
        print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} idx: {idx}")
        return v
  
    v       = operator_D(w, c)
    i       = n-1-r
    k       = K[i]
    p       = (c-1) * ((c ** (n-1)) - 1) + k
    allOnes = True

    for i in range(n-1):
        if v[i] != 1:
            allOnes = False
    
    print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} i: {i} k: {k} p: {p} allOnes: {allOnes}")

    if allOnes:
        lv = L[n-1-r][k]
        e  = lv + ((c-1) ** 2)
        j  = p  + mod(w[0]-e, c)
        print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} lv: {lv} e: {e} j: {j} allOnes: true")

    else:
        f = createBruijnDecoder(T, K, L, v, c)

        if (f>k):
            f -= 1

        lv= L[n-1-r][f]
        e = lv
        j = f + (((c ** (n-1)) - 1) * mod(e-w[0], c))
        
        print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} lv: {lv} e: {e} j: {j} allOnes: false f: {f}")

        if (j<0) or (j>p-1):
            j = j+c

    print_debug(f"createBruijnDecoder :: w: {join_list(w)} c: {c} j: {j}")

    return j

def vocabToDeBruijn(vocab, n):
    print_debug(f"vocabToDeBruijn :: vocab: {vocab} n: {n}")

    c = len(vocab)
    T = [None] * (c ** 2)
    K = [None] * (n-1)
    L = [None] * (n-2)

    print_debug(f"vocabToDeBruijn :: vocab: {vocab} c: {c} n: {n} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)}")

    T, K, L, t = genDecodableDeBruijn(T, K, L, c, n)

    print_debug(f"vocabToDeBruijn :: vocab: {vocab} c: {c} n: {n} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)} t: {join_list(t)}")

    return T, K, L, t, c

def decodeDeBruijn(T, K, L, s, c, n):
    J = [None] * (c ** n)

    print_log(f"decodeDeBruijn ::\n\tT: {join_list(T)} ({len(T)})\n\tK: {join_list(K)} ({len(K)})\n\tL: {L}\n\tJ: {join_list(J)} ({len(J)})\n\ts: {join_list(s)} ({len(s)})\n\tc: {c}\n\tn: {n}")

    for i in range(c ** n):
        w  = [None] * n
        i_ = i
        print_debug(f"decodeDeBruijn :: i: {i} i_ {i_} w {w} ")
        
        for j in range(n-1, -1, -1):
            v = i_ % c
            print_debug(f"decodeDeBruijn :: i: {i} i_ {i_} w {w} j: {j} v: {v}")
            w[j] = v
            i_   = math.floor(i_ / c)

        j_ = findWord(s, w)
        d  = createBruijnDecoder(T, K, L, w, c)
        print_debug(f"decodeDeBruijn :: i: {i} i_ {i_} w {w} d: {d}")
        assert j_==d, f"j_ == d. j_: {j_} d: {d}"
        J[i] = d

    print_debug(f"decodeDeBruijn :: J: {J}")
    
    return J

def main(vocab, n):
    test()

    T, K, L, t, c = vocabToDeBruijn(vocab, n)
    J             = decodeDeBruijn(T, K, L, t, c, n)
    sstr          = "".join([vocab[l] for l in t])

    print_info(f" Vocabulary: {vocab} ({len(vocab)})")
    print_info(f" T: {T} ({len(T)})")
    print_info(f" K: {K} ({len(K)})")
    print_info(f" L: {L} ({len(L)})")
    for l in range(len(L)):
        print_info(f"  L[{l}]: {L[l]} ({len(L[l])})")
    print_info(f" t: {t} ({len(t)})")
    print_info(f" t: {sstr}")
    print_info(f" l: {c}")
    print_info(f" n: {n}")
    print_info(f" J: {J} ({len(J)})")

def test():
    vocab = "ACGT"
    n     = 3
    c_    = len(vocab)
    L_    = [[0,0,0,1,1,3,3,2,3,1,2,1,3,1,0]]
    K_    = [7,45]
    T_    = [0,2,4,15,1,7,9,6,3,8,12,11,5,10,13,14]
    t_    = [0,0,0,1,1,3,3,2,3,1,2,1,3,1,0,3,3,3,0,0,2,2,1,2,0,1,0,2,0,3,2,2,2,3,3,1,1,0,1,3,0,3,1,3,2,1,1,1,2,2,0,0,3,0,1,2,3,0,2,3,2,0,2,1]

    T, K, L, t, c = vocabToDeBruijn(vocab, n)

    assert c_ == c, f"{c_} != {c}"
    assert t_ == t, f"{L_} != {L}"
    assert T_ == T, f"{T_} != {T}"
    assert K_ == K, f"{K_} != {K}"
    assert L_ == L, f"{L_} != {L}"

    sstr = "".join([vocab[l] for l in t])

    print_info(f"test :: T: {T} ({len(T)})")
    print_info(f"test :: K: {K} ({len(K)})")
    print_info(f"test :: L: {L} ({len(L)})")
    for l in range(len(L)):
        print_info(f"test ::  L[{l}]: {L[l]} ({len(L[l])})")
    print_info(f"test :: t: {t} ({len(t)})")
    print_info(f"test :: t: {sstr}")
    print_info(f"test :: l: {c}")
    print_info(f"test :: n: {n}")

    J = decodeDeBruijn(T, K, L, t, c, n)

    print_info(f"test :: J {J} ({len(J)})")

    sstr_c = sstr + sstr[:n-1] # circular
    keys = {}
    for kpos in range(len(sstr)):
        kmer = sstr_c[kpos:kpos+n]
        karr = [vocab.index(p) for p in kmer]
        d    = createBruijnDecoder(T, K, L, karr, c)
        print(kpos, kmer, karr, d)
        keys[kmer] = d

    #TODO: generalize
    for st in range(len(vocab)):
        for nd in range(len(vocab)):
            for rd in range(len(vocab)):
                kmer = vocab[st] + vocab[nd] + vocab[rd]
                kid  = 4**2*st + 4**1*nd + 4**0*rd
                karr = [st, nd, rd]
                d    = createBruijnDecoder(T, K, L, karr, c)
                assert keys[kmer] == d
                print(kmer, kid, d)

    print_info("all tests passed")


if __name__ == '__main__':
    vocab     =     sys.argv[1]
    kmer_size = int(sys.argv[2])
    main(vocab, kmer_size)
