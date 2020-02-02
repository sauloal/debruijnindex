import sys
import math
# https://jgeisler0303.github.io/deBruijnDecode/#decoderTest

from logger import *

def join_list(lst, sep=" "):
    return sep.join(map(str,lst))

def mod(n, m):
    n = n % m
    
    if (n<0):
        n += m
    
    print_debug(f"mod :: n: {n} m: {m} r: {n}")

    return n

def gcd(a,b):
    ao = a
    bo = b

    if (a < 0): a = -a
    if (b < 0): b = -b
    if (b > a):
        temp = a
        a = b
        b = temp

    while True:
        a %= b
        if (a == 0):
            print_debug(f"gcd :: a: {ao} b: {bo} r: {b}")
            return b
        b %= a
        if (b == 0):
            print_debug(f"gcd :: a: {ao} b: {bo} r: {a}")
            return a

def diff(a):
    s = [None]*(len(a)-1)
    
    for i in range(1, len(a)):
        s[i-1] = a[i] - a[i-1]

    print_debug(f"diff :: a: {a} s: {s}")

    return s

def cumsum(a):
    ao = list(a)

    for i in range(1, len(ao)):
        a[i] = a[i-1] + a[i]

    print_debug(f"cumsum :: ao: {ao} a: {a}")

    return a

def wordIndex(w, c):
    print_debug(f"wordIndex :: w: {w} c: {c}")

    idx = 0
    b   = 1
    # w.reverse();

    for i in range(len(w)):
        print_debug(f"wordIndex :: w: {w} c: {c} i: {i} b: {b} w[i]: {w[i]} idx: {idx} idx+: {b * w[i]} idx: {idx + b * w[i]} b*c: {b*c}")
        idx += b * w[i]
        b   *= c

    return idx

def wordIndexOrig(w, c):
    print_debug(f"wordIndex :: w: {w} c: {c}")
    
    idx = 0
    b   = 1
    # w.reverse();

    for i in range(len(w)):
        idx += b * int(w[i])
        b   *= c

    return idx

def findWord(s,w):
    print_debug(f"findWord :: w: {w}\n\ts: {s}")

    ls = len(s)
    lw = len(w)
    s_ = s+s[0:lw]
    
    print_debug(f"findWord :: w: {w} ls: {ls} lw: {lw} s_: {s_} ({len(s_)})")

    k_ = -1
    for k in range(ls):
        print_debug(f"findWord :: w: {w} ls: {ls} lw: {lw} k: {k} k+lw: {k+lw} s_[k:k+lw]: {s_[k:k+lw]} s_[k:k+lw] == w: {s_[k:k+lw] == w}")
        if s_[k:k+lw] == w:
            if k > len(s):
                k_ = -2
            k_ = k
    
    print_debug(f"findWord :: s: {s} ls: {ls} lw: {lw} k_: {k_}")

    return k_

def findWordOrig(s, w):
    print_debug(f"findWord :: s: {s} w: {w}")
    
    n   = len(w)
    s_  = "".join([str(x) for x in s+s[0:n]])
    wi_ = 0
    k_  = 0
    
    print_debug(f"findWord :: s: {s} w: {w} n: {n} s_: {s_} wi_: {wi_} k_: {k_}")

    # for(k= 0; k<s.length && wi!=n; k++):
    for k in range(len(w)):
        print_debug(f"findWord :: s: {s} w: {w} n: {n} s_: {s_} wi_: {wi_} k_: {k_} k: {k}")
        print_debug(f"findWord :: k: {k}")
        k_ = k
        if wi_ == n: break
        # for(wi= 0; wi<n && s_[k+wi]==w[wi]; wi++)
        for wi in range(n):
            print_debug(f"findWord :: s: {s} w: {w} n: {n} s_: {s_} wi_: {wi_} k_: {k_} k: {k} wi: {wi}")
            wi_ = wi
            if s_[k+wi] == w[wi]: break
        
    k_ -= 1

    if k_ >= len(s):
        k_ = -1

    print_debug(f"findWord :: s: {s} w: {w} n: {n} s_: {s_} wi_: {wi_} k_: {k_}")

    return k_

def findOnes(s, n):
    w = [None]*n

    for i in range(n):
        w[i] = 1
    
    return findWord(s, w)

def operator_D_inv(s, b, c):
    print_debug(f"operator_D_inv :: s: {join_list(s)} b: {b} c: {c}")

    w = sum(s) % c

    print_debug(f"operator_D_inv :: s: {join_list(s)} b: {b} c: {c} w: {w}")

    if w != 0:
        s_ = list(s)
        mr = math.floor((c/gcd(c, w))-1.0)
        print_debug(f"operator_D_inv :: s: {join_list(s)} b: {b} c: {c} w: {w} mr: {mr}")
        for i in range(mr):
        	s += s_
    
    print_debug(f"operator_D_inv :: s: {join_list(s)} b: {b} c: {c} w: {w} len(s): {len(s)}")

    s = cumsum(s[0:-1])
    
    print_debug(f"operator_D_inv :: cumsum s: {join_list(s)} [{len(s)}]")

    s.insert(0,0)
    
    print_debug(f"operator_D_inv :: cumsum s0: {join_list(s)} [{len(s)}]")

    for i in range(len(s)):
        s[i] = mod(s[i]+b, c)
      
    print_debug(f"operator_D_inv :: cumsum s0 mod: {join_list(s)} [{len(s)}]")

    return s

def operator_D(w, c):
    dw = diff(w)
    for i in range(len(dw)):
        dw[i] = mod(dw[i], c)

    print_debug(f"operator_D :: w: {w} c: {c} dw: {dw}")
      
    return dw

def rho(c, p, e, s):
    print_debug(f"rho :: c: {c} p: {p} e: {e} s: {s} ({len(s)})")

    s_ = list(s)
    s  = s[0:p]
    e_ = [None] * c

    print_debug(f"rho :: s: {s} ({len(s)})")

    for i in range(e, e+c):
        e_[i-e] = i % c

    print_debug(f"rho :: e: {e} ({len(e_)})")

    s += e_

    print_debug(f"rho :: s: {s} ({len(s)})")

    s += s_[p:]

    print_debug(f"rho ::\n\tc: {c}\n\tp: {p}\n\te: {e}\n\ts: {s} ({len(s)})\n\ts_: {s_} ({len(s_)})\n\te_: {e_} ({len(e_)})")

    return s

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

def decodableDeBruijn(T, K, L, c, n):
    print_debug(f"decodableDeBruijn :: c: {c} n: {n}\n\tT: {join_list(T)}\n\tK: {join_list(K)}\n\tL: {join_list(L)}")

    if n <= 2:
        db = DeBruijn(c, 2)
        t  = db.s
        # t_ = "".join([str(x) for x in t[0:2]])
        t_ = t + t[0:2]

        print_debug(f"decodableDeBruijn :: c: {c} n: {n} t: {join_list(t)} t_: {join_list(t_)}")
        
        for i in range(c ** 2):
            k   = t_[i:i+2]
            idx = wordIndex(k, c)
            print_debug(f"decodableDeBruijn :: c: {c} n: {n} t: {join_list(t)} t_: {join_list(t_)} i: {i} k: {k} idx: {idx}")
            T[idx] = i

    else:
        n_      = n-1
        _,_,_,s = decodableDeBruijn(T, K, L, c, n_)
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} s    : {join_list(s)}")
        
        k = findOnes(s, n_)
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} k    : {k}")
        
        s_    = s[0:k] + s[k+1:]
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} k    : {k} s_   : {join_list(s_)} [{len(s_)}]")

        s_hat = operator_D_inv(s_, 0, c)
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} k    : {k} s_hat: {join_list(s_hat)}  [{len(s_hat)}]")

        p = (c-1) * ((c ** n_) - 1) + k
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} p    : {p}")
        
        e = s_hat[p]
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} e    : {e}")
        
        t = rho(c, p, e, s_hat)
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} t    : {join_list(t)} [{len(t)}]")

        lidx = n-3 
        tidx = ((c ** n_)-1)
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} lidx : {lidx}")
        print_debug(f"decodableDeBruijn :: c: {c} n: {n} tidx : {tidx}")
        L[lidx] = t[0:tidx]

    K[n-2] = findOnes(t, n)
    
    print_debug(f"decodableDeBruijn ::\n\tt: {t} ({len(t)})\n\tT: {T} ({len(T)})\n\tK: {K} ({len(K)})\n\tL: {L} ({len(L)})")

    return T, K, L, t

def DeBruijnDecoder(T, K, L, w, c):
    print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c}")

    r = 2
    n = len(w)

    print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n}")

    if (n==r):
        idx = wordIndex(w, c)
        v   = T[idx]
        print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} idx: {idx}")
        return v
  
    v       = operator_D(w, c)
    i       = n-1-r
    k       = K[i]
    p       = (c-1) * ((c ** (n-1)) - 1) + k
    allOnes = True

    for i in range(n-1):
        if v[i] != 1:
            allOnes = False
    
    print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} i: {i} k: {k} p: {p} allOnes: {allOnes}")

    if allOnes:
        lv = L[n-1-r][k]
        e  = lv + ((c-1) ** 2)
        j  = p  + mod(w[0]-e, c)
        print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} lv: {lv} e: {e} j: {j} allOnes: true")

    else:
        f = DeBruijnDecoder(T, K, L, v, c)

        if (f>k):
            f -= 1

        lv= L[n-1-r][f]
        e = lv
        j = f + (((c ** (n-1)) - 1) * mod(e-w[0], c))
        
        print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c} r: {r} n: {n} v: {v} lv: {lv} e: {e} j: {j} allOnes: false f: {f}")

        if (j<0) or (j>p-1):
            j = j+c

    print_debug(f"DeBruijnDecoder :: w: {join_list(w)} c: {c} j: {j}")

    return j

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
        d  = DeBruijnDecoder(T, K, L, w, c)
        print_debug(f"decodeDeBruijn :: i: {i} i_ {i_} w {w} d: {d}")
        assert j_==d, f"j_ == d. j_: {j_} d: {d}"
        J[i] = d

    print_debug(f"decodeDeBruijn :: J: {J}")
    
    return J

def encodeDeBruijn(voc, n):
    print_debug(f"encodeDeBruijn :: voc: {voc} n: {n}")

    c = len(voc)
    T = [None] * (c ** 2)
    K = [None] * (n-1)
    L = [None] * (n-2)

    print_debug(f"encodeDeBruijn :: voc: {voc} c: {c} n: {n} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)}")

    T, K, L, t = decodableDeBruijn(T, K, L, c, n)

    print_debug(f"encodeDeBruijn :: voc: {voc} c: {c} n: {n} T: {join_list(T)} K: {join_list(K)} L: {join_list(L)} t: {join_list(t)}")

    return T, K, L, t, c

def main(voc, n):
    test()

    T, K, L, t, c = encodeDeBruijn(voc, n)
    J             = decodeDeBruijn(T, K, L, t, c, n)
    sstr          = "".join([voc[l] for l in t])

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
    voc = "ACGT"
    n   = 3
    c_  = len(voc)
    L_  = [[0,0,0,1,1,3,3,2,3,1,2,1,3,1,0]]
    K_  = [7,45]
    T_  = [0,2,4,15,1,7,9,6,3,8,12,11,5,10,13,14]
    t_  = [0,0,0,1,1,3,3,2,3,1,2,1,3,1,0,3,3,3,0,0,2,2,1,2,0,1,0,2,0,3,2,2,2,3,3,1,1,0,1,3,0,3,1,3,2,1,1,1,2,2,0,0,3,0,1,2,3,0,2,3,2,0,2,1]

    T, K, L, t, c = encodeDeBruijn(voc, n)

    assert c_ == c, f"{c_} != {c}"
    assert t_ == t, f"{L_} != {L}"
    assert T_ == T, f"{T_} != {T}"
    assert K_ == K, f"{K_} != {K}"
    assert L_ == L, f"{L_} != {L}"

    sstr = "".join([voc[l] for l in t])

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
        karr = [voc.index(p) for p in kmer]
        d    = DeBruijnDecoder(T, K, L, karr, c)
        print(kpos, kmer, karr, d)
        keys[kmer] = d

    for st in range(len(voc)):
        for nd in range(len(voc)):
            for rd in range(len(voc)):
                kmer = voc[st] + voc[nd] + voc[rd]
                kid  = 4**2*st + 4**1*nd + 4**0*rd
                karr = [st, nd, rd]
                d    = DeBruijnDecoder(T, K, L, karr, c)
                assert keys[kmer] == d
                print(kmer, kid, d)

    print_info("all tests passed")


if __name__ == '__main__':
    voc =     sys.argv[1]
    n   = int(sys.argv[2])
    main(voc, n)
