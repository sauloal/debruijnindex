from logger import *
import math

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