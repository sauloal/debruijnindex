import sys
# https://jgeisler0303.github.io/deBruijnDecode/#decoderTest


class DeBruijn():
    def __init__(self, c, n):
        self.c = c
        self.n = n
        self.a = [0 for j in range(c*n)]
        self.s = []
        
        self.generate(1, 1)

    def generate(self, t, p):
        # print("t: ", t, ", p: ", p, ", c: ", self.c)
        
        if t > self.n:
            if ((self.n%p)==0):
                for j in range(p):
                    self.s.append(self.a[j+1])
                    # print(" ".join([str(x) for x in self.a]))

        else:
            self.a[t] = self.a[t-p]
            self.generate(t+1, p)

            for j in range(self.a[t-p]+1, self.c):
                # print("j: ", j, ", a: ", " ".join([str(x) for x in self.a]))
                self.a[t] = j
                self.generate(t+1, t)

def decodableDeBruijn(T, K, L, c, n):
    if n <= 2:
        db = DeBruijn(c, 2)
        t  = db.s
        # t_ = "".join([str(x) for x in t[0:2]])
        t_ = t + t[0:2]
        
        for i in range(c ** 2):
            T[wordIndex(t_[i:i+2], c)] = i

    else:
        n_      = n-1
        s,_,_,_ = decodableDeBruijn(T, K, L, c, n_)
        print("generating for n= ", n, ", n_= ", n_)
        
        k = findOnes(s, n_)
        print("k= ", k)
        
        s_    = s[0:k] + [s[k+1]]
        s_hat = operator_D_inv(s_, 0, c)
        print("s_hat= ", s_hat)

        p = (c-1) * ((c ** n_) - 1) + k
        print("(c-1)*(Math.pow(c, n_) - 1)= ", n_)
        print("p= ", p)
        
        e = s_hat[p]
        print("e= ", e)
        
        t = rho(c, p, e, s_hat)
        L[n-3] = t[0:((c ** n_)-1)]

    K[n-2] = findOnes(t, n)
    
    return t, T, K, L

def decodeDeBruijn(T, K, L, s, c, n):
    print("decodeDeBruijn")
    print(" T", T)
    print(" K", K)
    print(" L", L)
    print(" c", c)
    print(" n", n)
    print(" s", s)
    J = [None] * c ** n
    for i in range(c ** n):
        w  = [None] * n
        i_ = i
        
        for j in range(n-1, -1, -1):
            # print(" j", j)
            w[j] = i_ % c
            i_ = i_//c

        # j_ = findWord(s, w)
        j  = DeBruijnDecoder(T, K, L, w, c)
        J[i] = j
    print("J", J)
    return J

def DeBruijnDecoder(T, K, L, w, c):
    # print("DeBruijnDecoder")
    # print(" T", T)
    # print(" K", K)
    # print(" L", L)
    # print(" w", w)
    # print(" c", c)

    r = 2
    n = len(w)

    if (n==r):
        return T[wordIndex(w, c)]
  
    v = operator_D(w, c)
    i = n-r-1
    k = K[n-r-1]
    p = (c-1) * ((c ** (n-1)) - 1) + k
    allOnes = True

    for i in range(n-1):
        if v[i] != 1:
            allOnes = False
    
    if allOnes:
        e = L[n-r-1][k] + ((c-1) ** 2)
        j = p + mod(w[0]-e, c)

    else:
        f = DeBruijnDecoder(T, K, L, v, c)

        if (f>k):
            f -= 1

        e = L[n-r-1][f]
        j = f + ((c ** (n-1)) - 1) * mod(e-w[0], c)
        if (j<0) or (j>p-1):
            j = j+c

    return j

def mod(n, m):
    n = n%m
    if (n<0):
        n += m
    return n

def gcd(a,b):
    if (a < 0): a = -a
    if (b < 0): b = -b
    if (b > a):
        temp = a
        a = b
        b = temp

    while True:
        a %= b
        if (a == 0): return b
        b %= a
        if (b == 0): return a

def diff(a):
    # print("diff", a)
    s = [None]*(len(a)-1)
    
    for i in range(1, len(a)):
        s[i-1] = a[i] - a[i-1]

    return s

def cumsum(a):
    for i in range(1, len(a)):
        a[i] = a[i-1] + a[i]

    return a

def wordIndex(w, c):
    # print("wordIndex", w, c)
    idx = 0
    b   = 1
    # w.reverse();

    for i in range(len(w)):
        # print( " idx", idx, "b", b, "i", i, "b * w[i]", b * w[i])
        idx += b * w[i]
        b   *= c

    return idx

def wordIndexOrig(w, c):
    print("wordIndex", w, c)
    idx = 0
    b   = 1
    # w.reverse();

    for i in range(len(w)):
        idx += b * int(w[i])
        b   *= c

    return idx

def findWord(s,w):
    # print("findWord", s, w)
    ls = len(s)
    lw = len(w)
    s_ = s+s[0:lw]
    for k in range(ls - lw):
        # print(" ", k, k+lw, s_[k:k+lw], w, s_[k:k+lw] == w)
        if s[k:k+lw] == w:
            if k > len(s):
                return -1
            return k
    return -1

def findWordOrig(s, w):
    print("findWord", s, w)
    n  = len(w)
    s_ = "".join([str(x) for x in s+s[0:n]])
    wi_ = 0
    k_ = 0
    print("S_", s_)

    # for(k= 0; k<s.length && wi!=n; k++):
    for k in range(len(w)):
        print(" k", k)
        k_ = k
        if wi_ == n: break
        # for(wi= 0; wi<n && s_[k+wi]==w[wi]; wi++)
        for wi in range(n):
            print("  wi", wi, k+wi)
            wi_ = wi
            if s_[k+wi] == w[wi]: break
        
    k_ -= 1

    if k_ >= len(s):
        k_ = -1

    return k_

def findOnes(s, n):
    w = [None]*n
    for i in range(n):
        w[i]= 1
    return findWord(s, w)

def operator_D_inv(s, b, c):
    w = sum(s) % c
    # print("w= ", w)

    if (w!=0):
        s_ = s
        for i in range(c//gcd(c, w)-1):
        	s+= s_
    
    # print("s= ", s)
    s = cumsum(s[0:-1])
    
    # print("s= ", s)
    s.insert(0,0)
    
    # print("s= ", s)
    for i in range(len(s)):
        s[i] = mod(s[i]+b, c)
      
    return s

def operator_D(w, c):
    # print("operator_D", w, c)
    dw = diff(w)
    for i in range(len(dw)):
        dw[i] = mod(dw[i], c)
      
    return dw

def rho(c, p, e, s):
    # print("rho c", c, "p", p, "e", e)
    s_ = s
    s = s[0:p]
    # print("s= ", s, len(s))
    e_ = [None] * c

    for i in range(e, e+c):
        e_[i-e] = i%c
        # print("e_= ", e_)

    s += e_
    # print("s_.slice(p)= ", s_[p])

    s += [s_[p]]
    # print("rho s", s)

    return s

def encodeDeBruijn(voc, n):
    c = len(voc)
    T = [None] * (c ** 2)
    K = [None] * (n-1)
    L = [None] * (n-2)

    s = decodableDeBruijn(T, K, L, c, n)

    return T, K, L, s, c

def main(voc, n):
    test()
    T, K, L, s, c = encodeDeBruijn(voc, n)
    J = decodeDeBruijn(T, K, L, s, c, n)

def test():
    voc = "ACGT"
    n   = 3
    L0  = [0,0,0,1,1,3,3,2,3,1,2,1,3,1,0]
    K   = [7,45]
    T   = [0,2,4,15,1,7,9,6,3,8,12,11,5,10,13,14]

    T_, K_, L_, s_, c_ = encodeDeBruijn(voc, n)

    assert c_    == len(voc), f"{c_} != {len(voc)}"
    assert L_[0] == L0      , f"{L_[0]} != {L0}, {L_}"
    assert K_    == K       , f"{K_} != {K}"
    assert T_    == T       , f"{T_} != {T}"


if __name__ == '__main__':
    voc =     sys.argv[1]
    n   = int(sys.argv[2])
    main(voc, n)






"""
function testDecoder(viewContainer) {
    clearDecoderTest(viewContainer);
    
    testTable = document.createElement("table");
    var testRow;
    var testCell;
    var testStr;
    testRow= document.createElement("tr");
    
    testCell= document.createElement("th");
    testCell.appendChild(document.createTextNode("No"));
    testRow.appendChild(testCell);
    
    testCell= document.createElement("th");
    testCell.appendChild(document.createTextNode("Pattern"));
    testRow.appendChild(testCell);
    
    testCell= document.createElement("th");
    testCell.appendChild(document.createTextNode("Found at"));
    testRow.appendChild(testCell);
    
    testCell= document.createElement("th");
    testCell.appendChild(document.createTextNode("Decoded at"));
    testRow.appendChild(testCell);
    
    testCell= document.createElement("th");
    testCell.appendChild(document.createTextNode("Correct"));
    testRow.appendChild(testCell);
    
    testTable.appendChild(testRow);
    
    for(var i= 0; i<Math.pow(c, n); i++) {
      var w= new Array(n);
      var i_= i;
      for(j= n-1; j>=0; j--) {
          w[j]= i_%c;
          i_= Math.floor(i_/c);
      }
      var j_= findWord(s, w);
      var j= decodeDeBruijn(w, c);
      
      testRow= document.createElement("tr");
  
      testCell= document.createElement("td");
      testCell.style.textAlign= "center";
      testCell.appendChild(document.createTextNode(i));
      testRow.appendChild(testCell);
      testCell= document.createElement("td");
      testCell.style.textAlign= "center";
      // testCell.appendChild(document.createTextNode(w));
      change(testCell, w);
      testRow.appendChild(testCell);
      testCell= document.createElement("td");
      testCell.style.textAlign= "center";
      testCell.appendChild(document.createTextNode(j_));
      testRow.appendChild(testCell);
      testCell= document.createElement("td");
      testCell.style.textAlign= "center";
      testCell.appendChild(document.createTextNode(j));
      testRow.appendChild(testCell);
      testCell= document.createElement("td");
      testCell.style.textAlign= "center";
      if(j==j_) {
        testCell.style.color= "green";
        testCell.appendChild(document.createTextNode("\u2714"));
      } else {
        testCell.style.color= "red";
        testCell.appendChild(document.createTextNode("\u2718"));
      }
      testRow.appendChild(testCell);
      testTable.appendChild(testRow);
    }
    viewContainer.appendChild(testTable);
  }
  
  function displayDecoderData(viewContainer, numberContainer) {
    var data= viewContainer.getElementsByTagName("code");
    while(data.length>0)
      viewContainer.removeChild(data[0]);
  
    data = document.createElement("code");
    var lenL= 0;
    for(var i= 0; i<L.length; i++) {
      data.appendChild(document.createTextNode("L[" + i + "]= {" + L[i] + "};\n"));
      lenL+= L[i].length;
    }
    data.appendChild(document.createTextNode("K= {" + K + "};\n"));
    data.appendChild(document.createTextNode("T= {" + T + "};"));
    
    viewContainer.appendChild(data);
    var dataBitsT= Math.log(Math.pow(c, 2))/Math.log(2) * T.length;
    var dataBitsK= Math.log(Math.pow(c, n))/Math.log(2) * K.length;
    var dataBitsL= Math.log(c)/Math.log(2) * lenL;
    
    numberContainer.replaceChild(document.createTextNode(dataBitsT+dataBitsK+dataBitsL), numberContainer.firstChild);
  }
  
  function displayCDecoderData(viewContainer, numberContainer) {
    var data= viewContainer.getElementsByTagName("code");
    while(data.length>0)
      viewContainer.removeChild(data[0]);
  
    data = document.createElement("code");
  
    data.appendChild(document.createTextNode("#define N " + n + "\n"));
    data.appendChild(document.createTextNode("#define C " + c + "\n"));
    data.appendChild(document.createTextNode("#define powC_2 " + ((c-1)*(c-1)) + "\n"));
  
    var powC_N= Array();
    for(var i= 2; i<=n; i++)
      powC_N[i-2]= (c-1)*(Math.pow(c, i) - 1);
    data.appendChild(document.createTextNode("const uint32_t powC_N[]= {" + powC_N + "};\n"));
    data.appendChild(document.createTextNode("const uint8_t L[]= {" + L + "};\n"));
    var idxL= Array();
    idxL[0]= 0;
    for(var i= 1; i<L.length; i++)
      idxL[i]= idxL[i-1] + L[i-1].length;
    data.appendChild(document.createTextNode("const uint16_t idxL[]= {" + idxL + "};\n"));
    data.appendChild(document.createTextNode("const uint8_t T[]= {" + T + "};\n"));
  
    var w= new Array(n);
    var locBD= new Array(Math.pow(c, n));
    for(var i= 0; i<Math.pow(c, n); i++) {
      var i_= i;
      for(j= n-1; j>=0; j--) {
      w[j]= i_%c;
      i_= Math.floor(i_/c);
      locBD[i]= decodeDeBruijn(w, c);
      }
    }
    data.appendChild(document.createTextNode("const uint32_t locDB[]= {" + locBD + "};\n"));
    
    viewContainer.appendChild(data);
    
  
    var lenL= 0;
    for(var i= 0; i<L.length; i++) {
      lenL+= L[i].length;
    }  
    var dataBitsT= Math.log(Math.pow(c, 2))/Math.log(2) * T.length;
    var dataBitsK= Math.log(Math.pow(c, n))/Math.log(2) * K.length;
    var dataBitsL= Math.log(c)/Math.log(2) * lenL;
    
    numberContainer.replaceChild(document.createTextNode(dataBitsT+dataBitsK+dataBitsL), numberContainer.firstChild);
  }
  """