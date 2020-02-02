var s;
var n;
var c;
var T;
var K;
var L;


function debugLog() {
  console.log(...arguments);
}

function mod(n, m) {
  n= n%m;
  if(n<0) n+= m;
  return n;
}

function gcd(a,b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (b > a) {var temp = a; a = b; b = temp;}
    while (true) {
        a %= b;
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
    }
}

function sum(a) {
  var s= 0;
  for(var i= 0; i<a.length; i++)
    s+= a[i];
  return s;
}

function diff(a) {
  var s= new Array(a.length-1);
  for(var i= 1; i<a.length; i++)
    s[i-1]= a[i]-a[i-1];
    
  return s;
}

function cumsum(a) {
  for(var i= 1; i<a.length; i++)
    a[i]= a[i-1]+a[i];
    
  return a;
}

function wordIndex(w, c) {
  var idx = 0;
  var b   = 1;
  // w.reverse();
  for(var i= 0; i<w.length; i++) {
    idx += b * w[i];
    b   *= c;
  }

  return idx;
}

function findWord(s, w) {
  var k;
  var n= w.length;
  var s_= s.concat(s.slice(0, n));
  var wi= 0;
  for(k= 0; k<s.length && wi!=n; k++)
    for(wi= 0; wi<n && s_[k+wi]==w[wi]; wi++);
  
  k--;
  if(k>=s.length)
    k= -1;
  return k;
}

function findOnes(s, n) {
  var w= new Array(n);
  for(var i= 0; i<n; i++)
    w[i]= 1;
  return findWord(s, w);
}

function operator_D_inv(s, b, c) {
  w = sum(s) % c;

  // debugLog("operator_D_inv :: s: " + s.join(" ") + " b: " + b + " c: " + c + " w: " + w);

  if (w!=0) {
    var s_= s;
    var mr = ((c/gcd(c, w))-1);
    // debugLog("operator_D_inv :: s: " + s.join(" ") + " b: " + b + " c: " + c + " w: " + w + " mr: " + mr);

    for(var i= 0; i<mr; i++) {
      s= s.concat(s_);
    }
  }

  // debugLog("operator_D_inv :: s             : " + s.join(" ") + " b: " + b + " c: " + c + " w: " + w + " s.length: " + s.length);
  
  s= cumsum(s.slice(0, -1));

  // debugLog("operator_D_inv :: s cumsum      : " + s.join(" ") + " b: " + b + " c: " + c + " w: " + w + " s.length: " + s.length);

  s.unshift(0);

  // debugLog("operator_D_inv :: s cumsum 0    : " + s.join(" ") + " b: " + b + " c: " + c + " w: " + w + " s.length: " + s.length);

  for(var i= 0; i<s.length; i++) {
    s[i]= mod(s[i]+b, c);
  }
    
  // debugLog("operator_D_inv :: s cumsum 0 mod: " + s.join(" ") + " b: " + b + " c: " + c + " w: " + w + " s.length: " + s.length);

  return s;
}

function rho(c, p, e, s) {
  var s_= s;
  
  // debugLog("rho :: s= ", s);

  s= s.slice(0, p);
  // debugLog("rho :: s= ", s);
  
  var e_= new Array();
  for(var i= e; i<e+c; i++) {
    e_[i-e]= i%c;
  }
  // debugLog("rho :: e_= ", e_);

  s= s.concat(e_);
  // debugLog("rho :: s=", s);
  
  s= s.concat(s_.slice(p));
  // debugLog("rho :: s=", s);
  
  return s;
}

function operator_D(w, c) {
  var dw= diff(w);
  for(var i= 0; i<dw.length; i++)
    dw[i]= mod(dw[i], c);
    
  return dw;
}




function DeBruijn(c, n) {
  this.c= c;
  this.n= n;
  this.a= new Array(c*n);
  this.s= new Array();
  debugLog("DeBruijn. :: c: "+ c + ", n: " + n + ", a: " + this.a.join(" ") + ", s: " + this.s.join(" "));
  for(var j= 0; j<c*n; j++) {
    this.a[j]= 0;
  }
  this.generate(1, 1);
  debugLog("DeBruijn. :: c: "+ c + ", n: " + n + ", a: " + this.a.join(" ") + ", s: " + this.s.join(" "));
}


DeBruijn.prototype.generate= function(t, p) {
  // debugLog("DeBruijn.generate :: c: " + this.c + ", n: " + this.n + ", t: "+ t + ", p: " + p);
  if(t>this.n) {
    if((this.n%p)==0) {
      for(var j= 0; j<p; j++){
	      this.s.push(this.a[j+1]);
        // debugLog("DeBruijn.generate :: c: " + this.c + ", n: " + this.n + ", t: "+ t + ", p: " + p + ", j: " + j + ", a: " + this.a.join(" ") + ", s: " + this.s.join(" "));
      }
    }
  } else {
    this.a[t]= this.a[t-p];
    this.generate(t+1, p);
    for(var j= this.a[t-p]+1; j<this.c; j++) {
      // debugLog("DeBruijn.generate :: c: " + this.c + ", n: " + this.n + ", t: "+ t + ", p: " + p + ", j: "+ j + ", a: " + this.a.join(" ") + ", s: " + this.s.join(" "));
      
      this.a[t]= j;
      this.generate(t+1, t);
    }
  }
}

function decodableDeBruijn(c, n) {  
  debugLog("decodableDeBruijn :: c: " + c + ", n: " + n);

  var t;
  if(n<=2) {
      db     = new DeBruijn(c, 2);
      t      = db.s;
      var t_ = t.concat(t.slice(0, 2));

      // debugLog("decodableDeBruijn :: c: " + c + ", n: " + n + ", t: " + t.join(" ") + ", t_: " + t_.join(" "));

      for(var i= 0; i<Math.pow(c, 2); i++) {
          k = t_.slice(i, i+2);
          idx = wordIndex(k, c);
          // debugLog("decodableDeBruijn :: c: " + c + ", n: " + n + ", t: " + t.join(" ") + ", t_: " + t_.join(" ") + ", i: " + i + ", k: " + k + ", idx: " + idx);
          T[idx]= i;
      }
  } else {
    var n_ = n-1;
    s = decodableDeBruijn(c, n_);
    debugLog("decodableDeBruijn :: c: " + c + ", n: " + n + ", s: " + s.join(" "));

    var k = findOnes(s, n_);
    debugLog("decodableDeBruijn :: c: " + c + ", n: " + n + ", k: " + k);
    
    var s_    = s.slice(0, k).concat(s.slice(k+1));
    var s_hat = operator_D_inv(s_, 0, c);
    debugLog("decodableDeBruijn :: c: " + c + ", n: " + n + ", s_   : " + s_   .join(" ") + " ["+s_   .length+"]");
    debugLog("decodableDeBruijn :: c: " + c + ", n: " + n + ", s_hat: " + s_hat.join(" ") + " ["+s_hat.length+"]");

    p= (c-1)*(Math.pow(c, n_) - 1) + k;
    console.log("(c-1)*(Math.pow(c, n_) - 1)= " + n_);
    console.log("p= " + p);
    e= s_hat[p];
    console.log("e= " + e);
    t= rho(c, p, e, s_hat);
    L[n-3]= t.slice(0, Math.pow(c, n_)-1);
  }    

  fo = findOnes(t, n)
  K[n-2] = fo;

  debugLog("decodableDeBruijn ::\nt: ", t, "\nT: ", T, "\nK: ", K, "\nL: ", L);

  return t;
}

function decodeDeBruijn(w, c) {  
  // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c);

  var r= 2;
  var n= w.length;

  // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c + ", r: " + r + ", n: " + n);

  if(n==r) {
    idx = wordIndex(w, c);
    v   = T[idx];
    // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c + ", r: " + r + ", n: " + n + ", v: " + v + ", idx: " + idx);
    return v;
  }
  
  var v= operator_D(w, c);
  var i= n-r-1;
  var k= K[n-r-1];
  var p= (c-1)*(Math.pow(c, n-1) - 1) + k;
  var allOnes= true;
  for(var ib= 0; ib<n-1; ib++) {
    if(v[ib]!=1){
      allOnes= false;
    }
  }

  // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c + ", r: " + r + ", n: " + n + ", v: " + v + ", i: " + i + ", k: " + k + ", p: " + p + ", allOnes: " + allOnes);
    
  var e;
  var j;
  if(allOnes) {
    lv= L[n-r-1][k];
    e= lv + Math.pow(c-1, 2);
    j= p + mod(w[0]-e, c);
    // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c + ", r: " + r + ", n: " + n + ", v: " + v + ", i: " + i + ", k: " + k + ", p: " + p + ", allOnes: " + allOnes + ", e: " + e + ", j: " + j + ", lv: " + lv);

  } else {
    var f= decodeDeBruijn(v, c);

    if(f>k) {
      f= f-1;
    }

    lv = L[n-r-1][f];
    e= lv;
    j= f + (Math.pow(c, n-1) - 1) * mod(e-w[0], c);

    // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c + ", r: " + r + ", n: " + n + ", v: " + v + ", i: " + i + ", k: " + k + ", p: " + p + ", allOnes: " + allOnes + ", e: " + e + ", j: " + j + ", lv: " + lv + ", f: " + f);

    if(j<0 || j>p-1) {
      j= j+c;
    }
  }

  // debugLog("decodeDeBruijn :: w: " + w + ", c: " + c + ", j: " + j);

  return j;
}