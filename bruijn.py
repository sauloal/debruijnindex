# https://en.wikipedia.org/wiki/De_Bruijn_sequence
def de_bruijn(alphabet, n: int) -> str:
    """
    de Bruijn sequence for alphabet k
    and subsequences of length n.
    """
    k = len(alphabet)

    a = [0] * k * n
    sequence = []

    def db(t, p):
        if t > n:
            if n % p == 0:
                sequence.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)
    db(1, 1)
    return "".join(alphabet[i] for i in sequence)

print(de_bruijn("abcd", 2))

#print(de_bruijn(2, 3))
b=de_bruijn("ACGT", 3)
print("{:3d} {}".format(len(b), b))
for s in range(len(b)-2):
    print("{:3d}{} {}".format(s, " "*s, b[s:s+3]))
