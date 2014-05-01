from fractions import Fraction

## STRUCTURES
class Z:
    one = 1
    zero = 0
    @staticmethod
    def divide(a, b):
        """
        Euclidian domain division:
        returns (q, r) s.t. a = qb + r with |r| < |b|
        """
        return (a // b), (a % b)

# TODO: fraction fields in general
class Q:
    """ Fraction field of Z, elements are Fractions """
    one = Fraction(1, 1)
    zero = Fraction(0, 1)


# TODO: use dicts instead for cleanliness
class Polynomial:
    def __init__(self, R, coeffs):
        """ A polynomial in R[x], for given ring R
        p(x) = c[0] + c[1]x + c[2]x^2 + ...
        """
        self.R = R

        # trim trailing zero coefficients
        max_nonzero = len(coeffs)-1
        while coeffs[max_nonzero] == R.zero and max_nonzero > 0:
            max_nonzero -= 1
        self.coeffs = coeffs[:max_nonzero+1] # represents the zero-polynomial as [0]

    def degree(self):
        return len(self.coeffs) - 1

    def __add__(self, other):
        ncoeffs = []
        for i in xrange(1 + max(self.degree(), other.degree())):
            c = self.R.zero
            if i < len(self.coeffs):
               c += self.coeffs[i]
            if i < len(other.coeffs):
               c += other.coeffs[i]
            ncoeffs.append(c)
        return Polynomial(self.R, ncoeffs)

    def __mul__(self, other):
        # TODO: FFT convolution ?!?!
        d1, d2 = self.degree(), other.degree()
        maxN = d1 + d2
        pc = []
        for n in xrange(maxN + 1):
            c = self.R.zero # coeff of x^n in result
            for i in xrange(max(0, n - d2),  min(d1, n) + 1):
                c += self.coeffs[i] * other.coeffs[n - i]
            pc.append(c)
        return Polynomial(self.R, pc)

    def __sub__(self, other):
        return self + (-other)
    def __neg__(self):
        return Polynomial(self.R, [-c for c in self.coeffs])
    def __eq__(self, other):
        return self.coeffs == other.coeffs
    def __ne__(self, other):
        return not self.coeffs == other.coeffs
    def __str__(self):
        return " + ".join([ "%sx^%d" % (str(c), n) for n, c in
            enumerate(self.coeffs)])
    def __repr__(self):
        return self.__str__()

class ConstPoly(Polynomial):
    def __init__(self, R, c):
        Polynomial.__init__(self, R, [c])

class MonicPoly(Polynomial):
    """ x^n in R[x] """
    def __init__(self, R, n):
        Polynomial.__init__(self, R, [R.one if i == n else R.zero for i in
            xrange(n+1)]) # use dicts...


class Fx:
    """ The Euclidian domain F[x], constructed from field F """
    def __init__(self, F):
        self.F = F
        self.one = Polynomial(F, [F.one])
        self.zero = Polynomial(F, [F.zero])

    def divide(self, a, b):
        """
        Euclidian domain division:
        returns (q, r) s.t. a = qb + r with deg r < deg b
        """
        dd = a.degree() - b.degree()
        if dd < 0 or a == self.zero:
            return self.zero, a
        q1 = ConstPoly(self.F, (a.coeffs[-1] / b.coeffs[-1])) \
                * MonicPoly(self.F, dd)
        r1 = a - q1*b
        qn, rn = self.divide(r1, b)
        return (q1 + qn, rn)



## ALGORITHMS
def egcd(R, a, b):
    """
    extended-gcd of (a, b) over Euclidian Domain R
    returns (d, x, y) s.t. d = xa + yb
    """
    q, r = R.divide(a, b)  # a = qb + r
    if r == R.zero:
        return (b, R.zero, R.one)
    d, j, k = egcd(R, b, r) # d = jb + kr
    return (d, k, j-k*q)



## EXAMPLES

F = Fx(Q)

"""
z1 = Polynomial(Z, [1, 2, 1])
z2 = Polynomial(Z, [5, 2])
print z1 * z2

f1 = Polynomial(Q, map(Fraction, [0, -1, 2]))   # 2x^2 - x = (2x+1)(x-1) + 1
p1 = Polynomial(Q, map(Fraction, [-1, 1]))      # x-1
p2 = Polynomial(Q, map(Fraction, [-5, 2]))      # 2x-5
print F.divide(f1, p1)

fA = p1*p2*p2
fB = p1*p1*p1
d, x, y = egcd(F, fA, fB)
print "gcd (%s, %s) = %s" % (fA, fB, d)
assert d == fA*x + fB*y
"""

#D&F 9.2.9
a = Polynomial(Q, map(Fraction, [1, 1, 1, 2, 0, 1]))
b = Polynomial(Q, map(Fraction, [1, 2, 2, 2, 1, 1]))
d, x, y = egcd(F, a, b)
assert d == a*x + b*y
print "gcd (a = %s, b = %s) = %s = (%s)*a + (%s)*b" % (a, b, d, x, y)

#D&F 13.1.2
a = Polynomial(Q, map(Fraction, [1, 1, 1]))
b = Polynomial(Q, map(Fraction, [-2, -2, 0, 1]))
d, x, y = egcd(F, a, b)
assert d == a*x + b*y
print "gcd (a = %s, b = %s) = %s = (%s)*a + (%s)*b" % (a, b, d, x, y)
