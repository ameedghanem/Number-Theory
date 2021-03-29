import sys
import logging
from math import floor

"""
@Author: Ameed S. Ghanem


Number Theory Algorithms.
I implemented in this module various algorithms taught in the number theory course.

1) GCD - the euclidean algorithm
2) XGCD - the extended euclidean algorithm
3) Diophantine equation
4) Modular Exponentiation
5) Legendre symbol
6) Jacobi symbol
7) inverse modulo p
8) finding square root modulo p
9) check if a quadratic equation is solvable
10) chinese remainder theorem
11) the cycle of fraction of 1/p where p is a prime number
12) the euler totient function. phi
13) convert a continued fraction into a rational value
14) convert a rational value into a continued fraction
15) Pell equation
16) finding premitive root
17) lcm
18) sieve of eratosthenes
19) factoring a number to it's prime factors
20) Primality testing using Fermat's little theorem

"""


def gcd(a, b):
	""" computes the gcd of a and b """
	if b == 0:
		return a
	return gcd(b, a%b)


def lcm(a, b):
    """ returns leam common multiplicity (lcm) of a and b """
    return a*b//gcd(a,b)


def xgcd(a, b):
	""" computes the gcd and finds x,y s.t. aX + bY = gcd(a,b) """
	if b == 0:
		return a, 1, 0
	g, x, y = xgcd(b, a%b)
	new_x = y
	new_y = x-y*(a//b)
	return g, new_x, new_y


def d_equation(a, b, c):
	""" solves the diaphontine equation: aX + bY = c """
	g, x, y = xgcd(a, b)
	if c%g == 0:
		factor = c//g
		return factor*x, factor*y
	return None

def my_pow(a, b, p):
	""" computes a**b % p in O(log(p)) """
	
	def my_pow_rec(a, b, p):
		if b == 0:
			return 1
		x = my_pow_rec(a, b//2, p) % p
		x = (x*x) % p
		if b%2 == 0:
			return x
		return (x*a) % p
	if not ((1).__class__ == type(a) == type(b) == type(p)):
		raise("a, b and p must be integer numbers")
	return my_pow_rec(a,b,p)


def legendre(a, p):
	""" computes the legendre symbol (a/p) """
	if p == 2:
		print("undefined")
		return None
	if a%p == 0:
		return 0
	sign = my_pow(a, (p-1)//2, p)
	if sign == p-1:
		return -1
	return sign


def jacobi(a, p):
	""" computes the jacobi symbol (a/p) """
	if a%p == 0:
		return 0
	if p%2 == 0:
		print("p must be an odd")
		return
	sign = 1
	factors = factorize(p)
	for q, r in factors:
		sign *= legendre(a, q)
	return sign


def sq_root(a, p):
	if legendre(a, p) == -1:
		return None
	if p%4 == 3:
		root = my_pow(a, (p+1)//4, p)
		return root, p-root
	return naive_sq_root(a, p)


def naive_sq_root(a, p):
	""" compues the square root of a modulo p """
	root1, root2 = 0, 0
	if a == 0:
		return 0
	for i in range(1, p//2 + 1):
		if i**2 % p == a:
			root1, root2 = i, p-i
	return root1, root2


def inverse(a, p):
	""" comutes a**(-1) modulo p """
	sols = d_equation(a, p, 1)
	return sols[0]%p


def q_equation(a, b, c, p):
	""" aX^2 + bX + c = 0 mod p """
	sols = []
	desc = b**2 - 4*a*c
	num_of_sols = legendre(desc, p)+1
	if num_of_sols == 0:
		return False
	return True


def crt(a, b):
	""" for each i:  x = a[i] mod b[i]  """

	assert len(a) == len(b)
	if a == []:
		return None
	l = len(a)
	mod = 1
	for num in b:
		mod *= num
	m = [mod//b[i] for i in range(l)]
	m_gag = [inverse(m[i], b[i]) for i in range(l)]
	sol = sum([m[i]*m_gag[i]*a[i]for i in range(l)]) % mod
	return sol, mod


def find_cycle_of_1_over_prime(p):
	""" find the cycle of the fraction 1/p where p is a prime number """
	ord = 1
	for i in range(1,phi(p)+1):
		if 10**i % p == 1:
			ord = i
			break
	x = (10**ord - 1)/p
	x = 1.0*x
	cycle = x/(10**ord)
	print("cycle length equals %d" % ord)
	return cycle


def phi(n):
	"""
	Euler totient function. computes how many numbers realtivly prime to n.
	Which means it computes #{x | gcd(x,n) = 1}
	"""
	res = n
	p = 2
	while p*p <= n:
		if n%p == 0:
			while n%p == 0:
				n /= p
			res -= res//p
		p += 1
	if n > 1:
		res -= res//n
	return res
	

def get_premitive_roots(n):
	""" finds all premitive roots  of n """
	phi_n = phi(n)
	phi_phi_n = phi(phi_n)


def get_rational_val(continued_fractions):
	"""given a continued fraction. return the relevant ratioal number"""
	a, a_next = 1, continued_fractions[0]
	b, b_next = 0, 1
	a_res, b_res = 0, 0
	for i in range(2, len(continued_fractions)+1):
		ci = continued_fractions[i-1]
		a_res = ci*a_next + a
		b_res = ci*b_next + b
		a, a_next = a_next, a_res
		b, b_next = b_next, b_res
	return a_next, b_next
	
	
def get_continued_fraction(d):
	"""given an integer d, returns the continued fraction of d**0.5"""
	
	d_tag = 1.0*int(d**0.5)
	d = d**0.5
	d0 = d_tag
	r = -1.0
	continued_fractions = [d0]
	while d_tag != 2*d0:
		r = 1/(d-d_tag)
		d_tag = floor(r)
		continued_fractions.append(d_tag)
		d = r
	return [int(x) for x in continued_fractions]


def pell_equation(d, a=1):
	""" X^2 - d*Y^2 = a.   a = 1 or -1 """
	
	def get_sol_from_odd_eq(x, y):
		ret_x = x**2+d*y**2
		ret_y = 2*x*y
		return int(ret_x), int(ret_y)

	if d < 0:
		print("Illegal d value")
		return None	
	root = d**0.5
	if root == int(root):
		print("Illegal d value")
		return None
	continued_fractions = get_continued_fraction(d)
	cycle_len = len(continued_fractions)-1
	x, y = get_rational_val(continued_fractions[:-1])
	x, y = int(x), int(y)
	if a == -1:
		if cycle_len%2 == 1:
			return x, y
		print("No solution for the odd pell equation.")
		return None
	if a == 1:
		if cycle_len % 2 == 0:
			return x,y
		
	x, y = get_sol_from_odd_eq(x, y)
	return x, y


def sieve(n):
    """
    Sieve of Eratosthenes: returns a boolean list
    list[i]=True means i+1 is a prime number
    lit[i]=False means i+1 is a composite number
    """
    isPrime = [True for i in range(n+1)]
    isPrime[0], isPrime[1] = False, False
    i = 2
    while i**2 <= n:
        if isPrime[i]:
            j = i**2
            while j <= n:
                isPrime[j] = False
                j += i
        i += 1
    return isPrime[1:]


def factorize(n, to_print=False):
    """ return all prime factors of n. if to_print=True it prints n's all prime factors """
    factors = []
    pwr = 0 # this will hold how many each prime appears in n
    while n%2 == 0:
        if to_print:
            print(2, end=' ')
        pwr += 1
        n /= 2
    if pwr > 0:
        factors.append((2, pwr))
    j = 3
    while j <= int(n**0.5)+1:
        pwr = 0
        while n%j == 0:
            if to_print:
                print(j, end=' ')
            pwr += 1
            n /= j
        if pwr > 0:
            factors.append((j, pwr))
        j += 2
    if n > 2:
        factors.append((int(n), 1))
    if to_print:
        print(int(n))       
    return factors


def is_prime(p, rounds=100):
	""" checks whether p is a prime number or not """
	for i in range(100):
		a = random.randint(1,p-1)
		if pow2(a, p-1, p) != 1:
			return False
	return True
