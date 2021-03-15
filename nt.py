import os, sys
import logging
from math import floor

"""

Number Theory Algorithms.
I implemented in this module various algorithms taught in the number theory course.

1) GCD - the euclidean algorithm
2) XGCD - the extended euclidean algorithm
3) Diophantine equation
4) Modular Exponentiation
5) Legendre symbol
6) square root modulo p
7) quadratic modular equation
8) chinese remainder theorem
9) the cycle of fraction of 1/p where p is a prime number
10) the euler totient function. phi
11) convert a continued fraction into a rational value
12) convert a rational value into a continued fraction
13) Pell equation
14) finding inverse modulu p
15) finding premitive root

"""


def gcd(a, b):
	if b == 0:
		return a
	return gcd(b, a%b)
	

def xgcd(a, b):
	if b == 0:
		return a, 1, 0
	g, x, y = xgcd(b, a%b)
	new_x = y
	new_y = x-y*(a//b)
	return g, new_x, new_y


def d_equation(a, b, c):
	g, x, y = xgcd(a, b)
	if c%g == 0:
		factor = c//g
		return factor*x, factor*y
	return None

def my_pow(a, b, p):
	
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
	sign = my_pow(a, (p-1)//2, p)
	if sign == p-1:
		return -1
	return sign


def sq_root(a, p):
	if legendre(a, p) == -1:
		return None
	if p%4 == 3:
		return my_pow(a, (p+1)//4, p)
	return None


def naive_sq_root(a, p):
	sols = []
	for i in range(p):
		if i**2 % p == a:
			sols.append(i)
	return sols


def inverse(a, p):
	sols = d_equation(a, p, 1)
	return sols[0]%p


def q_equation(a, b, c, p):
	""" aX^2 + bX + c = 0 mod p """
	sols = []
	desc = b**2 - 4*a*c
	num_of_sols = legendre(desc, p)+1
	if num_of_sols == 0:
		return None
	root_desc = naive_sq_root(desc, p)
	sols.append(((-b + root_desc)*inverse(2*a, p)) % p)
	sols.append(((-b - root_desc)*inverse(2*a, p)) % p)
	return sols[0], sols[1]


def crt(a, b):
	""" for each i:  x = a[i] mod b[i]  """

	assert len(a) == len(b)
	l = len(a)
	mod = 1
	for num in b:
		mod *= num
	m = [mod//b[i] for i in range(l)]
	m_gag = [inverse(m[i], b[i]) for i in range(l)]
	sol = sum([m[i]*m_gag[i]*a[i]for i in range(l)]) % mod
	return sol, mod


def find_cycle_of_1_over_prime(p):
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
