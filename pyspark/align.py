#!/usr/bin/env python

import sys
import operator

S	= -11.
E	= -2.
MAT = 8.
MIS = -3.
MIN = -float("inf")
B = 4


def _match(x, y):
	if x == y:
		return MAT
	else:
		return MIS

def initMat(M, I, D):
	nI = len(M)
	nJ = len(M[0])

	M[0][0] = 0
	for j in xrange(1, nJ):
		M[0][j] = MIN
	for i in xrange(1, nI):
		M[i][0] = MIN

	for j in xrange(nJ):
		I[0][j] = MIN
	for i in xrange(1, nI):
		I[i][0] = S + E*(i-1)

	for i in xrange(nI):
		D[i][0] = MIN
	for j in xrange(1, nJ):
		D[0][j] = S + E*(j-1)


def align_affine(q, r, bFwd, M, I, D):
	nI = len(q) + 1
	nJ = len(r) + 1

	nMax = MIN
	iMax = 0
	jMax = 0
	mMax = 0 # 0:M, 1:I, 2:D
	for i in xrange(1, nI):
		#for j in xrange(1, nJ):
		for j in xrange(max(1,i-B), min(i+B+1,nJ)):
			#print i, j
			I[i][j] = max(M[i-1][j]+S, I[i-1][j]+E)
			if nMax < I[i][j]:
				nMax = I[i][j]
				iMax = i
				jMax = j
				mMax = 1
			
			D[i][j] = max(M[i][j-1]+S, D[i][j-1]+E)
			if nMax < D[i][j]:
				nMax = D[i][j]
				iMax = i
				jMax = j
				mMax = 2

			M[i][j] = max(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + _match(q[i-1],r[j-1])
			if nMax < M[i][j]:
				nMax = M[i][j]
				iMax = i
				jMax = j
				mMax = 0

			if nMax < -60:
				return ('', 0, 0)

	#print 'M'
	#print np.array(M)[0:nI, 0:nJ]
	#print '*'*15
	#print 'I'
	#print np.array(I)[0:nI, 0:nJ]
	#print '*'*15
	#print 'D'
	#print np.array(D)[0:nI, 0:nJ]
	#print nMax, mMax, iMax, jMax

	if nMax < 0:
		return ('', 0, 0)

	lq = []
	lr = []
	i = iMax
	j = jMax
	m = mMax
	while i>0 and j>0:
		if m == 0:
			#print m, i, j
			lq.append(q[i-1])
			lr.append(r[j-1])
			if M[i][j] == M[i-1][j-1] + _match(q[i-1],r[j-1]):
				m = 0
			elif M[i][j] == I[i-1][j-1] + _match(q[i-1],r[j-1]):
				m = 1
			elif M[i][j] == D[i-1][j-1] + _match(q[i-1],r[j-1]):
				m = 2
			i -= 1
			j -= 1
		elif m == 1:
			#print m, i, j
			lq.append(q[i-1])
			lr.append('-'*len(q[i-1]))
			if I[i][j] == M[i-1][j] + S:
				m = 0
			elif I[i][j] == I[i-1][j] + E:
				m = 1
			i -= 1
		elif m == 2:
			#print m, i, j
			lq.append('-'*len(r[j-1]))
			lr.append(r[j-1])
			if D[i][j] == M[i][j-1] + S:
				m = 0
			elif D[i][j] == D[i][j-1] + E:
				m = 2
			j -= 1

	s = ''
	if True == bFwd:
		s += ' '.join(lq[::-1])
		s += '\n'
		s += ' '.join(lr[::-1])
		s += '\n'
		s += ' '.join(['|'*len(c1) if c1==c2 else '*'*len(c1) for c1,c2 in zip(lq[::-1],lr[::-1])])
	else:
		s = ' '.join(lq)
		s += '\n'
		s = ' '.join(lr)
		s += '\n'
		s = ' '.join(['|'*len(c1) if c1==c2 else '*'*len(c1) for c1,c2 in zip(lq,lr)])
	n = sum(map(operator.eq, lq[::-1], lr[::-1]))

	return (s, n, len(lq))


def extend(q, r, bFwd):
	M = [[MIN for j in xrange(200)] for i in xrange(200)]
	I = [[MIN for j in xrange(200)] for i in xrange(200)]
	D = [[MIN for j in xrange(200)] for i in xrange(200)]

	initMat(M, I, D)
	s, n, m = align_affine(q, r, bFwd, M, I, D)
	print s
	print n
	print m


if __name__ == '__main__':
	qr = ['GCCCC', 'CCCCG', 'CCCGC', 'CCGCC', 'CGCCC', 'GCCCA', 'CCCAG', 'CCAGG', 'CAGGG', 'AGGGC', 'GGGCC', 'GGCCA', 'GCCAC', 'CCACC', 'CACCA', 'ACCAC', 'CCACG', 'CACGA', 'ACGAC', 'CGACA', 'GACAA', 'ACAAC', 'CAACT', 'AACTG', 'ACTGC', 'CTGCA', 'TGCAA', 'GCAAA', 'CAAAC', 'AAACA', 'AACAC', 'ACACC', 'CACCG', 'ACCGC', 'CCGCA', 'CGCAC', 'GCACC', 'CACCT', 'ACCTA', 'CCTAC', 'CTACC', 'TACCT', 'ACCTA', 'CCTAA', 'CTAAC', 'TAACA', 'AACAG', 'ACAGC', 'CAGCA', 'AGCAA', 'GCAAC', 'CAACC', 'AACCA', 'ACCAG', 'CCAGC', 'CAGCT', 'AGCTG', 'GCTGC', 'CTGCC']
	rr = ['GCCCC', 'CCCCG', 'CCCGC', 'CCGCC', 'CGCCC', 'GCCCA', 'CCCAG', 'CCAGG', 'CAGGG', 'AGGGC', 'GGGCC', 'GGCCA', 'GCCAC', 'CCACC', 'CACCA', 'ACCAC', 'CCACG', 'CACGA', 'ACGAC', 'CGACC', 'GACCA', 'ACCAC', 'CCACT', 'CACTG', 'ACTGC', 'CTGCA', 'TGCAA', 'GCAAA', 'CAAAC', 'AAACA', 'AACAC', 'ACACC', 'CACCG', 'ACCGC', 'CCGCA', 'CGCAC', 'GCACC', 'CACCT', 'ACCTG', 'CCTGC', 'CTGCC', 'TGCCT', 'GCCTA', 'CCTAC', 'CTACC', 'TACCA', 'ACCAG', 'CCAGC', 'CAGCC', 'AGCCA', 'GCCAC', 'CCACC', 'CACCA', 'ACCAG', 'CCAGC', 'CAGCT', 'AGCTG', 'GCTGC', 'CTGCC', 'TGCCT', 'GCCTA', 'CCTAC', 'CTACA', 'TACAA', 'ACAAA', 'CAAAG', 'AAAGA', 'AAGAC', 'AGACA']
	#extend(ql, rl, False)
	extend(qr, rr, True)
