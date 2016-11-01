#!/usr/bin/env python
## extend


import sys
import math
import bisect as bs
import subprocess as sp
import os
import os.path as op
import cPickle as pk
import itertools as it

from pyspark import SparkContext, SparkConf, SparkFiles

from align import align_affine, initMat

N = 100000
MIN = -float("inf")

sHdfsDir = 'hdfs:/user/yongzhao'


def extend(idx, iters):
	M = [[MIN for j in xrange(200)] for i in xrange(200)]
	X = [[MIN for j in xrange(200)] for i in xrange(200)]
	Y = [[MIN for j in xrange(200)] for i in xrange(200)]
	initMat(X, Y, M)

	for t in iters:
		sId,sLeft,sRight,sRl,sRr = t.split('$')
		left = sLeft.split()
		right = sRight.split()
		rl = sRl.split()
		rr = sRr.split()

		s1,n1,m1 = align_affine(left, rl[:(len(left)+10)], False, X, Y, M) # left
		s2,n2,m2 = align_affine(right, rr[:len(right)+10], True, X, Y, M) # right

		sRes = ''
		#if n1+n2 >= m1+m2-20:
		if n1+n2 >= 40:
			#sRes += sId + ' aligned\n'
			sRes += str(n1) + '\t' + str(m1) + '\t' + str(n2) + '\t' + str(m2) + '\n'
			sRes += ' '.join(left) + '\n'
			sRes += ' '.join(rl[:(len(left)+10)]) + '\n'
			sRes += ' '.join(right) + '\n'
			sRes += ' '.join(rr[:len(right)+10]) + '\n'
			sRes += s1 + '\n'
			sRes += s2 + '\n\n'
			yield (sId, sRes)


if __name__ == '__main__':

	sApp = 'spark'
	nPart = 38*14
	#sInput = op.join(sHdfsDir, 'extend.mer')
	sInput = op.join(sHdfsDir, 'half.extend.mer')

	# print default SparkConf
	sf = SparkConf()
	print sf.toDebugString()
	sc = SparkContext(appName=sApp)

	sc.addPyFile('align.py')

	rdd = sc.textFile(sInput, nPart)

	rdd = rdd.mapPartitionsWithIndex(extend)
	rdd = rdd.reduceByKey()
	rdd.saveAsText(op.join(sHdfsDir,'extend_res'))
	lRes = rdd.collect()
	#with open('extend.txt', 'w') as f:
		#for x in lRes:
			#print >> f, x
	#print len(lRes)

	sc.stop()
