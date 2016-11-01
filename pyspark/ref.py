#!/usr/bin/env python
## preprocess references


import sys
import math
import bisect as bs
import subprocess as sp
import os
import os.path as op
import cPickle as pk
import shutil

from pyspark import SparkContext, SparkConf


sHdfsDir = 'hdfs:/user/yongzhao'


def split2KV(x):
	r = x.strip().split('$', 1)
	return (r[0], r[1])


def partFunc(x):
	idx = bs.bisect_right(ptter.value, x)
	return idx


def genPtter(rdd, fFrac, nPart):
# partition pilots
	rddSampled = rdd.sample(False, fFrac)
	rddSampled = rddSampled.map(lambda x: x[0])
	lSampled = rddSampled.collect()
	lSampled.sort()

	nSampled = len(lSampled)
	nStep = int(math.floor(nSampled/nPart)) + 1
	lPtter = [lSampled[x] for x in range(nStep, nSampled, nStep)]

	return lPtter


def count(idx, iters):
	n = 0
	for it in iters:
		n += 1
	return (idx, n)


if __name__ == "__main__":

	sApp = 'spark'
	nPart = 38*14*4*4
	sDir = op.join(sHdfsDir, 'hg38.fa.nb.enc.gzip')
	sPtter = op.join(sHdfsDir, 'ptter')
	codec = "org.apache.hadoop.io.compress.GzipCodec"

	# print default SparkConf
	sf = SparkConf()
	print sf.toDebugString()
	sc = SparkContext(appName=sApp)

	rdd = sc.textFile(sDir, use_unicode=False)
	rdd = rdd.map(split2KV)

	#lPtter = genPtter(rdd, 0.001, nPart)
	#sc.parallelize(lPtter).saveAsTextFile(sPtter)
	ptter = sc.broadcast(sc.textFile(sPtter, use_unicode=False).collect())

	nTime = 4
	nOne = nPart / nTime 
	lIndex = [i*nOne for i in xrange(1,nTime)]
	s0 = ptter.value[lIndex[0]]
	s1 = ptter.value[lIndex[1]]
	s2 = ptter.value[lIndex[2]]

	#print ptter.value[lIndex[0]], ptter.value[:lIndex[0]]
	#print ptter.value[lIndex[0]], ptter.value[lIndex[1]], ptter.value[lIndex[0]:lIndex[1]]
	#print ptter.value[lIndex[1]], ptter.value[lIndex[2]], ptter.value[lIndex[1]:lIndex[2]]
	#print ptter.value[lIndex[2]], ptter.value[lIndex[2]:]

	for i in xrange(4):
		sp.call('hdfs dfs -rm -r '+op.join(sHdfsDir, 'nb.'+str(i)), shell=True)

	rddPart = rdd.filter(lambda x: x[0]<s0)
	rddPart = rddPart.repartitionAndSortWithinPartitions(nOne, partFunc)
	rddPart = rddPart.map(lambda x: x[0]+'$'+x[1])
	rddPart = rddPart.saveAsTextFile(op.join(sHdfsDir, 'nb.0'), codec)

	rddPart = rdd.filter(lambda x: x[0]>=s0 and x[0]<s1)
	rddPart = rddPart.repartitionAndSortWithinPartitions(nOne, partFunc)
	rddPart = rddPart.map(lambda x: x[0]+'$'+x[1])
	rddPart = rddPart.saveAsTextFile(op.join(sHdfsDir, 'nb.1'), codec)

	rddPart = rdd.filter(lambda x: x[0]>=s1 and x[0]<s2)
	rddPart = rddPart.repartitionAndSortWithinPartitions(nOne, partFunc)
	rddPart = rddPart.map(lambda x: x[0]+'$'+x[1])
	rddPart = rddPart.saveAsTextFile(op.join(sHdfsDir, 'nb.2'), codec)

	rddPart = rdd.filter(lambda x: x[0]>=s2)
	rddPart = rddPart.repartitionAndSortWithinPartitions(nOne, partFunc)
	rddPart = rddPart.map(lambda x: x[0]+'$'+x[1])
	rddPart = rddPart.saveAsTextFile(op.join(sHdfsDir, 'nb.3'), codec)

	sc.stop()
