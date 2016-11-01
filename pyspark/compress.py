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


def count(idx, iters):
	n = 0
	for it in iters:
		n += 1
	return (idx, n)


if __name__ == "__main__":

	sApp = 'spark'
	#sName = 'chr21.fa.nb.enc'
	sName = 'hg38.fa.nb.enc'
	sNbInput = op.join(sHdfsDir, sName)
	sNbComp = op.join(sHdfsDir, sName+'.gzip')

	# print default SparkConf
	sf = SparkConf()
	print sf.toDebugString()
	sc = SparkContext(appName=sApp)

	sp.call('hdfs dfs -rm -r '+sNbComp, shell=True)
	rdd = sc.textFile(sNbInput, use_unicode=False).coalesce(3000)
	codec = "org.apache.hadoop.io.compress.GzipCodec"
	rdd.saveAsTextFile(sNbComp, codec)

	sc.stop()
