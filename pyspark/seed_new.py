#!/usr/bin/env python


import sys
import os.path as op
import bisect as bs
import operator
import marshal as pk
import math
import hashlib
import bitarray as ba
from struct import unpack, pack, calcsize
import subprocess as sp

from pyspark import SparkContext, SparkConf, SparkFiles



sHdfsDir = 'hdfs:/user/yongzhao'


def calcBfPar(capacity, error_rate):
	num_slices = int(math.ceil(math.log(1 / error_rate, 2)))
	bits_per_slice = int(math.ceil(
		(2 * capacity * abs(math.log(error_rate))) /
		(num_slices * (math.log(2) ** 2))))
	num_bits = num_slices * bits_per_slice

	if bits_per_slice >= (1 << 31):
		fmt_code, chunk_size = 'Q', 8
	elif bits_per_slice >= (1 << 15):
		fmt_code, chunk_size = 'I', 4
	else:
		fmt_code, chunk_size = 'H', 2
	total_hash_bits = 8 * num_slices * chunk_size
	if total_hash_bits > 384:
		hashfn = hashlib.sha512
	elif total_hash_bits > 256:
		hashfn = hashlib.sha384
	elif total_hash_bits > 160:
		hashfn = hashlib.sha256
	elif total_hash_bits > 128:
		hashfn = hashlib.sha1
	else:
		hashfn = hashlib.md5

	fmt = fmt_code * (hashfn().digest_size // chunk_size)
	num_salts, extra = divmod(num_slices, len(fmt))
	if extra:
		num_salts += 1
	seeds = []
	for i in xrange(num_salts):
		seeds.append(hashfn(pack('I', i)).digest())

	bf = ba.bitarray(num_bits, endian='little')
	bf.setall(False)

	return bf, num_slices, bits_per_slice, hashfn, seeds, fmt


def add(bf, num_slices, bits_per_slice, hashfn, seeds, fmt, key):
	off = []
	for x in seeds:
		y = hashfn(x)
		y.update(key.encode('utf-8'))
		for unit in unpack(fmt, y.digest()):
			off.append(unit % bits_per_slice)

	nOffset = 0
	for k in off[:num_slices]:
		bf[nOffset + k] = True
		nOffset += bits_per_slice


def contain(bf, num_slices, bits_per_slice, hashfn, seeds, fmt, key):
	off = []
	for x in seeds:
		y = hashfn(x)
		y.update(key.encode('utf-8'))
		for unit in unpack(fmt, y.digest()):
			off.append(unit % bits_per_slice)

	nOffset = 0
	for k in off[:num_slices]:
		if not bf[nOffset + k]:
			return False
		nOffset += bits_per_slice
	return True


def count(idx, iters):
	bitarray = bcBitarray.value
	bits_per_slice = bcBitsPerSlice
	hashes = bcHashes(key)
	n = 0
	for it in iters:
		n += 1
	return (idx, n)


def build(idx, iters):
	capacity = 661512672
	error_rate = 0.10
	bf, num_slices, bits_per_slice, hashfn, seeds, fmt = calcBfPar(capacity, error_rate)

	#for t in iters:
		#add(bf, num_slices, bits_per_slice, hashfn, seeds, fmt, t[0])

	yield (None, bf.tostring())


def merge(bf1, bf2):
	x = ba.bitarray(endian='little')
	x.fromstring(bf1)
	y = ba.bitarray(endian='little')
	y.fromstring(bf2)
	return (x|y).tostring()


def seed(idx, iters):
	bf = ba.bitarray(endian='little')
	bf.fromstring(bcBitarray.value)
	num_slices = bcNumSlices.value
	bits_per_slice = bcBitsPerSlice.value
	hashfn = bcHashfn.value
	seeds = bcSeeds.value
	fmt = bcFmt.value

	for t in iters:
		off = []
		for x in seeds:
			y = hashfn(x)
			y.update(t[0].encode('utf-8'))
			for unit in unpack(fmt, y.digest()):
				off.append(unit % bits_per_slice)

		bHit = True
		nOffset = 0
		for k in off[:num_slices]:
			if not bf[nOffset + k]:
				bHit = False
				break
			nOffset += bits_per_slice
		if bHit == True:
			yield t

	


if __name__ == '__main__':

	sApp = 'spark'
	nPart = 14*10
	#sRef = op.join(sHdfsDir, 'hg38.fa.nb.enc.gzip')
	#sRef = op.join(sHdfsDir, 'chr21.fa.nb.enc.gzip')
	sInput = op.join(sHdfsDir, 'half.enc')
	sSeeds = op.join(sHdfsDir, 'seed.enc')

	# print default SparkConf
	sf = SparkConf()
	print sf.toDebugString()
	sc = SparkContext(appName=sApp)

	#capacity = 661512672
	#error_rate = 0.1
	#bf, num_slices, bits_per_slice, hashfn, seeds, fmt = calcBfPar(capacity, error_rate)

	#with open('/hddscratch/half.enc.bf', 'r') as fOut:
		#bf = ba.bitarray(endian='little')
		#bf.fromstring(pk.load(fOut))

	#bcBitarray = sc.broadcast(bf.tostring())
	#bcNumSlices = sc.broadcast(num_slices)
	#bcBitsPerSlice = sc.broadcast(bits_per_slice)
	#bcHashfn = sc.broadcast(hashfn)
	#bcSeeds = sc.broadcast(seeds)
	#bcFmt = sc.broadcast(fmt)

	#rddRef = sc.textFile(sRef, use_unicode=False).map(lambda x: x.split('$',1))
	#rddPart = rddRef.mapPartitionsWithIndex(seed, preservesPartitioning=True)

	#rddInput = sc.textFile(sInput, use_unicode=False).map(lambda x: x.split('$',1))

	#rddRes = rddPart.join(rddInput)

	#rddRes = rddRes.map(lambda x: (int(x[1][1].split('$')[1]), x[0]+'$'+x[1][0]+'$'+x[1][1]))
	#rddRes = rddRes.groupByKey().mapValues(list)

	#rddRes = rddRes.filter(lambda x: len(x[1])>40).sortByKey()
	#rddRes = rddRes.map(lambda x: '$'.join(x[1]))

	#sp.call('hdfs dfs -rm -r '+sSeeds, shell=True)
	##rddRes.saveAsTextFile(sSeeds) # use this for performance measurement
	#codec = "org.apache.hadoop.io.compress.GzipCodec"
	#rddRes.coalesce(1, shuffle=True).saveAsTextFile(sSeeds, codec)

	rddInput = sc.textFile(sInput, use_unicode=False).coalesce(nPart).map(lambda x: x.split('$',1))
	rddBf = rddInput.mapPartitionsWithIndex(build, preservesPartitioning=True)
	rddBf = rddBf.reduceByKey(merge)
	print rddBf.count()
	bf = rddBf.collect()
	print len(bf[0][1])

	sc.stop()
