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
	nPart = 38*14*4
	#sRef = op.join(sHdfsDir, 'hg38.fa.nb.enc.gzip')
	sRef = op.join(sHdfsDir, 'chr21.fa.nb.enc.gzip')
	sInput = op.join(sHdfsDir, 'first1M.fa.nb.enc')
	sSeeds = op.join(sHdfsDir, 'seed.enc')

	# print default SparkConf
	sf = SparkConf()
	print sf.toDebugString()
	sc = SparkContext(appName=sApp)

	rdd = sc.textFile(op.join(sHdfsDir,'half.enc'), use_unicode=False)
	nTotal = rdd.count()

	sc.stop()

	print nTotal
