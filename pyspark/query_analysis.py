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
#from pybloom import BloomFilter




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
	#bf = ba.BitArray(length=num_bits)
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

	print off[:num_slices]
	nOffset = 0
	for k in off[:num_slices]:
		if not bf[nOffset + k]:
			return False
		nOffset += bits_per_slice
	return True


if __name__ == '__main__':
	#capacity = 64000000
	capacity = 1323025344/4
	error_rate = 0.05

	bf, num_slices, bits_per_slice, hashfn, seeds, fmt = calcBfPar(capacity, error_rate)
	print len(bf.tostring())
	print hashfn
	sys.exit()

	#with open('first1M.fa.nb.enc', 'r') as f:
		#n = 0
		#for l in f:
			#sKey = l.split('$',1)[0]
			#add(bf, num_slices, bits_per_slice, hashfn, seeds, fmt, sKey)
			#n += 1
		#with open('bf.pk', 'w') as fOut:
			#pk.dump(bf.tostring(), fOut)

	with open('bf.pk', 'r') as fOut:
		bf = ba.bitarray(endian='little')
		bf.fromstring(pk.load(fOut))
		sKey = 'C7E54606DFADBB5F832549617FE6358A'
		print contain(bf, num_slices, bits_per_slice, hashfn, seeds, fmt, sKey)
		#with open('first1M.fa.nb.enc', 'r') as f:
			#for l in f:
				#sKey = l.split('$',1)[0]
				#assert True == contain(bf, num_slices, bits_per_slice, hashfn, seeds, fmt, sKey)
