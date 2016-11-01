#!/usr/bin/env python

# test seeding code


import sys
import operator


def correct():
	with open('chr21.fa.nb', 'r') as fNb:
		with open('first1M.fa.ker', 'r') as fRead:
			lRes = []
			mNb = {}
			n = 0
			for l in fNb:
				r = l.strip().split('$', 1)
				mNb[r[0]] = r[1]

			for l in fRead:
				r = l.strip().split('$')
				for x in r[1:]:
					if x in mNb:
						n += 1
						#print x, r[0], mNb[x]
						lRes.append((x, r[0], mNb[x]))

			lRes.sort(key=operator.itemgetter(0))
			for x in lRes:
				print x
			print n


def test():
	with open('chr21.fa.nb', 'r') as fNb:
		r = fNb.next().split()
		with open('first1M.fa.ker', 'r') as fRead:
			for l in fRead:
				k = l.strip().split()


if __name__ == '__main__':
