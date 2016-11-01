#!/usr/bin/env python


if __name__ == '__main__':
	sBwa = set()
	#with open('matched', 'r') as f:
	with open('eval.sort', 'r') as f:
		for l in f:
			l = l.strip()
			sBwa.add(l.split()[0])

	sMine = set()
	with open('test.id.new', 'r') as f:
		for l in f:
			l = l.strip()
			sMine.add(l.split()[0][1:])

	print len(sBwa)
	print len(sMine)
	print len(sMine & sBwa)
	print sBwa - sMine
