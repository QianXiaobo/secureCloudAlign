#ifndef __util_hpp__
#define __util_hpp__

// utility functions
// yongzhao


#include <vector>
#include <string>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>



const char aRC[] = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
	'T', '0', 'G', '0', '0', '0', 'C', '0', '0', '0', '0', '0', '0',
	'N', '0', '0', '0', '0', '0', 'A', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', 'N', '0', '0', '0', '0', '0', '0',
	'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};

const uint32_t aC2I[] = {9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  0, 9,  1,
	9, 9, 9,  2, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  3,
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
	9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
	9, 9, 9, 9, 9, 9, 9, 9, 9};

const char aI2C[] = {'A', 'C', 'G', 'T'};

const uint64_t OFFSHIFT = 32;
const uint64_t SELFSHIFT = 31;
const uint64_t MERMASK = ((uint64_t)1<<30) - 1;


typedef std::vector<uint32_t> ITEM;
typedef std::vector<ITEM> REF;

void ToRC(std::string& s)
{
	uint32_t i = 0;
	uint32_t j = s.size() - 1;
	while (i < j)
	{
		std::swap(s[i], s[j]);
		s[i] = aRC[s[i]];
		s[j] = aRC[s[j]];

		++i;
		--j;
	}

	if (i == j)
	{
		s[i] = aRC[s[i]];
	}
}


// ref kmer 2 uint64_t
void TransRef(const std::string& sFw, uint32_t uOff, uint32_t uK1, uint32_t uK, uint64_t& uIdx, uint64_t& uMer)
{
	std::string sRC = sFw;
	ToRC(sRC);
	bool bSelf = sFw < sRC ? true : false;
	const std::string& s = bSelf ? sFw : sRC;

	// kmer to int
	uIdx = 0;
	for (uint32_t i = 0; i < uK1; ++i)
	{
		uIdx = (uIdx << 2) | aC2I[s[i]];
	}

	uMer = 0;
	for (uint32_t i = uK1; i < s.size(); ++i)
	{
		uMer = (uMer << 2) | aC2I[s[i]];
	}
	uMer |= (uint64_t)uOff << OFFSHIFT;
	if (true == bSelf)
	{
		uMer |= (uint64_t)1 << SELFSHIFT;
	}
}


// uint64_t 2 ref kmer
std::tuple<std::string, uint32_t> TransRef(uint64_t uIdx, uint64_t uMer, uint32_t uK1, uint32_t uK)
{
	std::string s;

	for (int uShift = (uK1-1)*2; uShift >= 0; uShift -= 2)
	{
		char c = aI2C[(uIdx>>uShift)&3];
		s += c;
	}

	for (int uShift = (uK-uK1-1)*2; uShift >= 0; uShift -= 2)
	{
		char c = aI2C[(uMer>>uShift)&3];
		s += c;
	}

	if (0 == ((uMer>>SELFSHIFT)&1))
	{
		ToRC(s);
	}

	return make_tuple(s, uMer>>OFFSHIFT);
}



// qry kmer 2 uint64_t
// assume uId is not greater than 2^22
// assume uOff is not greater than 2^10
void TransQry(bool bSelf, uint32_t uId, uint32_t uOff, uint64_t& uIdx, uint64_t& uMer)
{
	uMer |= (uint64_t)uId << (OFFSHIFT+10);
	uMer |= (uint64_t)uOff << OFFSHIFT;
	if (true == bSelf)
	{
		uMer |= (uint64_t)1 << SELFSHIFT;
	}
}


// uint64_t 2 qry kmer
// assume uId is not greater than 2^22
std::tuple<std::string, uint32_t, uint32_t> TransQry(uint64_t uIdx, uint64_t uMer, uint32_t uK1, uint32_t uK)
{
	std::string s;

	for (int uShift = (uK1-1)*2; uShift >= 0; uShift -= 2)
	{
		char c = aI2C[(uIdx>>uShift)&3];
		s += c;
	}

	for (int uShift = (uK-uK1-1)*2; uShift >= 0; uShift -= 2)
	{
		char c = aI2C[(uMer>>uShift)&3];
		s += c;
	}

	if (0 == ((uMer>>SELFSHIFT)&1))
	{
		ToRC(s);
	}

	return make_tuple(s, uMer>>(OFFSHIFT+10), (uMer>>OFFSHIFT)&((1<<10)-1));
}


void Kmer(const std::string& s, uint32_t uOff, uint32_t uK1, uint32_t uK, uint64_t& uIdx, uint64_t& uMer)
{
	// kmer to int
	uIdx = 0;
	for (uint32_t i = uOff; i < uOff+uK1; ++i)
	{
		uIdx = (uIdx << 2) | aC2I[s[i]];
	}

	uMer = 0;
	for (uint32_t i = uOff+uK1; i < uOff+uK; ++i)
	{
		uMer = (uMer << 2) | aC2I[s[i]];
	}
}

// append new char
void KmerF(const char c, uint32_t uK1, uint32_t uK, uint64_t& uIdx, uint64_t& uMer)
{
	uIdx = ((uIdx<<2)&((1<<(2*uK1))-1)) + ((uMer>>(2*(uK-uK1-1)))&3);
	uMer = ((uMer<<2)&((1<<(2*(uK-uK1)))-1)) + aC2I[c];
}

// insert new char
void KmerB(const char c, uint32_t uK1, uint32_t uK, uint64_t& uIdx, uint64_t& uMer)
{
	uMer = ((uMer>>2)|((uIdx&3)<<(2*(uK-uK1-1))));
	uIdx = ((uIdx>>2)|(aC2I[c]<<(2*(uK1-1))));
}


std::tuple<std::string, uint32_t, uint32_t> TransQry1(uint64_t uIdx, uint64_t uMer, uint32_t uK1, uint32_t uK)
{
	std::string s;

	for (int uShift = (uK1-1)*2; uShift >= 0; uShift -= 2)
	{
		char c = aI2C[(uIdx>>uShift)&3];
		s += c;
	}

	for (int uShift = (uK-uK1-1)*2; uShift >= 0; uShift -= 2)
	{
		char c = aI2C[(uMer>>uShift)&3];
		s += c;
	}

	return make_tuple(s, uMer>>(OFFSHIFT+10), (uMer>>OFFSHIFT)&((1<<10)-1));
}


void VerifyKmer1(uint64_t uIdx, uint64_t uMer, uint32_t uK1, uint32_t uK, std::string& s, uint32_t u)
{
	std::string sBack;
	uint32_t u1, u2;

	tie(sBack,u1,u2) = TransQry1(uIdx, uMer, uK1, uK);
	if (sBack != s.substr(u, uK))
	{
		std::cout << s << std::endl;
		std::cout << s.substr(u, uK) << std::endl;
		std::cout << sBack << std::endl;
		std::cerr << "wrong wrong wrong STR" << std::endl;
		exit(-2);
	}
}


void VerifyKmer(uint64_t uIdx, uint64_t uMer, uint32_t uK1, uint32_t uK, std::string& s, uint32_t uId, uint32_t i)
{
	std::string sBack;
	uint32_t u1, u2;

	tie(sBack,u1,u2) = TransQry(uIdx, uMer, uK1, uK);
	if (sBack != s.substr(i,uK))
	{
		std::cout << s << std::endl;
		std::cout << s.substr(i, uK) << std::endl;
		std::cout << sBack << std::endl;
		std::cout << uId << "\t" << i << std::endl;
		std::cout << ((uMer>>SELFSHIFT)&1) << std::endl;
		std::cerr << "wrong wrong wrong STR" << std::endl;
		exit(-2);
	}
	if (uId != u1)
	{
		std::cout << uId << "\t" << i << std::endl;
		std::cerr << uId << std::endl;
		std::cerr << u1 << std::endl;
		std::cerr << "wrong wrong wrong IDX" << std::endl;
		exit(-2);
	}
	if (i != u2)
	{
		std::cout << uId << "\t" << i << std::endl;
		std::cerr << i << std::endl;
		std::cerr << u2 << std::endl;
		std::cerr << "wrong wrong wrong POS" << std::endl;
		exit(-2);
	}
}


void CompKmer(uint64_t uI1, uint64_t uM1, uint64_t uI2, uint64_t uM2, uint64_t& uIdx, uint64_t& uMer, bool& bSelf)
{
	if (uI1 < uI2)
	{
		uIdx = uI1;
		uMer = uM1;
		bSelf = true;
	}
	else if (uI1 > uI2)
	{
		uIdx = uI2;
		uMer = uM2;
		bSelf = false;
	}
	else
	{
		if (uM1 <= uM2)
		{
			uIdx = uI1;
			uMer = uM1;
			bSelf = true;
		}
		else
		{
			uIdx = uI2;
			uMer = uM2;
			bSelf = false;
		}
	}
}
#endif
