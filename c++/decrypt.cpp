// encrypt kmers and its sequences
// yongzhao


#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <algorithm>
#include <bitset>
#include <unordered_set>

#include <boost/serialization/vector.hpp>
#include <boost/algorithm/string.hpp>

#include "cryptopp/aes.h"
#include "cryptopp/sha.h"
#include "cryptopp/hex.h"
#include "cryptopp/vmac.h"
#include "cryptopp/osrng.h"
#include "cryptopp/filters.h"
#include "cryptopp/modes.h"
#include "cryptopp/files.h"

#include "util.hpp"

using namespace std;
using namespace boost;
using namespace CryptoPP;



void PrintInfo()
{
	// print AES info
	cerr << "AES parameters:" << endl;
	cerr << "Algorithm name:\t" << AES::StaticAlgorithmName() << endl;

	cerr << "Block size:\t" << AES::BLOCKSIZE*8 << "bits,\t" << AES::BLOCKSIZE << "bytes" << endl;
	cerr << "Mini key length:\t" << AES::MIN_KEYLENGTH*8 << "bits,\t" << AES::MIN_KEYLENGTH << "bytes" << endl;
	cerr << "Max key length:\t" << AES::MAX_KEYLENGTH*8 << "bits,\t" << AES::MAX_KEYLENGTH << "bytes" << endl;
}


void CreateKey(string sKey)
{
	AutoSeededRandomPool rng;
	byte aKey[AES::DEFAULT_KEYLENGTH];
	rng.GenerateBlock(aKey, sizeof(aKey));
	ArraySource(aKey, sizeof(aKey), true, new FileSink(sKey.c_str()));
}


void Decrypt(string sKey, string sInput)
{
	char aH2C[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
	char aC2H[256] = {};
	fill(aC2H, aC2H+256, 0);
	aC2H['0'] = 0;
	aC2H['1'] = 1;
	aC2H['2'] = 2;
	aC2H['3'] = 3;
	aC2H['4'] = 4;
	aC2H['5'] = 5;
	aC2H['6'] = 6;
	aC2H['7'] = 7;
	aC2H['8'] = 8;
	aC2H['9'] = 9;
	aC2H['A'] = 10;
	aC2H['B'] = 11;
	aC2H['C'] = 12;
	aC2H['D'] = 13;
	aC2H['E'] = 14;
	aC2H['F'] = 15;

	// aes
	byte aKeyEnc[AES::DEFAULT_KEYLENGTH];
	FileSource((sKey+"/encrypt").c_str(), true, new ArraySink(aKeyEnc, sizeof(aKeyEnc)));
	byte aIv[AES::BLOCKSIZE];
	fill(aIv, aIv+AES::BLOCKSIZE, 0);

	CTR_Mode<AES>::Decryption dec;
	dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);

	ifstream iff(sInput.c_str());
	if (!iff.good())
	{
		cerr << "can not open the input file" << endl;
		exit(-2);
	}

	int nPrint = 0;
	string sPrint; sPrint.reserve(1024*1024*256);
	string sLine;
	byte a1[AES::BLOCKSIZE*2048];
	byte a2[AES::BLOCKSIZE*2048];
	string sOut, sHash;
	string sCipher, sRecovered;
	while(getline(iff, sLine))
	{
		if (sLine.empty()) continue;
		nPrint += 1;

		int nMostLeftOff = 1000;
		string sMostLeftOff;
		string sId;
		vector<string> vSplits;
		split(vSplits, sLine, is_space());
		for (auto x : vSplits)
		{
			size_t n1 = x.find('$');
			size_t n2 = x.find('$', n1+1);
			size_t n3 = x.find('$', n2+1);

			string s1 = x.substr(n2+1, n3-n2-1);
			string sDecode;
			StringSource(s1, true,
					new HexDecoder(
						new StringSink(sDecode)));
			dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
			dec.ProcessData(a2, (byte*)sDecode.data(), sDecode.size());
			sRecovered.assign(a2, a2+sDecode.size());

			for (size_t i = 0; i < s1.size(); i += 2)
			{
				a1[i/2] = (aC2H[s1[i]]<<4) | (aC2H[s1[i+1]]);
			}
			dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
			dec.ProcessData(a2, a1, s1.size()/2);
			sOut.assign(a2, a2+sDecode.size());
			//sPrint.append(sOut);
			//sPrint.push_back('\n');

			assert (sRecovered==sOut);

			size_t m = sOut.find('|');
			int n = stoi(sOut.substr(m+1, sOut.size()-m-1));
			if (n < nMostLeftOff)
			{
				nMostLeftOff = n;
				sMostLeftOff = x.substr(n1+1, n2-n1-1);
				sId = sOut.substr(0, m);
			}
		}

		for (size_t i = 0; i < sMostLeftOff.size(); i += 2)
		{
			a1[i/2] = (aC2H[sMostLeftOff[i]]<<4) | (aC2H[sMostLeftOff[i+1]]);
		}
		dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
		dec.ProcessData(a2, a1, sMostLeftOff.size()/2);
		sOut.assign(a2, a2+sMostLeftOff.size()/2);
		sPrint.append(sId);
		sPrint.append(" ");
		sPrint.append(to_string(nMostLeftOff));
		sPrint.append(" ");
		sPrint.append(sOut);
		sPrint.push_back('\n');

		if (sPrint.size() >= sPrint.capacity()*0.8)
		{
			cout << sPrint;
			sPrint.clear();
		}

		if (nPrint%1000000 == 0)
		{
			cerr << nPrint << endl;
		}
	}

	if (sPrint.size() > 0)
	{
		cout << sPrint;
		sPrint.clear();
	}

	iff.close();
}


int main(int argc, char** argv)
{
	if (argc == 3)
	{
		// argv[1]: key directory
		// argv[2]: seq file
		Decrypt(argv[1], argv[2]);
	}
	else
	{
		cerr << "./prog key-dir seed-file" << endl;
		exit(-1);
	}

	return 0;
}
