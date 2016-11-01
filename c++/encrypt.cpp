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


int nK = 27;

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


void EncryptRef(string sKey, string sInput)
{
	char aH2C[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};

	// hash
	byte aKeyHash[AES::DEFAULT_KEYLENGTH];
	FileSource((sKey+"/hash").c_str(), true, new ArraySink(aKeyHash, sizeof(aKeyHash)));

	VMAC<AES> hash;
	cerr << hash.StaticAlgorithmName() << endl;
	cerr << "Digest Size: " << hash.DigestSize() << endl;
	cerr << "AES Size: " << AES::BLOCKSIZE << endl;

	byte aIv[AES::BLOCKSIZE];
	fill(aIv, aIv+AES::BLOCKSIZE, 0);
	cerr << "Hash IV: " << string(aIv, aIv+sizeof(aIv)) << endl;

	hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);

	// aes
	byte aKeyEnc[AES::DEFAULT_KEYLENGTH];
	FileSource((sKey+"/encrypt").c_str(), true, new ArraySink(aKeyEnc, sizeof(aKeyEnc)));

	CTR_Mode<AES>::Encryption enc;
	enc.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
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
	byte a1[AES::BLOCKSIZE];
	byte a2[AES::BLOCKSIZE];
	byte a3[AES::BLOCKSIZE*2048];
	byte a4[AES::BLOCKSIZE*2048];
	string sOut, sHash;
	string sCipher, sRecovered;
	while(getline(iff, sLine))
	{
		if (sLine.empty()) continue;
		nPrint += 1;

		hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);
		hash.CalculateDigest(a1, (byte*)sLine.c_str(), nK);
		sHash.clear();
		for (int i = 0; i < AES::BLOCKSIZE; ++i)
		{
			sPrint.push_back(aH2C[a1[i]>>4]);
			sPrint.push_back(aH2C[a1[i]&15]);
			sHash.push_back(aH2C[a1[i]>>4]);
			sHash.push_back(aH2C[a1[i]&15]);
		}
		sPrint.push_back('$');

		hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);
		hash.CalculateDigest(a2, (byte*)sLine.c_str(), nK);
		sOut.clear();
		StringSource(a2, sizeof(a2), true, 
				new HexEncoder(new StringSink(sOut)));
		assert (sHash == sOut);

		// encrypt
		enc.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
		enc.ProcessData(a3, (byte*)sLine.data(), sLine.size());
		for (int i = 0; i < sLine.size(); ++i)
		{
			sPrint.push_back(aH2C[a3[i]>>4]);
			sPrint.push_back(aH2C[a3[i]&15]);
		}
		sPrint.push_back('\n');

		// decrypt
		dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
		dec.ProcessData(a4, a3, sLine.size());
		string sRecovered(a4, a4+sLine.size());

		assert (sRecovered==sLine);

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


void EncryptSeq(string sKey, string sInput)
{
	char aH2C[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};

	// hash
	byte aKeyHash[AES::DEFAULT_KEYLENGTH];
	FileSource((sKey+"/hash").c_str(), true, new ArraySink(aKeyHash, sizeof(aKeyHash)));

	VMAC<AES> hash;
	cerr << hash.StaticAlgorithmName() << endl;
	cerr << "Digest Size: " << hash.DigestSize() << endl;
	cerr << "AES Size: " << AES::BLOCKSIZE << endl;

	byte aIv[AES::BLOCKSIZE];
	fill(aIv, aIv+AES::BLOCKSIZE, 0);
	cerr << "Hash IV: " << string(aIv, aIv+sizeof(aIv)) << endl;

	hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);

	// aes
	byte aKeyEnc[AES::DEFAULT_KEYLENGTH];
	FileSource((sKey+"/encrypt").c_str(), true, new ArraySink(aKeyEnc, sizeof(aKeyEnc)));

	CTR_Mode<AES>::Encryption enc;
	enc.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
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
	string sId, sSeq;
	byte a1[AES::BLOCKSIZE];
	byte a2[AES::BLOCKSIZE];
	byte a3[AES::BLOCKSIZE*2048];
	byte a4[AES::BLOCKSIZE*2048];
	string sOut, sHash;
	string sCipher, sRecovered;
	while(getline(iff, sId))
	{
		getline(iff, sSeq);
		nPrint += 1;

		for(size_t j = 0; j < sSeq.size()-nK+1; ++j)
		{
			hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);
			hash.CalculateDigest(a1, (byte*)sSeq.c_str()+j, nK);
			sHash.clear();
			for (int i = 0; i < AES::BLOCKSIZE; ++i)
			{
				sPrint.push_back(aH2C[a1[i]>>4]);
				sPrint.push_back(aH2C[a1[i]&15]);
				sHash.push_back(aH2C[a1[i]>>4]);
				sHash.push_back(aH2C[a1[i]&15]);
			}
			sPrint.push_back('$');

			hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);
			hash.CalculateDigest(a2, (byte*)sSeq.c_str()+j, nK);
			sOut.clear();
			StringSource(a2, sizeof(a2), true, 
					new HexEncoder(new StringSink(sOut)));
			assert (sHash == sOut);

			// encrypt
			string sInfo(sSeq.begin()+j, sSeq.begin()+j+nK);
			sInfo.push_back('$');
			sInfo.append(sId);
			sInfo.push_back('|');
			sInfo.append(to_string(j));

			enc.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
			enc.ProcessData(a3, (byte*)sInfo.data(), sInfo.size());
			for (int i = 0; i < sInfo.size(); ++i)
			{
				sPrint.push_back(aH2C[a3[i]>>4]);
				sPrint.push_back(aH2C[a3[i]&15]);
			}
			sPrint.push_back('$');
			sPrint.append(to_string(nPrint));
			sPrint.push_back('\n');

			// decrypt
			dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
			dec.ProcessData(a4, a3, sInfo.size());
			string sRecovered(a4, a4+sInfo.size());

			assert (sRecovered==sInfo);
		}

		if (sPrint.size() >= sPrint.capacity()*0.8)
		{
			cout << sPrint;
			sPrint.clear();
		}

		if (nPrint%100000 == 0)
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
	if (argc == 2)
	{
		// argv[1]: key file
		CreateKey(string(argv[1])+"/hash");
		CreateKey(string(argv[1])+"/encrypt");
	}
	else if (argc == 4)
	{
		// argv[1]: key directory
		// argv[2]: seq file
		if (string(argv[3]) == "ref")
		{
			EncryptRef(argv[1], argv[2]);
		}
		else if (string(argv[3]) == "seq")
		{
			EncryptSeq(argv[1], argv[2]);
		}
		else
		{
			cerr << "usage ./prog key-dir" << "\tor\t" << "./prog key-dir seq-file seq|ref" << endl;
			exit(-1);
		}
	}
	else
	{
		cerr << "usage ./prog key-dir" << "\tor\t" << "./prog key-dir seq-file seq|ref" << endl;
		exit(-1);
	}

	return 0;
}
