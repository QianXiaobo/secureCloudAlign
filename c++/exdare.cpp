// generate m-mer sequences for seeds
// yongzhao


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <chrono>

#include <boost/archive/binary_iarchive.hpp>
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
using namespace std::chrono;
using namespace boost;
using namespace CryptoPP;

typedef high_resolution_clock TClock;


int main(int argc, char** argv)
{
	// argv[1]: refs
	// argv[2]: seqs
	// argv[3]: seeds
	// argv[4]: keys
	if (argc != 5)
	{
		cerr << "./prog ref-file seq-file seed-file key-dir" << endl;
		exit(-1);
	}

	int nK = 27;
	int nM = 5;
	int nReadLen = 100;
	auto t1 = TClock::now();

	// read reference
	string sRefFile = argv[1];
	vector<string> vId;
	vector<uint64_t> vPos;
	string sFw, sBw;
	ifstream iffRef(sRefFile+".all");
	if (iffRef.good())
	{
		archive::binary_iarchive ia(iffRef);
		ia >> vId;
		ia >> vPos;
		ia >> sFw;
		ia >> sBw;
		iffRef.close();
	}
	iffRef.close();
	//cout << vId.size() << endl;
	//copy(vId.begin(), vId.end(), ostream_iterator<string>(cout, "\n"));
	unordered_map <string, size_t> mId2N;
	for (size_t i = 0; i < vId.size(); ++i) mId2N[vId[i]] = i;

	// hash
	byte aKeyHash[AES::DEFAULT_KEYLENGTH];
	FileSource((string(argv[4])+"/hash").c_str(), true, new ArraySink(aKeyHash, sizeof(aKeyHash)));

	VMAC<AES> hash;
	byte aIv[AES::BLOCKSIZE];
	fill(aIv, aIv+AES::BLOCKSIZE, 0);
	hash.SetKeyWithIV(aKeyHash, sizeof(aKeyHash), aIv);

	// aes
	byte aKeyEnc[AES::DEFAULT_KEYLENGTH];
	FileSource((string(argv[4])+"/encrypt").c_str(), true, new ArraySink(aKeyEnc, sizeof(aKeyEnc)));

	CTR_Mode<AES>::Encryption enc;
	enc.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);
	CTR_Mode<AES>::Decryption dec;
	dec.SetKeyWithIV(aKeyEnc, sizeof(aKeyEnc), aIv);

	string sPrint; sPrint.reserve(1024*1024*256);
	string sId, sSeq;
	string sSeed;
	ifstream iffSeq(argv[2]);
	ifstream iffSeed(argv[3]);
	while (getline(iffSeed, sSeed))
	{
		vector<string> vSplits;
		split(vSplits, sSeed, is_any_of(" $"));
		//assert (vSplits[0] == vSplits[3]);
		assert (vSplits[0] == vSplits[4]); // CAUTIONS: remove

		while(getline(iffSeq,sId)
				&& getline(iffSeq,sSeq)
				//&& sId != vSplits[1])
				&& sId.find(vSplits[1])==string::npos)  // CAUTIONS: remove
		{}

		//cout << sId << endl;
		//for (auto x : vSplits) cout << "\t" << x;
		//cout << endl;

		//size_t nQryOff = stoi(vSplits[2]);
		size_t nQryOff = stoi(vSplits[3]);
		vector<string> vSeed;
		//for (size_t i = 4; i < vSplits.size(); ++i)
		for (size_t i = 5; i < vSplits.size(); ++i)
		{
			split(vSeed, vSplits[i], is_any_of("|"));
			size_t j = mId2N[vSeed[0]];

			size_t nLeft = vPos[j-1]; 
			size_t nRight = vPos[j] - 1; 
			size_t nInOff = stol(vSeed[1]);
			size_t n1 = max(nLeft, nLeft+nInOff-nReadLen+nK);
			size_t n2 = min(nRight, nLeft+nInOff+nReadLen);

			size_t n3 = nLeft + nInOff;
			while (n3>n1 && sFw[n3-1]!='N') --n3;
			size_t n4 = nLeft + nInOff + nK;
			while (n4<n2 && sFw[n4]!='N') ++n4;
			n1 = n3;
			n2 = n4;
			if (n2-n1 < nReadLen) continue;

			string sss = sFw.substr(n1, n2-n1);
			assert (sss.find('N') == string::npos);
			int nS = nLeft + nInOff - n1;

			if ("0" != vSeed[2])
			{
				ToRC(sss);
				nS = n2 -n1- nS - nK;
			}

			size_t diff = 0;
			auto it1 = sSeq.begin() + nQryOff;
			auto it2 = sss.begin() + nS;
			for (; it1 < sSeq.begin()+nQryOff+nK; ++it1,++it2) if (*it1!=*it2) ++diff;
			assert (diff <= 1);

			sPrint.append(sId);
			sPrint.push_back('$');
			//sPrint.push_back('\n');
			//sPrint.append(sSeq);
			//sPrint.push_back('\n');
			//sPrint.append(sss);
			//sPrint.push_back('\n');
			for (int k = nQryOff-nM; k > -1; --k)
			{
				sPrint.append(sSeq.begin()+k, sSeq.begin()+k+nM);
				sPrint.push_back(' ');
			}
			sPrint.push_back('$');
			//sPrint.push_back('\n');
			for (int k = nQryOff+nK; k < sSeq.size()-nM+1; ++k)
			{
				sPrint.append(sSeq.begin()+k, sSeq.begin()+k+nM);
				sPrint.push_back(' ');
			}
			sPrint.push_back('$');
			//sPrint.push_back('\n');

			for (int k = nS-nM; k > -1; --k)
			{
				sPrint.append(sss.begin()+k, sss.begin()+k+nM);
				sPrint.push_back(' ');
			}
			sPrint.push_back('$');
			//sPrint.push_back('\n');
			for (int k = nS+nK; k < sss.size()-nM+1; ++k)
			{
				sPrint.append(sss.begin()+k, sss.begin()+k+nM);
				sPrint.push_back(' ');
			}
			sPrint.push_back('\n');
		}

		if (sPrint.size() >= sPrint.capacity()*0.8)
		{
			cout << sPrint;
			sPrint.clear();
		}
	}

	if (sPrint.size() > 0)
	{
		cout << sPrint;
		sPrint.clear();
	}

	iffSeq.close();
	iffSeed.close();

	return 0;
}
