// seed kmers of reference sequences
// yongzhao


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>
#include <algorithm>
#include <bitset>
#include <chrono>
#include <tuple>
#include <utility>
#include <iterator>

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/bit_vectors.hpp>

#include "util.hpp"

using namespace std;
using namespace std::chrono;
using namespace boost;
using namespace sdsl;

typedef high_resolution_clock TClock;
typedef std::tuple<size_t, size_t, size_t> TOcc;
typedef std::tuple<size_t, size_t> TRange;


boost::mutex MT;
boost::mutex MT_BV;
boost::mutex MT_OUT;
size_t IDX = 0;


int read(string sIn, vector<string>& vId, vector<uint64_t>& vPos, string& sFw, string& sBw)
{
	vector<char> vC(128, 0);
	vC['A'] = 'T';
	vC['C'] = 'G';
	vC['G'] = 'C';
	vC['T'] = 'A';
	vC['N'] = 'N';

	ifstream iff(sIn.c_str());
	if (!iff.good())
	{
		cerr << "ERROR: cannot open the file" << endl;
		return -1;
	}

	iff.seekg(0, iff.end);
	size_t nEst = iff.tellg();
	iff.seekg(0, iff.beg);
	sFw.clear(); sFw.reserve(nEst);
	sBw.clear(); sBw.reserve(nEst);
	vId.clear();
	vId.push_back("");
	vPos.clear();
	vPos.push_back(0);

	string id;
	string seq;
	string sLine;
	while (getline(iff, sLine))
	{
		if (sLine.empty())
		{
			continue;
		}
		else if ('>' == sLine[0])
		{
			if (!seq.empty())
			{
				to_upper(seq);
				sFw.append(seq.begin(), seq.end()); sFw.append("X");
				sBw.append(seq.rbegin(), seq.rend()); sBw.append("X");
				vId.push_back(id);
				vPos.push_back(sFw.length());

				ToRC(seq);
				sFw.append(seq.begin(), seq.end()); sFw.append("X");
				sBw.append(seq.rbegin(), seq.rend()); sBw.append("X");
				vId.push_back(id+"_rc");
				vPos.push_back(sFw.length());
			}
			id = sLine;
			seq.clear();
		}
		else
		{
			trim(sLine);
			seq += sLine;
		}
	}

	to_upper(seq);
	sFw.append(seq.begin(), seq.end()); sFw.append("X");
	sBw.append(seq.rbegin(), seq.rend()); sBw.append("X");
	vId.push_back(id);
	vPos.push_back(sFw.length());

	ToRC(seq);
	sFw.append(seq.begin(), seq.end()); sFw.append("X");
	sBw.append(seq.rbegin(), seq.rend()); sBw.append("X");
	vId.push_back(id+"_rc");
	vPos.push_back(sFw.length());
	iff.close();

	return 0;
}


void verify(const string& sSeq, size_t i, size_t j, size_t flag, int nK, int nE)
{
	string ss(sSeq.substr(j, nK));
	if (flag != 0)
	{
		ToRC(ss);
	}

	int nn = 0;
	for (size_t o = 0; o < nK; ++o)
	{
		if (sSeq[i+o]!=ss[o])
		{
			++nn;
		}
	}
	assert (nn == nE);
}


template<typename TStringIter>
int matche(const csa_wt<>& fm_index, size_t nL, size_t nR, TStringIter beg, TStringIter end, size_t nE, char* a, int nMaxOcc, vector<size_t>& vSeed)
{
	if (0 == nE)
	{
		backward_search(fm_index, nL, nR, beg, end, nL, nR);
		if (nL > nR) return INT_MAX;

		size_t nOcc = nR - nL + 1;
		if (vSeed.size()+nOcc > nMaxOcc) return vSeed.size()+nOcc;

		for (size_t l = 0; l < nOcc; ++l)
		{
			vSeed.push_back(fm_index[nL+l]);
		}
		return vSeed.size();
	}

	TStringIter it = end;
	while (beg < it)
	{
		--it;
		for (size_t j = 0; j < 4; ++j)
		{
			if (*it == a[j]) continue;

			size_t n1 = nL;
			size_t n2 = nR;
			backward_search(fm_index, n1, n2, a[j], n1, n2);
			if (n1 > n2) continue;
			int nRes = matche(fm_index, n1, n2, beg, it, nE-1, a, nMaxOcc, vSeed);
			//if (nRes > nMaxOcc) return nRes;
		}
		backward_search(fm_index, nL, nR, *it, nL, nR);
	}

	return vSeed.size();
}


size_t getIdx(const vector<uint64_t>& vPos, size_t nOcc)
{
	return distance(vPos.begin(), upper_bound(vPos.begin(), vPos.end(), nOcc));
}


TOcc getFwOcc(const vector<uint64_t>& vPos, size_t nOcc, int nK)
{
	size_t nIdx = getIdx(vPos, nOcc);

	if (1 == nIdx%2) // fw
	{
		nOcc -= vPos[nIdx-1];
		return TOcc(nIdx, nOcc, 0);
	}
	else // rv
	{
		nOcc = vPos[nIdx] - 1 - nOcc - nK;
		return TOcc(nIdx-1, nOcc, 16);
	}
}


TOcc getBwOcc(const vector<uint64_t>& vPos, size_t nOcc, int nK)
{
	size_t nIdx = getIdx(vPos, nOcc);

	if (1 == nIdx%2) // r-fw
	{
		nOcc = vPos[nIdx] - 1 - nOcc - nK;
		return TOcc(nIdx, nOcc, 0);
	}
	else // rv of r-fw
	{
		nOcc -= vPos[nIdx-1];
		return TOcc(nIdx-1, nOcc, 16);
	}
}


int buildOne(const vector<TRange>& vRange, bit_vector& bv, const vector<string>& vId, const vector<uint64_t>& vPos, const string& sWhole, const csa_wt<>& fw_fm_index, const csa_wt<>& bw_fm_index, int nK, int nM, int nE, int nMaxOcc, ostream& offNb, ostream& offMer)
{
	// k-mer & nb
	int nK1 = nK%2==0 ? nK/2 : (nK+1)/2;
	int nK2 = nK - nK1;
	char a[] = {'A', 'C', 'G', 'T'};
	vector<size_t> vSeed(128);
	vector<TOcc> vOcc(128);

	// m-mer sequences
	int nReadLen = 100;
	uint32_t nMark = (1<<((nM-1)<<1)) - 1;
	size_t nBuffLen = 1<<28;
	string ssMer;
	ssMer.reserve(nBuffLen);
	string ssNb;
	ssNb.reserve(nBuffLen);

	while (true)
	{
		MT.lock();
		size_t nIdx = IDX;
		++IDX;
		MT.unlock();

		if (nIdx >= vRange.size()) break;
		size_t nBeg = std::get<0>(vRange[nIdx]);
		size_t nCur = std::get<1>(vRange[nIdx]);

		MT_OUT.lock();
		cout << this_thread::get_id() << "\tworking on\tseg\t" << nIdx << "\t" << nBeg << "\t" << nCur << endl;
		MT_OUT.unlock();

		ssNb.clear();
		ssMer.clear();
		for (size_t k = nBeg; k < nCur-nK+1; ++k)
		{
			vSeed.clear();
			vOcc.clear();
			size_t nSame = 0;

			auto it1 = sWhole.begin() + k;
			auto it2 = it1 + nK;
			int nRes = matche(fw_fm_index, 0, fw_fm_index.size()-1, it1, it2, 0, a, nMaxOcc, vSeed);
			if (nRes > nMaxOcc) continue;

			bool b = false;
			MT_BV.lock();
			for (auto nSeed : vSeed)
			{
				if (bv[nSeed] == 1) 
				{
					b = true;
					break;
				}
				else
				{
					bv[nSeed] = 1;
				}
			}
			MT_BV.unlock();
			if (b == true) continue;

			for (auto seed : vSeed)
			{
				size_t nIdx;
				size_t nInOff;
				size_t nFlag;
				auto occ = getFwOcc(vPos, seed, nK);
				std::tie(nIdx,nInOff,nFlag) = occ;
				vOcc.push_back(occ);
				verify(sWhole, k, vPos[nIdx-1]+nInOff, nFlag, nK, 0);
			}
			nSame = vSeed.size();

			size_t nSeed = vSeed.size();
			size_t nL = 0;
			size_t nR = fw_fm_index.size() - 1;
			backward_search(fw_fm_index, nL, nR, it2-nK1, it2, nL, nR);
			nRes = matche(fw_fm_index, nL, nR, it1, it2-nK1, 1, a, nMaxOcc, vSeed);
			if (nRes > nMaxOcc) continue;
			for (auto nSeedIdx = nSeed; nSeedIdx < vSeed.size(); ++nSeedIdx)
			{
				size_t seed = vSeed[nSeedIdx];
				size_t nIdx;
				size_t nInOff;
				size_t nFlag;
				auto occ = getFwOcc(vPos, seed, nK);
				std::tie(nIdx,nInOff,nFlag) = occ;
				vOcc.push_back(occ);
				verify(sWhole, k, vPos[nIdx-1]+nInOff, nFlag, nK, 1);
			}

			nSeed = vSeed.size();
			nL = 0;
			nR = fw_fm_index.size() - 1;
			auto it3 = std::reverse_iterator<string::const_iterator>(it2);
			auto it4 = std::reverse_iterator<string::const_iterator>(it1);
			backward_search(bw_fm_index, nL, nR, it4-nK2, it4, nL, nR);
			nRes = matche(bw_fm_index, nL, nR, it3, it4-nK2, 1, a, nMaxOcc, vSeed);
			if (nRes > nMaxOcc) continue;
			for (auto nSeedIdx = nSeed; nSeedIdx < vSeed.size(); ++nSeedIdx)
			{
				size_t seed = vSeed[nSeedIdx];
				size_t nIdx;
				size_t nInOff;
				size_t nFlag;
				auto occ = getBwOcc(vPos, seed, nK);
				std::tie(nIdx,nInOff,nFlag) = occ;
				vOcc.push_back(occ);
				verify(sWhole, k, vPos[nIdx-1]+nInOff, nFlag, nK, 1);
			}

			// print k-mer & nb
			ssNb.append(sWhole.substr(k, nK));
			for (auto x : vOcc)
			{
				ssNb.append("$");
				ssNb.append(vId[get<0>(x)]);
				ssNb.append("|");
				ssNb.append(to_string(get<1>(x)));
				ssNb.append("|");
				ssNb.append(to_string(get<2>(x)));
				if ((get<2>(x)) == 0)
				{
					auto it1=sWhole.begin()+k;
					auto it2=sWhole.begin()+vPos[get<0>(x)-1]+get<1>(x); 
					size_t diff = 0;
					for (; it1 < sWhole.begin()+k+nK; ++it1,++it2) if (*it1!=*it2) ++diff;
					assert (diff<=1);
				}
				else
				{
					string sss;
					sss = sWhole.substr(vPos[get<0>(x)-1]+get<1>(x), nK);
					ToRC(sss);
					auto it1=sWhole.begin()+k;
					auto it2=sss.begin(); 
					size_t diff = 0;
					for (; it1 < sWhole.begin()+k+nK; ++it1,++it2) if (*it1!=*it2) ++diff;
					assert (diff<=1);
				}
			}
			ssNb.append("\n");

			//print m-mer sequence
			for (auto nSameIdx = 0; nSameIdx < nSame; ++nSameIdx)
			{
				auto x = vOcc[nSameIdx];
				size_t nLeft = vPos[get<0>(x)-1]; 
				size_t nRight = vPos[get<0>(x)] - 1; 
				size_t nInOff = get<1>(x);
				size_t n1 = max(nLeft, nLeft+nInOff-nReadLen+nK);
				size_t n2 = min(nRight, nLeft+nInOff+nReadLen);

				size_t n3 = nLeft + nInOff;
				while (n3>n1 && sWhole[n3-1]!='N') --n3;
				size_t n4 = nLeft + nInOff + nK;
				while (n4<n2 && sWhole[n4]!='N') ++n4;
				n1 = n3;
				n2 = n4;
				if (n2-n1 < nReadLen) continue;

				string sss = sWhole.substr(n1, n2-n1);
				assert (sss.find('N') == string::npos);
				int nS = nLeft + nInOff - n1;
				if (get<2>(x) != 0)
				{
					ToRC(sss);
					nS = n2 -n1- nS - nK;
				}
				assert (sWhole.substr(k,nK) == sss.substr(nS,nK));

				auto it1 = sss.begin() + nS;
				auto it2 = it1 + nK;

				// print k-mer
				ssMer.append(it1, it2);
				ssMer.append("|");
				ssMer.append(vId[get<0>(x)]);
				ssMer.append("|");
				ssMer.append(to_string(get<1>(x)));
				ssMer.append("|");
				ssMer.append(to_string(get<2>(x)));
				ssMer.push_back('$');

				// print left m-mers
				//for (auto itX = 0; itX < nS-nM+1; ++itX)
				for (auto itX = nS-nM; itX >= 0; --itX)
				{
					ssMer.append(sss.begin()+itX, sss.begin()+itX+nM);
					ssMer.push_back(' ');
				}
				ssMer.push_back('$');

				// print right m-mers
				for (auto itX = nS+nK; itX < sss.size()-nM+1; ++itX)
				{
					ssMer.append(sss.begin()+itX, sss.begin()+itX+nM);
					ssMer.push_back(' ');
				}
				ssMer.push_back('\n');
			}

			if (ssNb.size() > nBuffLen-1000000)
			{
				MT_OUT.lock();
				offNb << ssNb;
				ssNb.clear();
				MT_OUT.unlock();
			}

			if (ssMer.size() > nBuffLen-1000000)
			{
				MT_OUT.lock();
				offMer << ssMer;
				ssMer.clear();
				MT_OUT.unlock();
			}
		}

		MT_OUT.lock();
		offNb << ssNb;
		ssNb.clear();
		MT_OUT.unlock();

		MT_OUT.lock();
		offMer << ssMer;
		ssMer.clear();
		MT_OUT.unlock();
	}

	return 0;
}


int build(const vector<string>& vId, const vector<uint64_t>& vPos, const string& sWhole, const csa_wt<>& fw_fm_index, const csa_wt<>& bw_fm_index, int nK, int nE, int nM, int nMaxOcc, string sRefFile)
{
	bit_vector bv(fw_fm_index.size()-1, 0);
	ofstream offNb(sRefFile+".nb");
	ofstream offMer(sRefFile+".mer");

	vector<TRange> vRange;
	for (size_t i = 1; i < vId.size(); ++i)
	{
		size_t nBeg = 0;
		if (i > 0)
		{
			nBeg = vPos[i-1];
		}
		size_t nEnd = vPos[i] - 1; // since 'X' is added at the end
		size_t nLen = nEnd - nBeg;

		while (nBeg < nEnd)
		{
			while (nBeg<nEnd && sWhole[nBeg]=='N') ++nBeg;
			size_t nCur = nBeg;
			while (nCur<nEnd && sWhole[nCur]!='N') ++nCur;
			if (nCur==nBeg) break;

			string sX(sWhole.substr(nBeg,nCur-nBeg));
			assert (sX.find('N') == string::npos);

			vRange.push_back(TRange(nBeg, nCur));
			cout << nBeg << "\t" << nCur << "\t=\t" << nCur-nBeg << "\t" << sWhole.size() << endl;

			nBeg = nCur;
		}
	}
	sort(vRange.begin(), vRange.end(), [](const TRange& x1, const TRange& x2)
			{
				return (std::get<1>(x1)-std::get<0>(x1))>(std::get<1>(x2)-std::get<0>(x2));
			});

	//buildOne(vRange, bv, vId, vPos, sWhole, fw_fm_index, bw_fm_index, nK, nM, nE, nMaxOcc, offNb, offMer);
	size_t nThd = 20;
	thread_group gNb;
	for (size_t i = 0; i < nThd; ++i)
	{
		gNb.create_thread(bind(buildOne, boost::ref(vRange), boost::ref(bv), boost::ref(vId), boost::ref(vPos), boost::ref(sWhole), boost::ref(fw_fm_index), boost::ref(bw_fm_index), nK, nM, nE, nMaxOcc, boost::ref(offNb), boost::ref(offMer)));
	}
	gNb.join_all();

	offNb.close();
	offMer.close();

	return 0;
}


int main(int argc, char** argv)
{
	string sRefFile = argv[1];
	int nK = 27;
	int nE = 1;
	int nM = 5;
	int nMaxOcc = 30;
	auto t1 = TClock::now();

	// read reference
	vector<string> vId;
	vector<uint64_t> vPos;
	string sFw, sBw;
	ifstream iff(sRefFile+".all");
	if (iff.good())
	{
		archive::binary_iarchive ia(iff);
		ia >> vId;
		ia >> vPos;
		ia >> sFw;
		ia >> sBw;
		iff.close();
	}
	else
	{
		if (0 != read(sRefFile, vId, vPos, sFw, sBw))
		{
			return -1;
		}

		ofstream off(sRefFile+".all");
		archive::binary_oarchive oa(off);
		oa << vId;
		oa << vPos;
		oa << sFw;
		oa << sBw;
		off.close();
	}

	cout << "read ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;

	copy(vId.begin(), vId.end(), ostream_iterator<string>(cout, "\n"));
	cout << endl;

	copy(vPos.begin(), vPos.end(), ostream_iterator<uint64_t>(cout, "\n"));
	cout << endl;

	// build fm_index
	t1 = TClock::now();
	csa_wt<> fw_fm_index;
	string sRefFw(sRefFile+".fw.fm"); 
	iff.open(sRefFw);
	if (iff.good())
	{
		iff.close();
		load_from_file(fw_fm_index, sRefFw);
	}
	else
	{
		construct_im(fw_fm_index, sFw.c_str(), 1);
		cout << "build fw ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;

		t1 = TClock::now();
		store_to_file(fw_fm_index, sRefFw);
		cout << "serialize fw ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;
	}

	t1 = TClock::now();
	csa_wt<> bw_fm_index;
	string sRefBw(sRefFile+".bw.fm"); 
	iff.open(sRefBw);
	if (iff.good())
	{
		iff.close();
		load_from_file(bw_fm_index, sRefBw);
	}
	else
	{
		construct_im(bw_fm_index, sBw.c_str(), 1);
		cout << "build bw ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;

		t1 = TClock::now();
		store_to_file(bw_fm_index, sRefBw);
		cout << "serialize bw ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;
	}

	int nRes = build(vId, vPos, sFw, fw_fm_index, bw_fm_index, nK, nE, nM, nMaxOcc, sRefFile);
	if (0 != nRes)
	{
		return nRes;
	}

	return 0;
}
