#include <fstream>
#include <vector>
#include <chrono>
#include <exception>
#include <utility>

#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>

#include "util.hpp"

using namespace std;
using namespace std::chrono;
using namespace boost;
using namespace sdsl;

typedef high_resolution_clock TClock;
typedef pair<uint64_t,uint64_t> TSeed;


boost::mutex MT;
size_t IDX = 0;


int read(string sIn, vector<string>& vId, vector<uint64_t>& vPos, string& sWhole)
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
				sWhole += seq + "X";
				vId.push_back(id);
				vPos.push_back(sWhole.length());

				ToRC(seq);
				sWhole += seq + "X";
				vId.push_back(id+"_rc");
				vPos.push_back(sWhole.length());
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
	sWhole += seq + "X";
	vId.push_back(id);
	vPos.push_back(sWhole.length());

	ToRC(seq);
	sWhole += seq + "X";
	vId.push_back(id+"_rc");
	vPos.push_back(sWhole.length());
	iff.close();

	return 0;
}


int test(csa_wt<>& fm_index)
{
	ifstream iff("lala");
	if (!iff.good())
	{
		cout << "ERROR: cannot open the file" << endl;
		return -1;
	}

	string sLine;
	while (getline(iff, sLine))
	{
		int nOcc = stoi(sLine.substr(1));

		getline(iff, sLine);
		int n = count(fm_index, sLine);
		if (n != nOcc)
		{
			cout << nOcc << "\t" << n << "\t" << sLine << endl;
		}

		size_t nL = 0;
		size_t nR = fm_index.size() - 1;
		for (auto it = sLine.rbegin(); it != sLine.rend(); ++it)
		{
			backward_search(fm_index, nL, nR, *it, nL, nR);
		}
		if (nL <= nR)
		{
			assert (nR-nL+1 == n);
			for (size_t i = 0; i < nR-nL+1; ++i)
			{
				string sExtract = extract(fm_index, fm_index[nL+i], fm_index[nL+i]+sLine.size()-1);
				for (size_t j = 0; j < sLine.size(); ++j)
				{
					if (sExtract[j] != sLine[j])
					{
						cout << sExtract << "\t" << sLine << endl;
						cout << sExtract[j] << "\t" << sLine[j] << endl;
						return -3;
					}
				}
			}
		}
		else
		{
			cerr << "ERROR: search wrong" << endl;
			return -2;
		}

	}

	iff.close();

	return 0;
}


int align(string& sId, string& sSeq, csa_wt<>& fm_index, int nK, vector<TSeed>& vSeed, vector<string>& vRId, vector<uint64_t>& vPos, string& sRes)
{
	sRes = sId.substr(1,sId.find('/')-1) + "\t" + to_string(4) + "\t*\t0\n";
	vector<vector<int> > vv1(200, vector<int>(200, 0));
	string sRefRL; sRefRL.reserve(200);
	string sQryRL; sQryRL.reserve(200);
	string sRefRR; sRefRR.reserve(200);
	string sQryRR; sQryRR.reserve(200);
	string sDiff; sDiff.reserve(200);
	size_t nMinDiff = 4;
	int nPenalty = 6;
	int nMismatch = 3;
	for (auto p : vSeed)
	{
		size_t ref = fm_index[p.first];
		size_t qry = p.second;

		string sRef = extract(fm_index, ref-qry-5, ref+sSeq.size()-qry-1+5);

		int rl = 5 + qry - 1;
		int rr = 5 + qry + nK;
		int ql = qry - 1;
		int qr = qry + nK;

		// forward
		for (size_t i = qr; i < sSeq.size(); ++i)
		{
			if (sRef[rr] != sSeq[qr])
			{
				break;
			}
			++rr;
			++qr;
		}

		string sRefR = sRef.substr(5+qry, rr-qry-5);
		string sQryR = sSeq.substr(qry, qr-qry);

		if (qr < sSeq.size())
		{
			for (auto& v : vv1) v.assign(200,0);
			for (size_t i = 0; i < vv1[0].size(); ++i) vv1[0][i]=-i*nPenalty;
			for (size_t i = 0; i < vv1.size(); ++i) vv1[i][0]=-i*nPenalty;
			size_t nI = 0;
			size_t nJ = 0;
			int nV = INT_MIN;
			for (size_t i = 1; i < sSeq.size()-qr+1; ++i)
			{
				for (size_t j = 1; j < sRef.size()-rr+1; ++j)
				{
					vv1[i][j] = max({vv1[i][j-1]-nPenalty,vv1[i-1][j]-nPenalty,vv1[i-1][j-1]+(sSeq[qr+i-1]==sRef[rr+j-1]?1:-nMismatch)});
				}
			}
			for (size_t j = 1; j < sRef.size()-rr+1; ++j)
			{
				if (vv1[sSeq.size()-qr][j] > nV)
				{
					nV = vv1[sSeq.size()-qr][j];
					nI = sSeq.size()-qr;
					nJ = j;
				}
			}

			sRefRR.clear();
			sQryRR.clear();
			while (nI>0 && nJ >0)
			{
				if (vv1[nI][nJ] == vv1[nI-1][nJ-1]+(sSeq[qr+nI-1]==sRef[rr+nJ-1]?1:-nMismatch))
				{
					sQryRR.push_back(sSeq[qr+nI-1]);
					sRefRR.push_back(sRef[rr+nJ-1]);
					--nI; --nJ;
				}
				else if (vv1[nI][nJ] == vv1[nI-1][nJ] - nPenalty)
				{
					sQryRR.push_back(sSeq[qr+nI-1]);
					sRefRR.push_back('-');
					--nI;
				}
				else if (vv1[nI][nJ] == vv1[nI][nJ-1] - nPenalty)
				{
					sQryRR.push_back('-');
					sRefRR.push_back(sRef[rr+nJ-1]);
					--nJ;
				}
			}
			while (nI > 0)
			{
				sQryRR.push_back(sSeq[qr+nI-1]);
				sRefRR.push_back('-');
				--nI;
			}
			while (nJ > 0)
			{
				sQryRR.push_back('-');
				sRefRR.push_back(sRef[rr+nJ-1]);
				--nJ;
			}
			reverse(sRefRR.begin(), sRefRR.end());
			reverse(sQryRR.begin(), sQryRR.end());

			sRefR += sRefRR;
			sQryR += sQryRR;

			for (auto c : sRefRR)
				if ('-' != c) ++rr;
		}


		// backward
		// note that - before comparison
		for (size_t i = 0; i < qry; ++i)
		{
			if (sRef[rl] != sSeq[ql])
			{
				break;
			}
			--rl;
			--ql;
		}

		sRefR.insert(0, sRef.substr(rl+1, 5+qry-rl-1));
		sQryR.insert(0, sSeq.substr(ql+1, qry-ql-1));
		if (ql >= 0)
		{
			for (auto& v : vv1) v.assign(200,0);
			for (size_t i = 0; i < vv1[0].size(); ++i) vv1[0][i]=-i*nPenalty;
			for (size_t i = 0; i < vv1.size(); ++i) vv1[i][0]=-i*nPenalty;
			size_t nI = 0;
			size_t nJ = 0;
			int nV = INT_MIN;
			for (size_t i = 1; i < ql+1+1; ++i)
			{
				for (size_t j = 1; j < rl+1+1; ++j)
				{
					vv1[i][j] = max({vv1[i][j-1]-nPenalty,vv1[i-1][j]-nPenalty,vv1[i-1][j-1]+(sSeq[ql-i+1]==sRef[rl-j+1]?1:-nMismatch)});
				}
			}
			for (size_t j = 1; j < rl+1+1; ++j)
			{
				if (vv1[ql+1][j] > nV)
				{
					nV = vv1[ql+1][j];
					nI = ql+1;
					nJ = j;
				}
			}

			sRefRL.clear();
			sQryRL.clear();
			while (nI>0 && nJ>0)
			{
				if (vv1[nI][nJ] == vv1[nI-1][nJ-1]+(sSeq[ql-nI+1]==sRef[rl-nJ+1]?1:-nMismatch))
				{
					sQryRL.push_back(sSeq[ql-nI+1]);
					sRefRL.push_back(sRef[rl-nJ+1]);
					--nI; --nJ;
				}
				else if (vv1[nI][nJ] == vv1[nI-1][nJ] - nPenalty)
				{
					sQryRL.push_back(sSeq[ql-nI+1]);
					sRefRL.push_back('-');
					--nI;
				}
				else if (vv1[nI][nJ] == vv1[nI][nJ-1] - nPenalty)
				{
					sQryRL.push_back('-');
					sRefRL.push_back(sRef[rl-nJ+1]);
					--nJ;
				}
			}
			while (nI > 0)
			{
				sQryRL.push_back(sSeq[ql-nI+1]);
				sRefRL.push_back('-');
				--nI;
			}
			while (nJ > 0)
			{
				sQryRL.push_back('-');
				sRefRL.push_back(sRef[rl-nJ+1]);
				--nJ;
			}

			sRefR.insert(0, sRefRL);
			sQryR.insert(0, sQryRL);

			for (auto c : sRefRL)
				if ('-' != c) --rl;
		}

		sDiff.clear();
		int nDiff = 0;
		for (size_t i = 0; i < sQryR.size(); ++i)
		{
			if (sQryR[i] == sRefR[i])
			{
				sDiff.push_back('|');
				++nDiff;
			}
			else
			{
				sDiff.push_back('*');
			}
		}
		nDiff = sSeq.size() - nDiff;
		if (nDiff <= nMinDiff)
		{
			size_t nIdx = distance(vPos.begin(), upper_bound(vPos.begin(), vPos.end(), ref));
			string sRefId;
			size_t nRefOff = ref - (5+qry-(rl+1)); // rl is the most left different nt
			size_t nFlag;
			if (0 == nIdx%2)
			{
				nRefOff = ref - (5+qry-(rl+1));
				nFlag = 0;
				sRefId = vRId[nIdx];
			}
			else
			{
				nRefOff = ref + (rr-5-qry);
				nFlag = 16;
				sRefId = vRId[nIdx-1];
			}
			if (0 != nIdx)
			{
				nRefOff -= vPos[nIdx-1]; // 0 based
			}
			if (nIdx%2 == 1)
			{
				nRefOff = vPos[nIdx]-vPos[nIdx-1]-1-nRefOff;
			}
			++nRefOff;
			if (nDiff == nMinDiff)
			{
				sRes += sId.substr(1,sId.find('/')-1) + "\t" + to_string(nFlag) + "\t" + sRefId.substr(1) + "\t" + to_string(nRefOff) + "\n" + sRefR + "\n" + sDiff + "\n" + sQryR + "\n";
			}
			else
			{
				nMinDiff = nDiff;
				sRes = sId.substr(1,sId.find('/')-1) + "\t" + to_string(nFlag) + "\t" + sRefId.substr(1) + "\t" + to_string(nRefOff) + "\n" + sRefR + "\n" + sDiff + "\n" + sQryR + "\n";
			}
		}
	}

	return nMinDiff;
}


void verify(string& sSeq, int i, csa_wt<>& fm_index, size_t j, int nK, int nE)
{
	string sExtract = extract(fm_index, fm_index[j], fm_index[j]+nK-1);
	int nn = 0;
	string sNN;
	for (size_t o = 0; o < nK; ++o)
	{
		if (sSeq[i+1-nK+o]!=sExtract[o])
		{
			++nn;
			sNN.push_back('*');
		}
		else
		{
			sNN.push_back('|');
		}
	}
	assert (nn == nE);
}


int match1e(string& sSeq, int i, csa_wt<>& fm_index, char* a, int nK, size_t nMaxOcc, vector<TSeed>& vSeed)
{
	size_t nL = 0;
	size_t nR = fm_index.size() - 1;
	for (int j = i; j >= i+1-nK; --j)
	{
		for (int k = 0; k < 4; ++k)
		{
			if (sSeq[j] == a[k]) continue;

			size_t n1 = nL;
			size_t n2 = nR;
			backward_search(fm_index, n1, n2, a[k], n1, n2);
			if (n1 > n2) continue;
			backward_search(fm_index, n1, n2, &sSeq[i-nK+1], &sSeq[j], n1, n2);
			if (n1 > n2) continue;

			size_t nVarOcc = n2 - n1 + 1;
			if (vSeed.size()+nVarOcc > nMaxOcc)
			{
				return vSeed.size()+nVarOcc;
			}

			for (size_t l = 0; l < nVarOcc; ++l)
			{
				vSeed.push_back(TSeed(n1+l,i-nK+1));
				verify(sSeq, i, fm_index, n1+l, nK, 1);
			}
		}
		backward_search(fm_index, nL, nR, sSeq[j], nL, nR);
		if (nL > nR) break;
	}

	return vSeed.size();
}


int match2e(string& sSeq, int i, csa_wt<>& fm_index, char* a, int nK, size_t nMaxOcc, vector<TSeed>& vSeed)
{
	size_t nL = 0;
	size_t nR = fm_index.size() - 1;
	for (int j = i; j > i+1-nK; --j)
	{
		for (int k = 0; k < 4; ++k)
		{
			if (sSeq[j] == a[k]) continue;

			size_t n1 = nL;
			size_t n2 = nR;
			backward_search(fm_index, n1, n2, a[k], n1, n2);
			if (n1 > n2) continue;

			for (int l = j-1; l >= i+1-nK; --l)
			{
				for (int m = 0; m < 4; ++m)
				{
					if (sSeq[l] == a[m]) continue;

					size_t n3 = n1;
					size_t n4 = n2;
					backward_search(fm_index, n3, n4, a[m], n3, n4);
					if (n3 > n4) continue;
					backward_search(fm_index, n3, n4, &sSeq[i-nK+1], &sSeq[l], n3, n4);
					if (n3 > n4) continue;

					size_t nVarOcc = n4 - n3 + 1;
					if (vSeed.size()+nVarOcc > nMaxOcc)
					{
						return vSeed.size()+nVarOcc;
					}

					for (size_t n = 0; n < nVarOcc; ++n)
					{
						vSeed.push_back(TSeed(n3+n,i-nK+1));
						verify(sSeq, i, fm_index, n3+n, nK, 2);
					}
				}
				backward_search(fm_index, n1, n2, sSeq[l], n1, n2);
				if (n1 > n2) break;
			}
		}
		backward_search(fm_index, nL, nR, sSeq[j], nL, nR);
		if (nL > nR) break;
	}

	return vSeed.size();
}


int mapOne(string& sId, string& sSeq, csa_wt<>& fm_index, int nK, int nE, vector<string>& vRId, vector<uint64_t>& vPos, string& sResult)
{
	char a[] = {'A', 'C', 'G', 'T'};
	size_t nMaxOcc = 30;

	string sSeedNum;
	vector<TSeed> vSeed(128);
	vector<int> vSelf(sSeq.size()-nK+1, INT_MAX);
	for (int i = nK-1; i < sSeq.size(); ++i)
	{
		vSeed.clear();

		// the exact match
		size_t nL = 0;
		size_t nR = fm_index.size() - 1;
		backward_search(fm_index, nL, nR, &sSeq[i-nK+1], &sSeq[i+1], nL, nR);
		if (nL > nR) continue;
		size_t nOcc = nR - nL + 1;
		vSelf[i-nK+1] = nOcc;
		if (nOcc > nMaxOcc)
		{
			continue;
		}
		for (size_t j = 0; j < nOcc; ++j)
		{
			vSeed.push_back(TSeed(nL+j, i-nK+1));
		}

		// one error match
		size_t res =  match1e(sSeq, i, fm_index, a, nK, nMaxOcc, vSeed);
		if (res > nMaxOcc)
		{
			continue;
		}

		// two errors match
		if (2 <= nE)
		{
			size_t res =  match2e(sSeq, i, fm_index, a, nK, nMaxOcc, vSeed);
			if (res > nMaxOcc)
			{
				continue;
			}
		}

		sSeedNum += "\t" + to_string(i-nK+1) + ")" + to_string(nOcc) + "|" + to_string(vSeed.size());
		if (vSeed.size() > nMaxOcc)
		{
			continue;
		}

		string sRes;
		size_t nDiff = align(sId, sSeq, fm_index, nK, vSeed, vRId, vPos, sRes);
		sResult += sId + sSeedNum + "\n";
		sResult += sRes;

		return 0;
	}

	// if all kmer+nb have many occ, then only consider kmer
	sSeedNum.clear();
	for (size_t i = 0; i < vSelf.size(); ++i)
	{
		sSeedNum += "\t" + to_string(i) + ")" + to_string(vSelf[i]);
		if (vSelf[i] <= nMaxOcc)
		{
			vSeed.clear();
			size_t nL = 0;
			size_t nR = fm_index.size() - 1;
			backward_search(fm_index, nL, nR, &sSeq[i], &sSeq[i+nK], nL, nR);
			int nOcc = nR - nL + 1;
			for (size_t j = 0; j < nOcc; ++j)
			{
				vSeed.push_back(TSeed(nL+j, i));
			}
			string sRes;
			size_t nDiff = align(sId, sSeq, fm_index, nK, vSeed, vRId, vPos, sRes);
			sResult += sId + sSeedNum + "\n";
			sResult += sRes;

			return 0;
		}
	}

	sResult += sId + sSeedNum + "\n";
	sResult += sId.substr(1,sId.find('/')-1) + "\t" + to_string(4) + "\t*\t0\n";
	return 0;
}


void mapBatch(vector<string>& vQId, vector<string>& vQSeq, csa_wt<>& fm_index, int nK, int nE, vector<string>& vRId, vector<uint64_t>& vPos, int nThd, int nThdIdx, vector<string>& vRes)
{
	size_t nStep = 50;	// better than size/thread ?
	auto t1 = TClock::now();
	while (true)
	{
		MT.lock();
		size_t nCur = IDX;
		IDX += nStep;
		MT.unlock();

		if (nCur >= vQId.size())
		{
			break;
		}

		for (size_t i = nCur; i < nCur+nStep; ++i)
		{
			if (i < vQId.size())
			{
				mapOne(vQId[i], vQSeq[i], fm_index, nK, nE, vRId, vPos, vRes[i]);
			}
			else
			{
				break;
			}
		}
	}

	MT.lock();
	cout << "thread:\t" << this_thread::get_id() << "\t" << duration_cast<seconds>(TClock::now()-t1).count() << endl;
	MT.unlock();
}


int mymap(csa_wt<>& fm_index, string sIn, int nK, int nE, vector<string>& vRId, vector<uint64_t>& vPos)
{
	ifstream iff(sIn);
	if (!iff.good())
	{
		cout << "ERROR: cannot open the file" << endl;
		return -1;
	}

	vector<string> vQId;
	vector<string> vQSeq;

	string sId; // assume two lines for a sequence
	string sSeq;
	while (getline(iff, sId))
	{
		getline(iff, sSeq);
		vQId.push_back(sId);
		vQSeq.push_back(sSeq);
	}

	vector<string> vResult(vQId.size());

	//mapBatch(boost::ref(vQId), boost::ref(vQSeq), boost::ref(fm_index), nK, nE, boost::ref(vRId), boost::ref(vPos), 1, 0, vResult);

	size_t nThd = 35;
	thread_group g;
	for (size_t i = 0; i < nThd; ++i)
	{
		g.create_thread(bind(mapBatch, boost::ref(vQId), boost::ref(vQSeq), boost::ref(fm_index), nK, nE, boost::ref(vRId), boost::ref(vPos), nThd, i, boost::ref(vResult)));
	}
	g.join_all();

	for (size_t i = 0; i < vResult.size(); ++i)
	{
		cout << vResult[i] << endl;
	}
}


int main(int argc, char** argv)
{
	auto t1 = TClock::now();
	if (2 == argc)
	{
		string sIn(argv[1]);
		t1 = TClock::now();
		vector<string> vId;
		vector<uint64_t> vPos;
		string sWhole;
		if (0 != read(sIn, vId, vPos, sWhole))
		{
			return -1;
		}
		cout << "read ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;

		for (size_t i = 0; i < vId.size(); ++i)
		{
			cout << vId[i] << "\t" << vPos[i] << endl;
		}

		memory_monitor::start();
		t1 = TClock::now();
		csa_wt<> fm_index;
		construct_im(fm_index, sWhole.c_str(), 1);
		cout << "build ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;
		memory_monitor::stop();
		//memory_monitor::write_memory_log<HTML_FORMAT>(cout);

		t1 = TClock::now();
		store_to_file(fm_index, (sIn+".fm"));
		ofstream out(sIn+".fm.html");
		write_structure<HTML_FORMAT>(fm_index, out);
		cout << "serialize ref in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;
	}
	else if (5 == argc)
	{
		t1 = TClock::now();
		string sRef(argv[1]);
		vector<string> vId;
		vector<uint64_t> vPos;
		ifstream iff(sRef+".il");
		if (!iff.good())
		{
			string sWhole;
			if (0 != read(sRef, vId, vPos, sWhole))
			{
				return -1;
			}
			sWhole.clear();

			ofstream off(sRef+".il");
			archive::binary_oarchive oa(off);
			oa << vId;
			oa << vPos;
			off.close();
		}
		else
		{
			archive::binary_iarchive ia(iff);
			ia >> vId;
			ia >> vPos;
			iff.close();
		}

		string sQry(argv[2]);
		csa_wt<> fm_index;
		load_from_file(fm_index, sRef+".fm");
		cout << "load in " << duration_cast<seconds>(TClock::now()-t1).count() << endl;
		for (size_t i = 0; i < vId.size(); i+=2)
		{
			cout << "@\t" << vId[i] << endl;
		}

		//if (0 != test(fm_index))
		//{
			//cerr << "ERROR: assert wrong" << endl;
			//return -2;
		//}

		int nK = stoi(argv[3]);
		int nE = stoi(argv[4]);
		mymap(fm_index, sQry, nK, nE, vId, vPos);
	}
	else
	{
		cout << "USAGE: ./index ref [qry] nK nE" << endl;
	}

	return 0;
}
