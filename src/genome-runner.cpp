#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <list>
using namespace std;

/* Genomic features and features of interest are represented as a map from chromsome to 
 * (sorted) vector of Features (positions) on that chromosome -- aka, a FeatureSet.
 *
 * Throughout the code, 
 * A = Feature of Interest
 * B = Genomic Feature
 */

struct Feature {
	int start;
	int end;
	int length() {return end - start;}
	bool operator < (const Feature& other) const {
		return start==other.start ? end < other.end : start < other.start;
	}
	bool overlaps(const Feature& other) {
		//WARNING: assumes chromsomes equal. Must not compare Features outside of FeatureSets!
		return other.start <= end && start <= other.end;
	}
};

typedef unordered_map<string, vector<Feature> > FeatureSet;

struct Genome {
	uint32_t size;
	vector<uint32_t> lengths;
	vector<uint32_t> cumlens;
	vector<string> chromosomes;
	pair<string, Feature> random(int length) {
		/* Inspired by BEDTools shuffleBed (in particular, shuffleBed.cpp)
		 * Generates a random number between 0 and Genome.size, find the correct chrom/start.
		 * */

		//FIXME binary search is faster
		int chrom_size;
		Feature f;
		string chrom;
		do {
			uint32_t rnd = uint32_t((rand() << 31) | rand()) % size;
			vector<uint32_t>::const_iterator low = 
				lower_bound(cumlens.begin(), cumlens.end(), rnd+1);
			int i = int(low - cumlens.begin());
			chrom = chromosomes[i];
			f.start = cumlens[i] - rnd;
			f.end = f.start + length;
			chrom_size = lengths[i];
		} while (f.end >= chrom_size);
		// Instead of storing the chromosome along with the feature, return the chromosome and 
		// start/end as a pair. The chromosome will be used to place the Feature within the
		// FeatureSet, then discarded.
		pair<string, Feature> pr(chrom, f);
		return pr;
	}
};

void fs_sort(Genome& g, FeatureSet& fs) {
	for (string chrom : g.chromosomes) {
		sort(fs[chrom].begin(), fs[chrom].end());
	}
}

// I/O
FeatureSet read_features(string path) {
	fstream strm(path);
	string line, chrom;
	int start, end;	
	FeatureSet fs;
	while (getline(strm,line)) {
		istringstream ss(line);
		ss >> chrom >> start >> end;
		Feature p = {start,end};
		fs[chrom].push_back(p);
	}
	return fs;	
}

void write_bed(FeatureSet& features, string path) {
	ofstream out(path);
	for (auto kv : features) {
		string chrom = kv.first;
		for (Feature f : kv.second)
			out << chrom << "\t" << f.start << "\t" << f.end << endl;
	}
}
Genome read_genome(string path) {
	/* A "genome" file is a tab-delimited file of chromsome names and lengths, exactly as used
 	* in BEDTools. It can be created, for example, as follows (for hg19):
 	*
 	* mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
 	*         "select chrom, size from hg19.chromInfo"  > hg19.genome
 	*/
	fstream strm(path);
	string line, chrom;
	int length;
	Genome g;
	getline(strm, line); //ignore header

	uint32_t cumlen = 0;
	while (getline(strm,line)) {
		istringstream ss(line);
		ss >> chrom >> length;
		g.lengths.push_back(length);
		cumlen += length;	
		g.cumlens.push_back(cumlen);
		g.chromosomes.push_back(chrom);
	}
	g.size = cumlen;
	return g;
}

/* Overlaps are determined for each chromosome individually using the chromSweep algorithm
 * https://github.com/arq5x/chrom_sweep
 * Then the number of overlaps for all chromosomes are merged ("overlaps" function).
 */

inline bool after(Feature& a, Feature& b) {
	return a.start > b.end;
}

int overlaps_single_chromosome(vector<Feature>& A, vector<Feature>& B) {
	//Return #of As overlapping with any element in B
	int ai, bi, n;
	ai = bi = n = 0;

	vector<Feature> cache;
	vector<Feature>::iterator c;

	Feature b = B[0];
	while (ai<A.size()) {
		Feature a = A[ai];
		bool overlap = false;
		c = cache.begin();
		// Check the existing cache to see if the current FOI overlaps 
		// any GFs in it. Remove any GFs from the cache that are no longer needed.
		while (c != cache.end()) {
			if (!after(a, *c)) {
				if (a.overlaps(*c))
					overlap = true;
				++c;
			} else {
				cache.erase(c);
			}
		}
		
		// Take new GFs from the stack of GFs, recording overlaps as we go,
		// until the next GF is after the current FOI.
		while (bi!=B.size() && !after(b,a)) {
			if (!overlap && a.overlaps(b))
				overlap = true;
			cache.push_back(b);
			b = B[bi++];
		}
			
		if (overlap)
			n++;
		ai++;
	}
	return n;
}	

int overlaps(Genome& g, FeatureSet& A, FeatureSet& B) {
	int result = 0;
	for (string c : g.chromosomes) {
		vector<Feature> fA = A[c];
		vector<Feature> fB = B[c];
		if (fA.empty() || fB.empty())
			continue;
		result += overlaps_single_chromosome(fA,fB);
	}
	return result;
}

FeatureSet shuffle(Genome& genome, FeatureSet& features) {
	// Build a shuffled FeatureSet with the same distribution of lengths as the original.
	FeatureSet shuffled;
	for (auto kv : features) {
		vector<Feature> v = kv.second;
		for (Feature f : v) {
			pair<string, Feature> pr = genome.random(f.length());
			shuffled[pr.first].push_back(pr.second);
		}
	}
	fs_sort(genome, shuffled);	
	return shuffled;
}

vector<int>
overlap_distribution(Genome& genome, FeatureSet& A, FeatureSet& B, int niter) {
	// Perform a monte-carlo test to find the overlap distribution. 
	vector<int> noverlap (niter);
	#pragma omp parallel for
	for (int i=0; i<niter; i++) {
		FeatureSet Ashuffled = shuffle(genome, A);
		noverlap[i] = overlaps(genome, Ashuffled, B);
	}	
	return noverlap;
}

double pvalue(int observed, vector<int>& distribution) {
	sort(distribution.begin(), distribution.end());
	double i = double(lower_bound(distribution.begin(), distribution.end(), observed)
		- distribution.begin());
	double n = distribution.size();
	return min(i,n-i) / n;
}

double mean(vector<int>& v) {
	double result = 0;
	for (int n : v) {
		result += n;
	}
	result /= v.size();
	return result;
}

int fs_size(FeatureSet& fs) {
	int n = 0;
	for (auto kv : fs)
		n += kv.second.size();
	return n;
}

extern "C" {
	struct Enrichment {
		const char* A;
		const char* B;
		int nA;
		int nB;
		int observed;
		double expected;
		double p_value;
	};

	Enrichment enrichment(const char* g,
			const char* a, const char* b, int n) {
		Genome genome = read_genome(g);
		FeatureSet A = read_features(a);
		FeatureSet B = read_features(b);
		fs_sort(genome, A);
		fs_sort(genome, B);

		Enrichment result;
		int obs = overlaps(genome, A, B);
		vector<int> dist = overlap_distribution(genome,A,B,n);

		result.observed = obs;
		result.expected = mean(dist);
		result.p_value = pvalue(obs,dist);
		result.nA = fs_size(A);
		result.nB = fs_size(B);
		result.A = a;
		result.B = b;
		return result;
	}
}

/*
int main(int argc, char* argv[]) {
	if (argc != 4)
		cerr << "USAGE: ./genome-runner genome-file A.bed B.bed" << endl;
	Genome genome = read_genome(argv[1]);
	FeatureSet A = read_features(argv[2]);
	FeatureSet B = read_features(argv[3]);
	feature_set_sort(genome, A);
	feature_set_sort(genome, B);

	int obs = overlaps(genome, A, B);
	vector<int> dist = overlap_distribution(genome,A,B,1000); 
	//cout << "Observed\tMean\tP-Value\n";
	cout << obs << "\t" << mean(dist) << "\t" << pvalue(obs, dist) << endl;
}
*/
