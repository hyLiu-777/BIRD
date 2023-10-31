#ifndef ALIAS_H
#define ALIAS_H

#include <random>
#include <algorithm>
#include <stack>

using namespace std;

typedef unsigned int uint;

class Alias{
public:
	double* p;
	uint* h;
	uint* map1;
	uint n;
	Alias(vector<pair<int, double> > pi){
		double sum = 0;
		n = pi.size();
		stack<uint> small;
		stack<uint> big;
		p = new double[n];
		h = new uint[n];
		map1 = new uint[n];
		for(uint i = 0; i < n; i++){
			sum += pi[i].second;
			map1[i] = pi[i].first;
		}
		for(uint i = 0; i < n; i++){
			p[i] = pi[i].second * n / sum;
			if(p[i] > 1)
				big.push(i);
			else
				small.push(i);
		}
		while(!small.empty() && !big.empty()){
			uint smallId = small.top();
			small.pop();
			uint bigId = big.top();
			h[smallId] = bigId;
			p[bigId] -= (1-p[smallId]);
			if(p[bigId] < 1){
				small.push(bigId);
				big.pop();
			}
		}
	}

	~Alias(){
		delete[] p;
		delete[] h;
		delete[] map1;
	}
	uint generateRandom(){
		uint firstId = rand() % RAND_MAX / (double) (RAND_MAX) * n;
		uint answer = (rand() % RAND_MAX / (double) (RAND_MAX) * n) < p[firstId] ? map1[firstId] : map1[h[firstId]];
		return answer;
	}
};

#endif
