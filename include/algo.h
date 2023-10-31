/*************************************************************************
    > File Name: algo.h
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:27:17 PM
 ************************************************************************/

#include<vector>
#include <unordered_map>
#include <map>
#include "graph.h"
#include "alias.h"

void powerIter(int src, double alpha, uint n_iter, std::vector<double>& ppr, const Graph& graph);
void powerIterGt(int src, double alpha, uint n_iter, std::vector<double>& sppr, const Graph& graph);
void monteCarlo(int src, double alpha, int64 n_walk, std::vector<double>& ppr, const Graph& graph);
void monteCarloWeighted(int src, double alpha, int64 n_walk, std::vector<double>& sppr, const Graph& graph, vector<Alias>& U_alias, vector<Alias>& V_alias);
void reversePush(int tgt, double alpha, double eps, std::vector<double>& ppr, const Graph& graph);
void BiPartialPush(int src, double alpha, double eps, double feps, double beps, std::vector<double>& ppr, const Graph& graph);
void fwdPush(int src, double alpha, double eps, double maxpr, std::vector<double>& ppr, const Graph& graph);
void RoughBiPartialPush(int src, double alpha, double eps, double delta, double gamma, std::vector<double>& sppr, const Graph& graph);