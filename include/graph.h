/*************************************************************************
    > File Name: graph.h
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:20:29 PM
 ************************************************************************/

#include<iostream>

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <unordered_map>
#include "utils.h"

class Graph{
    public:
        std::vector<std::vector<std::pair<uint, double>>> m_uedges; // weighted edges of each node in U
		std::vector<std::vector<std::pair<uint, double>>> m_vedges; // weighted edges of each node in V
        std::vector<std::vector<std::pair<uint, double>>> m_uedges_f; // weighted edges of each node in U
		std::vector<std::vector<std::pair<uint, double>>> m_vedges_f; // weighted edges of each node in V
        std::vector<uint> m_udeg;     // degree of each node in U
		std::vector<uint> m_vdeg;     // degree of each node in V
        std::vector<double> m_uwsum;     // edge weight sum of each node in U
		std::vector<double> m_vwsum;     // edge weight sum of each node in V
    private:
        std::string m_folder;
        std::string m_graph;
        uint m_nu; // the number of nodes in U
        uint m_nv; // the number of nodes in V
        uint64 m_m; // the number of edges
        float m_pr; // the max pagerank*nu
        double m_muwsum;
        double m_mvwsum;

        void readNM();
        void readGraph();
        void addEdge(uint u, uint v, double w);
    public:
        uint getUDeg(uint u) const;
		uint getVDeg(uint v) const;
        double getUWsum(uint u) const;
		double getVWsum(uint v) const;
        uint64 getM() const;
        uint getNu() const;
        uint getNv() const;
        float getMaxPR() const;
        double getMaxUwsum() const;
        std::string getGraphFolder() const;
        // const std::vector<std::pair<uint, double>>& operator [] (uint u) const;

        Graph(const std::string& t_folder, const std::string& t_graph);
};
#endif
