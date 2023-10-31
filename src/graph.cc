/*************************************************************************
    > File Name: graph.cc
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:20:50 PM
 ************************************************************************/
#include "graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstdlib>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

// #include "mtwist.h"

#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}

using namespace std;

void handle_error(const char* msg) {
	perror(msg);
	exit(255);
}

bool maxScoreCmp(const pair<uint, double>& a, const pair<uint, double>& b){
	if(a.second==b.second){
		return a.first > b.first;
	}
	else{
		return a.second > b.second;
	}
}

Graph::Graph(const string& t_folder, const string& t_graph): 
            m_folder(t_folder), m_graph(t_graph) {
    // mt_goodseed();
    // xorshifinit();
    readNM();
    cout << "graph info: nu=" << this->m_nu << " nv=" << this->m_nv << " m=" << this->m_m << " max-pr=" << this->m_pr << endl;

    this->m_uedges = std::vector<std::vector<std::pair<uint, double>>>(this->m_nu, std::vector<std::pair<uint, double>>());
    this->m_uedges_f = std::vector<std::vector<std::pair<uint, double>>>(this->m_nu, std::vector<std::pair<uint, double>>());
    this->m_udeg = std::vector<uint>(m_nu, 0);
	this->m_uwsum = std::vector<double>(m_nu, 0);
    this->m_vedges = std::vector<std::vector<std::pair<uint, double>>>(this->m_nv, std::vector<std::pair<uint, double>>());
    this->m_vedges_f = std::vector<std::vector<std::pair<uint, double>>>(this->m_nv, std::vector<std::pair<uint, double>>());
    this->m_vdeg = std::vector<uint>(m_nv, 0);
	this->m_vwsum = std::vector<double>(m_nv, 0);

    cout << "loading graph..." << endl;
    readGraph();
    this->m_muwsum = *std::max_element(this->m_uwsum.begin(),this->m_uwsum.end());
    this->m_mvwsum = *std::max_element(this->m_vwsum.begin(),this->m_vwsum.end());
}

void Graph::readNM(){
    ifstream fin((m_folder + "/" + m_graph + "/stat.txt").c_str());
    string s;
    if (fin.is_open()){
        while (fin >> s){
            if (s.substr(0, 2) == "u="){
                this->m_nu = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "v="){
                this->m_nv = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "m="){
                this->m_m = atoi(s.substr(2).c_str());
                continue;
            }
	    if (s.substr(0, 2) == "p="){
		this->m_pr = atof(s.substr(2).c_str());
		continue;
	    }
        }
        fin.close();
    }
    else handle_error("Fail to open attribute file!");
}

void Graph::readGraph(){
    FILE *fin = fopen((m_folder + "/" + m_graph + "/graph.txt.new").c_str(), "r");
    uint64 readCnt = 0;
    uint u, v;
	double w;
    while (fscanf(fin, "%d%d%lf", &u, &v, &w) != EOF) {
        // cout << u << " " << v << " " << w << endl; ll db
        readCnt++;
        ASSERT( u < this->m_nu );
        ASSERT( v < this->m_nv );
        // if(isnan(w) || w<0) w = 1.0;
        addEdge(u, v, (double)w);
        // cout << u << " " << v << " " << w << endl;
    }
    fclose(fin);

    // ifstream infile(m_folder + "/" + m_graph + "/graph.txt");
    // uint u, v;
    // int w;
    // ASSERT(readCnt == this->m_m);
    cout << "read Graph Done!" << endl;

    /*
    for(uint i=0;i<m_nu;i++){
        for(uint j=0;j<m_udeg[i];j++){
            if(i==2){
                cout<<i<<" "<<j<<" "<<m_uedges[i][j].first<<" "<<m_uedges[i][j].second<<" "<<m_vwsum[m_uedges[i][j].first]<<" "<<m_uedges[i][j].second / m_vwsum[m_uedges[i][j].first]<<endl;
            }
        }
    }
    */
    
    cout<<"Backward: Sort the adjacency list of each node in U with the ascending order of the w_uv / wv ..."<<endl;
    for(uint i=0;i<m_nu;i++){
        pair<uint, double> *templist;
        templist=new pair<uint, double>[m_udeg[i]];
        for(uint j=0;j<m_udeg[i];j++){
            uint tmpNbr=m_uedges[i][j].first;
            templist[j]=MP(tmpNbr, m_uedges[i][j].second / m_vwsum[tmpNbr]);
        }
        sort(templist,(templist+m_udeg[i]), maxScoreCmp);
        for(uint j=0;j<m_udeg[i];j++){
            m_uedges[i][j]=MP(templist[j].first, m_vwsum[templist[j].first] * templist[j].second);
            //if(i==2){
            //    cout<<i<<" "<<j<<" "<<m_uedges[i][j].first<<" "<<m_uedges[i][j].second<<" "<<m_vwsum[m_uedges[i][j].first]<<" "<<m_uedges[i][j].second / m_vwsum[m_uedges[i][j].first]<<endl;
            //}
        }
    }

    cout<<"Backward: Sort the adjacency list of each node in V with the ascending order of the w_vu / wu ..."<<endl;
    for(uint i=0;i<m_nv;i++){
        pair<uint, double> *templist;
        templist=new pair<uint, double>[m_vdeg[i]];
        for(uint j=0;j<m_vdeg[i];j++){
            uint tmpNbr=m_vedges[i][j].first;
            templist[j]=MP(tmpNbr, m_vedges[i][j].second / m_uwsum[tmpNbr]);
        }
        sort(templist,(templist+m_vdeg[i]), maxScoreCmp);
        for(uint j=0;j<m_vdeg[i];j++){
            m_vedges[i][j]=MP(templist[j].first, m_uwsum[templist[j].first] * templist[j].second);
            //if(i==2){
            //    cout<<i<<" "<<j<<" "<<m_vedges[i][j].first<<" "<<m_vedges[i][j].second<<" "<<m_uwsum[m_vedges[i][j].first]<<" "<<m_vedges[i][j].second / m_uwsum[m_vedges[i][j].first]<<endl;
            //}
        }
    }

    cout<<"Forward: Sort the adjacency list of each node in U with the ascending order of the w_uv / wv ..."<<endl;
    for(uint i=0;i<m_nu;i++){
        pair<uint, double> *templist;
        templist=new pair<uint, double>[m_udeg[i]];
        for(uint j=0;j<m_udeg[i];j++){
            uint tmpNbr=m_uedges_f[i][j].first;
            templist[j]=MP(tmpNbr, m_uedges_f[i][j].second);
        }
        sort(templist,(templist+m_udeg[i]), maxScoreCmp);
        for(uint j=0;j<m_udeg[i];j++){
            m_uedges_f[i][j]=MP(templist[j].first, templist[j].second);
            //if(i==2){
            //    cout<<i<<" "<<j<<" "<<m_uedges[i][j].first<<" "<<m_uedges[i][j].second<<" "<<m_vwsum[m_uedges[i][j].first]<<" "<<m_uedges[i][j].second / m_vwsum[m_uedges[i][j].first]<<endl;
            //}
        }
    }

    cout<<"Forward: Sort the adjacency list of each node in V with the ascending order of the w_vu / wu ..."<<endl;
    for(uint i=0;i<m_nv;i++){
        pair<uint, double> *templist;
        templist=new pair<uint, double>[m_vdeg[i]];
        for(uint j=0;j<m_vdeg[i];j++){
            uint tmpNbr=m_vedges_f[i][j].first;
            templist[j]=MP(tmpNbr, m_vedges_f[i][j].second);
        }
        sort(templist,(templist+m_vdeg[i]), maxScoreCmp);
        for(uint j=0;j<m_vdeg[i];j++){
            m_vedges_f[i][j]=MP(templist[j].first, templist[j].second);
            //if(i==2){
            //    cout<<i<<" "<<j<<" "<<m_uedges[i][j].first<<" "<<m_uedges[i][j].second<<" "<<m_vwsum[m_uedges[i][j].first]<<" "<<m_uedges[i][j].second / m_vwsum[m_uedges[i][j].first]<<endl;
            //}
        }
    }
    
}

void Graph::addEdge(uint u, uint v, double w){
    // cout << u << " " << v << " " << w << endl;
    this->m_uedges[u].push_back(MP(v,w)); // to do : sort the adjacency list with weights.
    this->m_vedges[v].push_back(MP(u,w));
    this->m_uedges_f[u].push_back(MP(v,w)); // to do : sort the adjacency list with weights.
    this->m_vedges_f[v].push_back(MP(u,w));
    this->m_udeg[u]+=1;
    this->m_vdeg[v]+=1;
    this->m_uwsum[u]+=w;
    this->m_vwsum[v]+=w;
}

uint Graph::getUDeg(uint u) const{
    return this->m_udeg[u];
}

uint Graph::getVDeg(uint v) const{
    return this->m_vdeg[v];
}

double Graph::getUWsum(uint u) const{
    return this->m_uwsum[u];
}

double Graph::getVWsum(uint v) const{
    return this->m_vwsum[v];
}

uint Graph::getNu() const{
    return this->m_nu;
}

uint Graph::getNv() const{
    return this->m_nv;
}

float Graph::getMaxPR() const{
    return this->m_pr;
}

double Graph::getMaxUwsum() const{
    return this->m_muwsum;
}

uint64 Graph::getM() const{
    return this->m_m;
}

std::string Graph::getGraphFolder() const{
    std::string folder(m_folder + "/" + m_graph + "/");
    return folder;
}
