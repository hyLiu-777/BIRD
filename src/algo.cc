/*************************************************************************
    > File Name: algo.cc
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:26:32 PM
 ************************************************************************/

#include "algo.h"
#include <math.h>
#include <iostream>
#include <queue>
#include <set>
#include <list>
#include <algorithm>
#include <random>
#include <alias.h>

using namespace std;


void monteCarlo(int src, double alpha, int64 n_walk, std::vector<double>& sppr, const Graph& graph){
    vector<double> ppr(graph.getNu());

    for(int i=0;i<n_walk;i++){
    	int cur = src;
        while(true){
        	if(graph.m_udeg[cur]==0){
        		break;
        	}
        	double stop = (double)rand()/(double)RAND_MAX;
            if(stop<=alpha){
                break;
            }
            else{
            	int k = rand()%graph.m_udeg[cur];
            	cur = graph.m_uedges[cur][k].first;
            	k = rand()%graph.m_vdeg[cur];
            	cur = graph.m_vedges[cur][k].first;
            }
        }
        ppr[cur] += 1;
    }

	for(uint u_i = 0; u_i < graph.getNu(); u_i++ ){
		sppr[u_i] += (double)ppr[u_i] / (double)n_walk;
	}
}

void monteCarloWeighted(int src, double alpha, int64 n_walk, std::vector<double>& sppr, const Graph& graph, vector<Alias>& U_alias, vector<Alias>& V_alias){
    
	vector<double> ppr(graph.getNu());
    for(int i=0;i<n_walk;i++){
    	int cur = src;
        while(true){
        	if(graph.m_udeg[cur]==0){
        		break;
        	}
        	double stop = (double)rand()/(double)RAND_MAX;
            if(stop<=alpha){
                break;
            }
            else{
            	int k = U_alias[cur].generateRandom(); 
            	cur = k;
				k = V_alias[cur].generateRandom();
            	// k = rand()%graph.m_vdeg[cur];
            	cur = k;
            }
        }
        ppr[cur] += 1;
    }
	for(uint u_i = 0; u_i < graph.getNu(); u_i++ ){
		sppr[u_i] += (double)ppr[u_i] / (double)n_walk;
	}
}

void powerIter(int src, double alpha, uint n_iter, std::vector<double>& sppr, const Graph& graph){
    vector<double> vecUResidue(graph.getNu());
    vector<double> vecVResidue(graph.getNv());
    vector<double> ppr(graph.getNu());

    uint nu = graph.getNu();
    uint nv = graph.getNv();

    vecUResidue[src]= 1.0;
    set<uint> Qu;
    set<uint> Qv;
    Qu.insert(src);

    Timer tm(2, "pi");
    for(int i=0;i<n_iter;i++){
    	for(std::set<uint>::iterator it=Qu.begin(); it!=Qu.end(); ++it){
	        uint u = *it;

	        double residue = vecUResidue[u];
			vecUResidue[u]=0;
			if(graph.m_udeg[u]>0){
				ppr[u] += alpha*residue;
				residue = (1-alpha)*residue;
				for(const auto& p: graph.m_uedges[u]){
					const uint v = p.first;
					const double w = p.second;
					Qv.insert(v);
					vecVResidue[v] += residue*w/(double)graph.m_uwsum[u];
				}
			}
			else{
				ppr[u] += residue;
			}
    	}
    	Qu.clear();

    	for(std::set<uint>::iterator it=Qv.begin(); it!=Qv.end(); ++it){
	        uint v = *it;

	        double residue = vecVResidue[v];
			vecVResidue[v]=0;
			if(graph.m_vdeg[v]>0){
				for(const auto& p: graph.m_vedges[v]){
					const uint u = p.first;
					const double w = p.second;
					Qu.insert(u);
					vecUResidue[u] += residue*w/(double)graph.m_vwsum[v];
				}
			}
    	}
    	Qv.clear();
    }

	for(uint u_i = 0; u_i < graph.getNu(); u_i++ ){
		sppr[u_i] += ppr[u_i];
	}
}

void powerIterGt(int src, double alpha, uint n_iter, std::vector<double>& sppr, const Graph& graph){
    vector<double> vecUResidue(graph.getNu());
    vector<double> vecVResidue(graph.getNv());
    vector<double> ppr(graph.getNu());

    uint nu = graph.getNu();
    uint nv = graph.getNv();

    vecUResidue[src]= 1.0;
    set<uint> Qu;
    set<uint> Qv;
    Qu.insert(src);

    Timer tm(2, "pi");
    for(int i=0;i<n_iter;i++){
    	for(std::set<uint>::iterator it=Qu.begin(); it!=Qu.end(); ++it){
	        uint u = *it;

	        double residue = vecUResidue[u];
			vecUResidue[u]=0;
			if(graph.m_udeg[u]>0){
				ppr[u] += alpha*residue;
				residue = (1-alpha)*residue;
				for(const auto& p: graph.m_uedges[u]){
					const uint v = p.first;
					const double w = p.second;
					Qv.insert(v);
					vecVResidue[v] += residue*w/(double)graph.m_uwsum[u];
				}
			}
			else{
				ppr[u] += residue;
			}
    	}
    	Qu.clear();

    	for(std::set<uint>::iterator it=Qv.begin(); it!=Qv.end(); ++it){
	        uint v = *it;

	        double residue = vecVResidue[v];
			vecVResidue[v]=0;
			if(graph.m_vdeg[v]>0){
				for(const auto& p: graph.m_vedges[v]){
					const uint u = p.first;
					const double w = p.second;
					Qu.insert(u);
					vecUResidue[u] += residue*w/(double)graph.m_vwsum[v];
				}
			}
    	}
    	Qv.clear();
    }

	for(uint u_i = 0; u_i < graph.getNu(); u_i++ ){
		sppr[u_i] += ppr[u_i] * (1.0 + (double) graph.m_uwsum[src] / (double) graph.m_uwsum[u_i]); // only PI
	}
}


void fwdPush(int src, double alpha, double eps, double maxpr, std::vector<double>& sppr, const Graph& graph){
    vector<double> vecUResidue(graph.getNu());
    vector<double> vecVResidue(graph.getNv());
    vector<double> ppr(graph.getNu());

    uint nu = graph.getNu();
    uint nv = graph.getNv();

    vecUResidue[src]= 1.0;
    set<uint> Qu;
    set<uint> Qv;
    Qu.insert(src);
	double thresh = eps/maxpr;
{   
	Timer tm(4, "fwd");
	uint itr=0;
    while(Qu.size()>0 || Qv.size()>0){
        itr++;
    	for(std::set<uint>::iterator it=Qu.begin(); it!=Qu.end(); ++it){
	        uint u = *it;

			double residue = vecUResidue[u];
			vecUResidue[u]=0;
			if(graph.m_udeg[u]>0){
				ppr[u] += alpha*residue;
				residue = (1-alpha)*residue;
				for(const auto& p: graph.m_uedges[u]){
					const uint v = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_uwsum[u];
					Qv.insert(v);
					vecVResidue[v] += mass;
				}
			}
			else{
				ppr[u] += residue;
			}
    	}
    	Qu.clear();


    	for(std::set<uint>::iterator it=Qv.begin(); it!=Qv.end(); ++it){
	        uint v = *it;

	        double residue = vecVResidue[v];
			vecVResidue[v]=0;
			if(graph.m_vdeg[v]>0){
				for(const auto& p: graph.m_vedges[v]){
					const uint u = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_vwsum[v];
					if(vecUResidue[u]<thresh && vecUResidue[u]+mass>=thresh){
						Qu.insert(u);
					}
					vecUResidue[u] += mass;
				}
			}
    	}
    	Qv.clear();
    }

    double rss=0;
    double rmax=0;
    double rmaxw=0;
    double rsum=0;
    for(int u=0; u<graph.getNu(); u++){
    	if(vecUResidue[u]>0){
    		rss+=vecUResidue[u]*vecUResidue[u];
    		if(vecUResidue[u]>rmax){
    			rmax=vecUResidue[u];
    		}
    		rsum+=vecUResidue[u];
    		if(vecUResidue[u]/graph.getUWsum(u)>rmaxw){
    			rmaxw=vecUResidue[u]/graph.getUWsum(u);
    		}
    	}
    }
    double err1 = sqrt(maxpr*rss);
    double err2 = rmax*maxpr;
    double err3 = rsum;
    double err4 = rmaxw*graph.getMaxUwsum();
    cout << "residue max:" << rmax << " rmax*maxpr:" << " max-pr:" << maxpr << err2 << " epsilon:" << eps << endl;
    cout << "residue/wsum max:" << rmaxw << " maxusum:" << graph.getMaxUwsum() << " error:" << err4 << " epsilon:" << eps << endl;
}	
}


void reversePush(int tgt, double alpha, double eps, std::vector<double>& sppr, const Graph& graph){
    vector<double> vecUResidue(graph.getNu());
    vector<double> vecVResidue(graph.getNv());
    vector<double> ppr(graph.getNu());

    uint nu = graph.getNu();
    uint nv = graph.getNv();

    vecUResidue[tgt]= 1.0;
    set<uint> Qu;
    set<uint> Qv;
    Qu.insert(tgt);
	double thresh = eps;

	Timer tm(3, "bwd");
    while(Qu.size()>0 || Qv.size()>0){
    	for(std::set<uint>::iterator it=Qu.begin(); it!=Qu.end(); ++it){
	        uint u = *it;

	        double residue = vecUResidue[u];
			vecUResidue[u]=0;
			if(graph.m_udeg[u]>0){
				ppr[u] += alpha*residue;
				residue = (1-alpha)*residue;
				for(const auto& p: graph.m_uedges[u]){
					const uint v = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_vwsum[v];
					Qv.insert(v);
					vecVResidue[v] += mass;
				}
			}
			else{
				ppr[u] += residue;
			}
    	}
    	Qu.clear();

    	for(std::set<uint>::iterator it=Qv.begin(); it!=Qv.end(); ++it){
	        uint v = *it;

	        double residue = vecVResidue[v];
			vecVResidue[v]=0;
			if(graph.m_vdeg[v]>0){
				for(const auto& p: graph.m_vedges[v]){
					const uint u = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_uwsum[u];
					if(vecUResidue[u]<thresh && vecUResidue[u]+mass>=thresh){
						Qu.insert(u);
					}
					vecUResidue[u] += mass;
				}
			}
    	}
    	Qv.clear();
    }
	for(uint u_i = 0; u_i < graph.getNu(); u_i++ ){
		sppr[u_i] += ppr[u_i];
	}
}


void BiPartialPush(int src, double alpha, double eps, double feps, double beps, std::vector<double>& sppr, const Graph& graph){
    vector<double> vecUResidue(graph.getNu());
    vector<double> vecVResidue(graph.getNv());
    vector<double> ppr(graph.getNu());

    uint nu = graph.getNu();
    uint nv = graph.getNv();

    vecUResidue[src]= 1.0;
    set<uint> Qu;
    set<uint> Qv;
    Qu.insert(src);
	double thresh = beps;

{
	Timer tm(3, "bwd");
	uint iter=0;
    while(Qu.size()>0 || Qv.size()>0){
    	iter++;
    	for(std::set<uint>::iterator it=Qu.begin(); it!=Qu.end(); ++it){
	        uint u = *it;

			double residue = vecUResidue[u];
			vecUResidue[u]=0;
			if(graph.m_udeg[u]>0){
				ppr[u] += alpha*residue;
				residue = (1-alpha)*residue;
				for(const auto& p: graph.m_uedges[u]){
					const uint v = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_vwsum[v];
					Qv.insert(v);
					vecVResidue[v] += mass;
				}
			}
			else{
				ppr[u] += residue;
			}
    	}
    	Qu.clear();

    	for(std::set<uint>::iterator it=Qv.begin(); it!=Qv.end(); ++it){
	        uint v = *it;

	        double residue = vecVResidue[v];
			vecVResidue[v]=0;
			if(graph.m_vdeg[v]>0){
				for(const auto& p: graph.m_vedges[v]){
					const uint u = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_uwsum[u];
					if(vecUResidue[u]<thresh && vecUResidue[u]+mass>=thresh){
						Qu.insert(u);
					}
					vecUResidue[u] += mass;
				}
			}
    	}
    	Qv.clear();
    }
}

    thresh = feps;
    double frsum=0;
    double oeps=1.0;
    double frsummax = thresh/oeps;
{
	Timer tm(2, "bwd to fwd");
    for(int u=0; u<graph.getNu(); u++){
    	if(ppr[u]>0){
    		ppr[u]+=graph.getUWsum(u)*ppr[u]/graph.getUWsum(src);
    	}
    	if(vecUResidue[u]>0){
	    	vecUResidue[u]=graph.getUWsum(u)*vecUResidue[u]/graph.getUWsum(src);
	    	frsum+=vecUResidue[u];
	    	if(vecUResidue[u]>thresh){
	    		Qu.insert(u);
	    	}
    	}
    }
}

{   
	Timer tm(4, "fwd");
	uint itr=0;
    while(Qu.size()>0 || Qv.size()>0){
    	if(frsum<=frsummax){
    		break;
    	}
        itr++;
    	for(std::set<uint>::iterator it=Qu.begin(); it!=Qu.end(); ++it){
	        uint u = *it;

			double residue = vecUResidue[u];
			vecUResidue[u]=0;
			if(graph.m_udeg[u]>0){
				ppr[u] += alpha*residue;
				frsum -= alpha*residue;
				residue = (1-alpha)*residue;
				for(const auto& p: graph.m_uedges[u]){
					const uint v = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_uwsum[u];
					Qv.insert(v);
					vecVResidue[v] += mass;
				}
			}
			else{
				ppr[u] += residue;
				frsum -= residue;
			}
    	}
    	Qu.clear();


    	for(std::set<uint>::iterator it=Qv.begin(); it!=Qv.end(); ++it){
	        uint v = *it;

	        double residue = vecVResidue[v];
			vecVResidue[v]=0;
			if(graph.m_vdeg[v]>0){
				for(const auto& p: graph.m_vedges[v]){
					const uint u = p.first;
					const double w = p.second;
					double mass = residue*w/(double)graph.m_vwsum[v];
					if(vecUResidue[u]<thresh && vecUResidue[u]+mass>=thresh){
						Qu.insert(u);
					}
					vecUResidue[u] += mass;
				}
			}
    	}
    	Qv.clear();
    }
}
	for(uint u = 0; u < nu; u++ ){
		sppr[u] = ppr[u];
	}
}


void RoughBiPartialPush(int src, double alpha, double eps, double delta, double gamma, std::vector<double>& sppr, const Graph& graph){

    uint nu = graph.getNu();
    uint nv = graph.getNv();
	
	// backward variables //
    vector<vector<double>> vecUResidueBack(2, vector<double>(nu, 0));
	vector<vector<int>> candidateUSetBack(2, vector<int>(nu, 0));
	vector<vector<int>> flagUBack(2, vector<int>(nu,0));
	vector<int> candidateUCountBack(2, 0); 
	vector<double> vecVResidueBack(nv, 0);
	vector<int> candidateVSetBack(nv, 0);
	vector<int> flagVBack(nv, 0);
	uint candidateVCountBack = 0;

	// forward variables //
    vector<vector<double>> vecUResidueFor(2, vector<double>(nu, 0));
	vector<vector<int>> candidateUSetFor(2, vector<int>(nu, 0));
	vector<vector<int>> flagUFor(2, vector<int>(nu,0));
	vector<int> candidateUCountFor(2, 0); 
	vector<double> vecVResidueFor(nv, 0);
	vector<int> candidateVSetFor(nv, 0);
	vector<int> flagVFor(nv, 0);
	uint candidateVCountFor = 0;

	vector<int> U_gamma(nu, 0); // node set shares computation. safe node set.
	vector<int> V_I(nv, 0);
	vector<int> U_gamma_left(nu, 0); // node set needs to push in each iteration.
	uint nu_gamma_left = 0;
	uint nu_gamma = 0;
	uint nv_i = 0;

	// cout << graph.m_uedges[27736][0].first <<endl;
	// common variables //
	vector<double> finalReserve(nu, 0);
	double temp_thre = (double)graph.m_uwsum[src] * gamma;
	// cout << "weight threshod: " << temp_thre <<", source weight:" << (double)graph.m_uwsum[src] << endl;
	// initialize U_gamma, V_I.
	for(uint s=0; s<nu; s++){
		if((double)graph.m_uwsum[s] <= temp_thre){ // safe nodes.
			// cout << "current node: " << s << ", weight: " << (double)graph.m_uwsum[s] << endl;
			if (U_gamma[s] == 0)
			{
				nu_gamma += 1;
			}
			
			U_gamma[s] = 1;
			// cout << "current node: " << s << ", weight: " << (double)graph.m_uwsum[s] << endl;
		}else{ // unsafe ones need to push.
			for (const auto& p: graph.m_uedges[s])
			{
				const uint v_j = p.first;
				// cout << v_j << endl;
				if(V_I[v_j]==0){
					nv_i += 1;
				}
				V_I[v_j] = 1;
			}
		}
	}

	cout << "total safe nodes number: " << nu_gamma << endl;
	// cout << "total to push V nodes number: " << nv_i << endl;
	// initialize U_gamma_left.
	for(uint s=0; s<nv; s++){
		// cout << s << endl;
		if (V_I[s] == 0)
		{
			continue;
		}
		
		for (const auto& p: graph.m_vedges[s])
		{
			const uint u_i = p.first;
			if (U_gamma_left[u_i] == 0)
			{
				nu_gamma_left += 1;
			}
			U_gamma_left[u_i] = 1;
		}
	}

	cout << "total to push U nodes number: " << nu_gamma_left << endl;

	// initialization //
	uint tempLevel = 0; // iteration number.
	// uint L = 5;
	uint L = (uint)ceil(log(eps/(double)nu)/log(1-alpha))+1;
	cout<<"Ineration amount L: "<<L<<endl;
	// double delta = 1.0 / nu; // delta is set to 1/n.
	double theta = eps*eps*delta/L/6.0/gamma/gamma;
	// cout<<"theta: "<<theta<<endl;

    vecUResidueBack[0][src] = 1;
	candidateUSetBack[0][0] = src;
	candidateUCountBack[0] = 1;

	 
    vecUResidueFor[0][src] = 1;
	candidateUSetFor[0][0] = src;
	candidateUCountFor[0] = 1;
	

	// cout<<"start main iteration"<<endl;
	// Main algorithm iteration.
{	
	Timer tm(2, "bwd");
	while(tempLevel < L){

		uint tempLevelID=tempLevel%2; // l
		uint newLevelID=(tempLevel+1)%2; // l+1

		// ********************************************* //
		// backward start.
		uint candidateUCntBack=candidateUCountBack[tempLevelID];
        // cout<<"Iteration "<<tempLevel<<": candidateUCntBack="<<candidateUCntBack<<endl;
		if(candidateUCntBack==0){
			//cout<<"candidateCnt=0 tempLevel="<<tempLevel<<endl;
			break;
		}
		candidateUCountBack[tempLevelID]=0;
		
		// ********************************************* //
		// backward: push from U to V.
		for(uint j = 0; j < candidateUCntBack; j++){
			uint tempNode = candidateUSetBack[tempLevelID][j];
			double tempR = vecUResidueBack[tempLevelID][tempNode];

			if (U_gamma[tempNode]==1)
			{
				vecUResidueFor[tempLevelID][tempNode] = tempR * graph.m_uwsum[tempNode] / graph.m_uwsum[src];
				// cout << tempNode << endl;
				if (U_gamma_left[tempNode]==1){
					flagUFor[tempLevelID][tempNode] = 1;
					candidateUSetFor[tempLevelID][candidateUCountFor[tempLevelID]++] = tempNode;
				}
			}
			

			flagUBack[tempLevelID][tempNode] = 0;
			vecUResidueBack[tempLevelID][tempNode] = 0; //push
			finalReserve[tempNode] += alpha * tempR; // add the reserve to result
			
			if(tempLevel>=L){
				continue;
			}
			
			if(graph.m_udeg[tempNode]>0){
				double ran = (double)rand()/(double)RAND_MAX;
				tempR = (1-alpha)*tempR;
				for(const auto& p: graph.m_uedges[tempNode]){ 
					const uint v_j = p.first;
					const double w = p.second;
					double mass = tempR*w/(double)graph.m_vwsum[v_j];
					if (mass > theta)
					{
						vecVResidueBack[v_j] += mass;
					}else{
						if(mass >= ran*theta){
							vecVResidueBack[v_j] += theta;
						}else{
							break; // according to the pre-ordering.
						}
					}
					// cout << tempNode << ":" << v_j << "," << flagVBack[v_j]<< "," << vecVResidueBack[v_j] << "," << ((flagVBack[v_j] == 0) && (vecVResidueBack[v_j] > 0)) << endl;
					if((flagVBack[v_j] == 0) && (vecVResidueBack[v_j] > 0))
					{
						flagVBack[v_j] = 1;
						candidateVSetBack[candidateVCountBack++] = v_j; // avoid repetitive nodes.
						// cout<<v_j<<endl;
					}
				}
			}
		}


		// backward: push from V to U.
		uint candidateVCntBack=candidateVCountBack;
        // cout<<"Iteration "<<tempLevel<<": candidateVCntBack="<<candidateVCntBack<<endl;
		candidateVCountBack = 0;
		for(uint j = 0; j < candidateVCntBack; j++){
			uint tempNode = candidateVSetBack[j];
			double tempR = vecVResidueBack[tempNode];
			flagVBack[tempNode] = 0;
			vecVResidueBack[tempNode] = 0; //push
			
			if(graph.m_vdeg[tempNode]>0){ // is "IF" costs unavoidable time? 
				double ran = (double)rand()/(double)RAND_MAX;
				for(const auto& p: graph.m_vedges[tempNode]){ 
					const uint u_j = p.first;
					const double w = p.second;
					double mass = tempR*w/(double)graph.m_uwsum[u_j];
					if (mass > theta)
					{
						vecUResidueBack[newLevelID][u_j] += mass;
					}else{
						if(mass >= ran*theta){
							vecUResidueBack[newLevelID][u_j] += theta;
						}else{
							break; // according to the pre-ordering.
						}
					}
					if((flagUBack[newLevelID][u_j] == 0) && (vecUResidueBack[newLevelID][u_j]) > theta)
					{
						flagUBack[newLevelID][u_j] = 1;
						candidateUSetBack[newLevelID][candidateUCountBack[newLevelID]++] = u_j; // avoid repetitive nodes.
					}
				}
			}
		}


		
		// *********************************************//
		// forward: reuse from U to V. 
		
		uint candidateUCntFor=candidateUCountFor[tempLevelID];
		// cout << "	forward U->V: " << candidateUCntFor << endl;
		if(candidateUCntFor==0){
			// cout<<"candidateCntFor=0 tempLevel="<<tempLevel<<endl;
			tempLevel++;
			continue;
		}
		// cout << "1" << endl;
		candidateUCountFor[tempLevelID]=0;
		for(uint j = 0; j < candidateUCntFor; j++){
			uint tempNode = candidateUSetFor[tempLevelID][j];
			// cout << tempNode << ": " << graph.m_udeg[tempNode] << endl;

			double tempR = vecUResidueFor[tempLevelID][tempNode];
			// cout << "ori tempR:" << tempR << endl;
			flagUFor[tempLevelID][tempNode] = 0;
			vecUResidueFor[tempLevelID][tempNode] = 0; //push
			finalReserve[tempNode] += alpha * tempR; // add the reserve to result

			if(tempLevel>=L){
				continue;
			}

			if(graph.m_udeg[tempNode]>0){
				double ran = (double)rand()/(double)RAND_MAX; //
				tempR = (1-alpha)*tempR/(double)graph.m_uwsum[tempNode];
				// cout << "tempR: " << tempR << endl;
				for(const auto& p: graph.m_uedges[tempNode]){ 
					const uint v_j = p.first;
					// cout << v_j << endl;
					if(V_I[v_j] == 0){ // no need to push, not in target set.
					// cout << "no need to push" << endl;
						continue;
					}

					const double w = p.second;
					double mass = tempR*w;
					// cout << "mass: " << mass << endl;
					if (mass > theta)
					{
						// cout << "1: determine" << endl;
						vecVResidueFor[v_j] += mass;
					}else{
						if(mass >= ran*theta){
							// cout << "2: random" << endl;
							vecVResidueFor[v_j] += theta;
						}else{
							// cout << "3: nothing" << endl;
							break; // according to the pre-ordering.
						}
					}
					if((flagVFor[v_j] == 0) && (vecVResidueFor[v_j] > 0))
					{
						flagVFor[v_j] = 1;
						candidateVSetFor[candidateVCountFor++] = v_j; // avoid repetitive nodes.
						// cout<<v_j<<endl;
					}
				}
			}
		}
		

		// forward: reuse from V to U.
		uint candidateVCntFor=candidateVCountFor;
		candidateVCountFor = 0;
		// cout << "forward V->U: " << candidateVCntFor << endl;

		for(uint j = 0; j < candidateVCntFor; j++){
			// cout << j << endl;
			uint tempNode = candidateVSetFor[j];
			double tempR = vecVResidueFor[tempNode];
			flagVFor[tempNode] = 0;
			vecVResidueFor[tempNode] = 0; //push
			
			if(graph.m_vdeg[tempNode]>0){ // is "IF" costs unavoidable time? 
				double ran = (double)rand()/(double)RAND_MAX;
				tempR = tempR/(double)graph.m_vwsum[tempNode];

				for(const auto& p: graph.m_vedges[tempNode]){ 
					const uint u_j = p.first;
					
					if (U_gamma[u_j]==1) // safe nodes reuse the estimation, don't push.
					{
						continue;
					}
					
					const double w = p.second;
					double mass = tempR*w;
					if (mass > theta)
					{
						vecUResidueFor[newLevelID][u_j] += mass;
					}else{
						if(mass >= ran*theta){
							vecUResidueFor[newLevelID][u_j] += theta;
						}else{
							break; // according to the pre-ordering.
						}
					}
					if((flagUFor[newLevelID][u_j] == 0) && (vecUResidueFor[newLevelID][u_j]) > theta)
					{
						flagUFor[newLevelID][u_j] = 1;
						candidateUSetFor[newLevelID][candidateUCountFor[newLevelID]++] = u_j; // avoid repetitive nodes.
					}
				}
			}
		}
	// cout<<"Finished"<<endl;
	tempLevel++;
	}
cout<<"finished in iteration: "<< tempLevel << endl;
}
for(uint u = 0; u < nu; u++ ){
	sppr[u] = finalReserve[u];
}
}