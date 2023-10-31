/*************************************************************************
    > File Name: bppr.cc
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:19:39 PM
 ************************************************************************/

#include <iostream>
#include <boost/program_options.hpp>
#include <string>
#include "graph.h"
#include "algo.h"
#include <sys/stat.h>
#include <alias.h>


using namespace std;
namespace po = boost::program_options;

namespace { 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace

Config parseParams(int argc, char** argv){
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("data-folder,f", po::value<string>()->required(), "graph data folder")
        ("graph-name,g", po::value<string>()->required(), "graph file name")
        ("algo,a", po::value<string>()->required(), "algorithm name")
        ("epsilon,e", po::value<double>()->default_value(0), "epsilon")
        ("pfail,p", po::value<double>()->default_value(0), "failure probability")
        ("delta,d", po::value<double>()->default_value(0), "delta")
        ("gamma,ga", po::value<double>()->default_value(1.0), "gamma")
        ("querynum,qn", po::value<int64>()->default_value(10), "querynum")
        ("if_percentile, pen", po::value<uint>()->default_value(0), "if_percentile")
    ;

    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw 
    po::notify(vm);

    Config config;

    if (vm.count("help")){
        cout << desc << '\n';
        exit(0);
    }
    if (vm.count("data-folder")){
        config.strFolder = vm["data-folder"].as<string>();
    }
    if (vm.count("graph-name")){
        config.strGraph = vm["graph-name"].as<string>();
    }
    if (vm.count("algo")){
        config.strAlgo = vm["algo"].as<string>();
    }
    if (vm.count("epsilon")){
        config.epsilon = vm["epsilon"].as<double>();
    }
    if (vm.count("delta")){
        config.delta = vm["delta"].as<double>();
    }
    if (vm.count("gamma")){
        config.gamma = vm["gamma"].as<double>();
    }
    if (vm.count("querynum")){
        config.querynum = vm["querynum"].as<int64>();
    }
    if (vm.count("if_percentile")){
        config.if_percentile = vm["if_percentile"].as<uint>();
    }
    return config;
}

int mkpath(string s, mode_t mode=0755){
	size_t pre=0, pos;
	string dir;
	int mdret;
	if(s[s.size()-1]!='/'){
		s+='/';
	}
	while((pos=s.find_first_of('/',pre))!=string::npos){
		dir=s.substr(0,pos++);
		pre=pos;
		if(dir.size()==0){
			continue;
		}
		if((mdret=::mkdir(dir.c_str(),mode)) && errno!=EEXIST){
			return mdret;
		}
	}
	return mdret;
}

double getPercentile(uint source, const std::vector<double>& vec, double gamma) {

    std::vector<double> sortedVec;
    for(uint i=0; i<vec.size(); i++){
        if(vec[i] > vec[source]){
            sortedVec.push_back(vec[i]);
        }
    }
    std::sort(sortedVec.begin(), sortedVec.end()); // Sort the copy
    
    int pos = gamma * sortedVec.size();
    cout << "pos: " << pos << ", value: " << sortedVec[pos] << "; Max weight: " << sortedVec[sortedVec.size()-1] << endl;
    return sortedVec[pos];

}

vector<int> loadSeed(string folder, string file_name, int count){
    FILE *fin = fopen((folder + "/" + file_name + "/seeds.txt").c_str(), "r");
    int s;
    vector<int> seeds;
    int i=0;
    while (fscanf(fin, "%d", &s) != EOF) {
        seeds.push_back(s);
        i++;
        if(i>=count)
            break;
    }
    fclose(fin);
    cout << "read seed Done!" << endl;
    return seeds;
}

int main(int argc, char **argv){
    Config config;
    try{
        config = parseParams(argc, argv);
        config.check();
        config.setDefault();
    }
    catch (const exception &ex){
        cerr << ex.what() << '\n';
        return ERROR_IN_COMMAND_LINE;
    }

    config.display();

    int query_count = config.querynum;
    Graph graph(config.strFolder, config.strGraph);
    vector<int> seeds_candidates = loadSeed(config.strFolder, config.strGraph, query_count);
    vector<int> seeds;
    for(uint i = 0; i < seeds_candidates.size(); i++){
        if (graph.m_udeg[seeds_candidates[i]] > 6){
            seeds.push_back(seeds_candidates[i]);
        }
    }


    stringstream ss_dir;
    ofstream fout;

    cout << "if_percentile: " << config.if_percentile << endl;

    if (config.if_percentile == 1){ // gamma is between [0,1].
        cout << "using percentile: " << config.gamma << endl;
        config.delta = 1.0 / (double) graph.getNu();
        ss_dir << "result/varry_gamma/" << config.strGraph << "/" << config.strAlgo << "/" << config.gamma << "/";
        mkpath(ss_dir.str());
    }else{
        if (config.delta == 0) // set to default with 1/U.
        {
            config.delta = 1.0 / (double) graph.getNu();
            ss_dir << "result/relative/" << config.strGraph << "/" << config.strAlgo << "/" << config.epsilon << "/";
            mkpath(ss_dir.str());
        }else{
            ss_dir << "result/tradeoff_new/" << config.strGraph << "/" << config.strAlgo << "/" << config.delta << "/";
            mkpath(ss_dir.str());
        }
    }


    if(config.strAlgo==PI){
        cout << "start ground-truth bppr with power-method!" << endl;
        Timer tm(1, "bppr");
        double epsilon = config.epsilon*(double)config.delta/2.0;
        config.setIteration(epsilon);
        cout << "#iterations: " << config.iteration << endl;
        cout << "epsilon: " << epsilon << endl;
        for(const auto& u: seeds){
            std::vector<double> ppr(graph.getNu(), 0);
            double u_weight = (double) graph.m_uwsum[u] / (double) graph.getMaxUwsum();
            config.setIteration(epsilon * u_weight); // ensure the error.
            powerIterGt(u, config.alpha, config.iteration, ppr, graph);
            stringstream ss;
            ss << ss_dir.str() << u << ".txt";
            fout.open(ss.str());
            fout.setf(ios::fixed,ios::floatfield);
            fout.precision(15);
            if(!fout){
                cout<<"Fail to open the writed file"<<endl;
            }
            for(uint u_i=0; u_i < graph.getNu(); u_i++){
                if(ppr[u_i]>1e-8){
                    fout<<u_i<<" "<<ppr[u_i]<<endl;
                }
            }
            fout.close();
            
        }
    }
    else if(config.strAlgo==PISP){
        cout << "start bppr with power-method!" << endl;
        Timer tm(1, "bppr");
        double epsilon = config.epsilon*(double)config.delta/2.0;
        config.setIteration(epsilon);
        cout << "#iterations: " << config.iteration << endl;
        cout << "epsilon: " << epsilon << endl;
        for(const auto& u: seeds){
            std::vector<double> ppr(graph.getNu(), 0);
            powerIter(u, config.alpha, config.iteration, ppr, graph);
            reversePush(u, config.alpha, epsilon, ppr, graph);
            stringstream ss;
            ss << ss_dir.str() << u << ".txt";
            fout.open(ss.str());
            fout.setf(ios::fixed,ios::floatfield);
            fout.precision(15);
            if(!fout){
                cout<<"Fail to open the writed file"<<endl;
            }
            for(uint u_i=0; u_i < graph.getNu(); u_i++){
                if(ppr[u_i]>1e-8){
                    fout<<u_i<<" "<<ppr[u_i]<<endl;
                }
            }
            fout.close();
        
        }
    }
    else if(config.strAlgo==MCSP){
        cout << "start bppr with monte-carlo!" << endl;
        Timer tm(1, "bppr");

        // construct alias tables for sampling random walks.
        cout << "constructing U_alias" << endl;
        vector<Alias> U_alias;
        for(uint i=0; i<graph.getNu(); i++){
            // cout << "i: " << i << endl;
            vector<pair<int, double>> aliasP(graph.m_udeg[i]);
            // cout << "1" << endl;
            for(uint j=0;j<graph.m_udeg[i];j++){
                uint tmpNbr=graph.m_uedges[i][j].first;
                aliasP.push_back(MP(tmpNbr, (double)graph.m_uedges[i][j].second / (double) graph.m_uwsum[i]));
            }
            // cout << "2" << endl;
            Alias alias(aliasP);
            // cout << "3" << endl;
            U_alias.push_back(alias);
            // cout << "4" << endl;
        }
        cout << "constructing V_alias" << endl;
        vector<Alias> V_alias;
        for(uint i=0; i<graph.getNv(); i++){
            vector<pair<int, double>> aliasP(graph.m_vdeg[i]);
            for(uint j=0;j<graph.m_vdeg[i];j++){
                uint tmpNbr=graph.m_vedges[i][j].first;
                aliasP.push_back(MP(tmpNbr, (double)graph.m_vedges[i][j].second / (double) graph.m_vwsum[i]));
            }
            Alias alias(aliasP);
            V_alias.push_back(alias);
        }
        cout << "Finished" << endl;


        double epsilon = config.epsilon*(double)config.delta/2.0;
        config.setNumWalks(config.epsilon/2.0, config.delta);
        cout << "#randomwalks: " << config.numwalks << endl;
        cout << "epsilon: " << epsilon << endl;
        srand(time(NULL));
        for(const auto& u: seeds){
            //monteCarlo(u, config.alpha, config.numwalks, ppr, graph);
            std::vector<double> ppr(graph.getNu(), 0);
            monteCarloWeighted(u, config.alpha, config.numwalks, ppr, graph, U_alias, V_alias);
            reversePush(u, config.alpha, epsilon, ppr, graph);
            stringstream ss;
            ss << ss_dir.str() << u << ".txt";
            fout.open(ss.str());
            fout.setf(ios::fixed,ios::floatfield);
            fout.precision(15);
            if(!fout){
                cout<<"Fail to open the writed file"<<endl;
            }
            for(uint u_i=0; u_i < graph.getNu(); u_i++){
                if(ppr[u_i]>1e-8){
                    fout<<u_i<<" "<<ppr[u_i]<<endl;
                }
            }
            fout.close();

        }
    }
    else if(config.strAlgo==ABHPP){
        cout << "start bppr with abhpp!" << endl;
        Timer tm(1, "bppr");
        double epsilon = config.epsilon * (double)config.delta;
    	double bepsilon = epsilon /(graph.getMaxPR()+1);
    	double fepsilon = epsilon - bepsilon;
    	cout << "fwd rmax=" << fepsilon << endl;
    	cout << "bwd rmax=" << bepsilon << endl;
        for(const auto& u: seeds){
            std::vector<double> ppr(graph.getNu(), 0);
            BiPartialPush(u, config.alpha, config.epsilon, fepsilon, bepsilon, ppr, graph);
            stringstream ss;
            ss << ss_dir.str() << u << ".txt";
            fout.open(ss.str());
            fout.setf(ios::fixed,ios::floatfield);
            fout.precision(15);
            if(!fout){
                cout<<"Fail to open the writed file"<<endl;
            }
            for(uint u_i=0; u_i < graph.getNu(); u_i++){
                if(ppr[u_i]>1e-8){
                    fout<<u_i<<" "<<ppr[u_i]<<endl;
                }
            }
            fout.close();
        }
    }
    else if(config.strAlgo==RBHPP){
        cout << "start bppr with rbhpp!" << endl;
        Timer tm(1, "bppr");
        for(const auto& u: seeds){
            std::vector<double> ppr(graph.getNu(), 0);
            double gamma_abosolute = config.gamma;

            if(config.if_percentile){
                gamma_abosolute = getPercentile(u, graph.m_uwsum, config.gamma) / (double) graph.m_uwsum[u];
                if(gamma_abosolute < 1){
                    cout << "weight threshed:" << gamma_abosolute << "less than 1, replace with 1." << endl;
                    gamma_abosolute = 1;
                }
            }
            cout << "current node weight: " << graph.m_uwsum[u] << "; " << "gamma_abosolute: " << gamma_abosolute << endl;
            RoughBiPartialPush(u, config.alpha, config.epsilon, config.delta, gamma_abosolute, ppr, graph);
            stringstream ss;
            ss << ss_dir.str() << u << ".txt";
            fout.open(ss.str());
            fout.setf(ios::fixed,ios::floatfield);
            fout.precision(15);
            // cout << ss.str() << endl;
            if(!fout){
                cout<<"Fail to open the writed file"<<endl;
            }
            for(uint u_i=0; u_i < graph.getNu(); u_i++){
                if(ppr[u_i]>1e-8){
                    fout<<u_i<<" "<<ppr[u_i]<<endl;
                }
            }
            fout.close();
        }
    }
    
    cout << Timer::used(1)*1000/query_count << " milli-seconds per query" << endl;
    Timer::show();

    return SUCCESS;
}
