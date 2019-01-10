#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include <bitset>

//important for reading a directory
#include <dirent.h>
#include <sys/types.h>
#include <getopt.h> //einlesen der argumente

//neccessary for function strcmp
#include <stdio.h>
#include <string.h>

//own classes
#include "Matrix_new.hpp"

using namespace std;

//Global variables default values 
int NUMBER_ALL_TFs = 515; // number TFs that are consider in DBDs and sum_TFs fuer K652
int MAXIMUM_LENGTH_SEQUENCE = 1000; // number of the maximum lenght of the input sequnces

//functions
vector<string> readDirectory(const char *path);
double DomainInfo(string& TF, string& motif, map<string, pair<string, int>>& TF_DBDs, map<string, vector<double>>& DBDs, map<string, double>& sum_TFs);
map<string,pair<string, int>> readTF_DBDs(string& TF_to_DBD_file);
map<string, vector<double>> readDBDs(string& DBDs_file);
map<string, double> readSumTFs(string& sum_TFs_file);

//main script
int main(int argc, char *argv[]){

	//input
	cout << "input arguments: "<< endl;	
	int opt = 0;
	while ((opt = getopt(argc, argv, "a:")) != -1) {
       		switch (opt) {
		case 'a':
			NUMBER_ALL_TFs = stoi(optarg);
			cout << "-a number all TFs of evaluation cluster " << NUMBER_ALL_TFs << endl;
			break;
		default:
			throw invalid_argument("so nicht");
        	}
    	}
	if (optind + 6 > argc)  // there should be 3 more non-option arguments
		throw invalid_argument("missing motif_dir, seq_dir and/or output_file, sum_TFs.txt, DBDs.txt and/or TF_to_DBD.txt");

	string path_to_pwms = argv[optind++]; //read path to motifs
	cout <<"path to motif dir: " << path_to_pwms << endl;

	string path_to_pro_files = argv[optind++]; //read path to seq dir
	cout <<"path to fasta file dir: " << path_to_pro_files << endl;

	string output_file = argv[optind++];
	cout <<"output_file: " << output_file << endl;
	ofstream output;
	output.open(output_file);
	if(!output.is_open()){
		throw invalid_argument ("cannot open output file");
	}
	
	string sum_TFs_file = argv[optind++];
	map<string, double> sum_TFs = readSumTFs(sum_TFs_file);
	cout <<"sum_TFs: " << sum_TFs_file << endl;
	//check sum_TF
//	for (auto& i : sum_TFs)
//		cout << "( " << i.first << "  " << i.second << ") ";
//	cout << endl;  

	string DBDs_file = argv[optind++];
	cout <<"DBDs: " << DBDs_file << endl;
	map<string, vector<double>> DBDs = readDBDs(DBDs_file);
	//check DBDs
//	for (auto& i: DBDs){
//		cout << i.first << ": ";
//		for (auto& j : i.second)
//			cout << j << " ";
//		cout << endl;
//	}

	string TF_to_DBD_file = argv[optind++];
	cout <<"TF_to_DBD: " << TF_to_DBD_file << endl;
	map<string,pair<string, int>> TF_DBDs = readTF_DBDs(TF_to_DBD_file);
	//check TF_DBDs
//	for (auto& i : TF_DBDs){
//		cout << i.first << "( ";
//		pair<string, int> j  = i.second;
//		cout << j.first << " "  << j.second;
//		cout << ") ";
//	}
//	cout << endl;
//	cout << "---------------------------" << endl;


	//store all names of PWMs and seq_files in vectors
	vector<string> PWM_files = readDirectory(path_to_pwms.c_str());
//	cout << "number pwm files: " << PWM_files.size() << endl;
	vector<string> pro_files = readDirectory(path_to_pro_files.c_str());
//	cout << "number pro files: " << pro_files.size() << endl;

	//calculate for each motif the DomainInfo for 
	string TF = "";
	string motif = "";
	double  pvalue_domainInfo = 0.0;
	for(unsigned int i = 0; i < pro_files.size(); ++i){
	        //skip info.txt 
                if(pro_files[i] == "info.txt"){
                         continue;
                }
		TF = pro_files[i].substr(0 , pro_files[i].size() -3); // set actual TF
//		cout << "TF: " << TF << endl;
		//write parts of the output
		output << "TF:\t" + TF + "\n.\tpredicted Motif\tscore\n\n";
		vector<pair<double, string>> results;
		for(unsigned int j = 0; j < PWM_files.size(); ++j){ // iterate over all given motifs
			motif = PWM_files[j].substr(0 , PWM_files[j].size() -4); // set motif
			pvalue_domainInfo = DomainInfo(TF, motif, TF_DBDs, DBDs, sum_TFs); // calculate pvalue of DomainInfo
			results.push_back(make_pair(pvalue_domainInfo, motif));
		}
//		for(auto & i : results)	
//			cout << i.first << " " << i.second << endl;
		sort(results.begin(), results.end());
		int counter = 0;
		for(auto & i : results){
			counter++;
			output << counter << "\t" << i.second << "\t" << i.first << "\n";
		}
		output << endl;
	}
}
/*
*input: path to a directory
*output: vector that contains the names of all files of the directory
*WATCH OUT: it contains two files that are not neccessary for further calcultions: . and .. (aktuell rausgenommen)
*/
vector<string> readDirectory(const char *path){

	vector <string> result;
	//stores file names
  	dirent* de;

  	DIR* dp;

  	dp = opendir(path);
	if (dp == NULL){
		throw invalid_argument("sorry can not open directory");
 	}
	de = readdir(dp); 

    	while (de != NULL ){
		if(strcmp(de->d_name, ".") and strcmp(de->d_name, "..")){
      			result.push_back(de->d_name);
		}
      		de = readdir(dp);
      	}	
    	closedir(dp);
	return result;
}


// TF ist der gechipte TF von dem man das motif nicht kennt, und motif ist ein motif au sder datenbank
double DomainInfo(string& TF, string& motif, map<string, pair<string, int>>& TF_DBDs, map<string, vector<double>>& DBDs, map<string, double>& sum_TFs){
	string DBD = "";
	double sumTF = 0.0;
	double Smax = 0.0;
	double info = 0.0;
	
	DBD = TF_DBDs[TF].first;
	//cout << "DBD: " << DBD << endl;
	sumTF = sum_TFs[motif];
	//cout << "sumTFs: " << sumTF << endl;
	Smax = DBDs[DBD][TF_DBDs[motif].second];
	//cout << "Smax: " << Smax << endl;
	info = Smax / sumTF;
	//cout << "DomainInfo: " << info << endl;
	return (info);
}
/*
*sum TF file format
*TF		max Smax Wert über alle TF-Cluster
*AHR::ARNT      87.022467
*ALX1    	70.639915
*ALX3    	78.627349
*ALX4    	69.968588
*AR      	90.020155
*ARID3A  	85.31772
*/
map<string, double> readSumTFs(string& sum_TFs_file){

	map<string, double> sum_TFs;
	string TF = "";
	string s_max_sum = "";
	ifstream sum_TFs_file_(sum_TFs_file);
	if (sum_TFs_file_.is_open()){
		while(getline(sum_TFs_file_, TF, '\t')){ // read line until first /t
			getline(sum_TFs_file_, s_max_sum); // reads the rest of the line 
			sum_TFs[TF] = stod(s_max_sum);
		}
	}else{
		throw invalid_argument( "cannot open sumTFs.txt");
	}
	return sum_TFs;
}

/*
*DBDs file format
* .	TF1     TF2	TF3	...
*DBD1	1.2	2.3	3.5	...
*DBD2	....
*/
map<string, vector<double>> readDBDs(string& DBDs_file){

	map<string, vector<double>> DBDs;
	ifstream DBDs_file_(DBDs_file);//open file
	if (not DBDs_file_.is_open())
		throw invalid_argument("Cannot open DBDs.txt");
	
	string line = "";
	string key = "";
	getline(DBDs_file_, line);// omit header
	vector<double> helper;// stores the value of each key
	int count = 0;
	while(true){ 
		if (count == NUMBER_ALL_TFs){ // checks if we need to read a new line in the file
			getline(DBDs_file_, line);
			count = 0;
		}else{
			getline(DBDs_file_, line, '\t');
			count++;
		}
		if( line.empty()){ // if the line is empty we are at the end of the file and stop reading
			break;
		}
		if (line[0] == 'D'){
			if (key != ""){
				DBDs[key] = helper;
				helper.clear();
			}
			key = line;
		}else{
			helper.push_back(stod(line));
		}	
	}
	DBDs[key] = helper;
	return DBDs;
}

/*
*TF_to_DBD file format
*TF		DBD	spalte in DBDs
*AHR::ARNT      DBD_1   0   
*ALX1   	DBD_4   1   
*ALX3    	DBD_4   2   
*ALX4    	DBD_4   3   
*/
map<string, pair<string, int>> readTF_DBDs(string& TF_to_DBD_file){
	map<string,pair<string, int>> TF_DBDs;
	ifstream TF_to_DBD(TF_to_DBD_file);
	if (not TF_to_DBD.is_open())
		throw invalid_argument("Cannot open TF_to_DBD.txt");
	string key = "";
	string DBD = "";
	string help = "";

	while(getline(TF_to_DBD, key, '\t')){
		
		getline(TF_to_DBD, DBD, '\t');
		getline(TF_to_DBD, help);
		TF_DBDs[key] =  make_pair(DBD, stoi(help));
	}
	return TF_DBDs;
}
