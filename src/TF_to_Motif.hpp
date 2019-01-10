#ifndef TF_TO_MOTIF_HPP
#define TF_TO_MOTIF_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <algorithm> 
#include <omp.h>
#include <bitset>

//own classes
#include "pvalue_copy.hpp"
#include "Matrix_new.hpp"

//important for reading a directory
#include <dirent.h>
#include <sys/types.h>
#include <getopt.h> //einlesen der argumente

//neccessary for function strcmp
#include <stdio.h>
#include <string.h>

//for log
#include <math.h>

using namespace std;

//------------------------------
//Global variables default values 
int COMPLEMENT[] = {0,4,0,3,1,0,0,2,0,0,0,0,0,0,5};// considers also N
int POSITION[] = {0,1,0,2,4,0,0,3,0,0,0,0,0,0,5}; // considers also N 
double EPSILON = 0.001;
//------------------------------

class TF_to_Motif{

	public:
	//constructors
	TF_to_Motif();
	TF_to_Motif(bool domain_info_, bool position_info_, bool ranked_, double p_value_, int num_threads_, int number_all_TFs_, int maximum_sequence_length_,string frequence,  string path_to_position_info_, string path_to_ranked_seq_, string path_to_pwms_, string path_to_pro_files_, string output_file_, string sum_TFs_file_, string DBDs_file_, string TF_to_DBD_file_, string pvalues_DomainInfo);
	~TF_to_Motif(); // deconstrutor
	TF_to_Motif(const TF_to_Motif&); //copy constructor

	//functions
	vector<string> readDirectory(const char *path); // stores names of all files in a directory

	//read input files with data from cluster-family analysis ->precomputet to save runtime
	void readSumTFs(map<string, double>& sum_TFs, string& path); 
	void readDBDs(map<string, vector<double>>& DBDs, string& path);
	void readTF_DBDs(map<string,pair<string, int>>& TF_DBDs, string& path);
	void readPvaluesDomainInfo(map<string, double>& pvalues_domainInfo, string& path);
	vector<double> readFrequence(string& frequence);


	//functions reading input files into vectors	
//	void read_user_defined_pvalues(vector<double>& user_defined_pvalues, vector<string>& PWM_files);
	void read_position_info(string& path_to_position_info);
	void read_ranked_sequences(string& path_to_ranked_sequences);
	vector<Matrix<double>> PFMsToPWMs(vector<string>& PWM_files, vector<double>& freq);

	vector<double> determineFreq(vector<string> sequences);
	Matrix<double> determineTransitionMatrix(vector<string> sequences, Matrix<double>&  background);

	//determine domain, rank and position info
	bool DomainInfo(string& TF, string& motif);	
	double PositionInfo(int& pos);
	double RankInfo(string& TF, string& header);
	double RankInfo(string& TF, string& header, unordered_map<string, double>& actual_ranked_seq);

	//determine probability of a given sequence
	double probProSeq(Matrix<double>& PWM, string& line, double pvalue, vector<double>&  pvalues);

	//determine probability of a given kmer
	double probKmerComplement(Matrix<double>& PWM,string::iterator start, string::iterator end);
	double probKmer(Matrix<double>& PWM,string::iterator start, string::iterator end);
//	double maxKmer(Matrix<double>& PWM);

	void info(string& sum_TFs, string& DBD, string& TF_to_DBD, string& path_to_position_info, string& path_to_ranked_seq, string& output_file, string& pvalue_domainInfo);

	void write(string text);
	//getter

	bool getDomainInfo();
	bool getPositionInfo();
	bool getRankInfo();
	double getPvalue();
	int getNumThreads();
	int getNumAllTFs();
	int getMaxSeqLen();
	string getPathToPWMs();
	string getPathToSeq();

	vector<string> getPWMFiles();
	vector<string> getProFiles();
	map<string, double> getSumTFs();
	map<string, vector<double>> getDBDs();
	map<string, pair<string, int>> getTF_DBDs();
	vector<Matrix<double>> getPWMs();
	vector<double> getValuesPositionInfo();
	unordered_map<string, double> getSum();
	unordered_map<string, unordered_map<string,double>> getRankedSeq();
	map<string, double> getPvaluesDomainInfo();

	vector<Matrix<double>> PWMs;

//TODO: getter schreiben
	private:
	bool domain_info;
	bool position_info;
	bool ranked;
	double p_value;
	int num_threads; 
	int number_all_TFs; 
	int maximum_sequence_length;
	string path_to_pwms; //read path to motifs
	string path_to_pro_files; //read path to seq dir
	//TODO die fehlen bei copy constructor
	vector<string> PWM_files; //stores names of all pwm files
	vector<string> pro_files; //stores names of all pro_seq files
	map<string, double> sum_TFs;
	map<string, vector<double>> DBDs;
	map<string, pair<string, int>> TF_DBDs;
	ofstream output;
	vector<double> values_position_info;
	map<string, double> pvalues_domainInfo;
	unordered_map<string, unordered_map<string, double>> ranked_seq;
};
//------------------------------
//functions

//constructor
TF_to_Motif::TF_to_Motif()
:domain_info(false),
 position_info(false), 
 ranked(false),
 p_value(0.05), 
 num_threads(1),
 number_all_TFs(102),
 maximum_sequence_length(1000),
 path_to_pwms(""),
 path_to_pro_files("")
{
//	cout << "default constructor" << endl;
}

//constructor
TF_to_Motif::TF_to_Motif(bool domain_info_, bool position_info_, bool ranked_, double p_value_, int num_threads_, int number_all_TFs_, int maximum_sequence_length_,string freq_, string path_to_position_info_, string path_to_ranked_seq_, string path_to_pwms_, string path_to_pro_files_, string output_file_, string sum_TFs_file_, string DBDs_file_, string TF_to_DBD_file_, string pvalues_domainInfo_file_)
:domain_info(domain_info_),
 position_info(position_info_),
 ranked(ranked_),
 p_value(p_value_),
 num_threads(num_threads_),
 number_all_TFs(number_all_TFs_),
 maximum_sequence_length(maximum_sequence_length_),
 path_to_pwms(path_to_pwms_),
 path_to_pro_files(path_to_pro_files_)
{
//	cout << "own defined constructor" << endl;
	//open outputfile
	output.open(output_file_);
	if(!output.is_open()){
		throw invalid_argument ("cannot open output file");
	}
	//write running info
	info(sum_TFs_file_, DBDs_file_, TF_to_DBD_file_, path_to_position_info_, path_to_ranked_seq_, output_file_, pvalues_domainInfo_file_);
	
	//store names of the pwm_files and seq_files
	PWM_files = readDirectory(path_to_pwms_.c_str()); //stores names of all pwm files
	pro_files = readDirectory(path_to_pro_files_.c_str()); //stores names of all pro_seq files

	//werden direkt so gespeichert wie wirs für die pvalue berechnung benötigen	
	vector<double> freq = readFrequence(freq_);// freq wird nicht gespeichert da es nur einmal benötigt wird
	PWMs = PFMsToPWMs(PWM_files,freq);
	

	//TODO: rausfinden wie das mit freq und transition richtig ist
	//store informations that are necessary for domain information
	if (domain_info){
		//store sum_TFs in a map
		readSumTFs(sum_TFs, sum_TFs_file_);
		// store DBDs in a map, key is DBD and value is a list of Smax values for the TFs, the Smax value is the largest among all clusters of the DBD
		readDBDs(DBDs, DBDs_file_); // read DBDs 
		readTF_DBDs(TF_DBDs, TF_to_DBD_file_);
		readPvaluesDomainInfo(pvalues_domainInfo, pvalues_domainInfo_file_);
	}
	//stores position_info (already normalized) in values_position_info
	if(position_info){
		read_position_info(path_to_position_info_);
	}
	//stores ranked info in sum and ranked_seq
	if(ranked){
		read_ranked_sequences(path_to_ranked_seq_); 

	}
}

//deconstructor
TF_to_Motif::~TF_to_Motif()
{
	cout << "deconstructor" << endl;
}

//copy constructor 
TF_to_Motif::TF_to_Motif(const TF_to_Motif& t_)
:domain_info(t_.domain_info),
 position_info(t_.position_info),
 ranked(t_.ranked),
 p_value(t_.p_value),
 num_threads(t_.num_threads), 
 number_all_TFs(t_.number_all_TFs),
 maximum_sequence_length(t_.maximum_sequence_length),
 path_to_pwms(t_.path_to_pwms),
 path_to_pro_files(t_.path_to_pro_files)
{
//TODO: hier fehlen noch die ganzen komplizierten sachen  
}


/*
*input: path to a directory
*output: vector that contains the names of all files of the directory
*WATCH OUT: it contains two files that are not neccessary for further calcultions: . and .. (aktuell rausgenommen)
*/
vector<string> TF_to_Motif::readDirectory(const char *path){

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
void TF_to_Motif::readSumTFs(map<string, double>& sum_TFs, string& sum_TFs_file){

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
}

/*
*DBDs file format
* .	TF1     TF2	TF3	...
*DBD1	1.2	2.3	3.5	...
*DBD2	....
*/
void TF_to_Motif::readDBDs(map<string, vector<double>>& DBDs, string& DBDs_file){

	ifstream DBDs_file_(DBDs_file);//open file
	if (not DBDs_file_.is_open())
		throw invalid_argument("Cannot open DBDs.txt");
	
	string line = "";
	string key = "";
	getline(DBDs_file_, line);// omit header
	vector<double> helper;// stores the value of each key
	int count = 0;
	while(true){ 
		if (count == number_all_TFs){ // checks if we need to read a new line in the file
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
}

/*
*TF_to_DBD file format
*TF		DBD	spalte in DBDs
*AHR::ARNT      DBD_1   0   
*ALX1   	DBD_4   1   
*ALX3    	DBD_4   2   
*ALX4    	DBD_4   3   
*/
void TF_to_Motif::readTF_DBDs(map<string,pair<string, int>>& TF_DBDs, string& TF_to_DBD_file){

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
}

/*void TF_to_Motif::read_user_defined_pvalues(vector<double>& user_defined_pvalues, vector<string>& PWM_files){
	
	if (user_defined_pvalues_file == "")
		throw invalid_argument("no path to user defined pvalues file is set!");
	string key = "";
	string p = "";
	int index_name_PWM = 0;
	ifstream file(user_defined_pvalues_file);
	if (not file.is_open())
		throw invalid_argument("Cannot open user_defined_pvalues.txt");
	while(getline(file, key, '\t')){
		getline(file, p);
		//int counter = 0;
	//	for(auto i = PWM_files.begin(); i <= PWM_files.end(); i++){
	//		if ((*i) == key){
	//
	//		break;
	//		}
	//		counter++;
	//	}
		index_name_PWM = find(PWM_files.begin(), PWM_files.end(), key) - PWM_files.begin();
//		cout << "index name PWM: " << index_name_PWM << endl;
		user_defined_pvalues[index_name_PWM] = stod(p);
		//cout << "index name PWM: " << counter << endl;
		//user_defined_pvalues[counter] = stod(p);
	}
}
*/
vector<Matrix<double>> TF_to_Motif::PFMsToPWMs( vector<string>& PWM_files, vector<double>& freq){
	
	if(path_to_pwms == "")
		throw invalid_argument("no path to pwm files is set!");
	
	Matrix<double> PWM_matrix;
	vector<Matrix<double>> PWMs(PWM_files.size(), PWM_matrix);

	double rounder = 1/EPSILON;
	double value = 0;

	double freq_max = *(max_element(freq.begin(), freq.end()));
//	cout << "freq_max: " << freq_max << endl;

	for (int k = 0; k < PWM_files.size(); ++k){// for all PWMs
		ifstream PWM(path_to_pwms + PWM_files[k]);//open actual PWM file
		PWM >> PWM_matrix; //read PWM file as matrix
//		cout << PWM_matrix << endl;
		// determine PWM, round matrix and shift it in such a way that only values > 0 are included in the matrix -> important for pvalue calculation
		for (int i = 1; i<= PWM_matrix.ncol(); i++){
//			cout << "i:" << i << endl;
			for (int j = 1; j <= 4; j++){
//				cout << "j: " << j << endl;
				value = log10(PWM_matrix(j,i)/freq[j-1]) - log10(EPSILON/freq_max) + EPSILON;
				PWM_matrix(j,i) =  (round(value * rounder ) / rounder);
			}
		}
	//	cout << PWM_matrix << endl;
		PWMs[k] = PWM_matrix; // store matrix in map
	//	cout << PWMs[0] << endl;
	}
/*	cout << PWMs.size() << endl;
	for(int i = 0; i < PWMs.size(); ++i){
		cout << PWMs[i] << endl;
	}*/
	return PWMs;
}

void TF_to_Motif::read_position_info(string& path_to_position_info){

	cout << path_to_position_info << endl;	
	double sum_dis = 0.0;
	if(path_to_position_info == "")
		throw invalid_argument("path to position info file is not set");
		
	ifstream position_info_file(path_to_position_info);//open file
//	if (not position_info_file.is_open())
//		throw invalid_argument("Cannot open position_info.txt");
	string word = "";
	string word_2 = "";
	while (position_info_file >> word){
		position_info_file >> word_2; 
//		cout << "position: " << word << " value: " << word_2 << endl;
		values_position_info.push_back(stod(word_2));	
		sum_dis += stod(word_2);
//		cout << sum_dis << endl;
	}
	for(int i = 0; i <  values_position_info.size(); ++i)
		values_position_info[i] = values_position_info[i] / sum_dis;
//	cout << "sum_dis: " << sum_dis << endl;
//	for(auto& i : values_position_info)
//		cout << i << endl;
	return;
}

void TF_to_Motif::read_ranked_sequences(string& path_to_ranked_seq){

	vector<string> files = readDirectory(path_to_ranked_seq.c_str());
	for (auto& x : files){
		double sum = 0.0;		
		unordered_map<string, double> temp_map;
		ifstream ranked_seq_file(path_to_ranked_seq + x);//open file
		if (not ranked_seq_file.is_open())
			throw invalid_argument("Cannot open file with ranked sequences");
		string name = "";
		ranked_seq_file >> name;
		string word = "";
		string word_2 = "";
		while (ranked_seq_file >> word){
			ranked_seq_file >> word_2; 
			temp_map[word] = stod(word_2);	
			sum+=stod(word_2);
		}
		for (auto& i: temp_map)
			temp_map[i.first] = temp_map[i.first]/sum;

		ranked_seq[name] = temp_map;
	}
	return;
}

vector<double> TF_to_Motif::readFrequence(string& frequence){
	
	vector<double> freq;
	if(frequence == "")
		throw invalid_argument("path to frequence.txt is not set");
	
	ifstream file(frequence);
	if (not file.is_open())
		throw invalid_argument("Cannot open frequence.txt");
	string word = "";
	while (file >> word){
		freq.push_back(stod(word));
	}
	return freq;	


}


void TF_to_Motif::readPvaluesDomainInfo(map<string, double>& pvalues_domainInfo, string& path){
	
	string DBD = ""; 
	string value = ""; 

	ifstream pvalue_domainInfo_file(path);//open file

	if (not pvalue_domainInfo_file.is_open())
		throw invalid_argument("Cannot open domainInfo.txt");
	
	while (pvalue_domainInfo_file >> DBD){
		pvalue_domainInfo_file >> value;
		pvalues_domainInfo[DBD] = stod(value);
	}

}

vector<double> TF_to_Motif::determineFreq(vector<string> sequences){
	
	vector<double> freq(4, 0.0);
	int pos = 0;
	int counter = 0;
	for(auto& s : sequences){
		for(int j = 0; j < s.length(); ++j){
			pos = POSITION[(s[j] & 0x0f)]-1; //da es auf matrix angepasst ist und die bei 1 anfaengt zu zaehlen
			if (pos != 4){
				counter++;
				freq[pos]++;
			}
		}
	}
	for(int i = 0; i < freq.size(); i++){
		freq[i] = freq[i]/counter;
	}
	return freq;
}

Matrix<double> TF_to_Motif::determineTransitionMatrix(vector<string> sequences, Matrix<double>& background){

//	Matrix<double> background(4,4); //initialize matrix
	int pos = 0; 
	int col = 0;
	char previous_letter = '\0';
	for(auto& s: sequences){
		for(int j = 0; j < s.length(); ++j){
			if(previous_letter == '\0'){

				previous_letter = s[j];
				
				pos = POSITION[(previous_letter & 0x0f)]; //calculates in which row we need to store the event,
				if (pos == 5){
					previous_letter = '\0';
				}
				//previous letter row, actual letter col in transition matrix	
			}else{	
				col = POSITION[(s[j] & 0x0f)]; //calculates the col in which we need to store the event
//				cout << "col: " << col << endl;
				if (col != 5){
					background(col,pos)++;
					pos = col;
				}else{
					previous_letter = '\0';
				}

			}
		}
	}
//	cout << "background vor normierung: \n" << background << endl; 
	int count = 0;
	for(int k = 1; k < 5; ++k){
		//counts how many letters there are pro column
		for(int i = 1; i < 5; i++){
			count = count + background(i,k);

		} 
//		cout << "count: " << count << endl;
		//divides counts by this value
		for(int i = 1; i< 5; ++i){
			background(i,k) = background(i, k)/count;
		}
		count = 0;
	}
//	cout <<"backgrund: " << background << endl;
	return background;	
}

// TF ist der gechipte TF von dem man das motif nicht kennt, und motif ist ein motif au sder datenbank
bool TF_to_Motif::DomainInfo(string& TF, string& motif){

	string DBD = "";
	double sumTF = 0.0;
	double Smax = 0.0;
	double info = 0.0;

	if (domain_info == true){
	
		DBD = TF_DBDs[TF].first;
//		cout << "DBD: " << DBD << endl;
		sumTF = sum_TFs[motif];
//		cout << "sumTFs: " << sumTF << endl;
		Smax = DBDs[DBD][TF_DBDs[motif].second];
//		cout << "Smax: " << Smax << endl;
		info = Smax / sumTF;
	}	
	//TODO: hier noch pvalue bestimmen bzw auslesen
//	cout << "domainInfo: " << info << " pvalue threshold: " << pvalues_domainInfo[DBD] << endl;
	if (info >= pvalues_domainInfo[DBD]){
		return true;
	}else{
		return false;
	}
	//return (info);
}

double TF_to_Motif::PositionInfo(int& pos){

	return(values_position_info[pos]);
}

double TF_to_Motif::RankInfo(string& TF, string& header){


	unordered_map<string, double> actual_ranked_seq = ranked_seq[TF]; //contains rank of the seq
//	double actual_sum = sum[TF];// contains sum over all ranks for this TF
	if (actual_ranked_seq.find(header) == actual_ranked_seq.end()){
		cout << "header not found" << header << endl;
		throw invalid_argument("header not found in ranked sequence list");
	}
	return(actual_ranked_seq[header]);
	
}

double TF_to_Motif::RankInfo(string& TF, string& header, unordered_map<string, double>& actual_ranked_seq){

//	unordered_map<string, double> actual_ranked_seq = ranked_seq[TF]; //contains rank of the seq
//	double actual_sum = sum[TF];// contains sum over all ranks for this TF
	if (actual_ranked_seq.find(header) == actual_ranked_seq.end()){
		cout << "header not found" << header << endl;
		throw invalid_argument("header not found in ranked sequence list");
	}
	return(actual_ranked_seq[header]);
	
}

double TF_to_Motif::probProSeq(Matrix<double>& PWM, string& line, double pvalue, vector<double>& pvalues){

//	cout << "line: " << line << endl;
//	cout << "PWM:\n" << PWM << endl;
/*	cout << "pvalue: " << pvalue << endl;
	cout << "position info" << position_info << endl;
*/

	int length_motif = PWM.ncol();
        double prob = 0.0; // stores prob for each k_mer
        double prob_comp = 0.0; // stores prob for each k_mer of the complement
	double result = 0.0; //stores overall prob for the actual line
	//double result_comp = 0.0; //stores overall prob for the actual line complement 
	string::iterator start = line.begin() - length_motif + 1;
	string::iterator end = line.begin();

//	cout << "start: " << *(start) << " end: " << *(end) <<  endl; // muesste am anfang null sein

	double rounder = 1/EPSILON;
	int size_pvalues = pvalues.size() -1;

	int len_line = line.length();
	int half_len_line = len_line/2;
//	cout << "len line: " << len_line << endl;
	//store pos for vector that contains the pvalues
	int pos = 0; 
	int pos_comp = 0;
	int check = 0; // checks if there is a N in the sequence

	// for PositionInfo
	int pos_kmer = 0;
	int pos_best_kmer = 0;

        for(int i = 0; i < len_line ; ++i){ // da index  bei 0 anfaengt oder?			
		if (line[i] == 'N'){
			check = length_motif; 
		}
		if (line[i] != 'N' && check == 0 && start >= line.begin()){
			if (i  > (half_len_line)){
				pos_kmer =  abs((i - (half_len_line)) - length_motif/2 + 1);		
			}else{
				pos_kmer = (((half_len_line) - i) + (half_len_line)-i+length_motif-1) / 2;
			}
			prob = probKmer(PWM, start, end); //calculate probability of k_mer
//			cout << "prob kmer: "<< prob << endl;
			prob_comp = probKmerComplement(PWM,start, end); // determine prob for complement
//			cout << "prob comp kmer: " << prob_comp << endl;
			pos = round(size_pvalues -(rounder*prob)); //berechnet an welcher position im vector der pvalue für prob steht
			pos_comp = round(size_pvalues -(rounder*prob_comp)); //berechnet an welcher position im vector der pvalue für prob steht


			if (pvalues[pos] <= pvalue){ 
				if (result < (-log10(pvalues[pos]))){
					pos_best_kmer = pos_kmer;
					result =  -log10(pvalues[pos]);
		//			cout << "prob kmer: "<< prob << endl;
					
				}
			}
//			cout << "normal: " << result << endl;	

		//cout <<"vorwaerts:" <<  pos_best_kmer << endl;
			if(pvalues[pos_comp] <= pvalue){
				if (result < (-log10(pvalues[pos_comp]))){
					pos_best_kmer = pos_kmer;
					result =  -log10(pvalues[pos_comp]);
		//			cout << "prob comp kmer: " << prob_comp << endl;
				}
			}
//			cout << "comp: " << result << endl;	
		}
		//cout <<"rueckwaerts" <<  pos_best_kmer << endl;

		end++;// update last letter of k_mer
		start++;
		if (check > 0)
			check--;
	}
//	cout << pos_best_kmer << endl;
	if (position_info)
       	 	return (result*PositionInfo(pos_best_kmer));
	else{
		return(result);
	}
}


/*
*Input: PWM, iterators that point to the first letter of the k_mer and the last
*Output: probability of the actual k_mer
*/
double TF_to_Motif::probKmer(Matrix<double>& PWM,string::iterator start, string::iterator end){

	
        double prob = 0.0; //muss 0.0 sein bei plus
	int counter = 1;
	for(auto i = start; i != end+1; i++){ //iterators over k_mer
//		cout << *i;
		prob+= PWM(POSITION[((*i) & 0x0f)], counter); //determine which letter we do consider and look up the accodirng prob in the PWM
		counter ++;
	}
//	cout << endl;
//	cout << prob << endl;
        return prob;
}
//same as for probKMer, only different iterates in reverse order and determne complement for actual letter
double TF_to_Motif::probKmerComplement(Matrix<double>& PWM,string::iterator start, string::iterator end){

        double prob = 0.0; //muss 0.0 sein bei plus
	int counter = 1;
	for(auto i = end; i != start-1; i--){// iterator in reverse order
//		cout << COMPLEMENT[((*i) & 0x0f)];
		prob+= PWM(COMPLEMENT[((*i) & 0x0f)], counter); //determine which letter we do consider and look up the accodirng prob in the PWM
		counter ++;
	}
//	cout << endl;
//	cout << prob << endl;
        return prob;
}
void TF_to_Motif::info(string& sum_TFs, string& DBD, string& TF_to_DBD, string& path_to_position_info, string& path_to_ranked_seq, string& output_file, string& pvalue_domainInfo){

	output << "#\t-t number threads " << num_threads << "\n";
	output << "#\t-s maximal lenght sequences " << maximum_sequence_length << "\n";
	output << "#\t-a number all TFs of evaluation cluster " << number_all_TFs << "\n";
//	output << "#\t-p use pvalue method: " << p_value <<  endl;
//	output << "#\t-m user defined transition matrix: " << TRANSITION << "\n";
//	output << "#\t-u user defined pvalue list: " << user_defined_pvalues_file << "\n"; 

	if(position_info)
		output << "#\t-g position info: " << path_to_position_info << "\n";
	if(ranked)
		output << "#\t-r use ranked sequences: " << path_to_ranked_seq << "\n";
	if(!domain_info)
		output << "#\t-n no structure prior: " << "\n";
	if (domain_info){
		output <<"#\tsum_TFs: " << sum_TFs << '\n';
		output <<"#\tDBDs: " << DBD << '\n';
		output <<"#\tTF_to_DBD: " << TF_to_DBD << '\n';
		output <<"#\tpvalue_domainInfo: " <<pvalue_domainInfo << "\n";
	}
	output << "#\tpath to motif dir: " << path_to_pwms << "\n";
	output << "#\tpath to fasta file dir: " << path_to_pro_files << "\n";
	output << "#\toutput_file: " << output_file << "\n";

}
void TF_to_Motif::write(string text){
	output << text;
	return;
}

//getter

bool TF_to_Motif::getDomainInfo(){
	return(domain_info);
}

bool TF_to_Motif::getPositionInfo(){
	return(position_info);
}

bool TF_to_Motif::getRankInfo(){
	return(ranked);
}

double TF_to_Motif::getPvalue(){
	return(p_value);
}

int TF_to_Motif::getNumThreads(){
	return(num_threads);
}

int TF_to_Motif::getNumAllTFs(){
	return(number_all_TFs);
}

int TF_to_Motif::getMaxSeqLen(){
	return(maximum_sequence_length);
}

string TF_to_Motif::getPathToPWMs(){
	return(path_to_pwms);
}

string TF_to_Motif::getPathToSeq(){
	return(path_to_pro_files);
}

vector<string> TF_to_Motif::getPWMFiles(){
//	for(auto& i : PWM_files)	
//		cout << i << " ";

	return(PWM_files);
}

vector<string> TF_to_Motif::getProFiles(){
//	for(auto& i : pro_files)	
//		cout << i << " ";

	return(pro_files);
}

map<string, double> TF_to_Motif::getPvaluesDomainInfo(){
	return(pvalues_domainInfo);
}


map<string, double> TF_to_Motif::getSumTFs(){
	for(auto& i : sum_TFs)	
		cout << "(" << i.first << " , " << i.second << ") ";

	return(sum_TFs);
}


map<string, vector<double>> TF_to_Motif::getDBDs(){
	for(auto& i : DBDs){	
		cout  << i.first << " ( ";
		for(auto& j : i.second)
			cout << j << " ";
		cout << ") ";
	}
	return(DBDs);
}


map<string, pair<string, int>> TF_to_Motif::getTF_DBDs(){
	for(auto& i : TF_DBDs){	
		cout  << i.first << " ( ";
		pair<string, int> j  = i.second;
		cout << j.first << " "  << j.second;
		cout << ") ";
	}
	return(TF_DBDs);
}

vector<Matrix<double>> TF_to_Motif::getPWMs(){
//	for(auto& i : PWMs)	
//		cout << i << endl;

	return(PWMs);
}

vector<double> TF_to_Motif::getValuesPositionInfo(){
	for(auto& i : values_position_info)	
		cout << i << " ";

	return(values_position_info);
}

/*unordered_map<string, double> TF_to_Motif::getSum(){
	for(auto& i : sum){	
		cout  << i.first << " ( " << i.second << ") ";
	}
	return(sum);
}*/

unordered_map<string, unordered_map<string,double>> TF_to_Motif::getRankedSeq(){
/*	for(auto& i : ranked_seq){	
		cout  << i.first << " ( ";
		unordered_map<string, double> help = i.second;
		for(auto& j : help )
			cout << j.first << " " << j.second;
	}*/
	return(ranked_seq);
}
#endif/*TF_TO_MOTIF_HPP*/
