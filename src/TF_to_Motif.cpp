#include "TF_to_Motif.hpp"

int main(){


	bool domain_info = true;
	bool position_info = false;
	bool ranked = false;
	double p_value = 0.05;
	int num_threads = 1;
	int num_all_TFs = 515;
	int maximum_sequence_length = 1000;
//	string transition  = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/Phase2/transition_matrix.txt";
	string frequence = "../frequence.txt";
	string path_to_position_info = "../../gamma_distribution/gamma_distribution_lenght_1000_0.02.txt";
	string path_to_ranked_seq = "../../PASTAA/data/gene_list_K652_one_ChIP_seq_experiment_pro_TF_delete_N_seq_count_own_method/seq_len_100/"; 
//	string user_defined_pvalues_file = "pvalue_grid_e_0.5_markov_100_MASSIF.txt";
	string path_to_pwm = "../random_sequence_examples/all_PWMs_0.001/";
	string path_to_pro_files = "../random_sequence_examples/seq_all_PWMs_len_100/";
	string output_file = "test.txt";
	string sum_TFs_file = "../../mosta_src/JASPAR_only_homo_sapiens/clusterJASPAR/sum_TFs.txt";
	string DBD_file = "../../mosta_src/JASPAR_only_homo_sapiens/clusterJASPAR/DBDs.txt";
	string TF_to_DBD_file = "../../mosta_src/JASPAR_only_homo_sapiens/clusterJASPAR/TF_to_DBD.txt";
	string thresholds_domainInfo = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/mosta_src/RandomMotifs/pvalueThresholdDomainInfo.txt";



	//TODO: macht es nicht mehr sinn die vectoren als private zu speichern anstatt die pfade?!
	//test constructor	
	TF_to_Motif a;	
	TF_to_Motif b(domain_info, position_info, ranked, p_value, num_threads, num_all_TFs, maximum_sequence_length, frequence, path_to_position_info, path_to_ranked_seq, path_to_pwm, path_to_pro_files, output_file, sum_TFs_file, DBD_file, TF_to_DBD_file, thresholds_domainInfo);

	cout << "Domain info " << b.getDomainInfo() << endl;
/*	cout << "Position info " << b.getPositionInfo() << endl;
	cout << "Rank info " << b.getRankInfo() << endl;
	cout << "Pvalue " << b.getPvalue() << endl;
	cout << "Num threads " << b.getNumThreads() << endl;
	cout << "Num all TFs " << b.getNumAllTFs() << endl;
	cout << "max seq len " << b.getMaxSeqLen() << endl;
	cout << " path to pwms " << b.getPathToPWMs() << endl;
	cout <<"path to seq " << b.getPathToSeq() << endl;
	
	cout << "------------------------------------------" << endl;
	cout << "PWM files " << endl;
	vector<string> PWMs = b.getPWMFiles();
	cout << "------------------------------------------" << endl;
	cout << "pro files " << endl;
	vector<string> pro = b.getProFiles();
	cout << "------------------------------------------" << endl;
	cout << "sum TFs " << endl;
	map<string, double> sumTFs = b.getSumTFs();
	cout << "------------------------------------------" << endl;
	cout << "DBD " << endl;
	map<string, vector<double>> DBD = b.getDBDs();
	cout << "------------------------------------------" << endl;
	cout << "TF_to_DBD " << endl;
	map<string, pair<string, int>> TF_DBD =  b.getTF_DBDs();
	cout << "------------------------------------------" << endl;
	cout << "PWM " << endl;*/
//	vector<Matrix<double>> PWM = b.getPWMs();
//	for (auto& i: PWM){
//		cout << i << endl;

//	}
/*	cout << "------------------------------------------" << endl;
	cout << "Values position info " << endl;
	vector<double> values = b.getValuesPositionInfo();
	cout << "------------------------------------------" << endl;
	cout << "sum " << endl;
	unordered_map<string, double> sum = b.getSum();
	cout << "------------------------------------------" << endl;
	cout << "ranked_seq " << endl;
	unordered_map<string, unordered_map<string, double>> r_seq = b.getRankedSeq();
*/

	//TODO rankde sse q aund pos info is nicht getestet und proprKmer und das ganze
	string ARNT = "ARNT";
	string ALX3 = "ALX3";
	int w = 2;
	b.write("test test!" + to_string(w));

	cout << "Domain Info TF ARNT und motif ARNT " << b.DomainInfo(ARNT, ARNT ) << endl;
	cout << "Domain Info TF ARNT und motif ALX3 " << b.DomainInfo(ARNT, ALX3 ) << endl;

	//string line = "ACGTACGTAGGGATTATTA";
	string line = "ACTAATTAACATGGGATTATACTAATTACGGT";
	cout << line << endl;
//	Matrix<double> motif = b.PWMs[0];
//	cout << motif << endl;
	pvalue  pvalue_obj(10, EPSILON); // EPSILON entspricht accuracy !!!
	vector<string> sequences; 
	sequences.push_back(line); 
	vector<double> freq = b.determineFreq(sequences);
	cout << "freq: " << endl;
	for(auto& i: freq)
		cout << i << " ";
	Matrix<double> transition_matrix(4,4);
	b.determineTransitionMatrix(sequences, transition_matrix);
/*	freq.push_back(0.25);
	freq.push_back(0.25);
	freq.push_back(0.25);
	freq.push_back(0.25);
	Matrix<double> transition_matrix(4,4);
	for (int i = 1; i <=4; ++i){
		for (int j = 1; j <= 4; ++j){
			transition_matrix(i,j) = 0.25;
		}
	}*/
	cout << "transition_matrix " << transition_matrix;
	vector<double> pvalues = pvalue_obj.calculatePvalues(b.PWMs[2], freq, transition_matrix);

/*	cout << "pvalues: \n";
	for( auto& i : pvalues)
		cout << i << " ";
*/

	cout << "prob pro seq: " << b.probProSeq(b.PWMs[2], line , 0.05,  pvalues) << endl;
	map<string, double> thresholds = b.getPvaluesDomainInfo();
	for(auto& i : thresholds)
		cout << i.first << " " << i.second << endl;

		
}


