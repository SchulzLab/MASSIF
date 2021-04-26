//Helge Roider, Max Planck Institute for Molecular Genetics
//14.10.2008
#include<stdio.h>
#include<cmath>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string>
#include<iomanip>
#include<algorithm>
#include<vector>

using namespace std;


//Required file formats:
//1. AFFINITY DATA
//GENEID  CpG  TF1Affy  TF2Affy  TF3Affy ...
//...
//2. USER LIST
//ENSGID1
//ENSGID2


void q_sort(double numbers[], int indices[], int left, int right){
  double pivot, l_hold, r_hold;
  int pivot_i;

  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  pivot_i = indices[left];

  while (left < right){
    while ((numbers[right] >= pivot) && (left < right)){
      right--;
    }
    if (left != right){
      numbers[left] = numbers[right];
      indices[left] = indices[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right)){
      left++;
    }
    if (left != right){
      numbers[right] = numbers[left];
      indices[right] = indices[left];
      right--;
    }
  }
  numbers[left] = pivot;
  indices[left] = pivot_i;

  pivot = left;
  pivot_i = left;
  left = l_hold;
  right = r_hold;
  if(left < pivot){
    q_sort(numbers, indices, left, pivot-1);
  }
  if(right > pivot){
    q_sort(numbers, indices, pivot+1, right);
  }
}

void quickSort(double numbers[], int indices[], int array_size){
  q_sort(numbers, indices, 0, array_size - 1);
}



int main(int argc, char *argv[]){
  if(argc < 2){
    cout << "Input file missing!\n";
    exit(1);
  }

  int numOfGenes = 0;
  int numOfFactors = 0;

  //--------------------------------------------------------------------
  //DETERMINING NUMBER OF GENES and FACTORS
  //--------------------------------------------------------------------

  string row = ""; 
  string delimiters = " \t"; //word seperators in each line
  int start = 0, end = 0;
  
  int rows = -1;
  int col = 0;
  int header = 1;

  ifstream affinities(argv[1]);
  if(!affinities){
    cout << "file not opened\n";
    exit(1);
  }
  while(!affinities.eof()){
    getline(affinities,row);
//	if (row.size() == 0){
//	cerr << "rows: " << rows<< endl;
//	}
    if(row.substr(0,1) == ""){continue;}
    rows++;
    if(header == 1){
      start = row.find_first_not_of(delimiters);
      while(start != string::npos){
	end = row.find_first_of(delimiters,start + 1);
	if(end == string::npos){
	  end = row.length();
	}
	col++;
	start = row.find_first_not_of(delimiters,end+1);
      }
      header = 0;
    }
  }
  affinities.close();

  numOfGenes = rows;
  numOfFactors = col - 1;


  //----------------------------------------------------------------------------------------
  //READING IN AFFINITY DATA
  //----------------------------------------------------------------------------------------

  //string genes[numOfGenes]; 
// changed such that this string array is stored dynamically on the heap needs to be deleted!!
  string * genes = new string[numOfGenes];
  string factors[numOfFactors];
//	cerr << "genes: " << genes << " factors: " << factors << endl;
  double ** psi;
  psi = new double * [numOfGenes];
  for(int j = 0; j < numOfGenes; j++){
    psi[j] = new double[numOfFactors];
  }

//new
//	vector<vector<double>> psi(numOfGenes, vector<double>(numOfFactors, 0));

  ifstream raffinities(argv[1]);
  header = 1;
  rows = -1;
  while(!raffinities.eof()){
    getline(raffinities,row);
    if(row.substr(0,1) == ""){continue;}
    col = -1;
    start = row.find_first_not_of(delimiters);
    while(start != string::npos){
      col++;
      end = row.find_first_of(delimiters,start + 1);
      if(end == string::npos){
	end = row.length();
      }
      if(header == 0){
	if(col == 0){
	  genes[rows] = row.substr(start,end-start);
	}
        else{
          psi[rows][col - 1] = strtod( row.substr(start,end-start).c_str(),NULL );
        }
      }
      else{
	if(col > 0){
	  factors[col - 1] = row.substr(start,end-start);
	}
      }
      start = row.find_first_not_of(delimiters,end+1);
    }
    if(header == 1){
      header = 0;
    }
    rows++;
  }
  raffinities.close();

//	return 0;


  //----------------------------------------------------------------------------------------
  //READING IN USER LIST
  //----------------------------------------------------------------------------------------

  ifstream userlist(argv[2]);
  if(!userlist){
    cout << "user list not opened\n";
    exit(1);
  }
  rows = 0;
  while(!userlist.eof()){
    getline(userlist,row);
    if(row.substr(0,1) == ""){continue;}
    rows++;
  }
  userlist.close();

  int numOfUsergenes = rows;
//  string usergenes[numOfUsergenes];
// changed such that this string array is stored dynamically on the heap needs to be deleted!!
  string * usergenes = new string[numOfGenes];
// changed such that this double array is stored dynamically on the heap needs to be deleted!!
  double * usergeneranks = new double[numOfUsergenes];
  //double usergeneranks[numOfUsergenes];
  //----------------------------------------------------------------------------------------

  ifstream ruserlist(argv[2]);
  rows = -1;
  while(!ruserlist.eof()){
    getline(ruserlist,row);
    if(row.substr(0,1) == ""){continue;}
    rows++;
    col = -1;
    start = row.find_first_not_of(delimiters);
    while(start != string::npos){
      end = row.find_first_of(delimiters,start + 1);
      if(end == string::npos){
        end = row.length();
      }
      col++;
      if(col == 0){
	usergenes[rows] = row.substr(start,end-start);
      }
      if(col == 1){
	usergeneranks[rows] = strtod( row.substr(start,end-start).c_str(),NULL );
      }
      start = row.find_first_not_of(delimiters,end+1);
    }
  }
  ruserlist.close();

  int numOfMatchedUsergenes = 0;
//  double tissuerankvalues[numOfGenes];  
// changed such that this double array is stored dynamically on the heap needs to be deleted!!
  double * tissuerankvalues = new double[numOfGenes];
  for(int g = 0; g < numOfGenes; g++){
    tissuerankvalues[g] = 99999;
    for(int u = 0; u < numOfUsergenes; u++){
      if(usergenes[u].compare(genes[g]) == 0){
	tissuerankvalues[g] = usergeneranks[u];
	numOfMatchedUsergenes++;
	break;
      }
    }
  }

  numOfUsergenes = numOfMatchedUsergenes;
  cerr << "\n";
  cerr << numOfUsergenes <<" genes in input list were matched\n";
  cerr << numOfGenes <<" genes with annotated affinity\n";
  cerr << numOfFactors <<" factors (PFMs)\n__________________\n\n";


  //---------------------------------------------------------------------
  //SETTING PREDEFINED CUTOFFS on AFFINITIES AND TISSUERANKS
  //---------------------------------------------------------------------
  //////////
  //TISSUE//
  //////////
  int tissuecutoffs [] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,300,400,500,600,700,800,900,1000};
  int numOfTissueCutoffs = sizeof(tissuecutoffs)/sizeof(int);
  bool tissueflag[numOfTissueCutoffs];
  for(int tc = 0; tc < numOfTissueCutoffs; tc++){
    tissueflag[tc] = false; //no new genes fall into this cutoff range
  }

  int tissueindices[numOfGenes];
  for(int g = 0; g < numOfGenes; g++){
    tissueindices[g] = g;
  }
  quickSort(tissuerankvalues, tissueindices, numOfGenes);

  int genes_in_tissue[numOfTissueCutoffs];
  double topT[numOfTissueCutoffs];
  for(int tc = 0; tc < numOfTissueCutoffs; tc++){
    if(tissuecutoffs[tc] <= numOfUsergenes){
      topT[tc] = tissuerankvalues[tissuecutoffs[tc] - 1];
      genes_in_tissue[tc] = tissuecutoffs[tc];

      int n = 1;
      while((tissuecutoffs[tc] - 1 + n <= numOfUsergenes)&&(tissuerankvalues[tissuecutoffs[tc] - 1 + n] == topT[tc])){//tied ranks
	genes_in_tissue[tc]++;
	n++;
      }
    }
  }
  //////////
  //FACTORS//
  ///////////
  int affycutoffs [] = {25,50,75,100,125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000};
  int numOfAffyCutoffs = sizeof(affycutoffs)/sizeof(int);

  //LOOP OVER FACTORS//
  for(int f = 0; f < numOfFactors; f++){
//        cerr << factors[f] << "\n";
//	cerr << "number genes " << numOfGenes << endl;
   // double inverseaffy1[numOfGenes];
// changed such that this double array is stored dynamically on the heap needs to be deleted!!
  double * inverseaffy1 = new double[numOfGenes];
  double * inverseaffy2 = new double[numOfGenes];
   // double inverseaffy2[numOfGenes];
    for(int g = 0; g < numOfGenes; g++){
      inverseaffy1[g] = 1/(psi[g][f]+1); //reverse ranking
      inverseaffy2[g] = inverseaffy1[g];
    }
    //    quickSort(inverseaffy1, dummyindices, numOfGenes);
    sort(inverseaffy1, inverseaffy1 + numOfGenes);

    int all_targets[numOfAffyCutoffs];
    for(int ac = 0; ac < numOfAffyCutoffs; ac++){
      if(affycutoffs[ac] < numOfGenes){
	all_targets[ac] = affycutoffs[ac] - 1;
	int n = 0;
	while((affycutoffs[ac] - 1 + n < numOfGenes)&&(inverseaffy1[affycutoffs[ac] - 1] == inverseaffy1[affycutoffs[ac] - 1 + n])){//tied ranks; run at least once
	  all_targets[ac]++;
	  n++;
	}
      }
    }
    int targets_within_tissue[numOfTissueCutoffs][numOfAffyCutoffs];
    for(int ac = 0; ac < numOfAffyCutoffs; ac++){
      for(int tc = 0; tc < numOfTissueCutoffs; tc++){
	targets_within_tissue[tc][ac] = 0;
      }
    }

    for(int u = 0; u < numOfUsergenes; u++){
      for(int ac = 0; ac < numOfAffyCutoffs; ac++){
	if(affycutoffs[ac] < numOfGenes){
	  if(inverseaffy2[tissueindices[u]] <= inverseaffy1[affycutoffs[ac] - 1]){
	    for(int tc = 0; tc < numOfTissueCutoffs; tc++){
	      if(tissuecutoffs[tc] <= numOfUsergenes){
		if(tissuerankvalues[u] <= topT[tc]){
		  targets_within_tissue[tc][ac]++;
		}
	      }
	    }
	  }
	}
      }
    }
    
    //------------------------------------------------------------------------------
    //HYPERGEOMETRIC TESTS
    //------------------------------------------------------------------------------
    
    //precalculated logs
// changed such that this double array is stored dynamically on the heap needs to be deleted!!
  double * storedlog = new double[numOfGenes + 1];
    storedlog[0] = 1;
    for(int n = 1; n < numOfGenes + 1; n++){
      storedlog[n] = log(n);
    }
    int optimaltargetsintissue = -1;
    int optimalgenesintissue = -1;
    int optimalalltargets = -1;
    double mostsignificant = 9999;
  
    for(int tc = 0; tc < numOfTissueCutoffs; tc++){
      if(tissuecutoffs[tc] > numOfUsergenes){
	break;
      }

      int previous_targetsintissue = 0;

      for(int ac = 0; ac < numOfAffyCutoffs; ac++){
	if(affycutoffs[ac] > numOfGenes){
	  break;
	}
	
	int alltargets = all_targets[ac];
	int genesintissue = genes_in_tissue[tc];
	int targetsintissue = targets_within_tissue[tc][ac];

	if(targetsintissue <= previous_targetsintissue){//no better p-value can emerge -> skip 
	  continue;
	}
	else{
	  previous_targetsintissue = targetsintissue;
	}

	double expected = alltargets * static_cast<double>(genesintissue) / static_cast<double>(numOfGenes);
	if((targetsintissue <= expected)||(expected == 0)){//if less targets than expected skip immediately
	  continue;
	}
	  	      

	long double probability = 0L;
	long double r = 0L;


	//-----------------------------------------------------------


	//Probability of NO hit in tissue: CAN HAPPEN ONLY IF genesintissue + alltargets[fc] <= numOfGenes!!!!
	if(alltargets <= numOfGenes - genesintissue){
	  for(int n = 0; n < alltargets; n++){
	    r += storedlog[numOfGenes - genesintissue - n] - storedlog[numOfGenes - n];
	  }
	}
	probability = 1 - exp(r);
	
	//Probability of 1 to targetsintissue hits in tissue
	for(int m = 1; m < targetsintissue; m++){ //use < to get probability 1 - sum_0_x-1 instead of 1 - sum_0_x
	  r = 0L;
	  for(int n = 0; n < m; n++){
	    r += storedlog[genesintissue - n] + storedlog[alltargets - n] - storedlog[m - n] - storedlog[numOfGenes - alltargets + m - n];  
	  }
	  for(int n = 0; n < alltargets - m; n++){
	    r += storedlog[numOfGenes - genesintissue - n] - storedlog[numOfGenes - n];
	  }
	  
	  probability -= exp(r);
	}


	//WARNING FLOATING POINT PRECISION INSUFFICIENT; estimate tail probability instead
	if(probability < 1e-9){
	  probability = 0L;
	  
	  int upperlimit = targetsintissue + 20;
	  if(upperlimit > alltargets){upperlimit = alltargets;}
	  if(upperlimit > genesintissue){upperlimit = genesintissue;}
	  int lowerlimit = targetsintissue; //lower limit includes x --> sum_X_X+20
	  
	  //sum over tail probabilities 
	  for(int m = lowerlimit; m <= upperlimit; m++){
	    r = 0L;
	    for(int n = 0; n < m; n++){
	      r += storedlog[genesintissue - n] + storedlog[alltargets - n] - storedlog[m - n] - storedlog[numOfGenes - alltargets + m - n];
	    }
	    for(int n = 0; n < alltargets - m; n++){
	      r += storedlog[numOfGenes - genesintissue - n] - storedlog[numOfGenes - n];
	    }
	    
	    probability += exp(r);
	  }
	}
	
	//Store most significant p-values
	if(mostsignificant > probability){
	  optimalgenesintissue = genesintissue;
	  optimaltargetsintissue = targetsintissue;
	  optimalalltargets = alltargets;	  
	  mostsignificant = probability;
	}

      }//END AFFY CUTOFFS
      
    }//END TISSUE CUTOFFS
    

    cout.precision(4);
    string FACTOR = factors[f];
    cout << FACTOR << "\t" << scientific << mostsignificant << "\t" << optimaltargetsintissue << "\t" << optimalgenesintissue << "\t" << optimalalltargets << "\t" << numOfGenes << "\t" << numOfUsergenes << "\n";
    
  delete [] inverseaffy1;
  delete [] inverseaffy2;
  delete [] storedlog;
  }//END FACTORS
  

  //--------------------------------------------------------------
  //CLEAR DYNAMIC VARIABLES
  //--------------------------------------------------------------

	//delete heap stored elements
	delete [] tissuerankvalues;
  	delete [] genes;
  	delete [] usergenes;
	delete [] usergeneranks;
  for(int j = 0; j < numOfGenes; j++){
    delete [] psi[j];
  }
  delete [] psi;
}
