
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#define samplesneeded 7040
#define slider 80
#define total_no_of_frame 85
#define no_of_obs_seq 85
#define codebook_size 32
#define training_utterances 15
#define avg_itration 3
#define p 12
#define N 5

using namespace std;

vector<long double> sample_vectr;     //Keep preprocessed samples
vector<long double> word_vectr;		  //Keep 7040 samples near to peak sample

long double codebook[32][12];
long double tokhura_wt[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
long double raise_sine[12]={2.552914271,4,5.242640687,6.196152423,6.795554958,7,6.795554958,6.196152423,5.242640687,4,2.552914271,1};   //TO STORE RAISED SINE WINDOW

long double A[N][N];
long double B[N][codebook_size];
long double PI[N];
long double A_mean[N][N];
long double B_mean[N][codebook_size];
long double PI_mean[N];


long double alpha[no_of_obs_seq][N];
long double beta[no_of_obs_seq][N];
long double delta[no_of_obs_seq][N];
long double psi[no_of_obs_seq][N];
long double gamma[no_of_obs_seq][N];
long double zyi[no_of_obs_seq][N][N];


long double pstar[training_utterances];       //probability
int qstar[training_utterances][no_of_obs_seq];  //optimal state sequence


string digit_arr[10] ={"zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"};



//Function to find dc_shift  value.
long double find_dc_shift()
{
    long double sum_of_samples_value =0.0; 
	long long no_of_samples =0;
	long double sample_value;
	long double dc_shift_value;
	ifstream fsilence;
	//open the silence_recorded.txt file.
	
	fsilence.open("silence_recorded.txt");
	if(!fsilence)
	{
		cout<<"Silence file is not there\n";
		return 0;
	}
	/* Read this silence file till end. sum every sample value. */
	while(!fsilence.eof())
	{
		fsilence>>sample_value;
		sum_of_samples_value +=sample_value;
		no_of_samples++;
	}
	fsilence.close();
	//find mean of the all sample value.
	dc_shift_value =sum_of_samples_value/no_of_samples;
	return dc_shift_value;

} //end of function find_dc_shift()

//Function to apply hamming window
void apply_hamming_window(long double sample_vec[])
{
	int itr;
	int cont =0;
	for(itr =0; itr<320; itr++)
	{
		sample_vec[itr]*=  (0.54 - 0.46*cos((2*3.14*cont)/319));
		cont++;

	} // end of for-loop.
}

//Applying dc_shift and normalization and returning the sample with peak value.
int preprocessing(string filename, double dc_shift)
{
	ifstream frawsample;
	frawsample.open(filename);
	if(!frawsample)
	{
		cout<<"Not able to open "<<filename.substr(10)<<endl;
	}
	//performing dc-shift.
	
	while(!frawsample.eof())
	{
		long double sample_val;
		frawsample>>sample_val;
		sample_val -=dc_shift;
		sample_vectr.push_back(sample_val);
	}

	//Apply Normalization
	int sample_value=0, sample_value_signed, peak_indx;
	long double peak_sample = LONG_MIN;
	long double normalization_factr;
	for(int i =1; i<sample_vectr.size(); i++)
	{
		if(peak_sample < abs(sample_vectr[i]) )
		{
			peak_sample =abs(sample_vectr[i]);
			peak_indx =i;
		}	
	}
	/*cout<<"*..............................................................................*\n";
	cout<<peak_sample<<endl;
	cout<<"*...............................................................................*\n";  */
	normalization_factr =10000/peak_sample;
	for(int i =1; i<sample_vectr.size(); i++)
	{
		//Normalizing each sample value by multiplying to the normalization_factr.
		sample_vectr[i] *= normalization_factr;
	}
	
	return peak_indx;
}

//Function to calculate R[0] to R[12]
void find_r0_to_r12(vector<long double> &R, long double frame_arr[])
{
	int lag;
				
	//Iterating through for loop for lag =0 to 12 and calculating R[0] to R[12] for each frame.
	for(lag =0; lag<=12; lag++)
	{
		int itr1;
		//Iterating for 320 samples of one frame.
		for(itr1 =0; itr1 <=319-lag; itr1++)
		{
			R[lag]+= frame_arr[itr1] * frame_arr[itr1 +lag];
		}
						
	}//end of for loop
	
}

//function to find A[1] to A[12]	
void find_a1_to_a12(vector<long double> &arr, vector<long double> R)
{
	long double alpha1[13][13]= {0};
	long double E[13] ={0}, k;
	//memset(alpha, 0.0, 13);
	//memset(E, 0.0, 13);
	E[0] =R[0];
	for(int i =1; i<=12; i++)
	{
		long double val =0.0;
		for(int j =1; j<=i-1; j++)
		{
			val += alpha1[i-1][j] *R[i-j];
		}

		k = (R[i] -val)/E[i-1];
		alpha1[i][i] =k;

		for(int j=1; j<=i-1; j++)
		{
			alpha1[i][j] =alpha1[i-1][j] - k*alpha1[i-1][i-j];   
		}

		E[i] = (1-k*k)*E[i-1];

	} //for-loop end

	//update A[1] to A[12]
	for(int i=1; i<=12; i++)
	{
		arr[i] = alpha1[12][i];
	}
}

/*Function to calculate C[1] to C[12]. As, discussed with Sir in the classRoom
  since, without negating A[] i.e..  without changing  sign of A[] the cepstral 
  coefficients value is coming in correct range and vowel detection is also more
  accurate. so not negating the A[] before finding cepstral coefficients. */
void find_c1_to_c12(vector<long double> &cep, vector<long double> arr)
{
	long double sum_value;
	for (int m=1; m<=12; m++)
	{
		cep[m] =0;
		sum_value = 0;
		for (int k=1; k<=m-1; k++)
		{
			sum_value += (double(k)/double(m)) * cep[k] * arr[m-k];  
		}
		cep[m] = arr[m] + sum_value;
		//cout<< cep[m];
	}
	//cout<<"check "<<endl;
	for (int m=1; m<=12; m++)
	{
		
		cep[m] *=raise_sine[m-1];
		//cout<< cep[m]<<" ";
	}
	//cout<<endl;
}

//Extracing the samples of a particular frame no. and putting into frame_arr to find ci's.
void extract_frame(int frame_no, long double frame_arr[])
{
	int sample_count =320;
	
	if(frame_no ==0)
	{
		for(int i =0; i<320; i++)
			frame_arr[i] =word_vectr[i];
	}
	else
	{
		int startindx =80 *frame_no;
		int endindx =80 *frame_no +sample_count;
		int j=0;
		for(int i =startindx; i<endindx; i++)
			frame_arr[j++] =word_vectr[i];
	}
}


void initialize_avg_A_B()
{

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			A_mean[i][j] =0;
	}

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<codebook_size; j++)
			B_mean[i][j] =0;
	}

	for(int i=0; i<N; i++)
	{
		PI_mean[i] =0;
	}
}

/**********************************************************************************************/

//Get observation  sequence from ci's.
void get_obs_seq_frm_featurevectr(int obs_seq_arr[], int frm_no, vector<long double> C)
{
	//cout<<"\n   Inside obs_seq function"<<frm_no<<"  \n";
	long double diff, tokhura_dist, min_dist =LONG_MAX;
	long double dist_from_codebook[codebook_size];
	int min_dist_codebook_indx;
	
	//cout<<"*****************  Observation sequence  *******************\n";
	for(int i=0; i<codebook_size; i++)
	{
		tokhura_dist =0;
		for(int j=0; j<p; j++)
		{
			//cout<<C[j+1]<<" ";
			diff = C[j+1] -codebook[i][j];
			tokhura_dist +=tokhura_wt[j] *diff*diff;
		}
		//cout<<endl;
		dist_from_codebook[i] =tokhura_dist;
		//cout<<i<<" index of codebook "<<dist_from_codebook[i]<<endl;
		if(dist_from_codebook[i] < min_dist)
		{
			min_dist =dist_from_codebook[i];
			min_dist_codebook_indx =i;
		}
	}
	obs_seq_arr[frm_no] =min_dist_codebook_indx+1;

	//system("pause");
}


//Read initial Model
void read_initial_model(int avg_num, string digit)
{
	ifstream freadA;
	ifstream freadB;
	ifstream freadPI;

	if(avg_num ==0)
	{
		freadA.open("A_MATRIX.txt");
		freadB.open("B_MATRIX.txt");
		freadPI.open("PI_MATRIX.txt");
		
	}

	else
	{
		string strA = "A_prev_model_for_"  + digit + ".txt";
		string strB = "B_prev_model_for_"  + digit + ".txt";
		string strPI = "PI_prev_model_for_"  + digit + ".txt";
		freadA.open(strA);
		freadB.open(strB);
		freadPI.open(strPI);
	}

	while(!freadA.eof())
	{
		for(int i = 0; i< N; i++)
		{
			for(int j = 0; j< N; j++)
				freadA >>A[i][j];
		}
	}

	while(!freadB.eof())
	{
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < codebook_size; j++)
				freadB >>B[i][j];
		}
	}


	while(!freadPI.eof())
	{
		for(int i = 0; i < N; i++)
			freadPI>> PI[i];
	}

	freadA.close();
	freadB.close();
	freadPI.close();
}


//Averaging the model for all utterances of a digit.
void add_to_average()
{
	//Adding A
	for(int i =0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			A_mean[i][j] += A[i][j];
	}

	//Adding B
	for(int i =0; i<N; i++)
	{
		for(int j=0; j<codebook_size; j++)
			B_mean[i][j] += B[i][j];
	}

	//Adding PI
	for(int i =0; i<N; i++)
	{
		PI_mean[i]+= PI[i];
	}

}




/********************* writing the final Model   into text file               ********************/
void write_final_model(string digit)
{
	ofstream fwriteA;
	ofstream fwriteB;
	ofstream fwritePI;

	string strA, strB, strPI;
	strA = "A_prev_model_for_" + digit + ".txt";
	strB = "B_prev_model_for_" + digit + ".txt";
	strPI = "PI_prev_model_for_" + digit + ".txt";

	fwriteA.open(strA);
	fwriteB.open(strB);
	fwritePI.open(strPI);
	
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			fwriteA <<A_mean[i][j]<<" ";
		fwriteA <<endl;
	}
	cout<<"\nWriting final model "<<endl;
	long double sum_row_b =0;
	for(int i=0; i<N; i++)
	{
		sum_row_b =0;
		for(int j=0; j<codebook_size; j++)
		{
			cout<<B_mean[i][j]<<"  ";
			sum_row_b += B_mean[i][j];
			fwriteB <<B_mean[i][j]<<" ";
		}
		fwriteB<<endl;
		cout<<"\nAfter doing sum of a row of B "<<sum_row_b<<endl;
	}
	

	
	for(int i=0; i<N; i++)
		fwritePI <<PI_mean[i]<<" ";
	

	fwriteA.close();
	fwriteB.close();
	fwritePI.close();
}




/******************************* Solutions of HMM Model   ***********************************/

long double forward_procedure(int obs_seq_arr[], int turn=9999999 )
{
	long double prob_obs_seq_given_model=0;
	/*********************    Forward Procedure    *******************/

	//Initialization
	//cout<<"....................   alpha  ......................"<<turn<<endl;
	for(int i=0; i<N; i++)
	{
		//cout<<PI[i]<<" "<<B[i][obs_seq_arr[0] -1]<<endl;
		alpha[0][i] = PI[i] *B[i][obs_seq_arr[0] -1];
		//cout<<alpha[0][i]<<" ";
	}
	//cout<<endl;
	//Induction
	//cout<<"\n ..............  Printing A[i][j]   .......................\n";
	/*for(int c=0; c<N; c++)
	{
		for(int y=0; y<N; y++)
		{
			cout<<A[c][y]<<" ";
		}
		cout<<endl;
	}   */
	//cout<<"..................  Alpha .............................. \n";
	for(int t=1; t<no_of_obs_seq; t++)
	{
		for(int j =0; j<N; j++)
		{
			long double sum_val= 0;
			for(int i=0; i<N; i++)
			{

				sum_val += alpha[t-1][i] *A[i][j];
				//cout<<"(" <<alpha[t-1][i]<<"  "<<A[i][j]<<"  )";
			}
			//if(j==0)
				//cout<<"sum_val "<<sum_val<<endl;
			alpha[t][j] = sum_val * B[j][obs_seq_arr[t] -1];
			//cout<<alpha[t][j]<<" ";
		}
		//cout<<endl;
	}
	//Termination
	//cout<<"addition of alpha to get prob  \n";
	for(int i=0; i<N; i++)
	{
		//cout <<alpha[no_of_obs_seq -1][i]<<" ";
		prob_obs_seq_given_model +=  alpha[no_of_obs_seq -1][i];
	}
	//cout<<endl;
	return prob_obs_seq_given_model;
}


void backward_procedure(int obs_seq_arr[])
{
	
	/******************     Backward Procedure     ******************/
	//Initialization
	for(int i=0; i<N; i++)
		beta[no_of_obs_seq-1][i] =1;

	//Induction
	for(int t=no_of_obs_seq-2; t>=0; t--)
	{
		for(int i =0; i<N; i++)
		{
			long double sum_val =0;
			for(int j=0; j<N; j++)
				sum_val += A[i][j] *B[j][obs_seq_arr[t +1] -1] *beta[t+1][j];

			beta[t][i] =sum_val;
		}
	}

}


/****************** Sol2(Viterbi Algorithm)  *****************************/
void sol2(int obs_seq_arr[], int model_no )
{
	long double maxm_delta, delta_output;

	//Initialization
	for(int i =0; i<N; i++)
	{
		delta[0][i] = PI[i] *B[i][obs_seq_arr[0] -1];
		psi[0][i] =0;
	}
	//Recursion
	for(int t =1; t<no_of_obs_seq; t++)
	{
		for(int j=0; j<N; j++)
		{
			int i=0;
			maxm_delta =delta[t-1][i] *A[i][j];
			psi[t][j] =i;
			for(i =1; i<N; i++)
			{
				delta_output =delta[t-1][i] *A[i][j];
				if(delta_output >maxm_delta)
				{
					maxm_delta =delta_output;
					psi[t][j] =i;
				}
			}
			delta[t][j] = maxm_delta * B[j][obs_seq_arr[t] -1];
		}
	}

	//Termination
	int i =0;
	long double maxm_prob_at_T;
	maxm_prob_at_T =delta[no_of_obs_seq-1][i];
	qstar[model_no][no_of_obs_seq -1] =i;

	for(int i =1; i<N; i++)
	{
		if(maxm_prob_at_T < delta[no_of_obs_seq-1][i] )
		{
			maxm_prob_at_T =delta[no_of_obs_seq-1][i];
			qstar[model_no][no_of_obs_seq -1] =i;
		}
		
	}
	pstar[model_no] =maxm_prob_at_T;

	//Back-Tracking
	for(int t =no_of_obs_seq-2; t>=0; t--)
	{
		qstar[model_no][t] =psi[t+1][qstar[model_no][t+1]];
	}
}


void calculate_gamma()
{
	//cout<<"\n                          gamma                                      \n\n";
	for(int t=0; t<no_of_obs_seq; t++)
	{
		long double product_sum =0 ;
		//cout<<"   Alpha   Beta      \n";
		for(int i=0; i<N; i++)
		{
			//cout<<alpha[t][i] <<" "<<beta[t][i]<<"      ";
			gamma[t][i] = alpha[t][i]* beta[t][i];
			product_sum += gamma[t][i];
		}
		//cout<<endl<<"     Gamma Array     \n";
		for(int i=0; i<N; i++)
		{
			gamma[t][i] /= product_sum;
			//cout<<gamma[t][i]<<" ";
		}
		//cout<<endl;
	}
}

/**********************  Sol3 (re-estimation problem)   *******************************/
void sol3(int obs_seq_arr[])
{
	int max_in_b_row_index;
	long double threshold = pow(10.0, -30.0);
	//Calculate gamma array.
	calculate_gamma();

	//cout<<"  gamma   calculated  ........................... \n\n";
	//system("pause");
	//Calculate zyi array.
	for(int t =0; t< no_of_obs_seq-1; t++)
	{
		long double product_sum =0;
		for(int i =0; i<N; i++)
		{
			for(int j=0; j<N; j++)
			{
				zyi[t][i][j] = alpha[t][i] *A[i][j] *B[j][obs_seq_arr[t+1]] *beta[t+1][j];
				product_sum += zyi[t][i][j];
			}
		}

		for(int i =0; i<N; i++)
		{
			for(int j=0; j<N; j++)
				zyi[t][i][j] /= product_sum;
		}
	}


	//MAKE LAST ROW OF zyi 0
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			zyi[no_of_obs_seq -1][i][j]=0;
		}
	}

	
	/**************   Re-estimation Of Model Parameter   *****************/

	//TO CALCULATE PI 
	for(int i=0;i<N;i++)
	{
		PI[i]=gamma[0][i];
		//cout<<pi_bar[i]<<"  ";
	}
	
	//update A (Transition Matrix)
	for(int i=0; i<N-1; i++)
	{
		long double row_sum_A=0;
		long double denomintr_sum =0;
		for(int t_hat=0; t_hat<no_of_obs_seq -1; t_hat++)	   //0 to (no_of_obs_seq -2) times loop will run.
			denomintr_sum +=gamma[t_hat][i];
		
		
		for(int j=0; j<N; j++)
		{
			long double numratr_sum =0;
			for(int t=0; t<no_of_obs_seq -1; t++)
				numratr_sum += zyi[t][i][j];

			//update A[i][j]
			//cout<<numratr_sum<<" ...... "<<denomintr_sum<<"   ";
			A[i][j] = numratr_sum/denomintr_sum;
		}

		//If row sum in A Matrix  is not 1 make it 1.
		long double max_in_a_row =0;
		int idx=-1;
		for(int j=0; j<N; j++)
		{
			if(max_in_a_row < A[i][j])
			{
				max_in_a_row =A[i][j];
				idx =j;
			}
			row_sum_A += A[i][j];
		}
		
		if(row_sum_A <1)
		{
			long double diffr = 1- row_sum_A;
			 A[i][idx] += diffr;
		}
		//cout<<endl;
	}

	//Update B Matrix.
	//cout<<"\n\n.,.,.,.,.,.,,.,..     B Matrix        .,.,.,.,.,.,,.,.,.,.,.,.,.,.,.,."<<endl;


	//TO CALCULATE B_BAR
	//long double threshold;
	long double val1=10;
	long double val2=-30;
	threshold=pow(val1,val2);

	for(int j=0;j<N;j++)
	{
		int count=0;
		long double max=0;
		int index=0;
		for(int k=0; k<codebook_size; k++)
		{
			long double temp1=0;
			long double temp2=0;
			for(int t=0;t< no_of_obs_seq;t++)
			{
				if((obs_seq_arr[t]-1)==k)
				{
					temp1=temp1+gamma[t][j];
				}
				temp2=temp2+gamma[t][j];
			}
			//cout<<temp1<<" "<<temp2<<endl;
			B[j][k]=temp1/temp2;

			//cout<<B[j][k]<<"  ";

			if(B[j][k]>max)
			{
				max=B[j][k];
				index=k;
			}

			if(B[j][k]<threshold)
			{
				B[j][k]=threshold;
				count++;
			} 
		}
		//cout<<endl;

		B[j][index]=max-count*threshold;  
	}
	///////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	
	/*for(int j = 0; j < N; j++)
	{
		long denomtr_sum = 0;
		for(int t_hat = 0; t_hat < no_of_obs_seq; t_hat++)
			denomtr_sum += gamma[t_hat][j];		// denominator computation done		b[j][obs_seqs[seq_num][t+1]]

		int count_b_zeros = 0;
		for(int k = 0; k < codebook_size; k++)
		{
			long double numer_sum = 0;
			for(int t = 0; t < no_of_obs_seq; t++)
			{
				if( (obs_seq_arr[t] -1) ==k)
					numer_sum += gamma[t][j];
			}
			cout<<numer_sum<<" "<<denomtr_sum<<endl;
			B[j][k] = numer_sum / denomtr_sum;
		}

		max_in_b_row_index = find_B_rowmax(j);

		//adding threshold
		for(int k = 0; k < codebook_size; k++)
		{
			if(B[j][k] < threshold)
			{
				count_b_zeros++;
				B[j][k] = threshold;
			}
		}

		// subtracting value from max
		B[j][max_in_b_row_index] -= count_b_zeros * threshold;  
	}  */

	//check_rowsum();			// check rowsum of matrices A and B after re-estimation; currrently commented out inside the function
	/*long double sum_B_row;
	for(int w=0; w<5; w++)
	{
		sum_B_row =0.0;
		for(int v=0; v<32;v++)
		{
			B[w][v] = (B[w][v]+(long double)0.03)/((long double)1+(long double)0.07*(long double)32);
			sum_B_row +=B[w][v];
			//cout<<B[w][v]<<" ";
		}
		//cout<<"Row sum is "<<sum_B_row<<endl;
	}  */
	//system("pause");

}


void run_model(int obs_seq_arr[])
{
	int turn =0;
	while(turn !=20)
	{
		//Solution 1.
		forward_procedure(obs_seq_arr, turn);
		backward_procedure(obs_seq_arr);
		//Solution 2.
		sol2(obs_seq_arr, turn);

		//Solution 3.
		sol3(obs_seq_arr);

		turn++;
	}
}

void process_digit(int training, int utterance_num, int avg_itr, string digit)
{

	/**************** getting each frame one by one using sliding window    *****************/
	
	long double frame_arr[320];
	int obs_seq_arr[total_no_of_frame];

	for(int frm_no=0; frm_no <total_no_of_frame; frm_no++)
	{
		memset(frame_arr, 0, 320);
		vector<long double> R(13, 0);
		vector<long double> arr(13, 0);
		vector<long double> C(13, 0);

		/*
		cout<<"\n*************************************************************\n";
		cout<<frm_no<<endl;
		cout<<"*************************************************************\n";
		system("pause");   */
		
		extract_frame(frm_no, frame_arr);
		apply_hamming_window(frame_arr);
		find_r0_to_r12(R, frame_arr);
		find_a1_to_a12(arr, R);
		find_c1_to_c12(C, arr);
		//cout<<"\n\nAfter ci........................................\n";
		//cout<<endl;
		get_obs_seq_frm_featurevectr(obs_seq_arr, frm_no, C);

	}
	
	//////////////
	//Delete
	//cout<<"/////////////////////////////////////////////////////////////////////////";
	//cout<<"/*********************************************************"<<endl;
	/*for(int i=0; i<total_no_of_frame; i++)
		cout<<obs_seq_arr[i]<<" ";
	cout<<endl;
	system("pause");  */
	//cout<<"/**********************************************************"<<endl;   



	/****************************************************************************************/
	/************** **         Start process of HMM Model                 ******************/

	
	if(training)
	{
		//Read Initial HMM Model.
		read_initial_model(avg_itr, digit);
		//cout<<"/**********************.........................................*********************\n";
		/*for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
				cout<<A[i][j]<< " ";
			cout<<endl;
		}  */
		//cout<<"\n******************.................................................*******************\n";
		//Run HMM Model for 20 iteration for a particular file.
		run_model(obs_seq_arr);
		
		//cout<<"/.................********************************************..........................\n";
		//cout<<"After run_model (i.e.. 20 iteration)";
		/*for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
				cout<<A[i][j]<< " ";
			cout<<endl;
		}  */
		//cout<<"\n..........................*****************************************.....................\n";

		//After 20 iterations of a file whatever A, B you got add to the average one so, that after 
		//All recorded files for a digit ends then we can average them.
		add_to_average();
		/*for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
				cout<<A_mean[i][j]<< " ";
			cout<<endl;
		}  */
	}
	//Testing
	else
	{
		cout<<"Training Done  \n";
		//system("pause");
		long double prob_arr[10], maxm_prob=0;
		int recognized_digit= -1;
		for(int num=0; num<10; num++)              // Digit iteration:: for each digit from training I will  iterate to get maxm probability
		{
			string A_train = "A_prev_model_for_" + digit_arr[num] +".txt";
			string B_train = "B_prev_model_for_" + digit_arr[num] +".txt";
			string PI_train = "PI_prev_model_for_" + digit_arr[num] +".txt";
			
			ifstream f_A;
			ifstream f_B;
			ifstream f_PI;

			f_A.open(A_train);
			f_B.open(B_train);
			f_PI.open(PI_train);
			
			//cout<<A_train<<endl;
			int i=0, j=0;
			while(!f_A.eof())
			{
				f_A>> A[i][j++];
				//cout<<"j "<<j<<" "<<A[i][j++]<<" ";
				if(j==5) 
				{	i++;
					j=0;
					//cout<<endl; 
				}
			}
			
			//system("pause");
			i=0, j=0;
			long double b_row_sum =0;
			cout<<B_train<<endl;
			while(!f_B.eof())
			{
				long double val;
				f_B>> val;
				B[i][j++] =val;
				//cout<<val<<" ";
				b_row_sum +=val;
				
				if(j==32)
				{
					i++;  j=0; 
					//cout<<"\nb_row_sum "<<b_row_sum<<endl;
					b_row_sum =0;
				}
			}
			i=0;
			while(!f_PI.eof())
			{
				
				f_PI>> PI[i++];
				//cout<< PI[i++];
			}
			//cout<<".\n.................................................................^^^^^^^^^^^^^^^^^^^^^^^^^^..............\n";
			//cout<<"....................................^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n";
			//cout<<"Testing obs_seq_arr   \n";

			//for(int v=0; v<84; v++)
				//cout<<obs_seq_arr[v]<<" ";
			//cout<<endl;
			prob_arr[num] =forward_procedure(obs_seq_arr);
			/*for(int f=0; f<N; f++)
			{
				cout<<endl<<"alpha "<<alpha[no_of_obs_seq -1][f];
			}  */
			//cout<<"\nprob_arr["<<num<<"] = "<<prob_arr[num]<<endl;
			if(maxm_prob < prob_arr[num])
			{
				maxm_prob =prob_arr[num];
				recognized_digit =num;
			}
		}
	
	//cout<<"\n *************************************************************************************************\n"<<endl;
	if(maxm_prob >1)
		cout<<maxm_prob<<endl;
	cout<<"\n\n\nRecognised digit is "<<recognized_digit<<" ";
	}

}

/*********************   Average the model      **********************************/
/* Getting the average model for all utterances of a digit.   */

void get_average_model()
{
	cout<<"final A for Average model for a digit\n"; 
	//Averaging parameter A
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			A_mean[i][j] /= training_utterances;
			//cout<<A_mean[i][j]<<" ";
		}
		//cout<<endl;
	}

	cout<<"final B for Average model for a digit\n";
	
	//Averaging parameter B
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<codebook_size; j++)
		{
			B_mean[i][j] /= training_utterances;
			//cout<<B_mean[i][j]<<" ";
		}
		
		//cout<<endl;
	}

	cout<<"final PI for Average model for a digit\n";
	//Averaging parameter PI
	for(int i=0; i<N; i++)
	{
		PI_mean[i] /= training_utterances;
		//cout<<PI_mean[i]<<" ";
	}
}
 

int _tmain(int argc, _TCHAR* argv[])
{
	int training =1;
	// populating codebook array.
	ifstream fcodebook;
	fcodebook.open("codebook.txt");
	for(int i = 0; i < codebook_size; i++)
	{
		for(int j = 0; j <p; j++)
		{
			fcodebook >>codebook[i][j];
		}
	}
	fcodebook.close();


	//Finding Dc-Shift() value.
	long double dc_shift_value;
	dc_shift_value = find_dc_shift();
	
	/*******************************************************************************************************/
	/************************       Training to build the HMM Model      ***********************************/

	//cout<<"Training Started   \n";
	//system("pause");
	/*for(int digit = 0; digit < 10; digit++)
	{
		string digit_str = digit_arr[digit];
		for(int avg_itr = 0; avg_itr < avg_itration; avg_itr++)	// averaging out the models for this digit
		{
			initialize_avg_A_B();
			//cout<<"   \nyes\n";
			for(long long utterance_num=1; utterance_num <= training_utterances; utterance_num++)		// since filenames start with utterance no. 1 
			{
				
				string filename = "hmm_digit/" + digit_str + "_" + to_string(utterance_num) + ".txt";
				
				int peak_indx= preprocessing(filename, dc_shift_value);
				//cout<<peak_indx<<endl;
				cout<<filename<<endl;
				//system("pause");
				/****************  finding the 7040 samples which needs to be considered  ***************/

				//Extract the 7040 samples near to peak sample.
			/*	int start_marker = peak_indx - (samplesneeded/2);
				int end_marker = peak_indx + (samplesneeded/2);
				//If peak_indx is shifted too much at left side.
				if(peak_indx <(samplesneeded/2)-1)
				{
					start_marker =540;
					end_marker =7580;
				}
				//If peak_indx is too much shifted right side.
				if(end_marker > sample_vectr.size() -1)
				{
					end_marker =sample_vectr.size()-1;
					start_marker =sample_vectr.size() -samplesneeded -1;
				}
				cout<<start_marker<<" "<<end_marker<<endl;
				//copying all 7040 samples near to peak sample into word_vectr
				for(int i=start_marker; i<end_marker; i++)
					word_vectr.push_back(sample_vectr[i]);
				sample_vectr.clear();
				process_digit(training, utterance_num, avg_itr, digit_str);
				word_vectr.clear();
				cout<<"utterance number    "<<utterance_num<<endl;
				system("pause");
			}
			
			get_average_model();
			write_final_model(digit_str);
			//cout<<"For a digit all utterances done "<<avg_itr<<endl;
			//system("pause"); 
		}
		//system("pause");
	}   */

	/////////////////////////////////////////////////////////////////////////////////////////////



	/**************************  Testing()         *********************************************/

	for(int test_digit = 0; test_digit < 10; test_digit++)
	{
		//int test_digit;
		//cout<<"Enter the digit you want to test \n";
		//cin>>test_digit;
		for(long long utterance_num=1; utterance_num <= 10; utterance_num++)		// since filenames start with utterance no. 1 
		{
			long long test_utterance =	utterance_num+15;
			string test_digit_str =digit_arr[test_digit];
			string filename2 = "hmm_digit/" + digit_arr[test_digit] + "_" + to_string(test_utterance) + ".txt";

			int peak_indx= preprocessing(filename2, dc_shift_value);
			/****************  finding the 7040 samples which needs to be considered  ***************/

			//Extract the 7040 samples near to peak sample.
			int start_marker = peak_indx - (samplesneeded/2);
			int end_marker = peak_indx + (samplesneeded/2);
			//If peak_indx is shifted too much at left side.
			if(peak_indx <(samplesneeded/2)-1)
			{
				start_marker =0;
				end_marker =7040;
			}
			//If peak_indx is too much shifted right side.
			if(end_marker > sample_vectr.size() -1)
			{
				end_marker =sample_vectr.size() -1;
				start_marker =sample_vectr.size() -samplesneeded -1;
			}
			//copying all 7040 samples near to peak sample into word_vectr
			for(int i=start_marker; i<end_marker; i++)
				word_vectr.push_back(sample_vectr[i]);
			sample_vectr.clear();
			training =0;       //for testing training is 0.
			//In testing no use of avg_itr so, in place of that passing 0.
			process_digit(training, utterance_num, 0, test_digit_str);	
			cout<<"for "<<test_digit<<" utterance number "<<utterance_num<<endl;
			system("pause");
			word_vectr.clear();		

		}
	}
}

