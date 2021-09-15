

/*  Error - some times 'o' becomes more than 'x' i.e.. computer turn is always after my turn so, no. 
	of 'o'  should not be greater than no. of 'x'   */


#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <Windows.h>
#include <MMSystem.h>


#define  TRUE 1
#define samplesneeded 7040
#define slider 80
#define total_no_of_frame 85
#define no_of_obs_seq 85
#define codebook_size 32
#define training_utterances 25
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

char initial_board[3][3] = {'0', '1', '2', '3', '4', '5', '6', '7', '8'};
char board[3][3] ={'0', '1', '2', '3', '4', '5', '6', '7', '8'};
int winning_matrix[][3] ={
                            {0, 1, 2},
                            {3, 4, 5},
                            {6, 7, 8},
                            {0, 3, 6},
                            {1, 4, 7},
                            {2, 5, 8},
                            {0, 4, 8},
                            {2, 4, 6}
                        };



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



int process_digit_test()
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
		
		get_obs_seq_frm_featurevectr(obs_seq_arr, frm_no, C);

	}
	//Testing
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
	
	
	//cout<<"\n\n\nRecognised digit is "<<recognized_digit<<" ";
	//system("pause");
	return recognized_digit;
}


int enter_place()
{
	//Finding Dc-Shift() value.
	long double dc_shift_value;
	dc_shift_value = find_dc_shift();
	
	int place;		
	system("Recording_Module.exe 3 voice_reccorded.wav voice_recorded.txt");
	/* Please change the recorded voice file path before running the program and make sure
	file path should be correct .*/
	//finobj.open("voice_recorded.txt");
		
		
	int peak_indx= preprocessing("voice_recorded.txt", dc_shift_value);
	/****************  finding the 7040 samples which needs to be considered  ***************/
	
	//Extract the 7040 samples near to peak sample.
	int start_marker = peak_indx - (samplesneeded/2);
	int end_marker = peak_indx + (samplesneeded/2);
	//cout<<start_marker<<"  "<<end_marker<<" "<<peak_indx<<endl;
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
	//cout<<endl<<word_vectr[3516]<<" "<<word_vectr[3517]<<" "<<word_vectr[3518]<<" "<<word_vectr[3519]<<" "<<word_vectr[3520]<<endl;
	sample_vectr.clear();
	//In testing no use of avg_itr so, in place of that passing 0.
	place =process_digit_test();	
	word_vectr.clear();		

	return place;
}

void show_board()
{
    cout<<"\n\n\n\n\n\n\n\n";
    cout<<"\t\t\t     |     |     "<<endl;
    cout<<"\t\t\t  "<<board[0][0]<<"  |  "<<board[0][1]<<"  |  "<<board[0][2]<<"   "<<endl;
    cout<<"\t\t\t_____|_____|_____"<<endl;
    cout<<"\t\t\t     |     |     "<<endl;
    cout<<"\t\t\t  "<<board[1][0]<<"  |  "<<board[1][1]<<"  |  "<<board[1][2]<<"   "<<endl;
    cout<<"\t\t\t_____|_____|_____"<<endl;
    cout<<"\t\t\t     |     |     "<<endl;
    cout<<"\t\t\t  "<<board[2][0]<<"  |  "<<board[2][1]<<"  |  "<<board[2][2]<<"   "<<endl;
    cout<<"\t\t\t     |     |     "<<endl;

}

void delay(int number_of_seconds) 
{ 
    // Converting time into milli_seconds 
    int milli_seconds = 1000 * number_of_seconds; 
  
    // Storing start time 
    clock_t start_time = clock(); 
  
    // looping till required time is not achieved 
    while (clock() < start_time + milli_seconds) 
        ; 
} 


//void
int player_choice()
{
    int place =-1;
    int row, column, count =0;
    do
    {
		if(count !=0)
		{
			cout<<"\n\t\t !!kindly enter  empty slots !!!\n";
			delay(5);
		}
        cout<<"\n Your Turn.\n";
		place =enter_place();

        //cin>>place;
        row =place/3;
        column = place %3;
		count++;
    }
    while(board[row][column] =='x' || board[row][column] =='o' );

    board[row][column] ='x';
	return place;
}

void computer_choice()
{
    int place =-1;
    int row, column;
    do
    {
        srand(time(0));
        place = rand();
        (place>8)?(place %=9):place;
        //cout<<"computer opted place \n";
        //cout<<place<<endl;
        row =place/3;
        column = place %3;
    }
    while(board[row][column] =='x' || board[row][column] =='o' );

    board[row][column] ='o';
}

int count_place(char ch)
{
    int total =0;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            if(board[i][j] == ch)
                total +=1;
        }

    }

    return total;
}

char check_winner()
{
    int winning_count;
    winning_count = sizeof(winning_matrix)/(sizeof(int) * 3);
    for(int i=0; i<winning_count; i++)
    {
        int row_0, column_0, row_1, column_1, row_2, column_2;
        row_0 = winning_matrix[i][0]/3;
        column_0 = winning_matrix[i][0] %3;
        row_1 = winning_matrix[i][1]/3;
        column_1 = winning_matrix[i][1] %3;
        row_2 = winning_matrix[i][2]/3;
        column_2 = winning_matrix[i][2] %3;
        if (board[row_0][column_0] == board[row_1][column_1]  && board[row_0][column_0] == board[row_2][column_2] )
            return board[row_0][column_0];
    }

}


void play_tic_tac_toe()
{
	int place;
	cout<<"\n\n\t\t  !!!! HELLO, "; 
	bool played = PlaySound(TEXT("Welcome.wav"), NULL, SND_SYNC);
	cout<<"Welcome to Tic Tac Toe!!!!!  \n";
	show_board();
    while(TRUE)
    {
        system("cls");
		cout<<"\n\n\t\t\t !!!!!! WELCOME TO TIC TAC  GAME  !!!"<<endl;
        show_board();
        if(count_place('x') == count_place('o'))
        {

            place= player_choice();
            system("cls");
			cout<<"\n position u opted is "<<place<<endl;
            show_board();
			delay(3);
            computer_choice();
            system("cls");
            show_board();

        }
        if(check_winner() =='x')
        {
			bool pl = PlaySound(TEXT("Congratulation.wav"), NULL, SND_SYNC);
            cout<<"\n\n\t\t!!!!  Congratulations! YOU WON   !!! \n";
            break;
        }
        else if(check_winner() =='o')
        {
            cout<<"\n\n\t\t!!! computer won Try again !!!!\n";
            break;
        }
        else if( count_place('x') + count_place('o') == 9)
        {
            cout<<"\n\n\t\t!!!! It's a   Tie  !!!!!\n";
            break;
        }
    }
	system("pause");
} 



bool detect_yes_no()
{
	bool flag =0;
	ifstream finobj;
	
	system("Recording_Module.exe 3 yes_no.wav yes_no.txt");
	/* Please change the recorded voice file path before running the program and make sure
	   file path should be correct .*/
	//finobj.open("voice_recorded.txt");
	finobj.open("yes_no.txt");
	if(!finobj)
	{
		cout<<"Input voice file can't be opened\n";
		return 0;
	}
	else
	{
		long long ste, sum =0;						//ste is Short term Energy variable. sum will contain summation of square of sample values.
		int count =0, zcr=0;
		long long samplevalue, prevsamplval =0;
		vector<int> zcrvec ;                      //"zcrvec", is  vector to store  zcr values of each frame(chunk of 320 samples). 
	    vector<long long> stevec;				  //"stevec", is  vector to store  ste values of each frame(chunk of 320 samples).                 
		while(!finobj.eof())						//Start reading the input file till end.
		{

			finobj>>samplevalue;					//sending the value read by "finobj" (object of ifstream) in "samplevalue".

			/*checking whether waveform crosses x-axis, to calculate zcr.*/
			if((prevsamplval>0 && samplevalue <0) || (prevsamplval <0 && samplevalue >0))
						zcr++;
			sum +=(samplevalue *samplevalue);		//squaring the samplevalue and adding into sum.
			count++;								//increasing the count by 1, after reading each samplevalue.
			
			if(count ==320)							 //As standard, 320 sample is considered as 1 frame
			{
				ste =sum/320;						//As, count reaches to 320 calculating short-term-energy(STE).
				zcrvec.push_back(zcr);				//pusing "zcr" value for this frame in "zcrvec" vector.
				stevec.push_back(ste);				//pusing "ste" value for this frame in "stevec" vector
				count =0;							//Making count as 0. So, that again start counting "320 sample" for next-frame.
				sum =0;								//completed sum calculation for current frame, now start for next frame so, making "sum" 0.
				zcr =0;								//completed zcr calculation for current frame, now start for next frame so, making "zcr" 0.
			}   
			prevsamplval =samplevalue;			   //Assigning samplevalue to prevsamplval so that, sign change can be detected to calculate zcr.

		} //while-loop end
		
		int word_start_marker =-1, word_end_marker =-1;		//variables to show word-boundaries.
		//bool flag =0;										//flag will be 1 only if zcr will be more than 100. ie.. for "yes" case.

		//Iterating through Short term energy vector "stevec" to get word-boundaries.
		vector<long long>::iterator itr;
		//cout<<"The Uttered sequence of yes no is :-\n\n";
		for(itr= stevec.begin(); itr< stevec.end(); itr++)
		{
			/* when I get ste(short term energy) greater than "65000" (taken as threshold value)  then
			   will mark that index as  word  start boundary, only if this is the first index 
			   (i.e .. word_start_marker ==-1) with ste value greater than threshold  for a 
			   particular word (i.e.. after silence). If word_start_marker != -1 and ste value> 65000 
			   this shows we are in the middle of word, and word_start_marker is already assigned with value.  */

			if(*itr> 65000 &&(word_start_marker ==-1))			
				word_start_marker =itr - stevec.begin();

			/* when I get ste(short term energy) less than "15000" and if already got word_start_
			   marker (i.e.. word_start_marker!=-1) for that word then only I will mark this index 
			   as word_end_marker.  */
			if((*itr < 10000) && (word_start_marker!=-1))
				word_end_marker = itr - stevec.begin();

			/* if I have valid word_start_marker and word_end_marker i.e.. if both  have 
			   value  greater than  0 this means word has started and already marked 
			   both start and end boundaries. so, I will investigate zcr value for that word 
			   If I get any frame(chunk of 320 sample) with zcr value >125 then this shows 
			   it  is "yes", for no zcr is low.(it is between 20-45).     */
			if((word_start_marker >=0) && (word_end_marker >0))   
			{
				while(word_start_marker< word_end_marker)         //we will iterate until a word ends.
				{
					/* if zcr value is greater than 100 with in word boundaries that means it is "yes".
					   So, flag is set.As, last part of yes is giving high zcr value. For "no" zcr is
					   not having value greater than  45. */
					if(zcrvec[word_start_marker++] >100)          
						flag =1;
				}
				//if(flag)                                 //flag is set so, this word is yes.
					//cout<<"yes ";
				//else									//otherwise, this is no.
					//cout<<"no "; 
				//unset the flag so, that again start investigation for next word.
				//flag =0;
				/* will again set both word start and end marker as -1 so, that can again iterate for next word. */
				//word_start_marker =-1, word_end_marker =-1;      
			}
			
		}	
	}
	//system("pause");
	return flag;
}


void initialize_board()
{

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			board[i][j] = initial_board[i][j];
	}
}



int _tmain(int argc, _TCHAR* argv[])
{
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
	
	/*******************************************************************************************************/
	bool flag_yes_no=0;
	system("pause");
	do{
		initialize_board();
		play_tic_tac_toe();
		delay(2);
		cout<<"\n\nDo You Want to play again? If Yes, speak \"Yes\" otherwise  speak \"NO\\n"<<endl;
		flag_yes_no =detect_yes_no();
		if(flag_yes_no)	
			cout<<"yes";
		else
			cout<<"no";
		delay(2);
	}
	while(flag_yes_no);

}

