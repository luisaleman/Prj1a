#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "time.h"
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[])
{
    //To make sure the user writes the outputfile name while running the program
    char *outfilename;
    if( argc <=1 ){
        cout << "Bad usage " << argv[0] << " read also output file on the same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }

    int numberpoints, i;
    double coef, timediagonal, timetotal,timeinitial, maxerror;
    ofstream myfile;
    myfile.open(outfilename);

    cout << endl << "Please enter the number of grid points: ";
    cin >> numberpoints;
    clock_t start, mids, midf, finish, defs, deff;                       //We declare the values to measure time
    start = clock();
    double h= (double)1/(numberpoints+1);
    cout << "Step length: "<< setprecision (8) << h << endl;
    defs = clock();
    //We delcare our 3 needed vectors for diagonalization
    vector<double> maindiagonal (numberpoints,2);
    vector<double> diagonal (numberpoints,-1);
    vector<double> sourcevector (numberpoints,0);
    //We declare the other vectors needed to print the results
    vector<double> resultvector (numberpoints,0);
    vector<double> relativeerror(numberpoints,0);
    vector<double> y (numberpoints,0);

    //We enter the values of our sourcevector and analytic solution "y"
    for(i=0;i<numberpoints;i++){
        sourcevector[i]=pow(h,2)*100*exp(-10*((i+1)*h));
        y[i]= (double) (1-((1-exp(-10))*((i+1)*h))-exp(-10*(i+1)*h));
    }
    deff = clock();
    timeinitial=(double)(deff-defs)/CLOCKS_PER_SEC;
    mids = clock();
    //We start forward substitution
    for(i=0;i<numberpoints-1;i++){
        coef= diagonal[i]/maindiagonal[i];                          //1 flop
        maindiagonal[i+1]-= diagonal[i]*coef;                       //2 flops
        sourcevector[i+1]-=sourcevector[i]*coef;                    //2 flops
     }
    //We start backward sbustitution, ending the "diagonalization" process
    for(i=numberpoints-1;i>0; i--){
        coef=diagonal[i]/maindiagonal[i];                           //1 flop
        sourcevector[i-1]-=sourcevector[i]*coef;                    //2 flops
    }
    midf = clock();
    timediagonal= ((midf-mids)/(double)CLOCKS_PER_SEC);            //diagonalization time

    //We calculate the results and print them
    for(i=0;i<numberpoints;i++){
        resultvector[i]=sourcevector[i]/maindiagonal[i];            //1 flop
        relativeerror[i]=log10(fabs((resultvector[i]-y[i])/y[i]));
        myfile << setw(15) << setprecision(8) << (i+1)*h;
        myfile << setw(15) << setprecision(8) << resultvector[i];
        myfile << setw(15) << setprecision(8) << y[i];
        myfile << setw(15) << setprecision(8) << relativeerror[i]<<endl;
    }

    myfile.close();
    //We calculate the maximum relative error, we set a really low number as default
    //value of maxerror
    maxerror=(double) -100;
    for(unsigned int j=0;j<numberpoints;j++){
        if(relativeerror[j]>maxerror){
            maxerror=relativeerror[j];
        }
    }
    cout << "Max error: "<< setprecision(8) << maxerror << endl;
    finish = clock();                                               //total time
    timetotal= ((finish-start)/(double)CLOCKS_PER_SEC);
    cout << "Time used for initialization: "<< setprecision(8) << timeinitial<< " seg" <<endl;
    cout << "Time used for diagonalization: "<< setprecision(8) << timediagonal <<" seg"<<endl;
    cout << "Total time: " << setprecision(8)<<timetotal <<" seg" << endl;

    return 0;

}

//I commented in my project that i made two programs, the one below was the first one, which
//defines a whole tridiagonal matrix, but a I couldn't go further than 10^4 grid points
// I decided to create another one working only in vectors. I just include it because i mention
//it in my pdf and compare it with the second program (above)
/*

int main(int argc, char* argv[])
{
    char *outfilename;
    // Declare the number of steps and two indeces I'll use in "for" loops
    int number_points,dimension, i,j;
    double coef, calculation,initial, total, maxerror;

    if( argc <=1 ){
        cout << "Bad usage " << argv[0] << " read also output file on the same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }
    ofstream myfile;
    myfile.open (outfilename);

    cout << endl << "Please enter the number of grid points: ";
    cin >> number_points;
    clock_t start, mids, midf, finish, decs,decf;
    start = clock();
    double h= (double)1/(number_points+1);
    dimension=number_points;
    //for the dimension is the number of grid points, and # of steps is grid points+1

    //First, we will allocate memory for an array which contains a set of pointers.
    //Next, we will allocate memory for each array which is pointed by the pointers
    decs = clock();
    double **tridmatrix=0;
    tridmatrix= new double*[dimension];
    for(i=0;i<number_points;i++){
        tridmatrix[i]= new double[dimension];
        }

    //Then we start the arrays of our source term, our result, our relative error
    double *sourcevector=0;
    sourcevector= new double[dimension];
    double *resultvector=0;
    resultvector= new double[dimension];
    double *relativeerror=0;
    relativeerror= new double[dimension];
    double *y=0;
    y= new double[dimension];

    //We enter the values of our tridiagonal matrix, source vector and analyticsolution
    for(i=0;i<dimension-1;i++){
        tridmatrix[i][i]=2;
        tridmatrix[i][i+1]=-1;
        tridmatrix[i+1][i]=-1;
        sourcevector[i]=pow(h,2)*100*exp(-10*((i+1)*h));
        y[i]= 1-((1-exp(-10))*((i+1)*h))-exp(-10*(i+1)*h);
    }
    tridmatrix[dimension-1][dimension-1]=2;
    sourcevector[dimension-1]=pow(h,2)*100*exp(-10*((dimension*h)));
    y[dimension-1]= 1-((1-exp(-10))*((dimension)*h))-exp(-10*(dimension)*h);
    decf=clock();
    initial=(double)(decf-decs)/CLOCKS_PER_SEC;
    mids = clock();
    //Now we start our forward substitution
    for(i=0;i<dimension-1;i++){
        coef=(double)(tridmatrix[i+1][i]/tridmatrix[i][i]);                      //1 flop
        //tridmatrix[i+1][i]-= tridmatrix[i][i]*coef;                     //2 flops  **I dont need to take this into acount
        tridmatrix[i+1][i+1]-= tridmatrix[i][i+1]*coef;                 //2 flops
        sourcevector[i+1]-=sourcevector[i]*coef;                        //2 flops
    }

    //Now we start our backward substitution
    resultvector[dimension-1]=sourcevector[dimension-1]/tridmatrix[dimension-1][dimension-1];      //1 flop
    for(i=dimension-1;i>0;i--){
        coef=tridmatrix[i-1][i]/(tridmatrix[i][i]);                         //1 flop
        sourcevector[i-1]-=sourcevector[i]*coef;                            //2 flops
        resultvector[i-1]=sourcevector[i-1]/tridmatrix[i-1][i-1];           //1 flop
    }
    midf = clock();
    calculation=(double)((midf - start)/(double)CLOCKS_PER_SEC);

    //cout<< "The results are:" << endl;
    for(i=0;i<dimension-1;i++){
        relativeerror[i]=log10(fabs((resultvector[i]-y[i])/y[i]));
        //cout << setw (15) << setprecision(8) << (i+1)*h<< setw(15)<< resultvector[i] <<
               // setw(15) << y[i] << setw(15)<< relativeerror[i] <<endl;
        myfile << setw(15) << setprecision(8) << (i+1)*h;
        myfile << setw(15) << setprecision(8) << resultvector[i];
        myfile << setw(15) << setprecision(8) << y[i];
        myfile << setw(15) << setprecision(8) << relativeerror[i]<<endl;
    }
    myfile.close();
    maxerror=(double) -100;
    for(unsigned int j=0;j<dimension;j++){
        if(relativeerror[j]>maxerror){
            maxerror=relativeerror[j];
        }
    }
    finish = clock();
    total= ((finish-start)/(double)CLOCKS_PER_SEC);
    cout << "Time used for initialization: "<< setprecision(8) << initial<< " seg" <<endl;
    cout << "Time used for diagonalization: "<< setprecision(8) << calculation <<" seg"<<endl;
    cout << "Total time: " << setprecision(8)<<total <<" seg" << endl;
    return 0;

}
*/


