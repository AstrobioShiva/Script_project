
//#include <Riostream.h>
#include <TMath.h>

//auto pi = TMath::Pi();

// function code in C
/*
float f1(float x) {
return 2*x;
}
float f2(float x)
{
    return sin(x);
}

void user_defined(){
float a = 20.;
cout<<f1(a)<<endl;
cout<<f2(a)<<endl;
}
*/
void user_defined (){

std::vector<vector<int>> vect{ { 1, 2, 3 },
                              { 4, 5, 6 },
                              { 7, 8, 9 } };

// int vect[3][3];
// vect[0][0] = 11;
std::vector<vector<int>>vec(3, std::vector<int>(3,0));
// int vec[3][3];
 vec[0][0] = 10;
int vec1[3][3];
for (int i = 0; i<3; i++){
    for (int j = 0; j<3; j++){
        vec1[i][j] = vect[i][j] + vec[i][j];
        std::cout <<vec1[i][j]<< std::endl;
    }
}
}
