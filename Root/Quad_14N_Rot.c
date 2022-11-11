#include <iostream>
#include <TH1D.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TRandom.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TROOT.h>
#include<TPaveText.h>
#include <TMath.h>

#include <vector>

/// @brief Code to extract antisymmetric shielding tensor components
auto pi = TMath::Pi();

//  void Rabc( float alfa1,  float beta1, float gama1){
//           double U[1][1];
//          U[0][0] = cos(alfa1)*cos(beta1)*cos(gama1) - sin(alfa1)*sin(gama1);
//    //      U[0][2] = sin(alfa1)*cos(beta1)*cos(gama1) + cos(alfa1)*sin(gama1);
//        U(0,2) = -sin(beta1)*cos(gama1);
//          U(1,0) = -cos(alfa1)*cos(beta1)*sin(gama1) - sin(alfa1)*cos(gama1);
//          U(1,1) = -sin(alfa1)*cos(beta1)*sin(gama1) + cos(alfa1)*cos(gama1);
//          U(1,2) = sin(beta1)*sin(gama1);
//          U(2,0) = cos(alfa1)*sin(beta1);
//          U(2,1) = sin(alfa1)*sin(beta1);
//          U(2,2) = cos(beta1);
         
    //  }


std::vector<std::vector<double>> Rabc(float alfa1,  float beta1, float gama1){
    std::vector<std::vector<double>> U(3, std::vector<double>(3, 0.));
  //  U[0][0] = TMath::Cos(alfa1) * TMath::Cos(beta1) * TMath::Cos(gama1) - TMath::Sin(alfa1) * TMath::Sin(gama1);
  U[0][0] = cos(alfa1)*cos(beta1)*cos(gama1) - sin(alfa1)*sin(gama1);
  
    return U;
}
// to be returned as function 2D array is to be defined as vector

void Quad_14N_Rot(){
    std::vector<std::vector<double>> matrix = Rabc(pi, pi, pi);
    for (unsigned int i = 0; i < matrix.size(); i++) {
        for (unsigned int j = 0; j < matrix[0].size(); j++) {
            std::cout << matrix[i][j] << std::endl;
        }
    }
}

/*

void Quad_14N_Rot() {
//    TCanvas *c1 = new TCanvas();
    int Nptx = 10; double freq2D[Nptx][1]; double freqSUM[Nptx][1]; double freqDIFF[Nptx][1]; 
    
     float XX[Nptx][1]; 
    
     int w0 = 100; int wx = w0*pow(10,6);
     int ntheta = 500; int nphi = 500;

     float dangle = 2*pi/Nptx;

     float Ispin = 1.; //spin of nuclei

    float CQ = 3.;            // CQ (MHz)
    float Qeta = 0.3;         // Quadrupole eta term
    float Siso = 0.;          // isotropic chemical shift + offset (ppm) 
    float  delta = 200.;        // chemical shift anisotrophy (CSA) (ppm)
    float eta = 0.5;          // CSA eta term
    float Sxy = 100.; float Sxz = -100.; float Syz = 100.;  // antisymmetric components (ppm) 

     Siso=Siso*w0;
     delta = delta*w0;
    //  cout<<delta<<endl;
      CQ = CQ*pow(10,6)/(2*Ispin*(2*Ispin-1));
       Sxy=Sxy*w0;  Sxz=Sxz*w0;  Syz=Syz*w0;

     float a=0*pi/180; float b=0*pi/180; float c=0*pi/180;
     float zeta = 0*pi/180; float lamda = 0*pi/180; float nu = 0*pi/180;
     float alpha = 0*pi/180; float beta = 0*pi/180; float gama = 0*pi/180;

    double QPAS[3][3];  double CSA[3][3]; double ACS[3][3];  
     QPAS[0][0] = (Qeta -1.)*CQ/2.;     // Quadrupolar interaction
     QPAS[1][1] = -(1.+Qeta)*CQ/2.; 
     QPAS[2][2] = CQ; 
    // for (int i = 0; i <3; i++){
    //     for (int j = 0; j <3; j++) {
    //         cout<<i<<j<<QPAS[i][j]<<endl;
    //     }
        
    // }

     CSA[0][0] = (eta -1.)*delta/2.; // 2nd-rank chemical shift anisotropy (CSA)
     CSA[1][1] = -(1.+eta)*delta/2.; 
     CSA[2][2] = delta;  

     ACS[0][1]= Sxy; ACS[0][2]= Sxz;   // 1st-rank antisymmetric chemical shift (ACS)
     ACS[1][0]=-Sxy; ACS[1][2]= Syz;
     ACS[2][0]=-Sxz; ACS[2][1]=-Syz;
     double CSPAS[3][3];
     for (int i = 0; i <3; i++){
        for (int j = 0; j<3; j++){
            CSPAS[i][j] = CSA[i][j] + ACS[i][j];
          
        }
     }
   
cout<<Rabc(20,30,30)<<endl;
    
}

*/
