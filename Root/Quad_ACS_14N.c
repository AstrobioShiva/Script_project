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
//auto cos = TMath::Cos();
 //auto sin = TMath::Sin();

const TMatrixD R(float alfa1,  float beta1, float gama1){
    TMatrixD U(3,3);
    U(0,0) = TMath::Cos(alfa1)*TMath::Cos(beta1)*TMath::Cos(gama1) - TMath::Sin(alfa1)*TMath::Sin(gama1);
    U(0,1) = TMath::Sin(alfa1)*TMath::Cos(beta1)*TMath::Cos(gama1) + TMath::Cos(alfa1)*TMath::Sin(gama1);
    U(0,2) = -TMath::Sin(beta1)*TMath::Cos(gama1);
    U(1,0) = -TMath::Cos(alfa1)*TMath::Cos(beta1)*TMath::Sin(gama1) - TMath::Sin(alfa1)*TMath::Cos(gama1);
    U(1,1) = -TMath::Sin(alfa1)*TMath::Cos(beta1)*TMath::Sin(gama1) + TMath::Cos(alfa1)*TMath::Cos(gama1); 
    U(1,2) = TMath::Sin(beta1)*TMath::Sin(gama1);
    U(2,0) = TMath::Cos(alfa1)*TMath::Sin(beta1);
    U(2,1) = TMath::Sin(alfa1)*TMath::Sin(beta1);
    U(2,2) = TMath::Cos(beta1);
    return U;
    
}
double Quad(double QXG, double ct, double st, double s2t, double c2t, double cp, double sp, double c2p, double s2p){
   sq15=sqrt(1.5); sq6=sqrt(6);
  R20Q=0.5*sq15*(QXG(2,2)*(3*ct*ct-1)+2*QXG(0,2)*s2t*cp+2*QXG(1,2)*s2t*sp+(QXG(0,0)-QXG(1,1))*st*st*c2p+2*QXG(0,1)*st*st*s2p);
  R2q1=0.75*QXG(2,2)*s2t;
  R2q2r=QXG(0,2)*c2t*cp; R2q2i=QXG(0,2)*ct*sp;
  R2q3r=QXG(1,2)*c2t*sp; R2q3i=QXG(1,2)*ct*cp;
  R2q4r=0.25*(QXG(0,0)-QXG(1,1))*s2t*c2p; R2q4i=0.5*(QXG(0,0)-QXG(1,1))*st*s2p;
  R2q5r=0.5*QXG(0,1)*s2t*s2p; R2q5i=QXG(0,1)]*st*c2p;
  R2p1Q=R2q1-R2q2r+1i*R2q2i-R2q3r-1i*R2q3i-R2q4r+1i*R2q4i-R2q5r-1i*R2q5i;
  R2m1Q=-R2q1+R2q2r+1i*R2q2i+R2q3r-1i*R2q3i+R2q4r+1*R2q4i+R2q5r-1i*R2q5i;
  
  R4q1=0.75*QXG(2,2)*st*st;
  R4q2r=0.5*QXG(0,2)*s2t*cp; R4q2i=QXG(0,2)*st*sp;
  R4q3r=0.5*QXG(1,2)*s2t*sp; R4q3i=QXG(1,2)*st*cp;
  R4q4r=0.25*(QXG(0,0)-QXG(1,1))*(1+ct*ct)*c2p; R4q4i=0.5*(QXG(0,0)-QXG(1,1))*ct*s2p;
  R4q5r=0.5*QXG(0,1)*(1+ct*ct)*s2p; R4q5i=QXG(0,1)*ct*c2p;
  R2p2Q=R4q1-R4q2r+1i*R4q2i-R4q3r-1i*R4q3i+R4q4r-1i*R4q4i+R4q5r+1i*R4q5i;
  R2m2Q=R4q1-R4q2r-1i*R4q2i-R4q3r+1i*R4q3i+R4q4r+1i*R4q4i+R4q5r-1i*R4q5i;
  return (R2m1Q, R2p1Q, R20Q, R2m2Q, R2p2Q))

}




void Quad_ACS_14N() {
//    TCanvas *c1 = new TCanvas();
    int Nptx = 10; 
    //initialization of matrices
    TMatrixD freq1D(Nptx,1);
    TMatrixD freq2D(Nptx,1);
    TMatrixD freqSUM(Nptx,1);
    TMatrixD freqDIFF(Nptx,1);
    
    TMatrixD X(Nptx,1); 
    
    int w0 = 100; int wx = w0*pow(10,6);
    int ntheta = 500; int nphi = 500;

    double dangle = 2*pi/Nptx;

    double Ispin = 1.; //spin of nuclei

    double CQ = 3.;            // CQ (MHz)
    double Qeta = 0.3;         // Quadrupole eta term
    double Siso = 0.;          // isotropic chemical shift + offset (ppm) 
    double  delta = 200.;        // chemical shift anisotrophy (CSA) (ppm)
    double eta = 0.5;          // CSA eta term
    float Sxy = 100.; float Sxz = -100.; float Syz = 100.;  // antisymmetric components (ppm) 

    Siso=Siso*w0;
    delta = delta*w0;
    //  cout<<delta<<endl;
    CQ = CQ*pow(10,6)/(2*Ispin*(2*Ispin-1));
    Sxy=Sxy*w0;  Sxz=Sxz*w0;  Syz=Syz*w0;

    float a=0*pi/180; float b=0*pi/180; float c=0*pi/180;
    float zeta = 0*pi/180; float lambda = 0*pi/180; float nu = 0*pi/180;
    float alpha = 0*pi/180; float beta = 0*pi/180; float gama = 0*pi/180;


    TMatrixD QPAS(3,3);
    TMatrixD  CSA(3,3);
    TMatrixD ACS(3,3);
    QPAS(0,0) = (Qeta -1.)*CQ/2.;     // Quadrupolar interaction
 //    cout<<CQ<<endl;
    QPAS(1,1) = -(1.+Qeta)*CQ/2.; 
    QPAS(2,2) = CQ; 

    cout<<"QPAS"<<endl;
 //   QPAS.Print();    
    

     CSA(0,0) = (eta -1.)*delta/2.; // 2nd-rank chemical shift anisotropy (CSA)
     CSA(1,1) = -(1.+eta)*delta/2.; 
     CSA(2,2) = delta;
      cout<<"CSA"<<endl;  
  //  CSA.Print();
     ACS(0,1)= Sxy; ACS(0,2)= Sxz;   // 1st-rank antisymmetric chemical shift (ACS)
     ACS(1,0)=-Sxy; ACS(1,2)= Syz;
     ACS(2,0)=-Sxz; ACS(2,1)=-Syz;
     cout<<"ACS"<<endl; 
  ACS.Print();
 // Adding the CSA and ACS
    TMatrixD CSPAS(3,3);
    CSPAS = CSA + ACS;
    // cout<<"CSPAS"<<endl; 
    // CSPAS.Print();
    // cout<<"CSA"<<endl; 
    // CSA.Print();
    // cout<<"ACS"<<endl; 
    // ACS.Print();

 //transformation to other frames
 /*
 CSA, ACS ------->Quad------------------>X-tal---------------->GON------------------------->LAB
          (a,b,c)      (zeta, lambda, nu)       (alpha, beta, gamma)   (phi, theta, 0)
 */


//initialization of matrices
 
 TMatrixD U(3,3);

 TMatrixD CSAQ(3,3);
 TMatrixD ACSQ(3,3);
 TMatrixD CSQ(3,3);

 TMatrixD QX(3,3);
 TMatrixD CSAQX(3,3);
 TMatrixD ACSQX(3,3);
 TMatrixD CSQX(3,3);

 TMatrixD QXG(3,3);
 TMatrixD CSAQXG(3,3);
 TMatrixD ACSQXG(3,3);
 TMatrixD CSQXG(3,3);


 U = R(a,b,c); 
 TMatrixD Uxabc = U; //Invert changes the matrix, defining a new matrix
 TMatrixD U_inv1 = Uxabc.Invert();


 CSAQ = U*CSA*U_inv1;
 ACSQ = U*ACS*U_inv1;
 CSQ = U*CSPAS*U_inv1;

U = R(zeta, lambda, nu);
TMatrixD Uxzln = U;
TMatrixD U_inv2 = Uxzln.Invert();


 QX = U*QPAS*U_inv2;
 CSAQX = U*CSAQ*U_inv2;
 ACSQX = U*ACSQ*U_inv2;
 CSQX= U*CSQ*U_inv2;

U=R(alpha,beta,gama);
TMatrixD Uxabg = U;
TMatrixD U_inv3 =Uxabg.Invert();

QXG = U*QX*U_inv3;
CSAQXG = U*CSAQX*U_inv3;
ACSQXG = U*ACSQX*U_inv3;
CSQXG = U*CSQX*U_inv3;


}

    



