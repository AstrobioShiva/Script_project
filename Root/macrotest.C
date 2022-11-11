#include <TMath.h>
#include "Riostream.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TF1.h"



void macrotest(){
  typedef ROOT::Math::SMatrix<double,3>                                       SMatrix33; 
  TMatrixD m(2,2);
  m(0,0)  = 1;
  m(1,1) = 5;

m.Print();
// m.Invert().Print();
const TMatrixD m1 = m.Invert();
m1.Print();
//const TMatrixD m2;
const TMatrixD m2(m1,TMatrixD::kMult,m);
m2.Print();
//m2.Mult(m,m1);
//m2.Print();
//m*(m.Inverse())
//m*(m.Invert()).Print();

// //  Invert a NxN matrix. The inverted matrix replace the existing one and returns if the result is successful_
// bool ret = m.Invert();  
// // return the inverse matrix of m. If the inversion fails ifail is different than zero_     
// int ifail = 0;
//  mInv = m.Inverse(ifail); 
//     for (int i = 0; i<3; i++){
//     for (int j = 0; j<3; j++){
//       cout<<ret(i,j)<<endl;
//     }
//   }
 }
