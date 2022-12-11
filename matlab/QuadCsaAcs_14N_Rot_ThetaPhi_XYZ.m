% Spectrl simulation Program of quadrupolar nuclei, CSA, and ACS
% The central transition of half-integer quadrupolr nucleus of 3/2
% It considers Quadrupolar coupling constant and CSA, ACS, and their relative
% tensor orientations 
% Theta and phi angles were considered explicitly
% %%%% Exact diagonalization & Frequency domain simulation %%%%%%%
%                         02/06/2021                           Sungsool Wi                      
function QCSAf()
tic
 global QXG Siso;
 global ct st s2t c2t cp sp c2p s2p;
 global R20Q R2p1Q R2m1Q R2p2Q R2m2Q; 
 global Siso R20cs R2p1cs R2m1cs wX;  
 global R1p1acs R1m1acs; 
 global wX CsaQXG AcsQXG CSQXG;

Nptx = 1024; freq1D=zeros(Nptx+1,1); freq2D=zeros(Nptx+1,1); freqSUM=zeros(Nptx+1,1); 
XX=zeros(Nptx+1,1); freq=zeros(3,1); freqDIFF=zeros(Nptx+1,1); 

w0 = 100;                              % Larmor frequency;
wX = w0*10^6;      % w0 (MHz), wX (The actual frequency in radian) 

ntheta =500;
nphi = 500;

  dangle = 2*pi/Nptx;

%%%%%% Select the spin quantum number 
              Ispin = 1.;               % spin quantum number (3/2, 5/2, or 7/2)             
%########### Quadrupolar coupling tensor #########################
     CQ = 3.;       % CQ (MHz)
     Qeta = 0.3;               % eta of Q
%###### Symmetric 2nd-rank chemical shift anisotropy tensor ######              
     Siso = 0.;          % isotropic chemical shift + offset (ppm) 
     delta = 200.;        % chemical shift anisotrophy (CSA) (ppm)
     eta = 0.5;           % eta of CSA
%######### Antisymmetric 1st-rank chemical shift tensor ##########
     Sxy = 100; Sxz = -100; Syz = 100;  % antisymmetric components (ppm)
     
     Siso=Siso*w0;
     delta = delta*w0;
     
     CQ = CQ*10^6/(2*Ispin*(2*Ispin-1));
     Sxy=Sxy*w0; Sxz=Sxz*w0; Syz=Syz*w0;
     
% Relative tensor orientations               
%  *************************** input parameters *************************************
%      {a,b,c)       {zeta,lamda,nu}     {alpha,beta,gama}      {phi,theta, 0}
%  CSA=========>Quad=================>X-tal===============> GON ==============> LAB FRAME
%  ACS 

a=0*pi/180; b=0*pi/180; c=0*pi/180;
zeta = 0*pi/180; lamda = 0*pi/180; nu = 0*pi/180;
alpha = 0*pi/180; beta = 0*pi/180; gama = 0*pi/180;

%**** tensor parameter assignments at their PASs **************************************            
  QPAS = zeros(3); Csa = zeros(3); Acs=zeros(3);            
  QPAS(1,1) = (Qeta -1.)*CQ/2.;     % Quadrupolar interaction
  QPAS(2,2) = -(1.+Qeta)*CQ/2.; 
  QPAS(3,3) = CQ; 
   
  Csa(1,1) = (eta -1.)*delta/2.; % 2nd-rank chemical shift anisotropy (CSA)
  Csa(2,2) = -(1.+eta)*delta/2.; 
  Csa(3,3) = delta;   
  
  Acs(1,2)= Sxy; Acs(1,3)= Sxz;   % 1st-rank antisymmetric chemical shift (ACS)
  Acs(2,1)=-Sxy; Acs(2,3)= Syz;
  Acs(3,1)=-Sxz; Acs(3,2)=-Syz;
  
    CSPAS=Csa+Acs;         % combining CSA and ACS
 
% transformations into different coordinate frames
 U=Rabc(a,b,c); CsaQ=U*Csa*U'; AcsQ=U*Acs*U'; CSQ=U*CSPAS*U'; 
 U=Rabc(zeta,lamda,nu); QX=U*QPAS*U'; CsaQX=U*CsaQ*U'; AcsQX=U*AcsQ*U'; CSQX=U*CSQ*U';
 U=Rabc(alpha,beta,gama); QXG=U*QX*U'; CsaQXG=U*CsaQX*U'; AcsQXG=U*AcsQX*U'; CSQXG=U*CSQX*U';             
    
     angle = 0.;
 for ii=1:Nptx+1
  %  theta = -angle;         % -z rotatopn
 %   phi = pi/2;
 
     theta = angle;         % y rotation
     phi = 0;

     %  theta = -angle;     % -x rotation
     %  phi = pi/2; 
     
     ct = cos(theta); st = sin(theta); s2t=2*ct*st; c2t=ct*ct-st*st;
     cp=cos(phi); sp=sin(phi); s2p=2*cp*sp; c2p=cp*cp-sp*sp;  

     Quad()
     ChemShift()
     AntiShift()
     
                      qcsa2 = -0.5*(R2m1Q*R2p1cs+R2p1Q*R2m1cs)/wX;
                      qacs2 = 0.5*(R2m1Q*R1p1acs-R2p1Q*R1m1acs)/wX;
     freq1D(ii) = Siso+sqrt(2/3)*R20cs+sqrt(1.5)*R20Q+0.5*(R2m2Q*R2p2Q-R2m1Q*R2p1Q)/wX+3*(qcsa2+qacs2); % 1 <-> 0 
     freq2D(ii) = Siso+sqrt(2/3)*R20cs-sqrt(1.5)*R20Q+0.5*(R2m2Q*R2p2Q-R2m1Q*R2p1Q)/wX-3*(qcsa2+qacs2); % 0 <-> -1  
     freqSUM(ii) = freq1D(ii)+freq2D(ii);           %1 <-> -1 transition 
     freqDIFF(ii) = freq1D(ii)-freq2D(ii);           %1 <-> -1 transition 
     XX(ii)=angle*180/pi;
     
    % Frequency position encoding (Provide a wider SW if a given SW is too narrow)
  angle = angle+dangle;
 end
 
  % Plot the spectrum after the frequency scale calibration
    figure
    plot(XX,freq1D,'b-')
    hold on
    plot(XX,freq2D,'m-')
    figure
    plot(XX,freqSUM,'r')
    figure
    plot(XX,freqDIFF,'r-')
    toc
% end of the main function      

  function U= Rabc(alfa1,beta1,gama1)
  U = zeros(3);
  U(1,1) = cos(alfa1)*cos(beta1)*cos(gama1) - sin(alfa1)*sin(gama1);
  U(1,2) = sin(alfa1)*cos(beta1)*cos(gama1) + cos(alfa1)*sin(gama1);
  U(1,3) = -sin(beta1)*cos(gama1);
  U(2,1) = -cos(alfa1)*cos(beta1)*sin(gama1) - sin(alfa1)*cos(gama1);
  U(2,2) = -sin(alfa1)*cos(beta1)*sin(gama1) + cos(alfa1)*cos(gama1);
  U(2,3) = sin(beta1)*sin(gama1);
  U(3,1) = cos(alfa1)*sin(beta1);
  U(3,2) = sin(alfa1)*sin(beta1);
  U(3,3) = cos(beta1);
  return;

 function Quad()
 global QXG;
 global R20Q R2p1Q R2m1Q R2p2Q R2m2Q; 
 global ct st s2t c2t cp sp c2p s2p;
 sq15=sqrt(1.5); sq6=sqrt(6);
 R20Q=0.5*sq15*(QXG(3,3)*(3*ct*ct-1)+2*QXG(1,3)*s2t*cp+2*QXG(2,3)*s2t*sp+(QXG(1,1)-QXG(2,2))*st*st*c2p+2*QXG(1,2)*st*st*s2p);
% 
 R2q1=0.75*QXG(3,3)*s2t;
 R2q2r=QXG(1,3)*c2t*cp; R2q2i=QXG(1,3)*ct*sp;
 R2q3r=QXG(2,3)*c2t*sp; R2q3i=QXG(2,3)*ct*cp;
 R2q4r=0.25*(QXG(1,1)-QXG(2,2))*s2t*c2p; R2q4i=0.5*(QXG(1,1)-QXG(2,2))*st*s2p;
 R2q5r=0.5*QXG(1,2)*s2t*s2p; R2q5i=QXG(1,2)*st*c2p;
 R2p1Q=R2q1-R2q2r+i*R2q2i-R2q3r-i*R2q3i-R2q4r+i*R2q4i-R2q5r-i*R2q5i;
 R2m1Q=-R2q1+R2q2r+i*R2q2i+R2q3r-i*R2q3i+R2q4r+i*R2q4i+R2q5r-i*R2q5i;
%
 R4q1=0.75*QXG(3,3)*st*st;
 R4q2r=0.5*QXG(1,3)*s2t*cp; R4q2i=QXG(1,3)*st*sp;
 R4q3r=0.5*QXG(2,3)*s2t*sp; R4q3i=QXG(2,3)*st*cp;
 R4q4r=0.25*(QXG(1,1)-QXG(2,2))*(1+ct*ct)*c2p; R4q4i=0.5*(QXG(1,1)-QXG(2,2))*ct*s2p;
 R4q5r=0.5*QXG(1,2)*(1+ct*ct)*s2p; R4q5i=QXG(1,2)*ct*c2p;
 R2p2Q=R4q1-R4q2r+i*R4q2i-R4q3r-i*R4q3i+R4q4r-i*R4q4i+R4q5r+i*R4q5i;
 R2m2Q=R4q1-R4q2r-i*R4q2i-R4q3r+i*R4q3i+R4q4r+i*R4q4i+R4q5r-i*R4q5i;
return;


function ChemShift()
 global CsaQXG;
 global R20cs R2p1cs R2m1cs; 
 global ct st s2t c2t cp sp c2p s2p;
sq15=sqrt(1.5);
R20cs=0.5*sq15*((2*CsaQXG(3,3)-CsaQXG(1,1)-CsaQXG(2,2))*(ct*ct-1/3)+2*CsaQXG(1,3)*s2t*cp+2*CsaQXG(2,3)*s2t*sp+(CsaQXG(1,1)-CsaQXG(2,2))*st*st*c2p+2*CsaQXG(1,2)*st*st*s2p);
%
R2cs0=0.25*(2*CsaQXG(3,3)-CsaQXG(1,1)-CsaQXG(2,2))*s2t;
R2cs1r=CsaQXG(1,3)*c2t*cp; R2cs1i=CsaQXG(1,3)*ct*sp;
R2cs2r=CsaQXG(2,3)*c2t*sp; R2cs2i=CsaQXG(2,3)*ct*cp;
R2cs3r=0.25*(CsaQXG(1,1)-CsaQXG(2,2))*s2t*c2p; R2cs3i=0.5*(CsaQXG(1,1)-CsaQXG(2,2))*st*s2p;
R2cs4r=0.5*CsaQXG(1,2)*s2t*s2p; R2cs4i=CsaQXG(1,2)*st*c2p;
R2p1cs=R2cs0-R2cs1r+i*R2cs1i-R2cs2r-i*R2cs2i-R2cs3r+i*R2cs3i-R2cs4r-i*R2cs4i;
R2m1cs=-R2cs0+R2cs1r+i*R2cs1i+R2cs2r-i*R2cs2i+R2cs3r+i*R2cs3i+R2cs4r-i*R2cs4i;
return;


function  AntiShift()
 global AcsQXG;
 global R1p1acs R1m1acs; 
 global ct st cp sp;
 R11=0.5*(AcsQXG(1,2)-AcsQXG(2,1))*st;
 R12=0.5*(AcsQXG(1,3)-AcsQXG(3,1))*cp; R13=0.5*(AcsQXG(1,3)-AcsQXG(3,1))*ct*sp;
 R14=0.5*(AcsQXG(2,3)-AcsQXG(3,2))*sp; R15=0.5*(AcsQXG(2,3)-AcsQXG(3,2))*ct*cp;
 R1p1acs=-i*R11+R12-i*R13+R14+i*R15;
 R1m1acs=i*R11+R12+i*R13+R14-i*R15;
 return;


