
void Zcen2quad(DVector zx, DVector zy, DVector zz, double sfrq, double ival,

               double *cq, double *etaq, DMatrix rq, DVector zxq, DVector zyq,

               DVector zzq)

{

  int     i, j;

  IVector sort;

  double  k, a1m2, a2m3;

  DVector qpas, qgon13, qgon23, qgon11, qgon22, qgon33, a1m3, test;

  DMatrix qgon, temp;

 

  a1m3   = DV_alloc(2);

  test   = DV_alloc(2);

  qgon13 = DV_alloc(2);

  qgon23 = DV_alloc(2);

  qgon11 = DV_alloc(2);

  qgon22 = DV_alloc(2);

  qgon33 = DV_alloc(2);

  qpas   = DV_alloc(3);

  qgon   = DM_alloc(3,3);

  temp   = DM_alloc(3,3);

  sort   = IV_alloc(3);

 

  if (fabs(zx[4]) < -1e5 && fabs(zy[4]) < -1e5 && fabs(zz[4]) < -1e5 &&

      fabs(zx[5]) < -1e5 && fabs(zy[5]) < -1e5 && fabs(zz[5]) < -1e5) {

    printf("Not in the mood for determining quad and CSA from the central transition\n");

    return;

  }

 

  k = 16.0*sfrq/(9.0*(4.0*ival*(ival+1.0) - 3.0))/1.0E3;  /* [qgon] = MHz */

 

  qgon[1][2] = sqrt(k*(sqrt(zz[4]*zz[4] + zz[5]*zz[5]) - zz[4]));

  a1m2 = sqrt(4.0*k*(sqrt(zz[4]*zz[4] + zz[5]*zz[5]) + zz[4]));

  qgon13[1] = sqrt(k*(sqrt(zy[4]*zy[4] + zy[5]*zy[5]) - zy[4]));

  qgon13[2] =-sqrt(k*(sqrt(zy[4]*zy[4] + zy[5]*zy[5]) - zy[4]));

  a2m3 = sqrt(4.0*k*(sqrt(zx[4]*zx[4] + zx[5]*zx[5]) + zx[4]));

  for (i = 1; i <= 2; i++) {

    qgon23[i] = sqrt(k*(sqrt(zx[4]*zx[4] + zx[5]*zx[5]) - zx[4]));

    a1m3[i] = sqrt(4.0*k*(sqrt(zy[4]*zy[4] + zy[5]*zy[5]) + zy[4]));

  }

 

   

  /* Bestemmelse af fortegn. Vaelger qgon[1][2] > 0 */

  /* Eng: Determination of the sign. Choosing qgon[1][2] > 0 */

  a1m2 *= zz[5] >= 0.0 ? 1.0:-1.0;

 

  /* Proever med qgon[1][3] hhv  > 0 og < 0 */

  /* Eng: Trying with qgon[1][3] respectively >0 and <0 */

  for (i = 1; i <= 2; i++) {

    a1m3[i] *= zy[5]*qgon13[i] >= 0.0 ? 1.0:-1.0;

    test[i] = zx[5]/(zy[5]/qgon13[i] - zz[5]/qgon[1][2]);

    qgon23[i] *= test[i] < 0.0 ? -1.0:1.0;

 

    qgon11[i] = (a1m2 + a1m3[i])/3.0;

    qgon22[i] = (a1m3[i] - 2.0*a1m2)/3.0;

    qgon33[i] = (a1m2 - 2.0*a1m3[i])/3.0;

    test[i] = fabs(fabs(qgon22[i] - qgon33[i])/a2m3 - 1.0);

  }

 

  i = test[2] < test[1] ? 2:1;

   

  qgon[1][1] = qgon11[i];

  qgon[1][3] = qgon13[i];

  qgon[2][2] = qgon22[i];

  qgon[2][3] = qgon23[i];

  qgon[3][3] = qgon33[i];

 

  /* Fortegn OK. Optimering */

  /* Eng: Sign Okay. Optimizing */

 

  opt_cen(qgon, sfrq, ival, zx, zy, zz);

 

  /* Nu til diqgononaliseringen */

  /* Eng: Now for the diagonalization */

 

  DM_cp(temp, qgon);

 

  DJacobi(temp, qpas, rq);

  Prsort(qpas, sort);

 

  if (qpas[3] < 0.0) {

    for (i = 1; i <= 3; i++) {

      for (j = 1; j <= 3; j++) {

        temp[i][j] = - qgon[i][j];

        qgon[i][j] = - qgon[i][j];

      }

    }

    DJacobi(temp, qpas, rq);

    Prsort(qpas, sort);

  }

 

  Matsort(sort, rq);

  Righthand(rq);

  Def_interval(rq);

 

  *cq = 2.0*ival*(2.0*ival-1.0)*qpas[3];

  *etaq = (qpas[2]-qpas[1])/qpas[3];

 

  /* Beregning af de forskellige Z'er */

  /* Eng: Calculating the different Z's */

 

  Quad2zcen(qgon, ival, sfrq, 0.0, zxq, zyq, zzq);

 

  DV_free(a1m3);

  DV_free(test);

  DV_free(qgon13);

  DV_free(qgon23);

  DV_free(qgon11);

  DV_free(qgon22);

  DV_free(qgon33);

  DV_free(qpas);

  DM_free(qgon); 

  IV_free(sort);

}