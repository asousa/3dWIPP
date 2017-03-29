#include <wipp.h>

// Misc. methods from the original WIPP code, with any 
// c -> c++ modifications as needed.


/*
 * FUNCTION: ionoAbsorp
 * --------------------
 * This function just calculates the ionospheric absorption of a 
 * ray in moving through the ionosphere from ~100km to ~1000km.
 * Uses data from Helliwell 1965, fig. 3-35, and interpolates for 
 * latitudes between 0-90 degrees, then uses data at 2 kHz, and 
 * 20 kHz to interpolate/extrapolate for other frequencies.
 *
 */

double ionoAbsorp(float lat, long f)
{
  float db2i, db20i;
  double  db2iLog, db20iLog, m, c; 

  float lats[46] = {  0.0,  2.0,  4.0,  6.0,  8.0,
      10.0, 12.0, 14.0, 16.0, 18.0,
      20.0, 22.0, 24.0, 26.0, 28.0,
      30.0, 32.0, 34.0, 36.0, 38.0, 
      40.0, 42.0, 44.0, 46.0, 48.0, 
      50.0, 52.0, 54.0, 56.0, 58.0,
      60.0, 62.0, 64.0, 66.0, 68.0,
      70.0, 72.0, 74.0, 76.0, 78.0, 
      80.0, 82.0, 84.0, 86.0, 88.0, 
      90.0  };

  float db20kHzLog[46] = {
    2.72890330135884,
    2.61413531724660,
    2.49936733313435,
    2.38459934902211,
    2.26983136490987,
    2.15506338079763,
    2.04029539668538,
    1.92531644206607,
    1.79924052160795,
    1.65746834723986,
    1.49802217224091,
    1.36186709269283,
    1.25598388505325,
    1.15652386401487,
    1.07820410788996,
    1.01300984614119,
    0.95644081081742,
    0.90646976883951,
    0.85864979711507,
    0.81518987704864,
    0.77367088822273,
    0.73407595275903,
    0.70105485532855,
    0.67236287394032,
    0.64568469779523,
    0.62025316330165,
    0.59707886563282,
    0.57435687349111,
    0.55630869316853,
    0.53826051284596,
    0.52280663628542,
    0.50796595891253,
    0.49540084969300,
    0.48392405611603,
    0.47327480980675,
    0.46355655505289,
    0.45432489421971,
    0.44715189958245,
    0.43997890494519,
    0.43289667529356,
    0.42589049202632,
    0.41888430875909,
    0.41279384916331,
    0.40818263512987,
    0.40428972648146,
    0.40189872848716  };

  float db2kHzLog[46] = {
    2.20354296279991,
    2.10139102288977,
    1.99923908297962,
    1.89708714306948,
    1.79493520315934,
    1.69278326324919,
    1.59070536198610,
    1.47122984416933,
    1.32946932452408,
    1.19113861945499,
    1.04312210740998,
    0.92225656425647,
    0.81909999940519,
    0.73129314437689,
    0.65802549162668,
    0.59490689076827,
    0.54352690196227,
    0.49737249399863,
    0.45587327008181,
    0.41437404616500,
    0.38243818008768,
    0.35131375279102,
    0.32279494989472,
    0.29651211576080,
    0.27148893501543,
    0.24935600431553,
    0.22747724000215,
    0.20903314606947,
    0.19058905213679,
    0.17380919954012,
    0.15966172931239,
    0.14551425908466,
    0.13161845036475,
    0.11976152945062,
    0.10790460853649,
    0.09604768762236,
    0.08688693663841,
    0.07777734644718,
    0.06866775625595,
    0.05958751374305,
    0.05050956953555,
    0.04345615362238,
    0.03810140109787,
    0.03270477356885,
    0.02717155386097,
    0.02163833415309   };
    


  db2iLog  = interpPt(lats, db2kHzLog  , 46, lat);
  db20iLog = interpPt(lats, db20kHzLog , 46, lat);
  
  m = (db20iLog - db2iLog); // officially should be 
        // m'=m/( log10(20)-log10(2) )
        // but denominator = 1

  c = (db2iLog+db20iLog)/2 - m * 0.8010; // this is just
           // (log10(2)+log10(20))/2

  // now extrapolate to desired f value in log10 domain
  // get 10^result, and return it
  return ( pow(10.0 , ( m*log10( (f/1000.0) ) + c ) ) );
}




/*
 * FUNCTION: interpPt
 * ------------------
 * This function is designed to work pretty much the same as matlab's
 * interp1 function.  Specify the input x and y vectors, and the output
 * x point.
 *
 */

float interpPt(float *xI, float *yI, int n, float xO)
{
  int i, iHigh, iLow, iMid;
  float yO;

  // Check that xO is within bounds
  if( (xO < xI[0]) || xO > xI[n-1] ) {
    printf("\nPoint is out of bounds! %g, {%g, %g}\a\n",xO, xI[0], xI[n-1]);
    return 0.0;
  }
  
  // Do a binary search for the correct index 
  iHigh = n-1;
  iLow = 0;  
  while(1) {
    iMid = (iHigh+iLow)/2;
    if( (xO >= xI[iLow]) && (xO < xI[iMid]) ) {
      iHigh = iMid;
    } else {
      iLow = iMid;
    }
    if(xO==xI[n-1]){printf("\nin interpPt\n"); return(yI[n-1]);}
    if(iHigh==iLow) {
      printf("\nexiting from binary search in 1st condtion\n");
      break;
    }
    if( (xO>=xI[iMid]) && (xO<xI[iMid+1]) ) break;
    if( (xO>=xI[iMid-1]) && (xO<xI[iMid]) ) {
      iMid--;
      break;
    }
  }     
  yO = (  ( yI[iMid+1]-yI[iMid] ) / ( xI[iMid+1]-xI[iMid] )  )
    *( xO-xI[iMid] ) + yI[iMid];
  
  return(yO);
}



/*
 * FUNCTION: Fresnel
 * -----------------
 * This function calculates the sine and cose Fresnel integrals and 
 * returns the values of either.  The integrals are defined as:
 *
 *    /x           /x
 *    | sin( pi/2*t^2 ) dt           | cos( pi/2*t^2 ) dt
 *    /0           /0
 *
 * This uses the methodology of  Klaus D. Mielenz, "Computation of
 * Fresnel Integrals II" Journal of Research of the National 
 * Institute of Standards and Technology, 105, 589 (2000)
 *
 */
void Fresnel(double x0, double *FS, double *FC)
{

  double fn[12] = { 0.318309844,
        9.34626e-8 , 
        -0.09676631, 
        0.000606222, 
        0.325539361, 
        0.325206461, 
        -7.450551455,
        32.20380908, 
        -78.8035274, 
        118.5343352, 
        -102.4339798,
        39.06207702 } ;

  double gn[12] = { 0.0 ,
        0.101321519, 
        -4.07292e-5, 
        -0.152068115, 
        -0.046292605, 
        1.622793598, 
        -5.199186089, 
        7.477942354, 
        -0.695291507, 
        -15.10996796,
        22.28401942, 
        -10.89968491 };
  
  double cn[12] = { 1.0 ,
        -0.24674011002723,
        0.02818550087789 ,
        -0.00160488313564 ,
        5.407413381408390e-05 ,
        -1.200097255860028e-06,
        1.884349911527268e-08,
        -2.202276925445466e-10, 
        1.989685792418021e-12,
        -1.430918973171519e-14,
        8.384729705118549e-17,
        -4.079981449233875e-19 } ;

  double sn[12] = {    0.52359877559830,
           -0.09228058535804,
           0.00724478420420,
           -3.121169423545791e-04,
           8.444272883545251e-06,
           -1.564714450092211e-07,
           2.108212193321454e-09,
           -2.157430680584343e-11,
           1.733410208887483e-13,
           -1.122324478798395e-15,
           5.980053239210401e-18,
           -2.667871362841397e-20 };

  double xpow, x_sq, fx=0, gx=0, x;
  int n;

  x = fabs(x0);
  *FS = 0.0;
  *FC = 0.0;
  x_sq = x*x;



  if(x<=1.6) {

    *FS =   sn[0]*pow(x,3) +  // it takes longer to write this out
      sn[1]*pow(x,7) +        // but we save valuable CPU cycles!
      sn[2]*pow(x,11) + 
      sn[3]*pow(x,15) + 
      sn[4]*pow(x,19) + 
      sn[5]*pow(x,23) + 
      sn[6]*pow(x,27) + 
      sn[7]*pow(x,31) + 
      sn[8]*pow(x,35) + 
      sn[9]*pow(x,39) + 
      sn[10]*pow(x,43) + 
      sn[11]*pow(x,47)  ; 

    *FC =   cn[0]*x + 
      cn[1]*pow(x,5) +  
      cn[2]*pow(x,9) + 
      cn[3]*pow(x,13) + 
      cn[4]*pow(x,17) + 
      cn[5]*pow(x,21) + 
      cn[6]*pow(x,25) + 
      cn[7]*pow(x,29) + 
      cn[8]*pow(x,33) + 
      cn[9]*pow(x,37) + 
      cn[10]*pow(x,41) + 
      cn[11]*pow(x,45)  ; 

  } else {
      
    
    for(n=0; n<=11; n++) {
      xpow = pow(x, (-2*n-1) );
      fx += fn[n]*xpow;
      gx += gn[n]*xpow;
    }     
    *FC = 0.5 + fx*sin(PI/2*x_sq) - gx*cos(PI/2*x_sq);
    *FS = 0.5 - gx*sin(PI/2*x_sq) - fx*cos(PI/2*x_sq);      
  }
 
  if(x0<0) {
    *FC = -(*FC);
    *FS = -(*FS);
  }
}


/*
 * FUNCTION: calcRes
 * -----------------
 * This function will calculate the resonant pitch angle change of the 
 * particle using the lat_arr data  
 *
 */
// void calcRes(cellT *lat_arr[][NUM_TARGS], long lower_freq)
void calcRes(cellT cell, double da_N[NUM_E][NUM_TIMES], double da_S[NUM_E][NUM_TIMES])
{
  FILE *resPtrN, *resPtrS;
  int i, j=0, kk, mres, noutN, noutS, ei, ti, e_toti;
  cellT *next;
  // char *prefN="pN", *prefS="pS", suff[64];
  double lat, L, t, f, pwr, psi, mu, stixP, stixR, stixL, latk;
  double Bxw, Byw, Bzw, Exw, Eyw, Ezw, stixD, stixS, stixA;
  double stixB, stixX, n_x, n_z, k, kx, kz, rho1, rho2, Byw_sq;
  double flt_const_N, flt_const_S, flt_time, eta_dot;

  // double flt_const_N[EA_SPLIT], flt_const_S[EA_SPLIT], flt_time, eta_dot;
  double wh, dwh_ds, gamma, alpha1, alpha2, beta, v_para, v_perp;
  double spsi, cpsi, spsi_sq, cpsi_sq, mu_sq, w, R1, R2, w1, w2;
  double alpha_lc, alpha_eq, epsm, slat, clat, slat_term, ds;
  double t1, t2, t3, direction, v_para_res, v_tot, v_tot_res;
  double salph, calph, wtau_sq, Y, dv_para_ds, AA, BB, T1;
  double Farg, Farg0, Fs, Fc, Fs0, Fc0, dFs_sq, dFc_sq, dalpha;
  double alpha_eq_p, dalpha_eq, E_res, e_starti, e_endi;
  double v_tot_arr[NUM_E], E_tot_arr[NUM_E], v_para_star, v_para_star_sq;
  float *arr_N, *arr_S;
  time_t start, end;
  int timei;

  double temp;
 
  // // start = time(NULL);
  // time(&start);
  // arr_N = NULL;
  // arr_S = NULL;
  
  L = cell.Lsh;
  epsm = (1/L)*(R_E+H_IONO)/R_E;
  alpha_eq = asin(sqrt( pow(epsm,3)/sqrt(1+3*(1-epsm)) ));
  // sprintf(suff, "_%g", L);
  // resPtrN = writeFile(lower_freq, prefN, suff);
  // resPtrS = writeFile(lower_freq, prefS, suff);
  // printf("\nNow doing L: %g\n", L);
  
  //initialize the velocity and energy arrays
  for(i=0; i<NUM_E; i++) {
    E_tot_arr[i] = pow(10, (E_EXP_BOT+ DE_EXP*i) ); // energy in eV
    v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
    // printf("%e ",v_tot_arr[i]);
  }

  // //initialize the Alpha output arrays (prevents having empty pN, pS files)
  // addToAlphaArr(&arr_N, &arr_S, 0.0, 
  //       0.0, 0.0);




  // // Go through all latitudes
  // for(i=0; i<NUM_EA; i++) {
  //   next = lat_arr[i][j];
    
    // CALCULATE:
    // wh
    // dwh/ds
    // flight-time constant
    // alpha_lc
    
    lat = cell.lat;
    slat = sin( lat*D2R );
    clat = cos( lat*D2R );
    slat_term = sqrt(1+3*slat*slat);
    wh = 2*PI*880000/pow(L,3)*slat_term/pow(clat,6);
    dwh_ds = 3*wh/(L*R_E)*slat/slat_term*
      (1/(slat_term*slat_term) + 2/(clat*clat));




    // for(kk=0; kk < EA_SPLIT; kk++) {
    //   latk = (lat-.5*EAIncr+.5*EAIncr/EA_SPLIT)+kk*(EAIncr/EA_SPLIT);
    //   getFltConst(L,latk,alpha_eq,&(flt_const_N[kk]),&(flt_const_S[kk]));
    // }
    latk = (lat-.5*EAIncr+.5*EAIncr)+kk*(EAIncr);
    getFltConst(L,latk,alpha_eq,&(flt_const_N),&(flt_const_S));

    alpha_lc = asin(sqrt( slat_term/pow(clat,6) )*sin(alpha_eq));
    salph = sin(alpha_lc);
    calph = cos(alpha_lc);
    // ds = L*R_E* slat_term*clat*EAIncr/180*PI;    
    // ds = L*R_E* slat_term*clat*EAIncr/180*PI * MULT ; 
    ds = L*R_E* slat_term*clat*EAIncr*D2R; 

    dv_para_ds = -0.5*pow(salph,2)/calph/wh*dwh_ds;
    
    printf("------------- EA at lat: %2.2f ---------------\n", lat);
    //printf("lat: %g, wh: %g, dwh_ds: %g, flt_const_N: %g \n",
    //     lat, wh, dwh_ds, flt_const_N);
    //printf("flt_const_S: %g, alpha_lc: %g, \nds: %g, dv_para_ds: %g\n", 
    //     flt_const_S, alpha_lc, ds, dv_para_ds);
        printf("wh: %g dwh_ds: %g\n",wh, dwh_ds);

    
    // Go through the cells in each latitude
    // while(next != NULL) {
      
      t = cell.t + TIME_STEP/2;   // We want the time and freq to be in the 
      f = cell.f + FREQ_STEP_SIZE/2;   // center of the cell, so add DT/2 or DF/2
      pwr = cell.pwr*FREQ_STEP_SIZE/cell.num_rays;
      psi = D2R*cell.psi/cell.num_rays;

      // pwr = (next->pwr)/(next->num_rays)/(T_STEP); // <-I know this looks wrong 
      // psi = (next->psi)/(next->num_rays)*D2R;      // but it's right!
                                                   // Grrr... back off!
      // printf("\npwr: %e\n",pwr);
      // printf("t: %g psi: %g pwr: %g num_rays: %i\n",t, R2D*psi, pwr, next->num_rays);
      // printf("num_rays: %d\n",next->num_rays);
      // printf("DT: %2.3f\n",DT);
      if(pwr > 1.0e-50) {
      // if(pwr > 0.0) {


      mu = cell.mu/cell.num_rays;
      stixP = cell.stixP/cell.num_rays;
      stixR = cell.stixR/cell.num_rays;
      stixL = cell.stixL/cell.num_rays;

      // mu  = (next->mu)/(next->num_rays);
      // stixP =   (next->stixP)/(next->num_rays);
      // stixR = (next->stixR)/(next->num_rays);
      // stixL =   (next->stixL)/(next->num_rays);
      
      // printf("t: %g, f: %g, pwr: %g, psi: %g,\nmu: %g, stixP: %g, stixR: %g, stixL: %g \n\n", t, f, pwr, psi, mu, stixP, stixR, stixL);
      spsi = sin(psi);
      cpsi = cos(psi);
      spsi_sq = pow(spsi,2);
      cpsi_sq = pow(cpsi,2);
      n_x = mu*fabs(spsi);
      n_z = mu*cpsi;
      mu_sq = mu*mu;
      w = 2.0*PI*f;
      k = w*mu/C;
      kx = w*n_x/C;
      kz = w*n_z/C;
      Y = wh / w ;
      stixS = ( stixR + stixL ) /2.0;
      stixD = ( stixR - stixL ) /2.0;
      stixA = stixS + (stixP-stixS)*cpsi_sq;
      stixB = stixP*stixS+stixR*stixL+(stixP*stixS-stixR*stixL)*cpsi_sq;
      stixX = stixP/(stixP- mu_sq*spsi_sq);
      
      // printf("t: %g f: %g\n",t, w/(2*PI));

      printf("t: %g, f: %g, pwr: %g, psi: %g, Num_rays: %g\n",t,f,pwr, R2D*psi, cell.num_rays);
    printf("wh: %g dwh_ds: %g alpha_lc: %g alpha_eq: %g ds: %g dv_para_ds: %g\n",
            wh,    dwh_ds,    alpha_lc,    alpha_eq,    ds,    dv_para_ds);   
      rho1=((mu_sq-stixS)*mu_sq*spsi*cpsi)/(stixD*(mu_sq*spsi_sq-stixP));
      rho2 = (mu_sq - stixS) / stixD ;
      
      // OLD - INCORRECT WAY!
      //Byw_sq = 2.0*MU0/C* pwr *stixX*stixX *mu* fabs(cpsi)*
      // sqrt( pow( ((spsi/cpsi)-rho1*rho2) ,2) + 
      //     pow( (1+rho2*rho2*stixX)  ,2) );
      Byw_sq =  2.0*MU0/C*pwr*stixX*stixX*rho2*rho2*mu*fabs(cpsi)/
           sqrt(  pow((tan(psi)-rho1*rho2*stixX),2) + 
           pow( (1+rho2*rho2*stixX), 2 ) );


    // printf("\nt: %2.3f, f: %g, pwr: %g, Byw_sq: %e\n",t,f,pwr,Byw_sq);

 
      // get all RMS wave components
      
      Byw = sqrt(Byw_sq);
      Exw = fabs(C*Byw * (stixP - n_x*n_x)/(stixP*n_z)); 
      Eyw = fabs(Exw * stixD/(stixS-mu_sq));
      Ezw = fabs(Exw *n_x*n_z / (n_x*n_x - stixP));
      Bxw = fabs(Exw *stixD*n_z /C/ (stixS - mu_sq));
      Bzw = fabs((Exw *stixD *n_x) /(C*(stixX - mu_sq)));
      // printf("\nn_x: %g, n_z: %g, stixX: %g, stixD: %g, stixA: %g, stixB: %g, mu_sq: %g\n", n_x, n_z, stixX, stixD, stixA, stixB, mu_sq);
    printf("Byw: %g Exw: %g Eyw: %g Ezw: %g Bxw: %g Bzw: %g\n",
              Byw,    Exw,    Eyw,    Ezw,    Bxw,    Bzw);

      // Oblique integration quantities
      R1 = (Exw + Eyw)/(Bxw+Byw);
      R2 = (Exw - Eyw)/(Bxw-Byw);
      w1 = Q_EL/(2*M_EL)*(Bxw+Byw);
      w2 = Q_EL/(2*M_EL)*(Bxw-Byw);
      alpha1 = w2/w1;
          printf("R1: %g R2: %g w1: %g w2: %g alpha1: %g\n",
            R1,    R2,    w1,    w2,    alpha1);
      // printf("R1: %g, R2: %g, w1: %g, w2: %g\n", R1, R2, w1, w2);
      // printf("alpha_eq: %g alpha_lc: %g\n",R2D*alpha_eq, R2D*alpha_lc);
      //begin MRES loop here
      for(mres=-SCATTERING_RES_MODES; mres <= SCATTERING_RES_MODES; mres++) {
        
        // get parallel resonance velocity
        t1 = w*w*kz*kz;
        t2 = pow((mres*wh),2)-w*w;
        t3 = kz*kz + pow((mres*wh),2)/(pow(C*cos(alpha_lc),2));
        if(mres==0) {
          direction = -kz/fabs(kz);
        } else {
          direction = kz/fabs(kz) * mres/fabs(mres) ;
        }
        v_para_res = ( direction*sqrt(t1 + t2*t3) - w*kz ) / t3;
        v_tot_res = v_para_res / cos(alpha_lc); 
        E_res = E_EL*( 1.0/sqrt( 1.0-(v_tot_res*v_tot_res/(C*C)) ) -1.0 );

        printf("t1: %g t2: %g t3: %g v_para_res: %g v_tot_res: %g E_res: %g\n",
                t1,    t2,    t3,    v_para_res,    v_tot_res,    E_res);
        // get starting and ending indices, +-20% energy band
        // temp = floor((log10(E_res) - E_EXP_BOT - 0.3)/(DE_EXP));

        // printf("temp: %g\n", temp);
        e_starti = floor((log10(E_res) - E_EXP_BOT - 0.3)/(DE_EXP));
        e_endi   =  ceil((log10(E_res) - E_EXP_BOT + 0.3)/(DE_EXP));
        // printf("e_res: %g, e_starti: %g, e_endi: %g\n",E_res, e_starti, e_endi);

        if(e_endi>NUM_E) e_endi=NUM_E;
        if(e_starti>NUM_E) e_starti=NUM_E;
        if(e_endi<0) e_endi=0;
        if(e_starti<0) e_starti=0;
        
        // printf("dir: %g t1: %g t2: %g t3: %g\n",direction, t1, t2, t3);

        // begin V_TOT loop here
        for(e_toti=e_starti; e_toti <= e_endi; e_toti++) {

          // printf("E_EXP_BOT: %g, E_EXP_TOP: %g, DE_EXP: %g\n",E_EXP_BOT,E_EXP_TOP,DE_EXP);
          // printf("e_res: %g, e_starti: %g, e_endi: %g\n",E_res, e_starti, e_endi);
       
          v_tot = direction*v_tot_arr[e_toti];
          v_para = v_tot * calph;
          v_perp = fabs(v_tot * salph);
          
          gamma = 1.0 / sqrt(1 - pow((v_tot/C),2)); 
          alpha2 = Q_EL*Ezw /(M_EL*gamma*w1*v_perp);
          beta = kx*v_perp / wh ;
          wtau_sq = pow((-1),(mres-1)) * w1/gamma * 
            ( jn( (mres-1), beta ) - 
              alpha1*jn( (mres+1) , beta ) +
              gamma*alpha2*jn( mres , beta ) ); 
          T1 = -wtau_sq*(1+ ( (calph*calph) / (mres*Y-1) ) );
          
          // Now - start analytical evaluation!!!
          
          if( fabs(lat)< 1e-3) {
            
            eta_dot = mres*wh/gamma - w - kz*v_para;

            if(fabs(eta_dot)<10) {
              dalpha_eq = fabs(T1/v_para)*ds/sqrt(2); 
            } else {
              dalpha_eq = fabs(T1/eta_dot)*sqrt(1-cos(ds*eta_dot/v_para)); 
            }

          } else {  
            
            // AA = (mres/(2.0*v_para*gamma))*dwh_ds - (kz/2.0)*dv_para_ds;
            // doesn't the second term need  to get divided by v_para?


            //BB = mres/v_para/gamma*(wh-dwh_ds*ds/2) - 
            //  w/v_para - kz + kz/v_para*dv_para_ds*ds/2;

            v_para_star = v_para - dv_para_ds*ds/2.0;
            v_para_star_sq = v_para_star * v_para_star;

            AA = (mres/(2.0*v_para_star*gamma))*dwh_ds* 
              (1 + ds/(2.0*v_para_star)*dv_para_ds) - 
              mres/(2.0*v_para_star_sq*gamma)*wh*dv_para_ds + 
              w/(2.0*v_para_star_sq)*dv_para_ds ;


            BB = mres/(gamma*v_para_star)*wh - 
              mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) -
              w/v_para_star - kz;

            // // Bortnik A.18 -- part A0
            // BB =   mres*wh/(gamma*v_para_star)
            //      - mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) * (w/v_para_star)*kz;

            //      printf("\nAA: %g, BB: %g, v_para_star: %g\n", 
            //     AA, BB, v_para_star);

            // MATLAB code
            //aa1 = mres/2./v_para_star/gammaiA.*dwhi_ds.* ...
            //  (1+dsiA./2./v_para_star.*dv_para_ds) - ...
            //  mres/2./(v_para_star.^2)./gammaiA.*whiA.*dv_para_ds + ...
            //  wi./2./(v_para_star.^2).*dv_para_ds ;
      
            //  bb1 = mres./gammaiA./v_para_star.*whiA - ...
            // mres./gammaiA./v_para_star.*dwhi_ds.*dsiA/2 - ...
            // wi./v_para_star - kzi;


            
            Farg = (BB + 2*AA*ds) / sqrt(2*PI*fabs(AA));
            Farg0 = BB / sqrt(2*PI*fabs(AA));  
            
            Fresnel(Farg, &Fs, &Fc);
            Fresnel(Farg0, &Fs0, &Fc0);
            
            dFs_sq = pow((Fs - Fs0),2);
            dFc_sq = pow((Fc - Fc0),2);
            
            dalpha = sqrt(PI/4/fabs(AA))*fabs(T1/v_para)*sqrt(dFs_sq+dFc_sq);
            
            alpha_eq_p = asin( sin(alpha_lc+dalpha)*pow(clat,3) / 
                       sqrt(slat_term) );
            dalpha_eq = alpha_eq_p - alpha_eq;

          // printf("dalpha_eq: %e\n",dalpha_eq);
            
            if (isnan(dalpha_eq)) {
              printf("isnans: w1: %g\n",w1);
            }
          }
          // printf("dalpha_eq: %g\n",R2D*dalpha_eq);
          if(direction>0) {
            flt_time = fabs(flt_const_N/v_para);
          } else {
            flt_time = fabs(flt_const_S/v_para);
          }
         
         printf("flt_time: %g dalpha_eq: %g alpha_eq_p: %g alpha_eq: %g\n",
                        flt_time,    dalpha_eq,    alpha_eq_p,    alpha_eq);
          // Get time index into output array
          timei = floor((t + flt_time)/TIME_STEP);
          // Save it!
          if (direction > 0) {
              da_N[e_toti][timei] += dalpha_eq*dalpha_eq;
              // cout << timei << " " << da_N[e_toti][timei] << "\n";

          } else {
              da_S[e_toti][timei] += dalpha_eq*dalpha_eq;
              // cout << timei << " " << da_S[e_toti][timei] << "\n";
          }

        } // v_para
      } // mres loop
    } // if pwr > 0.0
}





/*
 * FUNCTION: getFltConst
 * ---------------------
 * This function calculates the 'flight-time-const' for a particular 
 * latitude,  and returns the constant for a particle flying adiabatically
 * from that latitude, to the Northern hemisphere, and a similar constant 
 * to the Southern hemisphere.  This constant is then multiplied by 
 * 1/v_tot for the particular particle, which gives the total flight
 * time in seconds to the appropriate hemisphere.
 *
 */
void getFltConst(double L, double lat, double alpha_eq, 
         double *flt_const_N, double *flt_const_S)
{
  double sin_alpha_eq_sq, x=0.0, dx, endx, walt_tau, xterm, I=0.0;
  int num_div=10000, n;

  sin_alpha_eq_sq = pow(sin(alpha_eq),2);
  endx = sin(lat*D2R);
  dx = endx/num_div;
  
  // Evaluate Walt's flight-time constant, from equator to mirror pt.
  walt_tau = 0.02925/R_E*C*(1 - 0.4635*pow(sin_alpha_eq_sq,0.375));

  // Evaluate flight-time integral from equator to lat  
  //for(x=0; x/endx < 1; x+=dx) {
  for(n=0; n<num_div; n++) {
    xterm = sqrt(1+3*x*x);
    I += dx*xterm/sqrt( 1 - (sin_alpha_eq_sq/pow((1-x*x),3))*xterm);
    x+=dx;
  }
  
  *flt_const_N = L*R_E*(walt_tau - I);
  *flt_const_S = L*R_E*(walt_tau + I);

  // printf("Flight constant: L: %g, lat: %g, alpha_eq: %g, fN: %g, fS: %g\n",L, lat, alpha_eq, *flt_const_N, *flt_const_S);

}









