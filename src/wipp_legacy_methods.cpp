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

