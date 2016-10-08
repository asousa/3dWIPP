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
