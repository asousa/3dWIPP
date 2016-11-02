/* -----------------------------------
 *       consts.h -- 3DWIPP        
 * -----------------------------------
 */
#ifndef consts_H
#define consts_H

// Fundamental constants:
// #define     PI          M_PI
#define		PI		    3.14159265358979311599796346854418516159
#define		D2R		    (PI/180.0)
#define		R2D		    (180.0/PI)
#define		Q_EL	    1.602E-19               // Electron charge
#define		M_EL	 	9.1E-31                 // Electron mass
#define		E_EL		5.105396765648739E5 
#define		MU0		    (PI*4E-7)		        
#define		EPS0		8.854E-12
#define		C		    2.997956376932163e+08
#define		R_E		    6378000.0               // Earth radius (meters)
#define		H_MAGNETO	1E6                     // Altitude of Magnetosphere (meters)
#define		H_IONO		1E5                     // Altitude of Ionosphere (meters)

// For input power scaling:
#define     Z0      377.0
#define		A		5E3	
#define		B		1E5
#define		H_E		5000.0

// Maximum number of species in the raytracer (to avoid using scalable arrays)
#define     NS_MAX  8

// Energy scaling (change these inputs):
//#define 	E_EXP_BOT	
#define 	E_MIN	1.0E1 	// ev
#define 	E_MAX	1.0E8 	// ev
#define 	NUM_E 	512 //32 //128 //1000 	// Number of energy bins
// #define 	SQUARE	1		// Square vs Sin particle distribution

// Calculated params
#define 	E_EXP_BOT	log10(E_MIN)
#define 	E_EXP_TOP 	log10(E_MAX)
#define 	DE_EXP	    ((E_EXP_TOP - E_EXP_BOT)/(NUM_E))

// // EA array grid settings:
#define		EALimS		-40.0
#define		EALimN	    40.0
#define     NUM_EA      11
#define		EAIncr	    ((EALimN - EALimS)/(NUM_EA - 1))	

#define     TRACER_MAX  20000    // Maximum steps in field line tracer
#define     TRACER_STEP 1e-3     // Step size in field line tracer (units of Earth radii)

// Radius around field line, in L-shells, in which to consider a crossing
#define		L_MARGIN    0.2

// Output time axis:
#define     TIME_MAX    10      // Seconds
#define     NUM_TIMES   1000    //ceil(TIME_MAX/TIME_STEP)

#define     TIME_STEP   (NUM_TIMES/TIME_MAX)     // Seconds

#define     SCATTERING_RES_MODES   5        

#endif