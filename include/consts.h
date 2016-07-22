/* -----------------------------------
 *             consts.h         
 * -----------------------------------
 */
#ifndef consts_H
#define consts_H

// Fundamental constants:
#define		PI		    3.14159265358979311599796346854418516159
#define		D2R		    (PI/180.0)
#define		R2D		    (180.0/PI)
#define		Q_EL	    1.602E-19
#define		M_EL	 	9.1E-31
#define		E_EL		5.105396765648739E5
#define		MU0		    (PI*4E-7)		
#define		EPS0		8.854E-12
#define		C		    2.997956376932163e+08
#define		R_E		    6378000.0
#define		H_MAGNETO	1E6
#define		H_IONO		1E5

// For input power scaling:
#define     Z0      377.0
#define		A		5E3	
#define		B		1E5
#define		H_E		5000.0

// Time axis for simulation (axis to be interpolated onto -- not raytracer!)
#define     T_MAX        60   // T_MAX
#define     NUM_STEPS    6000
#define     T_STEP       (1.0*((1.0*T_MAX)/NUM_STEPS))

// Width in degrees around center latitude for which we'll load rays
#define     LAT_SPREAD 20 

#define     RAYTRACER_STEPS 65535 // Maximum number of steps to read from the rayfiles
                                  // (i.e., length to allocate; real vector may be less)

#define     F_STEP   1   // Hz. Separation in frequency (i.e., do 1-hz interpolation)
#define     LAT_STEP 0.01  // Degrees interpolation between launch rays

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

// EA array grid settings:
#define		EALimS		-50.0
#define		EALimN		50.0
#define		EAIncr		0.1
#define     NUMLATS     ((EALimN - EALimS)/EAIncr + 1)

// Width around field line, in L-shells, in which to consider a crossing
#define		L_MARGIN    6E-4
// // Delete these pls
// #define		EA_SPLIT	1
// #define		MULT		1 //2.0

#define 	WAVE_PWR_THRESH		0 //1e-6	// minimum wave power (ratio of initial)
										// Once a ray is damped to this power, all
										// subsequent steps are zero.


// For flux calculation:
#define     ALPHA_DISTRIBUTION 0    // "SQUARE" in previous versions

#endif