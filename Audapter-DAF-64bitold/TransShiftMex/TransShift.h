/* 01/07/2007
TransShift.h

Defines the interface for the [a]-->[i] transition shifting algorithm

(c) 2007 Marc Boucek
speech communication group, RLE @MIT
*/

#pragma once

//#include "pvocWarpAtom.h"

typedef double mytype;

#define N_COEFFS_SRFILT 21
//#define DOWNSAMP_FACT 4
#define DOWNSAMP_FACT 3
#define MAX_FRAMELEN (384 / DOWNSAMP_FACT)
#define MAX_BUFLEN (15 * MAX_FRAMELEN)  // max ndelay< 8
#define MAX_NWIN   16
#define M_1_SQRTPI 0.564189583547756
#define M_PI       3.14159265358979323846
#define MAX_NFMTS 5
#define MAX_NLPC	20
#define MAX_NPOLES (MAX_NLPC / 2 + 2)
#define MAX_NLPC_SQR (MAX_NLPC * MAX_NLPC)
#define MAX_FMT_TRACK_JUMP MAX_NPOLES
#define MAX_NTRACKS 5
#define MAX_PITCHLEN 100	//100
#define MAX_NTRANS	10		//SC How many dipthong transitions are allowed at most? In practice, only one is used.
#define PF_NPOINTS	257
#define PF_NBIT	8

#define	MAX_ITER	10		//SC Maximum number of iteractions for the polynomial perpendicular foot finding

#define NFFT 1024			//SC(2008/05/06)
#define MAX_NFFT 4096		//SC(2012/03/05): for frequency/pitch shifting

#define MAX_DATA_VEC 2*MAX_NTRACKS+10+MAX_NLPC
#define MAX_REC_SIZE  120000 // 10 sec recording at p.sr = 48kHz/4 
#define MAX_DATA_SIZE 120000 // 10 sec data p.sr = 48kHz/4 and frameshift>=1

#define MAX_PB_SIZE   120000 // 2.5 sec, playback rate = 48 kHz.

#define MAX_DELAY_FRAMES  600 // The maximum number of added frame delays

// subfunction roots specific definitions
#define SIGN(a,b) ((b)>=0.0 ? fabs((a)) : -fabs((a)))
#define IMAX(k,j) ((k)<=(j) ? (j) : (k))
#define aMat(k,j) AHess[((j)-1)*nLPC+(k)-1]
#define ISABOVE(a,b)((a)>=(b) ? true : false)
#define ZEROIFABOVE(a,b) ((a)>=(b) ? 0 : (a))

int algoCallbackFunc(char *buffer, int buffer_size, void * data); // algorithm callback function 
int algoCallbackFuncSineGen(char *buffer, int buffer_size, void * data); //SC algorithm sine wave generator
int algoCallbackFuncWavePB(char *buffer, int buffer_size, void * data); //SC algorithm sine wave generator

// DSPLib routines - in case we have to go back to the TI board
void	DSPF_dp_blk_move(const mytype * x, mytype * r, const int nx);
mytype	DSPF_dp_vecsum_sq(const mytype *x,int n);                           
void	DSPF_dp_biquad(mytype * x, mytype * b, mytype * a, mytype * delay, mytype * r, int nx);
void	DSPF_dp_vecmul(const mytype * x, const mytype * y, mytype * r, int n);
void	DSPF_dp_autocor(mytype * r, mytype * x, int nx, int nr);

class pvocWarpAtom{
public:
	mytype tBegin;
	mytype rate1;	// Should aways be in the interval of (0, 1)
	mytype dur1;
	mytype durHold;
	mytype rate2;	// Should aways be > 1
	mytype dur2;

	pvocWarpAtom(){
		tBegin = 0.25;
		rate1 = 0.5;	
		dur1 = 0.25; 
		durHold = 1;
		rate2 = 2;
			
		dur2 = (1 - rate1) / (rate2 - 1) * dur1;
	}
	
	pvocWarpAtom(mytype t_tBegin, mytype t_rate1, mytype t_dur1, mytype t_durHold, mytype t_rate2){
		tBegin = t_tBegin;
		rate1 = t_rate1;
		dur1 = t_dur1;
		durHold = t_durHold;
		rate2 = t_rate2;

		dur2 = (1 - rate1) / (rate2 - 1) * dur1;
	}
};

// Class TransShift
// Implements the audio input - process - output queue
class TransShift{
public:
	// Constructor initializes all variables
	TransShift();

	// The Reset function reinitializes all internal buffers
	void reset();

	// The main function that controls the processing
	int myProcess(mytype *inFrame_ptr, mytype *outFrame_ptr, int frame_size);
	int myProcessSineGen(mytype *inFrame_ptr, mytype *outFrame_ptr, int frame_size);	//SC Since wave generator
	int myProcessWavePB(mytype *inFrame_ptr, mytype *outFrame_ptr, int frame_size);	//SC Wave playback

	// Allows external entities to set / get various parameters of the process
	int setparams(void * name, void * value);
	int getparams(void * name);
	const mytype* getsignal(int & size);
	const mytype* getdata(int & size,int & vecsize);	

	int pbCounter;	//SC The integer counter used in wave playback 

	pvocWarpAtom *warpCfg;

private:
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VARIABLES & COUNTERS  *****************************************************%%%%%%%		

	// Variables	
	int doEdge;
	bool   bTransReset;			// set true when resetted, set to false after first process

	bool	bFrameLenShown;

	// counters
	int frame_counter;		// frame counter : increments at every new frame from soundcard	(frame rate)						
	int data_counter;		// data counter  : increment rate is process rate ,i.e. frame rate * nWIn
	int circ_counter;		// circular data counter , loops over MAX_PITCHLEN

	// 
	mytype time_step;       // process time unit 
	mytype ma_rms1;			// movinga average rms 
	mytype ma_rms2;	

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BUFFERS  *****************************************************%%%%%%%%%%%	

	// buffers stores samples at original sample rate
	mytype inFrameBuf[MAX_FRAMELEN*DOWNSAMP_FACT];      // stores incoming  frame from soundcard @p.sr*downsamp_fact
	//mytype outFrameBuf[MAX_FRAMELEN*DOWNSAMP_FACT];     // stroes outgoing frame to soundcard @p.sr*downsamp_fact
	//SC(2012/02/28) DAF: expanded outFrameBuf for DAF
	mytype outFrameBuf[MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES]; // Shifting output buffer after formant shifting, before pitch shifting
	mytype outFrameBufPS[MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES]; // Shifting output buffer after pitch shifting

	int outFrameBuf_circPtr;

	mytype srfilt_buf[MAX_FRAMELEN*DOWNSAMP_FACT];      // multiple purpose before for down / upsampling

	// buffers stores samples at downsampled rate
	mytype filtbuf[MAX_FRAMELEN];                       // filter buffer (formant shift)
	mytype oBuf[MAX_BUFLEN];							// Buffer stores original input samples
	mytype pBuf[MAX_BUFLEN];							//SC-Mod(2008/05/15) Buffer stores preemphasized input samples
	mytype inBuf[MAX_FRAMELEN];							// Buffer stores input samples after downsampling
	mytype outBuf[MAX_FRAMELEN];							// Buffer stores processed output
	mytype fakeBuf[MAX_FRAMELEN];							// trash buffer to flush filter states
	mytype zeros[MAX_FRAMELEN];						    // utility buffer for filtering


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VARIOUS EXTRACTED DATA  *****************************************************%%%%%%%%%%%	


	// other buffers : process rate = downsampled rate * nWin
	mytype lpcAi[MAX_NLPC+1];						// lpc coeeficients	
	mytype hwin[MAX_BUFLEN];						// Hanning window
	mytype hwin2[MAX_BUFLEN];						// Hanning window for frequency/pitch shifting
	mytype realRoots[MAX_NLPC];						// real part of roots
	mytype imagRoots[MAX_NLPC];						// imag part of roots
	mytype Acompanion[(MAX_NLPC_SQR)];				// companion matrix for (eigenvalue ---> )roots calculation
	mytype amps[MAX_NPOLES];						// radius of poles
	mytype orgPhis[MAX_NPOLES];						// angle of original poles
	mytype bw[MAX_NPOLES];							// bandwith of poles

	mytype newPhis[MAX_NPOLES];						// new (shifted) pole angles 
	mytype gtot[MAX_NWIN];							// gain factor (for gain adaption)

	mytype weiMatPhi[MAX_NTRACKS][MAX_PITCHLEN];	// weithed matrix of past pole angles (for moving average formant smoothing)
	mytype weiMatBw[MAX_NTRACKS][MAX_PITCHLEN];		// weithed matrix of past bandwith (for moving average formant smoothing)
	mytype weiVec[MAX_PITCHLEN];					// weigthing vector of past weigths (default weigth : short time rms ... )
	mytype sumWeiPhi[MAX_NTRACKS];					// sum of past weights (rms) *  angles (over avglen )
	mytype sumWeiBw[MAX_NTRACKS];                   // sum of past weights (rms) *  angles (over avglen )
	mytype sumWei;									// sum of weigths (over pitchlen)

	mytype wmaPhis[MAX_NTRACKS];					// weighted moving average pole angles
	mytype wmaR[MAX_NTRACKS];						// weighted moving average pole radius
	mytype deltaPhi[MAX_NFMTS];						// derivative of angle
	mytype fmts[MAX_NTRACKS];						// estimated formant tracks
	mytype dFmts[MAX_NFMTS];						// formant derivatives
	mytype sFmts[MAX_NFMTS];						// shifted formants (=0 when no shift)

	// getDFmt
	mytype deltaFmt[2];								// formant derivatives
	mytype deltaMaFmt[2];							// moving average deltafmt
	mytype lastFmt[2];								// last formant estimate

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILTER COEEFS AND DELAYS  *****************************************************	

	// preemphasis 
	mytype a_preemp[2];								// denominator coefficients 
	mytype b_preemp[2];								// numerator coefficients
	mytype preemp_delay[1];							// filter delay

	// deemphasis
	mytype a_deemp[2];								// denominator coefficients 
	mytype b_deemp[2];								// numerator coefficients
	mytype deemp_delay[1];								// filter delay


	// first filter (f1 shift)
	mytype a_filt1[2];								// denominator coefficients 
	mytype b_filt1[3];								// numerator coefficients
	mytype filt_delay1[2];							// filter delays


	// second filter (f2 shift)
	mytype a_filt2[2];								// denominator coefficients 
	mytype b_filt2[3];								// numerator coefficients
	mytype filt_delay2[2];							// filter delays


	// sample rate conversion filter 
	mytype srfilt_a[N_COEFFS_SRFILT];				// denominator coefficients
	mytype srfilt_b[N_COEFFS_SRFILT];				// numerator coefficients
	mytype srfilt_delay_up[N_COEFFS_SRFILT-1];		// filter delays for upsampling
	mytype srfilt_delay_down[N_COEFFS_SRFILT-1];    // filter delays for downsampling

	//SC(2008/05/07)
	mytype ftBuf1[NFFT*2];
	mytype ftBuf2[NFFT*2];
	mytype ftBuf1ps[MAX_NFFT * 2]; // For frequency/pitch shifting
	mytype ftBuf2ps[MAX_NFFT * 2]; // For frequency/pitch shifting

	//SC(2012/03/05): Pitch Shifting
	mytype lastPhase[MAX_NFFT];
	mytype sumPhase[MAX_NFFT];
	mytype outputAccum[MAX_NFFT];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DATA   RECORDING  *****************************************************%%%%%%%%%%%	


	//recorded data is stored here
	mytype signal_recorder[2][MAX_REC_SIZE];            // downsampled input signal and downsampled output signal
	mytype data_recorder[MAX_DATA_VEC][MAX_DATA_SIZE];  // stores other data
	mytype transdata[2];								// only transition data

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DATA PLAYBACK     **************************************************************
	//SC-Mod(2007/12/28)
	mytype data_pb[MAX_PB_SIZE];

	//SC-Mod(2008/04/07) Vowel related
	mytype f1m;     // F1 (mel)
	mytype f2m;     // F2 (mel)
	mytype f2mp;	// F2 (mel) of the foot of perpendicular to the fit polynomial.

	//SC-Mod(2008/05/15) FFT related
	mytype fftc[NFFT*2];
	mytype fftc_ps0[MAX_NFFT * 2];
	mytype fftc_ps[MAX_NFFT * 2];	//SC(2012/03/05) For frequency/pitch shifting

	//SC(2012/02/29) For data writing:
	int dataFileCnt;

	//SC(2012/03/13) For PVOC time warping
	mytype pvocWarpCache[MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64][1024 * 2];

	//pvocWarpAtom warpCfg;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS  *****************************************************%%%%%%%%%%%

	// Externally adjustable parameters
	struct{

		// framing and processing 
		int	   frameLen;				// length of one frame ( framelen = nWin * frameshift)
		int	   nWin;					// number of processes per frame (of frameLen samples)
		int	   frameShift;				// number of samples shift between two processes ( = size of processed samples in 1 process)	
		int	   bufLen;					// main buffer length : buflen stores (2*nDelay -1)*frameLen samples
		int    nDelay;					// number of delayed framas (of size framelen) before an incoming frame is sent back
		// thus the overall process latency (without souncard) is :Tproc=nDelay*frameLen/sr
		int	   anaLen;					// size of lpc analysis (symmetric around window to be processed)
		int	   pvocFrameLen;			// size of analysis window for pitch shifting (aim for power of 2) 
		int	   pvocHop;					// step size of analysis window (must be multiples of frameLen)

		int    avgLen;				// length of smoothing ( should be approx one pitch period, 
		// can be greater /shorter if you want more / lesss moothing)
		// avgLen = 1 ---> no smoothing ( i.e. smoothing over one value)

		int    sr;						// internal samplerate

		int    nLPC;					// lpc order ... number of lpc coeefs= nLPC +1
		int    nFmts;					// originally the number of formants you want shift ( hardseted = 2)

		mytype dRMSThresh;				// RMS threshhold for voicing detection
		mytype dRMSRatioThresh;			// preemp / original RMS ratio threshhold, for fricative detection 
		mytype rmsFF;                   // rms fogetting factor for long time rms 


		mytype dPreemp;					// preemphasis factor
		mytype dScale;					// scaling the output (when upsampling) (does not affect internal signal

		// for transition detection
		mytype dFmtsFF;					// formant forgeeting factor for s
		mytype maxDelta;				// maximal allowed formant derivative 		
		mytype fmtDetectStart[2];		// formant frequencies at start of transition (goal region),i.e. [a]		
		int    minDetected;				// min transition detection in row before starting transition		


		// for formant tracking algorithm
		mytype trackFF;
		int	   nTracks;					// number of tracked formants 
		int    nCands;					// number of possible formant candiates  ( > ntracks     but  < p.nLPC/2!!!! (choose carefully : not idiot proofed!)

		mytype trackIntroTime;          // intro tracking time 

		//mytype aFact1;				// 
		mytype bFact;
		//mytype gFact1;
		mytype aFact;
		//mytype bFact2;
		mytype gFact;

		mytype fn1;		// Neutral F1: used in formant tracking
		mytype fn2;		// Neutral F2: used in formant tracking			

		// for shifting
		mytype dDev;					// deviation vector (in angle euclidien distance) @ max of hanning window

		// booleans						
		int    bRecord;					// record signal
		int    bTrack;					// use formant tracking algorithm
		int    bShift;					// do shifting
		int    bGainAdapt;				// use gain adaption
		int    bDetect;					// detect transition	
		int    bRelative;				// shift relative to actual formant point, (otherwise absolute coordinate)			
		int    bWeight;					// do weighted moving average formant smoothing (over pitchlen) , otherwise not weigthed (= simple moving average)				
		int	   bCepsLift;				//SC-Mod(2008/05/15) Whether the cepstral lifting is done before the autocorrelation

		int	   bRatioShift;				//SC(2009/01/20). 
		int	   bMelShift;				//SC(2009/01/20). 

		//SC(2012/03/05) Frequency/pitch shifting
		int		bPitchShift;
		mytype	pitchShiftRatio;
		
		//SC-Mod(2008/01/05). Arrays: for intensity correction during formant shifting (mainly for the use of parallel shifts)		
		
		mytype wgFreq;										//SC Wave generator frequency (Hz)
		mytype wgAmp;										//SC Wave generator amplitude (digitized peak amplitude)
		mytype wgTime;										//SC Wave generator current phase (rad)

		//SC(2008/04/03). Perturbation-related variables		
		int transCounter;
		mytype F2Min;
		mytype F2Max;
		//mytype triF2MinMin;
		//mytype triF2MaxMax;
		mytype F1Min;
		mytype F1Max;
		mytype LBk;
		mytype LBb;
		mytype pertF2[PF_NPOINTS];		
		mytype pertPhi[PF_NPOINTS];
		mytype pertAmp[PF_NPOINTS];
		mytype minVowelLen;
		
		bool transDone;

		//SC(2008/06/20) Feedback type: fb=0: silent; fb=1: voice only (pert. or unpert.); fb=2: noise only; fb=3: voice+noise;
		int fb;		
		
		//SC(2008/06/22)
		int bRepData;		// Replicate data or replay, in mode 1. 0: no; 1: yes.

		//SC(2008/05/15) Cepstral lifting related
		int	cepsWinWidth;

		//SC(2008/06/22)
		mytype trialLen;
		mytype rampLen;

		//SC(2012/02/28)
		int delayFrames;	// For DAF: the number of frames delayed. 0 correponds to no _added_ delay (~11 ms).

		int bDownSampFilt;


		
	} p;



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  *****************************************************%%%%%%%%%%%

	void     getAi(mytype* xx, mytype* aa, const int size, const int nlpc);
	mytype  getGain(mytype * r, mytype * ophi,mytype * sphi, int nfmts);
	int     hqr_roots (	mytype *c, 	mytype *wr, mytype *wi,	mytype *Acompanion, const int nLPC);
	void    getRPhiBw (mytype *wr,  mytype *wi, mytype *radius,  mytype *phi ,mytype *bandwith);
	void trackPhi(mytype *r_ptr,mytype *phi_ptr,mytype time);
	void myFilt (mytype *xin_ptr, mytype* xout_ptr,mytype *oldPhi_ptr,mytype *newPhi_ptr,mytype *r_ptr,const int size);
	mytype calcRMS1(const mytype *xin_ptr, int size);	
	mytype calcRMS2(const mytype *xin_ptr, int size);	
	int getWma(mytype *phi_ptr, mytype *bw_ptr , mytype * wmaPhi_ptr, mytype * wmaR_ptr);
	void    levinson(mytype *R, mytype* aa, int size);
	int gainAdapt(mytype *buffer,mytype *gtot_ptr,int framelen, int frameshift);
	int sign(mytype x);
	bool  detectTrans(mytype *fmt_ptr, mytype *dFmt_ptr,int datcnt, mytype time);
	int getDFmt(mytype *fmt_ptr,mytype *dFmt_ptr, mytype time);	
	void upSampSig (mytype *b, mytype *a, mytype *x, mytype *buffer, mytype *r,mytype  *d,const int nr, const int n_coeffs, const int upfact,const mytype scalefact);
	void downSampSig(mytype *b, mytype *a, mytype *x, mytype *buffer, mytype *r,mytype  *d,const int nr, const int n_coeffs, const int downfact);
	void downSampSig_noFilt(mytype *b, mytype *a, mytype *x, mytype *buffer, mytype *r,mytype  *d,const int nr, const int n_coeffs, const int downfact);
	void iir_filt (mytype *b, mytype *a,  mytype *x, mytype *r,mytype  *d,const int nr, const int n_coeffs,  mytype g);
			
	mytype TransShift::hz2mel(mytype hz);
	mytype TransShift::mel2hz(mytype hz);
	mytype TransShift::locateF2(mytype f2);

	void	DSPF_dp_cfftr2(int n, mytype * x, mytype * w, int n_min);
	void	DSPF_dp_icfftr2(int n, double * x, double * w, int n_min);
	void	gen_w_r2(double* w, int n);
	void	bit_rev(double* x, int n);

	void	smbFft(mytype *fftBuffer, double fftFrame_Size, int sign);
};

