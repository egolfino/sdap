/* 12/28/2008
TransShift.cpp

Vowel formant shifting algorithm
Can be used to shift formant frequencies of pure vowels or time-varying vowels (diphthongs or tripththongs) in 
	static or time-varying fashion. 

Incorporated 
	1) Cepstral liftering (optional)
	2) Dynamic programming style formant tracking (Xia and Espy-Wilson, 2000, ICSLP)
	3) Gain adaptation after shifting (optional)
	4) RMS-based formant smoothing (to reduce the influence of subglottal coupling on formant estimates) (optional)

(c) 2007 Marc Boucek
(c) 2008 Shanqing Cai (cais@mit.edu)

Speech Communication Group, RLE, MIT
*/

//#include <string.h> 
#include <stdio.h>

#define TRACE printf

#include "TransShift.h"
#include <math.h>

//#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define RSL(INTEGER,SHIFT) (int)( ( (unsigned)INTEGER ) >> SHIFT )
#define NPARAMS 29
char *intparamarray[NPARAMS]={"srate","framelen","ndelay","nwin","nlpc","nfmts","ntracks",
//								 1		   2		3		4				5		6		7								
	"avglen","cepswinwidth","fb","minvowellen", "delayframes", "bpitchshift", "pvocframelen", "pvochop", "bdownsampfilt"};
//		8			9		 10		11				12				13				14			15			16

char *doubleparamarray[NPARAMS]={"scale","preemp","rmsthr","rmsratio","rmsff","dfmtsff",
//									  1        2         3       4        5        6
 "wgfreq","wgamp","wgtime","datapb",
//	  7		  8		   9		10
 "f2min","f2max","pertf2","pertamp","pertphi","f1min","f1max","lbk","lbb","triallen","ramplen",
//   11		  12      13       14       15        16      17	  18   19     20          21
 "afact","bfact","gfact","fn1","fn2", "pitchshiftratio", "pvocwarp"};
//   22	     23      24     25    26		27				28

char *boolparamarray[NPARAMS]={"bgainadapt","bshift","btrack","bdetect","bweight","bcepslift","bratioshift","bmelshift"};
//									1          2         3         4          5           6			7			8

char *vectarray[NPARAMS]={"trajectory"};
//                             1 
#ifdef TIME_IT      
LARGE_INTEGER freq, time1, time2;
LONGLONG overhead;
#endif

// Non-class function that is called by the audioIO routines when data becomes available
int algoCallbackFunc(char *buffer, int buffer_size, void * data)	//SC The input is 8-bit, so char type is proper for buffer.
{
	TransShift *algo = (TransShift*)data;	//SC copy data to algo
#ifdef TIME_IT
	QueryPerformanceCounter(&time1);
#endif


	algo->myProcess((mytype*)buffer, (mytype*)buffer, buffer_size);		//SC(12/19/2007)

#ifdef TIME_IT
	QueryPerformanceCounter(&time2);
	TRACE("%.6f\n",double(time2.QuadPart - time1.QuadPart - overhead)/double(freq.QuadPart));
#endif
	return 0;
}

int algoCallbackFuncSineGen(char *buffer, int buffer_size, void * data)	//SC sine wave generator
{
	TransShift *algo = (TransShift*)data;	//SC copy data to algo
#ifdef TIME_IT
	QueryPerformanceCounter(&time1);
#endif
	//TRACE("%d\n",buffer_size);
	//DSPF_dp_blk_move((MY_TYPE*)buffer, (MY_TYPE*)buffer, buffer_size);

	algo->myProcessSineGen((mytype*)buffer, (mytype*)buffer, buffer_size);


#ifdef TIME_IT
	QueryPerformanceCounter(&time2);
	TRACE("%.6f\n",double(time2.QuadPart - time1.QuadPart - overhead)/double(freq.QuadPart));
#endif
	return 0;
}

int algoCallbackFuncWavePB(char *buffer, int buffer_size, void * data)	//SC sine wave generator
{
	TransShift *algo = (TransShift*)data;	//SC copy data to algo
#ifdef TIME_IT
	QueryPerformanceCounter(&time1);
#endif

	algo->myProcessWavePB((mytype*)buffer, (mytype*)buffer, buffer_size);


#ifdef TIME_IT
	QueryPerformanceCounter(&time2);
	TRACE("%.6f\n",double(time2.QuadPart - time1.QuadPart - overhead)/double(freq.QuadPart));
#endif
	return 0;
}


TransShift::TransShift()		//SC construction function
{//modifiable parameters ( most of them can be modified externally) 
		int		n;

		p.sr				= 48000 / DOWNSAMP_FACT;				// internal samplerate (souncard samplerate = p.sr*DOWNSAMP_FACT)
		p.nLPC				= 13;					// LPC order ... number of lpc coeefs= nLPC +1
		p.nFmts				= 2;					// originally the number of formants you want shift ( hardseted = 2)

		// framing and processing 
		p.frameLen			= 96 / DOWNSAMP_FACT;		// length of one internal frame ( framelen = nWin * frameshift) (souncard frame length = p.frameLen *DOWNDAMP_FACT)
		p.nDelay			= 7;					// number of delayed framas (of size framelen) before an incoming frame is sent back
													// thus the overall process latency (without souncard) is :Tproc=nDelay*frameLen/sr
		p.bufLen			= (2*p.nDelay-1)*p.frameLen;	// main buffer length : buflen stores (2*nDelay -1)*frameLen samples
		
		p.nWin				= 1;					// number of processes per frame (of frameLen samples)	
		p.frameShift		= p.frameLen/p.nWin;	// number of samples shift between two processes ( = size of processed samples in 1 process)	
		p.anaLen			= p.frameShift+2*(p.nDelay-1)*p.frameLen;// size of lpc analysis (symmetric around window to be processed)
		p.pvocFrameLen			= p.frameShift + (2 * p.nDelay - 3) * p.frameLen;// For frequency/pitch shifting: size of lpc analysis (symmetric around window to be processed)
		p.avgLen			= 10;				    // length of smoothing ( should be approx one pitch period, 
		// can be greater /shorter if you want more / lesss smoothing)
		// avgLen = 1 ---> no smoothing ( i.e. smoothing over one value)

		// RMS
		p.dRMSThresh		= 0.02;	// RMS threshhold for voicing detection
		p.dRMSRatioThresh	= 1.3;	// preemp / original RMS ratio threshhold, for fricative detection 
		p.rmsFF				= 0.9;  // rms fogetting factor for long time rms 

		p.dPreemp			= .98;	// preemphasis factor
		p.dScale			= 1;	// scaling the output (when upsampling) (does not affect internal signal

		// for transition detection
		p.dFmtsFF			= 0;	// formant forgeeting factor for s
		p.maxDelta			= 40;	// maximal allowed formant derivative 		
		p.fmtDetectStart[0] = 800;	// formant frequencies at start of transition (goal region),i.e. [a]
		p.fmtDetectStart[1] = 1600;	// formant frequencies at start of transition (goal region),i.e. [a]		
		p.minDetected		= 10;	// min transition detection in row before starting transition		

		// for formant tracking algorithm
		p.trackIntroTime	= 100;	// time in ms
				
		p.aFact			= 1;	// end value ....		
		p.bFact			= 0.8;	// end value ....
		p.gFact			= 1;	// end value ....

		p.fn1			= 633;	// A priori expectation of F1 (Hz)
		p.fn2			= 1333; // A priori expectation of F2 (Hz)

		p.nTracks			= 4;	// number of tracked formants 
		p.nCands			= 6;	// number of possible formant candiates  ( > ntracks     but  < p.nLPC/2!!!! (choose carefully : not idiot proofed!)
		p.trackFF			= 0.95;		

		// booleans						
		p.bRecord			= 1;	// record signal, should almost always be set to 1. 
		p.bTrack			= 1;	// use formant tracking algorithm
		p.bShift			= 1;	// do shifting	//SC-Mod(09/23/2007): changed from 0 to 1
		p.bGainAdapt		= 1;	// use gain adaption
		p.bDetect			= 0;	// detect transition
		p.bRelative			= 1;	// shift relative to actual formant point, (otherwise absolute coordinate)			
		p.bWeight			= 1;	// do weighted moving average formant smoothing (over avglen) , otherwise not weigthed (= simple moving average)				
		p.bCepsLift			= 0;	//SC-Mod(2008/05/15) Do cepstral lifting by default
		
		p.bRatioShift		= 0;	//SC(2009/01/20)
		p.bMelShift			= 1;	//SC(2009/01/20)

		//SC(2012/03/05)
		p.bPitchShift		= 0;
		p.pitchShiftRatio   = 1.;	// 1. = no shift.
		
		//SC Initialize the playback data and the counter
		for(n=0;n<MAX_PB_SIZE;n++){
			data_pb[n]      = 0;
		}
		pbCounter			= 0;

		p.wgFreq			= 1000.;
		p.wgAmp				= 0.1;
		p.wgTime			= 0.;
		
		p.F2Min		= 0;		// Lower boundary of the perturbation field (mel)
		p.F2Max		= 0;		// Upper boundary of the perturbation field (mel)
		p.F1Min		= 0;		// Left boundary of the perturbation field (mel)
		p.F1Max		= 0;		// Right boundary of the perturbation field (mel)
		p.LBk		= 0;		// The slope of a tilted boundary: F2 = p.LBk * F1 + p.LBb. (mel/mel)
		p.LBb		= 0;		// The intercept of a tilted boundary (mel)

		for(n=0;n<PF_NPOINTS;n++){
			p.pertF2[n]=0;			// Independent variable of the perturbation vectors
			p.pertAmp[n]=0;			// Magnitude of the perturbation vectors (mel)
			p.pertPhi[n]=0;			// Angle of the perturbation vectors (rad). 0 corresponds to the x+ axis. Increase in the countetclockwise direction. 
		}

		p.minVowelLen=60;
		p.transDone=false;
		p.transCounter=0;

		//SC(2008/05/15)
		p.cepsWinWidth=30;

		//SC(2008/06/20)
		p.fb=1;		

		//SC(2008/06/22)
		p.trialLen = 9;	//sec
		p.rampLen=0.05;	//sec

		//SC(2012/02/28) DAF
		p.delayFrames = 0; // Unit: # of frames (no delay by default)

		//SC(2012/02/29) For data saving
		dataFileCnt = 0;

		//SC(2012/03/09) Phase vocoder (pvoc) pitch shifting
		p.pvocFrameLen = 1024;
		p.pvocHop = 256;		// 256 = 16 * 16

		p.bDownSampFilt = 1;

		//SC(2012/03/13) PVOC Time warp
		warpCfg = new pvocWarpAtom();
//************************************** Initialize filter coefs **************************************	



	// preemphasis filter
	a_preemp[0] = 1;
	a_preemp[1] = 0;
	b_preemp[0] = 1;
	b_preemp[1] = -p.dPreemp;

	// deemphasis filter
	a_deemp[0] = 1;
	a_deemp[1] = -p.dPreemp;
	b_deemp[0] = 1;
	b_deemp[1] = 0;


    // shift filter
	b_filt1[0] = 1;
	b_filt2[0] = 1;

	/*
	srfilt_b[0] = 0.003207504864776675100000;
	srfilt_b[1] = -0.015051748179863501000000;
	srfilt_b[2] = 0.045149574027758405000000;
	srfilt_b[3] = -0.090740842850500311000000;
	srfilt_b[4] = 0.144017161473200290000000;
	srfilt_b[5] = -0.177763639685871390000000;
	srfilt_b[6] = 0.185429987362228570000000;
	srfilt_b[7] = -0.156965404667578690000000;
	srfilt_b[8] = 0.120750787331069550000000;
	srfilt_b[9] = -0.082508106409673099000000;
	srfilt_b[10] = 0.073479882151955889000000;
	srfilt_b[11] = -0.082508106409664439000000;
	srfilt_b[12] = 0.120750787331029840000000;
	srfilt_b[13] = -0.156965404667645190000000;
	srfilt_b[14] = 0.185429987362568880000000;
	srfilt_b[15] = -0.177763639686368770000000;
	srfilt_b[16] = 0.144017161473438790000000;
	srfilt_b[17] = -0.090740842850170200000000;
	srfilt_b[18] = 0.045149574027470510000000;
	srfilt_b[19] = -0.015051748179931381000000;
	srfilt_b[20] = 0.003207504864907421600000;

	srfilt_a[0] = 1.000000000000000000000000;
	srfilt_a[1] = -8.054528775835033000000000;
	srfilt_a[2] = 33.086573089411893000000000;
	srfilt_a[3] = -90.193446136595213000000000;
	srfilt_a[4] = 181.256203524881440000000000;
	srfilt_a[5] = -283.571375474201600000000000;
	srfilt_a[6] = 356.755781218958760000000000;
	srfilt_a[7] = -368.305783722966910000000000;
	srfilt_a[8] = 316.054159733307870000000000;
	srfilt_a[9] = -227.162501157152780000000000;
	srfilt_a[10] = 137.260014043998750000000000;
	srfilt_a[11] = -69.738296640133896000000000;
	srfilt_a[12] = 29.705121931131252000000000;
	srfilt_a[13] = -10.537214215529730000000000;
	srfilt_a[14] = 3.080153549165016300000000;
	srfilt_a[15] = -0.729595107204623840000000;
	srfilt_a[16] = 0.136881161697956440000000;
	srfilt_a[17] = -0.019547624983513916000000;
	srfilt_a[18] = 0.002060528139489559900000;
	srfilt_a[19] = -0.000143860211730633450000;
	srfilt_a[20] = 0.000014362805789346628000;
	*/

	srfilt_b[0] = 0.005985366448016847400000;
	srfilt_b[1] = -0.000000000068663473436596;
	srfilt_b[2] = 0.029926833561855812000000;
	srfilt_b[3] = 0.014963399494903253000000;
	srfilt_b[4] = 0.072946803942075492000000;
	srfilt_b[5] = 0.066399110082245749000000;
	srfilt_b[6] = 0.128831523706446540000000;
	srfilt_b[7] = 0.141307195958322970000000;
	srfilt_b[8] = 0.183515087779119460000000;
	srfilt_b[9] = 0.196038702055692930000000;
	srfilt_b[10] = 0.207578586483177310000000;
	srfilt_b[11] = 0.196038702055692630000000;
	srfilt_b[12] = 0.183515087779119760000000;
	srfilt_b[13] = 0.141307195958322280000000;
	srfilt_b[14] = 0.128831523706446790000000;
	srfilt_b[15] = 0.066399110082245277000000;
	srfilt_b[16] = 0.072946803942075533000000;
	srfilt_b[17] = 0.014963399494903067000000;
	srfilt_b[18] = 0.029926833561855826000000;
	srfilt_b[19] = -0.000000000068663525932820;
	srfilt_b[20] = 0.005985366448016846500000;

	srfilt_a[0] = 1.000000000000000000000000;
	srfilt_a[1] = -4.137689759094149300000000;
	srfilt_a[2] = 11.417342955970334000000000;
	srfilt_a[3] = -21.230389508442666000000000;
	srfilt_a[4] = 31.507204607241498000000000;
	srfilt_a[5] = -36.677292780605917000000000;
	srfilt_a[6] = 36.042584528469732000000000;
	srfilt_a[7] = -28.996821243768743000000000;
	srfilt_a[8] = 20.262367357856544000000000;
	srfilt_a[9] = -11.637468104552259000000000;
	srfilt_a[10] = 5.968975493498319000000000;
	srfilt_a[11] = -2.417954280896708500000000;
	srfilt_a[12] = 0.941027354810217260000000;
	srfilt_a[13] = -0.241109659478893040000000;
	srfilt_a[14] = 0.083935453370180629000000;
	srfilt_a[15] = -0.005511361553189712100000;
	srfilt_a[16] = 0.006142808678570149300000;
	srfilt_a[17] = 0.001292100725808184000000;
	srfilt_a[18] = 0.000588047191250507470000;
	srfilt_a[19] = 0.000146757274221299580000;
	srfilt_a[20] = 0.000035865709068928935000;

	//SC-Mod(2008/05/15) FFT related
	gen_w_r2(fftc, NFFT);
	gen_w_r2(fftc_ps, MAX_NFFT);

	// Pitch and frequency shifting-related
	for (int i0 = 0; i0 < MAX_NFFT; i0++){
		fftc_ps0[i0] = 0;
		fftc_ps[i0] = 0;
	}

	reset();
}

void TransShift::reset()
{// resets all 
	int i0,j0;
	bTransReset				= true;

	bFrameLenShown			= false;

//*****************************************************  BUFFERS   *****************************************************
	// Initialize input, output and filter buffers (at original sample rate!!!)
	for(i0=0;i0<MAX_FRAMELEN*DOWNSAMP_FACT;i0++)
	{
		inFrameBuf[i0]=0;		
		srfilt_buf[i0]=0;
	}

	for (i0 = 0; i0 < MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES; i0 ++){
		outFrameBuf[i0] = 0;
		outFrameBufPS[i0] = 0;
	}
	outFrameBuf_circPtr = 0;

	// Initialize internal input, output buffers  (at downsampled  rate !!!)
	for(i0=0;i0<MAX_BUFLEN;i0++)
	{
		inBuf[i0] = 0;
		outBuf[i0] = 0;
		zeros[i0]  = 0;
	}

//*****************************************************  FILTER STATES  *****************************************************

	// reinitialize formant shift filter states
	for(i0=0;i0<2;i0++)
	{
		filt_delay1[i0]=0;
		filt_delay2[i0]=0;
	}

	// reinitialize preempahsis and deemphasis filter states
		deemp_delay[0]=0;
		preemp_delay[0]=0;


	// reinitialize down - and upsampling filter states
	for(i0=0;i0<N_COEFFS_SRFILT-1;i0++)
	{
		srfilt_delay_up[i0]=0;
		srfilt_delay_down[i0]=0;
	}
	
	

//*****************************************************  RECORDING  *****************************************************

	// Initialize signal recorder
	for(i0=0;i0<2;i0++)
	{
		for(j0=0;j0<MAX_REC_SIZE;j0++)
		{
			signal_recorder[i0][j0]=0;
		}
	}


	// Initialize data recorder
	for(i0=0;i0<MAX_DATA_VEC;i0++)
	{
		for(j0=0;j0<MAX_DATA_SIZE;j0++)
		{
			data_recorder[i0][j0]=0;
		}
	}

	// Initialize variables used in the pitch shifting algorithm 
	// SC(2012/03/05)
	for (i0 = 0; i0 < NFFT; i0++){
		lastPhase[i0] = 0;
		sumPhase[i0] = 0;
		outputAccum[i0] = 0;
	}

	// initialize recording counters
	frame_counter=0;
	data_counter=0;
	circ_counter=0;

// %%%%%%%%%%%%%%%%%%%%%%%%%%              INITIALYZE FUNCTIONS VARIABLES         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

//*****************************************************  getDfmt   *****************************************************


	for(i0=0;i0<2;i0++)
	{
		deltaFmt[i0]=0;
		deltaMaFmt[i0]=0;
		lastFmt[i0]=0;
	}


//*****************************************************  calcRms    *****************************************************	
		
	//initialize  moving average
	ma_rms1=0;
	ma_rms2=0;	

//*****************************************************  getWma    *****************************************************

	// initialize weight matrix and moving sum
	for(j0=0;j0<MAX_NTRACKS;j0++)
	{
		sumWeiPhi[j0]=0;
		sumWeiBw[j0]=0;
		for(i0=0;i0<MAX_PITCHLEN;i0++)
		{
			weiVec[i0]=0;
			weiMatPhi[j0][i0]=0;
			weiMatBw[j0][i0]=0;
		}
	}
	sumWei=0;

//*****************************************************  getAI    *****************************************************

	// Initialize hanning window
	for(i0=0;i0<p.anaLen/2;i0++){
		hwin[i0] = 0.5*cos(mytype(2*M_PI*(i0+1))/mytype(p.anaLen+1)); 
		hwin[i0] = 0.5 - hwin[i0];
		hwin[p.anaLen-i0-1] = hwin[i0];
	}

	// Initialize hanning window (for frequency / pitch shifting) SC(2012/03/05)
	for(i0=0;i0<p.pvocFrameLen;i0++){
		hwin2[i0] = 0.5*cos(mytype(2*M_PI*i0)/mytype(p.pvocFrameLen)); 
		hwin2[i0] = 0.5 - hwin2[i0];
		//hwin2[p.pvocFrameLen-i0-1] = hwin2[i0];
	}

//*****************************************************  hqr_Roots    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	//  The following forms a template for the initial companion matrix
    for (i0=(p.nLPC*p.nLPC) -1; i0>=0; i0--) 
		Acompanion[i0] = 0.0F;
 	for (i0=(p.nLPC-2);  i0>=0; i0--) 
		Acompanion[(p.nLPC+1)*i0+1] = 1.0F;	
	
	p.transDone=false;
	p.transCounter=0;

	//SC(2012/02/29) For data saving (Writing to disk)
	dataFileCnt = 0;

	//SC(2012/03/13) PVOC time warp
	for (int i0 = 0; i0 < MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64; i0++){
		for (int i1 = 0; i0 < 1024 * 2; i0++){
			pvocWarpCache[i0][i1] = 0;
		}
	}
	
}

int sw(char *table[], const char *str, const int &entries)
{
	char **bp;
	int i;
	for(i=0,bp=&table[0]; i<entries; ++i,++bp) // scan all test strings
	{
		if(*str != **bp)continue; // 1st char's must match
		if(strcmp(str,*bp)==0)return i+1; // match whole string
	}
return 0;
}


int TransShift::setparams(void * name, void * value){
	char *arg = (char *)name;
	int k,n;
	bool param_set=true;

//	char *intparamarray[NPARAMS]={"srate","framelen","ndelay","nwin","nlpc","nfmts","ntracks",
//									1		   2		3		4		5		6		7								
//	"avglen","cepswinwidth","fb","minvowellen", "delayframes", "bpitchshift", "pvocframelen", "pvochop", "bdownsampfilt"};
//		8			9		 10		11             12				13				14			15				16
	k=sw(intparamarray, arg, 16);
	//TRACE("ALGO: setting param: %s\n", arg);
	switch (k){
	case 1: // sample rate
		p.sr				= (int)(*(mytype *)value);
		break;
	case 2: // frame length = ASIO IO buffer 
		p.frameLen			= (int)(*(mytype *)value);
		p.anaLen			= p.frameShift + 2 * (p.nDelay - 1) * p.frameLen;
		//p.pvocFrameLen			= p.frameShift + (2 * p.nDelay - 3) * p.frameLen;
		//for (int i0 = 0; i0 < MAX_NFFT; i0++){
		//	fftc_ps0[i0] = 0;
		//}
		//gen_w_r2(fftc_ps0, p.pvocFrameLen);
		break;
	case 3: // number of delayed buffers : total procesing delay = frameLen*nDelay/samprate
		p.nDelay			= (int)(*(mytype *)value);
		p.anaLen			= p.frameShift + 2 * (p.nDelay - 1) * p.frameLen;
		//p.pvocFrameLen			= p.frameShift + (2 * p.nDelay - 3) * p.frameLen;
		//for (int i0 = 0; i0 < MAX_NFFT; i0++){
		//	fftc_ps0[i0] = 0;
		//}
		//gen_w_r2(fftc_ps0, p.pvocFrameLen);
		break;
	case 4: // each frame can be divided in nWin windows (=mini_frames), on which the whole process will be applied 
		p.nWin				= (int)(*(mytype *)value);
		break;
	case 5: // LPC order ( number of coeffs = nLPC +1)
		p.nLPC			    = (int)(*(mytype *)value);
		break;
	case 6: // number of formants to be shifted
		p.nFmts				= (int)(*(mytype *)value);
		break;
	case 7: // number of formants the tracking algorithm will track
		p.nTracks			= (int)(*(mytype *)value);
		break;
	case 8: // weithed moving average len (should be approx one pitch period)
		p.avgLen			= (int)(*(mytype *)value);
		break;
	case 9:	// Window width of low-pass liftering
		p.cepsWinWidth		= (int)(*(mytype *)value);
		break;
	case 10: // Auditory feedback mode
		p.fb				= (int)(*(mytype *)value);
		break;	
	case 11: // Minimum allowed vowel length
		p.minVowelLen		= (int)(*(mytype *)value);
		break;
	case 12:
		p.delayFrames		= (int)(*(mytype *)value);
		if (p.delayFrames < 0){
			TRACE("WARNING: delayFrames < 0. Set to 0 automatically.\n");
			p.delayFrames = 0;
		}
		if (p.delayFrames > MAX_DELAY_FRAMES){
			TRACE("WARNING: delayFrames > %d. Set to %d automatically.\n", MAX_DELAY_FRAMES, MAX_DELAY_FRAMES);
			p.delayFrames = MAX_DELAY_FRAMES;
		}
		break;
	case 13:
		p.bPitchShift		= (int)(*(mytype *)value);
		break;
	case 14:
		p.pvocFrameLen		= (int)(*(mytype *)value);
		for(int i0=0;i0<p.pvocFrameLen;i0++){
			hwin2[i0] = 0.5*cos(mytype(2*M_PI*i0)/mytype(p.pvocFrameLen)); 
			hwin2[i0] = 0.5 - hwin2[i0];
			//hwin2[p.pvocFrameLen-i0-1] = hwin2[i0];
		}
		for (int i0 = 0; i0 < MAX_NFFT; i0++){
			fftc_ps0[i0] = 0;
		}
		gen_w_r2(fftc_ps0, p.pvocFrameLen);
		break;
	case 15:
		p.pvocHop			= (int)(*(mytype *)value);
		break;
	case 16:
		p.bDownSampFilt		= (int)(*(mytype *)value);
		break;
	default:
		param_set=false;
	}

	if (param_set)
	{
		p.frameShift	= p.frameLen/p.nWin;	
		p.anaLen		= p.frameShift+2*(p.nDelay-1)*p.frameLen;
		time_step		= (mytype)p.frameShift*1000/p.sr;	//SC Unit: ms
		p.minDetected	= p.nWin;
		return 0;
	}
	else
		param_set=true;

//char *doubleparamarray[NPARAMS]={"scale","preemp","rmsthr","rmsratio","rmsff","dfmtsff",
//									  1        2         3       4        5        6        7        8         9         10       11     12    13    14       15    16
// "wgfreq","wgamp","wgtime","datapb",
//	  7		  8		   9		10
// "f2min","f2max","pertf2","pertamp","pertphi","f1min","f1max","lbk","lbb","triallen","ramplen",
//   11		  12      13       14       15        16      17	  18   19     20          21
// "afact","bfact","gfact","fn1","fn2", "pitchshiftratio"};
//   22	     23      24     25    26		27
	k=sw(doubleparamarray, arg, 28);
	switch (k){
	case 1: // scaling factor for IO output ( 1 = 0db)
		p.dScale			= *(mytype *)value;
		break;
	case 2: // preemphasize value
		p.dPreemp			= *(mytype *)value;
		break;
	case 3: // rms threshold value ( to detect voiced regions)
		p.dRMSThresh		= *(mytype *)value;
		break;
	case 4: // ratio threshold between preemphasized and original frame ( to detect fricatives)
		p.dRMSRatioThresh 	= *(mytype *)value;
		break;
	case 5: // rms forgetting factor (for smoothing)
		p.rmsFF 			= *(mytype *)value;
		break;
	case 6: // fmts forgetting factor (for smoothing in trajectory detection)
		p.dFmtsFF			= *(mytype *)value;
		break;
	case 7: // Frequency of the sine wave (Hz)
		p.wgFreq			= *(mytype *)value;
		break;
	case 8: // Amplitude of the sine wave (wav amp)
		p.wgAmp				= *(mytype *)value;
		break;
	case 9: // Initial time of the sine wave (equivalent to the initial phase)
		p.wgTime			= *(mytype *)value;
		break;
	case 10: // Wave data for playback
		for (n=0;n<MAX_PB_SIZE;n++){
			data_pb[n] = *((mytype *)value+n);
		}
		pbCounter=0;
		break;
	case 11:
		p.F2Min			= *(mytype *)value;
		break;
	case 12:
		p.F2Max			= *(mytype *)value;
		break;
	case 13:
		for (n=0;n<PF_NPOINTS;n++){
			p.pertF2[n] = *((mytype *)value+n);
		}
		break;	
	case 14:
		for (n=0;n<PF_NPOINTS;n++){
			p.pertAmp[n] = *((mytype *)value+n);
		}
		break;
	case 15:
		for (n=0;n<PF_NPOINTS;n++){
			p.pertPhi[n] = *((mytype *)value+n);
		}
		break;	
	case 16:
		p.F1Min			= *(mytype *)value;
		break;
	case 17:		
		p.F1Max			= *(mytype *)value;
		break;
	case 18:
		p.LBk			= *(mytype *)value;
		break;
	case 19:
		p.LBb			= *(mytype *)value;
		break;	
	case 20:
		p.trialLen			= *(mytype *)value;
		break;
	case 21:
		p.rampLen			= *(mytype *)value;
		break;
	case 22:
		p.aFact				= *(mytype *)value;
		break;	
	case 23:
		p.bFact				= *(mytype *)value;
		break;
	case 24:
		p.gFact				= *(mytype *)value;
		break;
	case 25:
		p.fn1				= *(mytype *)value;
		break;
	case 26:
		p.fn2				= *(mytype *)value;
		break;
	case 27:
		p.pitchShiftRatio   = *(mytype *)value;
		break;
	case 28:
		delete(warpCfg);
		warpCfg = new pvocWarpAtom(*((mytype *)value), 
								   *((mytype *)value + 1), 
								   *((mytype *)value + 2), 
								   *((mytype *)value + 3), 
								   *((mytype *)value + 4));
		break;
	default:
		param_set=false;
	}

	if (param_set)
		return 0;
	else
		param_set=true;

//char *boolparamarray[NPARAMS]={"bgainadapt","bshift","btrack","bdetect","bweight","bcepslift","bratioshift","bmelshift"};
//									1          2         3         4          5           6			7				8
	k=sw(boolparamarray, arg, 8);
	
	switch (k){
	case 1: // 1 = do gainadaption  
		p.bGainAdapt		=  (int)(*(mytype *)value);
		break;
	case 2: // If do shifting: set to 1
		p.bShift			=  (int)(*(mytype *)value);
		break;
	case 3: // Should almost always be set to 1
		p.bTrack			=  (int)(*(mytype *)value);
		break;
	case 4: // If do shifting: set to 1
		p.bDetect			=  (int)(*(mytype *)value);
		break;	
	case 5: // 1 = do temporal smoothing of formant frequencies based on RMS weighted averaging
		p.bWeight		=  (int)(*(mytype *)value);
		break;
	case 6:	// 1 = do low-pass liftering
		p.bCepsLift			=  (int)(*(mytype *)value);
		break;
	case 7:	
		p.bRatioShift       =  (int)(*(mytype *)value);
		break;
	case 8:
		p.bMelShift			=  (int)(*(mytype *)value);
		break;
	default:
		param_set=false;
	}

	if (param_set)
		return 0;
	else
		return 1;

}


int TransShift::getparams(void * name){
	char *arg = (char *)name;
	int k;
//char *intparamarray[NPARAMS]={"srate","framelen","ndelay","nwin","nlpc","nfmts","ntracks",
//								 1		   2		3		4		 5		 6		 7								
//	"avglen","cepswinwidth","fb"};
//		8			9		 10
	k=sw(intparamarray, arg, 10);
	TRACE("ALGO: Getting param: %s\n", arg);
	switch (k){
		case 1: // sample rate
			return p.sr;	
		case 2: //
			return p.frameLen;
		case 3: // 
			return p.nDelay;
		case 4: // 
			return p.nWin;
		case 5: // LPC order ( number of coeffs = nLPC +1)
			return p.nLPC;
		case 6: // number of formants to be shifted
			return p.nFmts;
		case 7: // 
			return p.nTracks;
		case 8: // transition length	
			return p.avgLen;			
		case 9: // avglen
			return p.cepsWinWidth;			
		case 10: // downsample factor
			return  p.fb;
		default:
			TRACE("ALGO: Unknown parameter: %s\n", arg);
			return 0;
	}
}

const mytype* TransShift::getsignal(int & size)	//SC 
{
	size = frame_counter*p.frameLen ;	//SC Update size
	return signal_recorder[0];			//SC return the
}

const mytype* TransShift::getdata(int & size, int & vecsize)
{
	size = data_counter;
	vecsize = 4 +2*p.nTracks +2*2+p.nLPC+1;
	return data_recorder[0];
}

// Assumes: MAX_BUFLEN == 3*p.frameLen 
//		    
int TransShift::myProcess(mytype *inFrame_ptr, mytype *outFrame_ptr, int frame_size)	// Sine wave generator
{
	static bool during_trans=false;
	static bool maintain_trans=false;
	bool above_rms=false;
	mytype sf1m,sf2m,loc,locfrac,mphi,mamp;
	int locint,n;

	static mytype time_elapsed=0;

	int fi=0,si=0,i0=0,offs=0,quit_hqr=0,indx=0,nZC=0;
	mytype rms_o,rms_p,rms_s,rms_ratio,wei;
	int optr; //SC(2012/02/28) DAF

	/* Temporary buffer for holding the output before duplexing into stereo (two channels) */
	mytype outputBuf[MAX_FRAMELEN];

	// ====== Variables for frequency/pitch shifting (smbPitchShift) ======
	// SC(2012/03/05)
	mytype xFrameW[MAX_BUFLEN+MAX_NLPC];
	mytype X_magn[MAX_NFFT];
	mytype X_phase[MAX_NFFT];
	mytype anaMagn[MAX_NFFT], anaFreq[MAX_NFFT];
	mytype synMagn[MAX_NFFT], synFreq[MAX_NFFT];

	/*mytype X_magn[NFFT];
	mytype X_phase[NFFT];
	mytype anaMagn[NFFT], anaFreq[NFFT];
	mytype synMagn[NFFT], synFreq[NFFT];*/
	mytype p_tmp, magn, phase;
	mytype expct, osamp, freqPerBin;
	int qpd, index;
	// ====== ~Variables for frequency/pitch shifting (smbPitchShift) ======


	if (!bFrameLenShown){
		printf("frame_size (w/o downsampling) = %d\n", frame_size);
		bFrameLenShown = true;
	}

	if (frame_size!=DOWNSAMP_FACT*p.frameLen)	//SC This should be satisfied. Just for safeguard.
		return 1;	

	for (i0=0;i0<p.nTracks;i0++)	//SC Initialize the frequencies and amplitudes of the formant tracks
	{
		fmts[i0]=0;
		amps[i0]=0;
	}

	for (i0=0;i0<2;i0++)			//SC Initialize the time derivative of formants (dFmts) and the shifted formants (sFmts)
	{
		dFmts[i0]=0;
		sFmts[i0]=0;
	}

	for (i0=0;i0<p.nLPC+1;i0++)		//SC Initialize the order LPC coefficients
	{
		lpcAi[i0]=0;
	}

	// downsample signal provided by soundcard
	//SC Notice that inFrame_ptr is the pointer to input, and inFrameBuf is the pointer to the output.
	if (p.bDownSampFilt == 1)
		downSampSig(&srfilt_b[0], &srfilt_a[0],inFrame_ptr,&srfilt_buf[0], &inFrameBuf[0],&srfilt_delay_down[0], p.frameLen , N_COEFFS_SRFILT , DOWNSAMP_FACT);
	else
		downSampSig_noFilt(&srfilt_b[0], &srfilt_a[0],inFrame_ptr,&srfilt_buf[0], &inFrameBuf[0],&srfilt_delay_down[0], p.frameLen , N_COEFFS_SRFILT , DOWNSAMP_FACT);
	//SC Output: algo.inFrame_ptr (length: p.frameLen)	

	if (p.bRecord)		//SC Record the downsampled original signal
	{// recording signal in (downsampled)
		DSPF_dp_blk_move(&inFrameBuf[0],&signal_recorder[0][frame_counter*p.frameLen],p.frameLen);	//SC Copying
	}

	// push samples in oBuf and pBuf
	//SC algo.oBuf: (downsampled) input, algo.pBuf: preemphasized (downsampled) input
	DSPF_dp_blk_move(&oBuf[p.frameLen],&oBuf[0],2*(p.nDelay-1)*p.frameLen);	
	DSPF_dp_blk_move(&pBuf[p.frameLen],&pBuf[0],2*(p.nDelay-1)*p.frameLen);	//SC Pre-emphasized buffer shift to left

	// move inFrame into oBuf
	DSPF_dp_blk_move(&inFrameBuf[0],&oBuf[2*(p.nDelay-1)*p.frameLen],p.frameLen);

	// preemphasize inFrame and move to pBuf
	//SC Preemphasis amounts to an high-pass iir filtering, the output is pBuf
	//SC(2008/05/07)
	iir_filt(&b_preemp[0], &a_preemp[0], &inFrameBuf[0],&pBuf[2*(p.nDelay-1)*p.frameLen], &preemp_delay[0],p.frameLen, 2,1);	

	// load inBuf with p.frameLen samples of pBuf
	//SC Copy pBuf to inBuf. pBuf is the signal based on which LPC and other
	// analysis will be done. inBuf is the signal that will be shifted (if in shifting mode).
	DSPF_dp_blk_move(&pBuf[(p.nDelay-1)*p.frameLen],&inBuf[0],p.frameLen);

	// move inBuf to outBuf  ( will be overwritten when filtering)
	DSPF_dp_blk_move(&inBuf[0],&outBuf[0],p.frameLen);

	for(fi=0;fi<p.nWin;fi++)// each incoming frame is divided in nwin windows
	{
		during_trans=false;
		gtot[fi]=1;			// initialize gain factor				
		si=fi*p.frameShift;// sample index
		rms_s=sqrt(DSPF_dp_vecsum_sq(&oBuf[(p.nDelay-1)*p.frameLen+si],p.frameShift)/((mytype)p.frameShift)); //short time rms of current window
		rms_o=calcRMS1(&oBuf[(p.nDelay-1)*p.frameLen+si],p.frameShift); // Smoothed RMS of original signal 
		rms_p=calcRMS2(&pBuf[(p.nDelay-1)*p.frameLen+si],p.frameShift); // RMS preemphasized (i.e., high-pass filtered) signal 	
		rms_ratio=rms_o/rms_p; // rmsratio indicates if there is a fricative around here...	

		//SC-Mod(2008/01/11)		
		//SC Notice that the identification of a voiced frame requires, 1) RMS of the orignal signal is large enough,
		//SC	2) RMS ratio between the orignal and preemphasized signals is large enough
		if (rms_o>=p.dRMSThresh*2){			
			above_rms=(ISABOVE(rms_o,p.dRMSThresh) && ISABOVE(rms_ratio,p.dRMSRatioThresh/1.3));
		}
		else{
			above_rms=(ISABOVE(rms_o,p.dRMSThresh) && ISABOVE(rms_ratio,p.dRMSRatioThresh));
		}


		if(above_rms)
		{
			if (p.bWeight)	//SC bWeight: weighted moving averaging of the formants
				wei=rms_s; // weighted moving average over short time rms
			else
				wei=1; // simple moving average

			weiVec[circ_counter]=wei; // weighting vector
			//SC the core code (getAi) is here
			getAi(&pBuf[si],&lpcAi[0],p.anaLen,p.nLPC); // get LPC coefficients
			quit_hqr=hqr_roots (&lpcAi[0],&realRoots[0],&imagRoots[0],&Acompanion[0],p.nLPC); // find the roots of polynomial
			getRPhiBw(&realRoots[0],&imagRoots[0],&amps[0],&orgPhis[0],&bw[0]); // get radius and angle of sorted! roots  
			if (p.bTrack)
				trackPhi(&amps[0],&orgPhis[0],time_elapsed); // formant tracking routine 

			getWma(&orgPhis[0],&bw[0],&wmaPhis[0],&wmaR[0]); // weighted moving average

			for (i0=0;i0<p.nTracks;i0++)						
				fmts[i0] = wmaPhis[i0]*p.sr/(2*M_PI);

			//SC Convert to mel. The perturbation field variables (F2Min, F2Max, pertF2, etc.) are all in mel. 			
			if (p.bMelShift){
				f1m=hz2mel(fmts[0]);
				f2m=hz2mel(fmts[1]);
			}
			else{
				f1m=fmts[0];
				f2m=fmts[1];
			}

			getDFmt(&fmts[0],&dFmts[0],time_elapsed); // formant derivatives (note utilized in the current version, but may be useful in the future)

			if (p.bDetect)	//SC bDetect: to detect transition?
				during_trans=detectTrans(&fmts[0],&dFmts[0],data_counter,time_elapsed); // [a] [i] transition detection
			else
				during_trans=false;

			time_elapsed+=time_step;
		}
		else
		{
			time_elapsed=0;
			weiVec[circ_counter]=0;
			for (i0=0;i0<p.nTracks;i0++)		//SC Put zeros in formant frequency and amplitude for unvoiced portion
			{
				weiMatPhi[i0][circ_counter]=0;
				weiMatBw[i0][circ_counter]=0;
			}
		}
		
		if (p.bShift)
		{					
			if(during_trans && above_rms) // Determine whether the current point in perturbation field
			{// yes : windowed deviation over x coordinate
				loc=locateF2(f2mp);	// Interpolation (linear)								
				locint=(int)floor(loc);
				locfrac=loc-locint;
				
				mamp=p.pertAmp[locint]+locfrac*(p.pertAmp[locint+1]-p.pertAmp[locint]);	// Interpolaton (linear)
				mphi=p.pertPhi[locint]+locfrac*(p.pertPhi[locint+1]-p.pertPhi[locint]);
				if (!p.bRatioShift){	// Absoluate shift					
					sf1m=f1m+mamp*cos(mphi);	// Shifting imposed
					sf2m=f2m+mamp*sin(mphi);
				}
				else{	// Ratio shift					
					mamp=p.pertAmp[locint]+locfrac*(p.pertAmp[locint+1]-p.pertAmp[locint]);
					sf1m=f1m*(1+mamp*cos(mphi));
					sf2m=f2m*(1+mamp*sin(mphi));
				}

				if (p.bMelShift){
					newPhis[0]=mel2hz(sf1m)/p.sr*2*M_PI;	// Convert back to Hz
					newPhis[1]=mel2hz(sf2m)/p.sr*2*M_PI;
				}
				else{
					newPhis[0]=sf1m/p.sr*2*M_PI;
					newPhis[1]=sf2m/p.sr*2*M_PI;
				}
			}
			else
			{// no : no force applied
				maintain_trans=false;
				during_trans=false;	
				newPhis[0]=wmaPhis[0];	// No shifting
				newPhis[1]=wmaPhis[1];
			}

			if(during_trans || maintain_trans)// directly after transition
			{
				for (i0=0;i0<2;i0++)
				{
					sFmts[i0]=newPhis[i0]*p.sr/(2*M_PI); // shifted fotmants for recording
				}

				myFilt(&inBuf[si],&outBuf[si],&wmaPhis[0],&newPhis[0],&amps[0],p.frameShift); // f1 f2 filtering 
				gtot[fi]=getGain(&amps[0],&wmaPhis[0],&newPhis[0],p.nFmts); // gain factor calculation

			}
			else// no transition
			{
	 
				if(above_rms)				
					myFilt(&inBuf[si],&fakeBuf[0],&wmaPhis[0],&wmaPhis[0],&amps[0],p.frameShift);
			}

		}
		// data recording
		data_recorder[0][data_counter]=frame_counter*p.frameLen+si+1;// matlab intervals
		data_recorder[1][data_counter]=rms_o;
		data_recorder[2][data_counter]=rms_p;
		data_recorder[3][data_counter]=rms_s;
		offs=4;

		//SC Write formant frequencies and amplitudes to data_recorder
		for (i0=0;i0<p.nTracks;i0++)
		{
			data_recorder[i0+offs][data_counter]         = fmts[i0];
			data_recorder[i0+offs+p.nTracks][data_counter] = amps[i0];
		}
		offs+=2*p.nTracks;

		//SC Write time derivative of formants and shifted formants to data_recorder
		for (i0=0;i0<2;i0++)
		{
			data_recorder[i0+offs][data_counter]		= dFmts[i0];
			data_recorder[i0+offs+2][data_counter]		= sFmts[i0];
		}
		offs+=4;

		//SC Write the LPC coefficients to data_recorder
		for (i0=0;i0<p.nLPC+1;i0++)
		{
			data_recorder[i0+offs][data_counter]		= lpcAi[i0];
		}

		data_counter++;
		circ_counter= data_counter % MAX_PITCHLEN;
	}


	// gain adaption: optional
	//SC Notice that gain adaptation is done after formant shifts
	if (p.bGainAdapt) 
		nZC=gainAdapt(&outBuf[0],&gtot[0],p.frameLen,p.frameShift); // apply gain adaption between zerocrossings
	
	// deemphasize last processed frame and send to outframe buffer
	//SC(2008/05/07)
	iir_filt(&b_deemp[0], &a_deemp[0], &outBuf[0],&outFrameBuf[outFrameBuf_circPtr], &deemp_delay[0], p.frameLen , 2, 1);
	//DSPF_dp_blk_move(&outBuf[0],&outFrameBuf[0],p.frameLen);

	// === Frequency/pitch shifting code ===
	mytype xBuf[MAX_NFFT];

	if (p.bPitchShift == 1){		// PVOC Pitch shifting
		//////////////////////////////////////////////////////////////////////
		if ((frame_counter - (p.nDelay - 1)) % (p.pvocHop / p.frameLen) != 0){
			
		}
		else{
			if (frame_counter - (p.nDelay - 1) >= p.pvocFrameLen / p.frameLen){
				expct = 2.* M_PI * (double)p.pvocHop / (double)p.pvocFrameLen;
				osamp = p.pvocFrameLen / p.pvocHop;
				freqPerBin = p.sr / (double)p.pvocFrameLen;

				if (outFrameBuf_circPtr - p.pvocFrameLen >= 0)
					DSPF_dp_blk_move(&outFrameBuf[outFrameBuf_circPtr - p.pvocFrameLen], xBuf, p.pvocFrameLen);
				else{ // Take care of wrapping-around
					DSPF_dp_blk_move(&outFrameBuf[0], &xBuf[p.pvocFrameLen - outFrameBuf_circPtr], outFrameBuf_circPtr);
					DSPF_dp_blk_move(&outFrameBuf[MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES - p.pvocFrameLen + outFrameBuf_circPtr], 
								     &xBuf[0], 
									 p.pvocFrameLen - outFrameBuf_circPtr);
				}
				DSPF_dp_vecmul(xBuf, hwin2, xFrameW, p.pvocFrameLen);

				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					ftBuf1ps[i0 * 2] = xFrameW[i0];
					ftBuf1ps[i0 * 2 + 1] = 0;
				}

				DSPF_dp_cfftr2(p.pvocFrameLen, ftBuf1ps, fftc_ps0, 1);
				bit_rev(ftBuf1ps, p.pvocFrameLen);

				for (i0=0; i0 < p.pvocFrameLen; i0++){
					X_magn[i0] = 2. * sqrt(ftBuf1ps[i0*2] * ftBuf1ps[i0*2] + ftBuf1ps[i0*2+1] * ftBuf1ps[i0*2+1]);
					X_phase[i0] = atan2(ftBuf1ps[i0*2+1], ftBuf1ps[i0*2]);
				}
				
				for (i0=0; i0 <= p.pvocFrameLen / 2; i0++){
					p_tmp = X_phase[i0] - lastPhase[i0];

					lastPhase[i0] = X_phase[i0];

					p_tmp -= (mytype)i0 * expct;

					qpd = (int)(p_tmp / M_PI);

					if (qpd >= 0) qpd += qpd&1;
					else qpd -= qpd&1;
					p_tmp -= M_PI * (mytype)qpd;

					p_tmp = osamp * p_tmp / (2. * M_PI);
					p_tmp = (mytype)i0 * freqPerBin + p_tmp * freqPerBin;

					anaMagn[i0] = X_magn[i0];
					anaFreq[i0] = p_tmp;
				}

				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					synMagn[i0] = 0.;
					synFreq[i0] = 0.;
				}

				for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++){
					index = (int)(i0 / p.pitchShiftRatio);

					if (index <= p.pvocFrameLen /2) {
						synMagn[i0] += anaMagn[index];
						synFreq[i0] = anaFreq[index] * p.pitchShiftRatio;
					}
				}

				// --- Synthesis ---
				for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++) {
					magn = synMagn[i0]; // get magnitude and true frequency from synthesis arrays
					p_tmp = synFreq[i0];
					
					p_tmp -= (double)i0 * freqPerBin;	// subtract bin mid frequency		
					p_tmp /= freqPerBin;	// get bin deviation from freq deviation
					p_tmp = 2. * M_PI * p_tmp / osamp;	// take osamp into account
					p_tmp += (double)i0 * expct;		// add the overlap phase advance back in
					
					sumPhase[i0] += p_tmp;		// accumulate delta phase to get bin phase
					phase = sumPhase[i0];

					/* get real and imag part and re-interleave */
					ftBuf2ps[2 * i0] = magn * cos(phase);
					ftBuf2ps[2 * i0 + 1] = magn * sin(phase);			// What causes the sign reversal here?
				}

				for (i0 = p.pvocFrameLen + 1; i0 < p.pvocFrameLen * 2; i0++){
					ftBuf2ps[i0] = 0.;
				}
				
				DSPF_dp_icfftr2(p.pvocFrameLen, ftBuf2ps, fftc_ps0, 1);
				bit_rev(ftBuf2ps, p.pvocFrameLen);
				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					ftBuf2ps[i0 * 2] /= (p.pvocFrameLen);
					ftBuf2ps[i0 * 2 + 1] /= (p.pvocFrameLen);
				}

				// --- Accumulate to buffer ---
				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					outFrameBufPS[(outFrameBuf_circPtr + i0) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES)] = 
						outFrameBufPS[(outFrameBuf_circPtr + i0) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES)] + 
						2 * ftBuf2ps[2 * i0] * hwin2[i0] / (osamp / 2);
				}
				for (i0 = 1; i0 <= p.pvocHop; i0++){
					outFrameBufPS[(outFrameBuf_circPtr - p.delayFrames * p.frameLen - i0) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES)] = 0;
				}

			}
		}
		//////////////////////////////////////////////////////////////////////
	}
	else if (p.bPitchShift == 2){	// PVOC Time warp
		//////////////////////////////////////////////////////////////////////
		if ((frame_counter - (p.nDelay - 1)) % (p.pvocHop / p.frameLen) != 0){
			
		}
		else{
			if (frame_counter - (p.nDelay - 1) >= p.pvocFrameLen / p.frameLen){
				expct = 2.* M_PI * (double)p.pvocHop / (double)p.pvocFrameLen;
				osamp = p.pvocFrameLen / p.pvocHop;
				freqPerBin = p.sr / (double)p.pvocFrameLen;
				mytype dp;

				if (outFrameBuf_circPtr - p.pvocFrameLen >= 0)
					DSPF_dp_blk_move(&outFrameBuf[outFrameBuf_circPtr - p.pvocFrameLen], xBuf, p.pvocFrameLen);
				else{ // Take care of wrapping-around
					DSPF_dp_blk_move(&outFrameBuf[0], &xBuf[p.pvocFrameLen - outFrameBuf_circPtr], outFrameBuf_circPtr);
					DSPF_dp_blk_move(&outFrameBuf[MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES - p.pvocFrameLen + outFrameBuf_circPtr], 
								     &xBuf[0], 
									 p.pvocFrameLen - outFrameBuf_circPtr);
				}
				DSPF_dp_vecmul(xBuf, hwin2, xFrameW, p.pvocFrameLen);

				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					ftBuf1ps[i0 * 2] = xFrameW[i0];
					ftBuf1ps[i0 * 2 + 1] = 0;
				}

				DSPF_dp_cfftr2(p.pvocFrameLen, ftBuf1ps, fftc_ps0, 1);
				bit_rev(ftBuf1ps, p.pvocFrameLen);

				for (i0=0; i0 < p.pvocFrameLen; i0++){
					X_magn[i0] = 2. * sqrt(ftBuf1ps[i0*2] * ftBuf1ps[i0*2] + ftBuf1ps[i0*2+1] * ftBuf1ps[i0*2+1]);
					X_phase[i0] = atan2(ftBuf1ps[i0*2+1], ftBuf1ps[i0*2]);
				}
				
				// --- Time warping ---
				mytype t0 = (mytype)(frame_counter - (p.nDelay - 1)) * p.frameLen / p.sr;
				mytype t1;
				int cidx0 = (frame_counter - (p.nDelay - 1)) / (p.pvocHop / p.frameLen);
				mytype cidx1_d, cidx1_f;
				int cidx1;

				
				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					pvocWarpCache[cidx0 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0] = X_magn[i0];
					pvocWarpCache[cidx0 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0 + p.pvocFrameLen] = X_phase[i0];
				}

				if (frame_counter - (p.nDelay - 1) == p.pvocFrameLen / p.frameLen){
					for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++){
						lastPhase[i0] = X_phase[i0];
					}
				}

				if (t0 < warpCfg->tBegin || t0 >= warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold + warpCfg->dur2){
					cidx1 = cidx0;

					for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++){
						magn = pvocWarpCache[cidx1 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0];
						phase = pvocWarpCache[cidx1 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0 + p.pvocFrameLen];

						ftBuf2ps[2 * i0] = magn * cos(phase);
						ftBuf2ps[2 * i0 + 1] = magn * sin(phase);

						lastPhase[i0] = phase;						
					}
				}
				else{
					if (t0 < warpCfg->tBegin + warpCfg->dur1){
						t1 = (t0 - warpCfg->tBegin) * warpCfg->rate1 + warpCfg->tBegin;
					}
					else if (t0 < warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold){
						t1 = warpCfg->rate1 * warpCfg->dur1 - warpCfg->dur1 + t0;
					}
					else if (t0 < warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold + warpCfg->dur2){
						//t1 = (t0 - (warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold)) * warpCfg->rate2 + warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold;
						t1 = warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold + warpCfg->dur2;
						t1 -= (warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold + warpCfg->dur2 - t0) * warpCfg->rate2;
					}

					cidx1_d = t1 * (mytype)p.sr / (mytype)p.pvocHop;
					cidx1 = (int)floor(cidx1_d);
					cidx1_f = cidx1_d - (mytype)cidx1;

					// [Pure interpolation]
					/* for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++){
						magn = pvocWarpCache[cidx1 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0] * (1 - cidx1_f) + 
							   pvocWarpCache[(cidx1 + 1) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0] * cidx1_f;
						phase = pvocWarpCache[cidx1 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0 + p.pvocFrameLen] * (1 - cidx1_f) + 
								pvocWarpCache[(cidx1 + 1) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0 + p.pvocFrameLen] * cidx1_f;

						ftBuf2ps[2 * i0] = magn * cos(phase);
						ftBuf2ps[2 * i0 + 1] = magn * sin(phase);
					} */

					for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++){
						magn = pvocWarpCache[cidx1 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0] * (1 - cidx1_f) + 
							   pvocWarpCache[(cidx1 + 1) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0] * cidx1_f;

						dp = pvocWarpCache[(cidx1 + 1) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0 + p.pvocFrameLen] - 
							 pvocWarpCache[cidx1 % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES / 64)][i0 + p.pvocFrameLen];
						dp -= (mytype)i0 * expct;
						
						qpd = (int)(dp / M_PI);
						if (qpd >= 0) qpd += qpd & 1;
						else qpd -= qpd & 1;
						dp -= M_PI * (mytype)qpd;

						ftBuf2ps[2 * i0] = magn * cos(lastPhase[i0]);
						ftBuf2ps[2 * i0 + 1] = magn * sin(lastPhase[i0]);

						lastPhase[i0] += (mytype)i0 * expct + dp;
					}
				}

				// --- Synthesis ---
				/* for (i0 = 0; i0 <= p.pvocFrameLen / 2; i0++){
					ftBuf2ps[2 * i0] = magn * cos(lastPhase[i0]);
					ftBuf2ps[2 * i0 + 1] = magn * sin(lastPhase[i0]);
				} */
				for (i0 = p.pvocFrameLen + 1; i0 < p.pvocFrameLen * 2; i0++){
					ftBuf2ps[i0] = 0.;
				}
				
				DSPF_dp_icfftr2(p.pvocFrameLen, ftBuf2ps, fftc_ps0, 1);
				bit_rev(ftBuf2ps, p.pvocFrameLen);
				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					ftBuf2ps[i0 * 2] /= (p.pvocFrameLen);
					ftBuf2ps[i0 * 2 + 1] /= (p.pvocFrameLen);
				}

				// --- Accumulate to buffer ---
				for (i0 = 0; i0 < p.pvocFrameLen; i0++){
					outFrameBufPS[(outFrameBuf_circPtr + i0) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES)] = 
						outFrameBufPS[(outFrameBuf_circPtr + i0) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES)] + 
						2 * ftBuf2ps[2 * i0] * hwin2[i0] / (osamp / 2);
				}
				for (i0 = 1; i0 <= p.pvocHop; i0++){
					outFrameBufPS[(outFrameBuf_circPtr - p.delayFrames * p.frameLen - i0) % (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES)] = 0;
				}
				
				

			}
		}
		//////////////////////////////////////////////////////////////////////
	}
	else{
		DSPF_dp_blk_move(&outFrameBuf[outFrameBuf_circPtr], &outFrameBufPS[outFrameBuf_circPtr], p.frameLen);
	}

	// === ~Frequency/pitch shifting code ===

	//SC(2012/02/28) For DAF: shift the data in outFrameBuf
	if (outFrameBuf_circPtr - p.delayFrames * p.frameLen > 0){
		optr = outFrameBuf_circPtr - p.delayFrames * p.frameLen;
	}
	else{
		optr = outFrameBuf_circPtr - p.delayFrames * p.frameLen + (MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES);
	}

	if (p.bRecord)
	{// recording signal output
		//DSPF_dp_blk_move(&outFrameBuf[0],&signal_recorder[1][frame_counter*p.frameLen],p.frameLen); 
		//SC(2012/02/28) For DAF: use a past frame for playback, if p.delayFrames > 0
		//DSPF_dp_blk_move(&outFrameBuf[p.delayFrames * p.frameLen], &signal_recorder[1][frame_counter*p.frameLen], p.frameLen);
		//DSPF_dp_blk_move(&outFrameBuf[optr], &signal_recorder[1][frame_counter*p.frameLen], p.frameLen);
		DSPF_dp_blk_move(&outFrameBufPS[optr], &signal_recorder[1][frame_counter*p.frameLen], p.frameLen);
	}

	// Upsample signal, scale and send to sound card buffer
	//upSampSig(&srfilt_b[0], &srfilt_a[0],&outFrameBuf[0],&srfilt_buf[0], &outFrame_ptr[0],&srfilt_delay_up[0], p.frameLen * DOWNSAMP_FACT, N_COEFFS_SRFILT ,DOWNSAMP_FACT,p.dScale);
	//SC(2012/02/28) For DAF: use a past frame for playback, if p.delayFrames > 0
	//upSampSig(&srfilt_b[0], &srfilt_a[0], &outFrameBuf[optr], 
	//	      &srfilt_buf[0], &outFrame_ptr[0], &srfilt_delay_up[0], p.frameLen * DOWNSAMP_FACT, 
	//		  N_COEFFS_SRFILT ,DOWNSAMP_FACT, p.dScale);
	/*upSampSig(&srfilt_b[0], &srfilt_a[0], &outFrameBufPS[optr], 
		      &srfilt_buf[0], &outFrame_ptr[0], &srfilt_delay_up[0], p.frameLen * DOWNSAMP_FACT, 
			  N_COEFFS_SRFILT ,DOWNSAMP_FACT, p.dScale);*/
	upSampSig(&srfilt_b[0], &srfilt_a[0], &outFrameBufPS[optr], 
		      &srfilt_buf[0], outputBuf, &srfilt_delay_up[0], p.frameLen * DOWNSAMP_FACT, 
			  N_COEFFS_SRFILT ,DOWNSAMP_FACT, p.dScale);


	outFrameBuf_circPtr += p.frameLen;
	if (outFrameBuf_circPtr >= MAX_FRAMELEN * DOWNSAMP_FACT * MAX_DELAY_FRAMES){
		outFrameBuf_circPtr = 0;
	}
	
	//SC(2008/06/20) Blend with noise, if fb=2 or 3; Mute, if fb=0;
	if (p.fb==3){	// voice + noise
		for(n=0;n<frame_size;n++){
			/*outFrame_ptr[n]=outFrame_ptr[n]+data_pb[pbCounter];*/
			outputBuf[n] = outputBuf[n] + data_pb[pbCounter];
			pbCounter=pbCounter+1;
			if (pbCounter==MAX_PB_SIZE){
				pbCounter=0;
			}
		}
	}
	else if (p.fb==2){	// noise only
		for(n=0;n<frame_size;n++){
			/*outFrame_ptr[n]=data_pb[pbCounter];*/
			outputBuf[n] = data_pb[pbCounter];
			pbCounter=pbCounter+1;
			if (pbCounter==MAX_PB_SIZE){
				pbCounter=0;
			}
		}
	}
	else if (p.fb==0){	// Mute
		for(n=0;n<frame_size;n++){
			/*outFrame_ptr[n]=0;*/
			outputBuf[n] = 0;
		}
	}

	/* Duplex into stereo */
	for (n = 0; n < frame_size; n++)
		outFrame_ptr[n * 2] = outFrame_ptr[n * 2 + 1] = outputBuf[n];

	//SC(2008/06/22) Impose the onset and offset ramps, mainly to avoid the unpleasant "clicks" at the beginning and end
	/* if ((mytype)frame_counter*(mytype)p.frameLen/(mytype)p.sr>p.trialLen){
		for (n=0;n<frame_size;n++)
			outFrame_ptr[n]=0;
	}
	if ((mytype)frame_counter*(mytype)p.frameLen/(mytype)p.sr<=p.rampLen){
		for (n=0;n<frame_size;n++)
			outFrame_ptr[n]=outFrame_ptr[n]/p.rampLen*((mytype)((frame_counter-1)*frame_size+n)/(mytype)p.sr/DOWNSAMP_FACT);
	}
	if ((mytype)frame_counter*(mytype)p.frameLen/(mytype)p.sr>=p.trialLen-p.rampLen){
		for (n=0;n<frame_size;n++)
	 		outFrame_ptr[n]=outFrame_ptr[n]/p.rampLen*(p.trialLen-(mytype)((frame_counter-1)*frame_size+n)/(mytype)p.sr/DOWNSAMP_FACT);
	} */

	frame_counter++;

	if (((frame_counter+1)*p.frameLen>=MAX_REC_SIZE) || ((data_counter+p.nWin)>=MAX_DATA_SIZE)) //avoid segmentation violation
	{
		//SC(2012/02/29) Save signal to file
		/* char dataFileName[200];
		sprintf(dataFileName, "data%.3d.bin", dataFileCnt);
		printf("Writing data to file %s: size = %d (%d bytes)\n", dataFileName, MAX_REC_SIZE, MAX_REC_SIZE * (sizeof mytype));

		ofstream dataFileCl(dataFileName);
		dataFileCl.close();

		fstream dataFile(dataFileName, ios::binary | ios::out);
		if (!dataFile){
			printf("WARNING: Cannot open file %s\n", dataFileName);
		}

		dataFile.write((char *) signal_recorder, MAX_REC_SIZE * (sizeof mytype));
		dataFile.close(); */
		
		frame_counter=0;
		data_counter=0;
		dataFileCnt ++;


		//SC(2012/02/29) Writing to binary files
		
	}

	return 0;
}

int TransShift::myProcessSineGen(mytype *inFrame_ptr, mytype *outFrame_ptr, int frame_size)	// Sine wave (pure tone) generator
{
	int n;
	double dt;

	dt=0.00002083333333333333;//((double)p.sr)*((double)DOWNSAMP_FACT);

	for(n=0;n<frame_size;n++){
		outFrame_ptr[n * 2] = outFrame_ptr[n * 2 + 1] = p.wgAmp*sin(2*M_PI*p.wgFreq*p.wgTime);
		p.wgTime=p.wgTime+dt;
	}

	return 0;
}

int TransShift::myProcessWavePB(mytype *inFrame_ptr, mytype *outFrame_ptr, int frame_size)	// Wave playback
{	// Wave playback
	int n;

	for(n=0;n<frame_size;n++){
		outFrame_ptr[n* 2] = outFrame_ptr[n * 2 + 1] = data_pb[pbCounter];
		pbCounter=pbCounter+1;
		if (pbCounter==MAX_PB_SIZE)
			pbCounter=0;
		
	}

	return 0;
}


int TransShift::getWma(mytype *phi_ptr, mytype *bw_ptr , mytype * wmaPhi_ptr, mytype * wmaR_ptr)
{// computational efficient weithed moving average 

	int i0=0;
	int circ_indx_sub=0;


	circ_indx_sub=(data_counter-p.avgLen) % MAX_PITCHLEN; // points to the data to withdraw from sum 

	sumWei+=weiVec[circ_counter]-weiVec[circ_indx_sub];  // update weighting sum


	for (i0=0;i0<p.nTracks;i0++)
		{						
			weiMatPhi[i0][circ_counter]=weiVec[circ_counter]*phi_ptr[i0];
			weiMatBw[i0][circ_counter]=weiVec[circ_counter]*bw_ptr[i0];

			sumWeiPhi[i0]+=weiMatPhi[i0][circ_counter]-weiMatPhi[i0][circ_indx_sub];

			sumWeiBw[i0]+=weiMatBw[i0][circ_counter]-weiMatBw[i0][circ_indx_sub];

			if (sumWei>0.0000001)
			{
				wmaPhi_ptr[i0]=sumWeiPhi[i0]/sumWei;
				wmaR_ptr[i0]=pow(10,(-(sumWeiBw[i0]/sumWei)*M_PI/(mytype(p.sr))));
			}
			else
				return 1;

		}	

	return 0;
}
int TransShift::getDFmt(mytype *fmt_ptr,mytype *dFmt_ptr, mytype time)
{// calculates the formant derivatives [Hz/ms] smoothed with forgetting factor 
	int i0=0;

	for(i0=0;i0<2;i0++)
	{

		if (time< 20)
			deltaFmt[i0]=0;
		else
			deltaFmt[i0]=(fmt_ptr[i0]-lastFmt[i0])/time_step;//0.5625;

		if (fabs(deltaFmt[i0])>p.maxDelta)
			deltaFmt[i0]=p.maxDelta*(mytype)sign(deltaFmt[i0]);

		deltaMaFmt[i0]=(1-p.dFmtsFF)*deltaFmt[i0]+p.dFmtsFF*deltaMaFmt[i0];
		dFmt_ptr[i0]=deltaMaFmt[i0];
		lastFmt[i0]=fmt_ptr[i0];
	
	}

	return 0; 
}


bool TransShift::detectTrans(mytype *fmt_ptr, mytype *dFmt_ptr,int datcnt, mytype time)
{ // detects a-i transition using f1 and f2 derivatives in time	
	static bool btransition=false;		//SC Note these are static addresses. Output argument: btransition.
	//static bool bIsTrans=false;			//SC bIsTrans is used as an internal state variable.
	//static int detect_counter = 0;			

	if (p.bShift){
		if (!p.transDone){
			f2mp=f2m;
			if (f2mp>=p.F2Min && f2mp<=p.F2Max && f1m>=p.F1Min && f1m<=p.F1Max 
				&& ((p.LBk<=0 && f1m*p.LBk+p.LBb<=f2m) || (p.LBk>0 && f1m*p.LBk+p.LBb>=f2m))){
				btransition=true;
				p.transCounter++;
			}
			else{
				btransition=false;				
				if (p.transCounter>=p.minVowelLen){
					p.transDone=true;
				}
				p.transCounter=0;
			}			
		}
		else{
			btransition=false;
			p.transCounter=0;
		}
	}
	else{
		btransition=false;
		p.transCounter=0;
	}	
	return btransition;
}



void TransShift::getAi(mytype* xx, mytype* aa,const int size,const int nlpc)	// Autocorrelation-based LPC analysis
//SC Input arguments
//SC	xx: input buffer pointer
//SC	aa: LPC coefficients pointer
//SC	size: 
//SC	nlpc: LPC order, i.e., number of LPC coefficients
//SC(2008/05/06): Incoroporating cepstral lifting
{// performs the LPC analysis on a given frame... returns the lpc coefficients

	int i0, nlpcplus1;
	int	size1;	//Debug
	mytype temp_frame[MAX_BUFLEN+MAX_NLPC];	// Utility buffer for various filtering operations
	mytype R[MAX_NLPC];					// lpc fit Autocorrelation estimate
	
	mytype a[16]={0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0};

	nlpcplus1 = nlpc+1;		//SC nlpc: order of LPC

	for(i0=0;i0<nlpcplus1;i0++)
		temp_frame[i0] = 0;
	// Window input
	//SC vecmu1: vector multiply
	//SC	hwin: a Hanning window
	DSPF_dp_vecmul(xx,hwin,&temp_frame[nlpcplus1],size);	//SC Apply a Hanning window

	////////////////////////////////////////////////////
	//SC(2008/05/07) ----- Cepstral lifting -----
	if (p.bCepsLift){
		for (i0=0;i0<NFFT;i0++){
			if (i0<size){
				ftBuf1[i0*2]=temp_frame[nlpcplus1+i0];
				ftBuf1[i0*2+1]=0;
			}
			else{
				ftBuf1[i0*2]=0;
				ftBuf1[i0*2+1]=0;
			}
		}
		DSPF_dp_cfftr2(NFFT,ftBuf1,fftc,1);	
		bit_rev(ftBuf1,NFFT);
		// Now ftBuf1 is X
		for (i0=0;i0<NFFT;i0++){
			if (i0<=NFFT/2){
				ftBuf2[i0*2]=log(sqrt(ftBuf1[i0*2]*ftBuf1[i0*2]+ftBuf1[i0*2+1]*ftBuf1[i0*2+1]));	// Optimize
				ftBuf2[i0*2+1]=0;
			}
			else{
				ftBuf2[i0*2]=ftBuf2[(NFFT-i0)*2];
				ftBuf2[i0*2+1]=0;
			}
		}
		DSPF_dp_icfftr2(NFFT,ftBuf2,fftc,1);
		bit_rev(ftBuf2,NFFT);
		// Now ftBuf2 is Xceps: the cepstrum
		for (i0=0;i0<NFFT;i0++){
			if (i0<p.cepsWinWidth || i0>NFFT-p.cepsWinWidth){			// Adjust! 
				ftBuf1[i0*2]=ftBuf2[i0*2]/NFFT;		// Normlize the result of the previous IFFT
				ftBuf1[i0*2+1]=0;
			}
			else{
				ftBuf1[i0*2]=0;
				ftBuf1[i0*2+1]=0;
			}
		}
		// Now ftBuf1 is Xcepw: the windowed cepstrum
		DSPF_dp_cfftr2(NFFT,ftBuf1,fftc,1);
		bit_rev(ftBuf1,NFFT);
		for (i0=0;i0<NFFT;i0++){
			if (i0<=NFFT/2){
				ftBuf2[i0*2]=exp(ftBuf1[i0*2]);
				ftBuf2[i0*2+1]=0;
			}
			else{
				ftBuf2[i0*2]=ftBuf2[(NFFT-i0)*2];
				ftBuf2[i0*2+1]=0;
			}
		}
		DSPF_dp_icfftr2(NFFT,ftBuf2,fftc,1);	// Need normalization
		bit_rev(ftBuf2,NFFT);
		
		size1=size;
		for (i0=0;i0<size1/2;i0++){
			temp_frame[nlpcplus1+size1/2+i0]=ftBuf2[i0*2]/NFFT;
		}
		for (i0=1;i0<size1/2;i0++){
			temp_frame[nlpcplus1+size1/2-i0]=ftBuf2[i0*2]/NFFT;
		}
		temp_frame[nlpcplus1]=0;
	}
	//~SC(2008/05/07) ----- Cepstral lifting -----
	///////////////////////////////////////////////////////////

	// Find autocorrelation values
	DSPF_dp_autocor(R, temp_frame, size, nlpcplus1);		//SC Get LPC coefficients by autocorrelation

	// Get unbiased autocorrelation
	for(i0=0;i0<nlpcplus1;i0++)
		R[i0] /= size;

	// levinson recursion
	levinson(R, aa, nlpcplus1);
}


mytype TransShift::getGain(mytype * r, mytype * ophi,mytype * sphi, int nfmts)
{// this routine calculates the gain factor used by the routine gainAdapt to compensate the formant shift
//SC-Mod(2008/01/05) arguments xiff0 and xiff1: for intensity correction during formant shifts
	int i0;
	mytype prod=1;
	
	for (i0 = 0; i0 < nfmts; i0++)
		prod  *= (1-2*r[i0]*cos(sphi[i0])+pow(r[i0],2.0))/(1-2*r[i0]*cos(ophi[i0])+pow(r[i0],2.0));
	
	return prod;
}

int TransShift::gainAdapt(mytype *buffer,mytype *gtot_ptr,int framelen, int frameshift)
{// this routine applies the gain factors collected in the vector gTot to the signal
 // for each window (nwin windows in one frame) the function finds the first zerocrossing
 // and updates the gain factor used to scale this window
 // gainadaption is done before deemphasis for two reasons:
	//1: more zerocrossings in the preemphasized signal
	//2: abrupt gain changings are more likely to introduce high frequency noise which will be filtered by the deemphasis filter
 // if no zerocrossing is found in the current window the gain factor will not be updated in this frame
 // assumtpion: zerocrossing is present in nearly each window (can be verified : the value returned by this fuction should equal nWin)
 // gain factor only changes slowely
	int i,index=0,updated=0;
	static mytype gain=1;// 1=no gain adaption =0dB
	static mytype lastSample=buffer[0];// used for continuous zero-crossing finding

	bool armed=true;

	if (lastSample * buffer[0] <= 0)// zero crossing at first sample ?
		{
			gain=gtot_ptr[0]; // update gain
			armed=false; // zero finding disabled untill next win
			updated++; // number of times gain has been updated (should be equal to nwin)
		}


	for (i=0;i<framelen-1;i++)// scan complete frame 
	{
		
		buffer[i]=gain * buffer[i];// apply gain 

		if (i>=(index+1)*frameshift-1) // next win ?
		{
			index++; // next win
			armed=true; // enable next zero search
		}


		if (armed && (buffer[i] * buffer[i+1] <= 0))// zero crossing in win ? 
		{
			gain=gtot_ptr[index]; // update gain
			armed=false; // zero finding disabled untill next win
			updated++; // number of times gain has been updated (should be equal to nwin)
		}

	
	}
	lastSample=buffer[framelen-1]; //store last sample for next function call
	buffer[framelen-1]=gain * buffer[framelen-1];//last sample
	return updated;

	
}

void TransShift::levinson(mytype *R, mytype* aa, int size)
{// Levinson recursion, used in the LPC analysis
	mytype ki, t;
    mytype E = R[0]; 
    int   i0, j0; 

	if (R[0] == 0.0)
	{
		for(i0=1; i0<size; i0++)  
			aa[i0] = 0.0;
		aa[0] = 1;
		return;
	}

    for(i0=1; i0<size; i0++)  
    { 
        ki = R[i0]; 
      
        // Update reflection coefficient: 
        for (j0=1; j0<i0; j0++) 
			ki += aa[j0] * R[i0-j0]; 
      
        ki   /= -E; 
        E    *= (1 - ki*ki); 
        
        // Update polynomial: 
        for (j0=1; j0<=RSL(i0-1,1); j0++)
        { 
            t = aa[j0]; 
            aa[j0]    += ki * aa[i0-j0]; 
            aa[i0-j0] += ki * t; 
        } 
  
        if (i0%2 == 0) aa[RSL(i0,1)] *= 1+ki; 
  
        // Record reflection coefficient
        aa[i0] = ki; 

    } // end of for loop
    aa[0] = 1.0;
} 



/*  The following takes in a polynomial stored in *c, and yields the roots of
	this polynomial (*wr stores the real comp, and *wi stores the imag comp)
	It forms a companion matrix, then uses the hqr algorithm 
	VV 19 June 2003 */

int TransShift::hqr_roots (	mytype *c, 	mytype *wr, mytype *wi,	mytype *Acompanion, const int nLPC)
{
	int nn,m,l,k,j,i,its,mmin,nLPC_SQR=nLPC*nLPC;
    mytype AHess[MAX_NLPC_SQR];   
    mytype z,y,x,w,v,u,t,s,r,q,p,anorm = 0.0F;
        
/*  generate companion matrix, starting off with an intialized version */
//	memcpy (&AHess[0],&Acompanion[0], (LPC_ORDER_SQR)*sizeof(mytype)); 
	
	DSPF_dp_blk_move(&Acompanion[0],&AHess[0],nLPC_SQR);


	for (i=0; i<nLPC; i++)  AHess[nLPC*i]= -c[i+1]; 

    /* end of companion matrix generation  */

    /* the following performs the hqr algoritm  */
    /* NOTE:  This was taken from Numerical Recipes in C, with modification
       Specifically, the wr and wi arrays were assumed to number from 1..n
       in the book.  I modified calls to these arrays so that they number 0..n-1
       Additionally, n (the order of the polynomial) is hardset to be 8 
       VV 19 June 2003 */

    for (i=1;i<(nLPC+1);i++)
        for (j=IMAX(i-1,1);j<(nLPC+1);j++)
            anorm += fabs(aMat(i,j));
    nn = nLPC;
    t=0.0;
    while (nn >= 1) {
        its=0;
        do {
           for (l=nn;l>=2;l--) {
               s=fabs(aMat(l-1,l-1))+fabs(aMat(l,l));
               if (s == 0.0) s=anorm;
               if ((mytype)(fabs(aMat(l,l-1)) + s) == s) break;
           }
         x=aMat(nn,nn);
         if (l == nn) {
                wr[(-1)+ nn]=x+t;
                wi[(-1) + nn--]=0.0;
             } else {
                    y=aMat(nn-1,nn-1);
                    w=aMat(nn,nn-1)*aMat(nn-1,nn);
                    if (l == (nn-1)) {
                        p=0.5*(y-x);
                        q=p*p+w;
                        z=sqrt(fabs(q));
                        x += t;
                        if (q >= 0.0) {
                              z=p+SIGN(z,p);
                              wr[(-1) + nn-1]=wr[(-1)+ nn]=x+z;
                              if (z) wr[(-1) + nn]=x-w/z;
                              wi[(-1) + nn-1]=wi[(-1) + nn]=0.0;
                           } else {
                                 wr[(-1) + nn-1]=wr[(-1) + nn]=x+p;
                                wi[(-1) + nn-1]= -(wi[(-1) + nn]=z);
                           }
                        nn -= 2;
                    } else {
                        if (its == 10 || its == 20) {
                             t += x;
                             for (i=1;i<=nn;i++) aMat(i,i) -= x;
                             s=fabs(aMat(nn,nn-1))+fabs(aMat(nn-1,nn-2));
                             y=x=0.75*s;
                             w = -0.4375*s*s;
                        }
                        ++its;
                        for (m=(nn-2);m>=l;m--) {
                             z=aMat(m,m);
                             r=x-z;
                             s=y-z;
                             p=(r*s-w)/aMat(m+1,m)+aMat(m,m+1);
                             q=aMat(m+1,m+1)-z-r-s;
                             r=aMat(m+2,m+1);
                             s=fabs(p)+fabs(q)+fabs(r);
                             p /= s;
                             q /= s;
                             r /= s;
                             if (m == l) break;
                             u=fabs(aMat(m,m-1))*(fabs(q)+fabs(r));
                             v=fabs(p)*(fabs(aMat(m-1,m-1))+fabs(z)+fabs(aMat(m+1,m+1)));
                             if ((mytype)(u+v) == v) break;
                        }
                        for (i=m+2;i<=nn;i++) {
                             aMat(i,i-2)=0.0F;
                             if (i != (m+2)) aMat(i,i-3)=0.0F;
                        }
                        for (k=m;k<=nn-1;k++) {
                             if (k != m) {
                                p=aMat(k,k-1);
                                q=aMat(k+1,k-1);
                                r=0.0F;
                                if (k != (nn-1)) r=aMat(k+2,k-1);
                                if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
                                   p /= x;
                                   q /= x;
                                   r /= x;
                                }
                             }
                             if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
                                 if (k == m) {
                                     if (l != m)
                                     aMat(k,k-1) = -aMat(k,k-1);
                                 } else
                                     aMat(k,k-1) = -s*x;
                                 p += s;
                                 x=p/s;
                                 y=q/s;
                                 z=r/s;
                                 q /= p;
                                 r /= p;
                                 for (j=k;j<=nn;j++) {
                                     p=aMat(k,j)+q*aMat(k+1,j);
                                     if (k != (nn-1)) {
                                           p += r*aMat(k+2,j);
                                           aMat(k+2,j) -= p*z;
                                     }
                                     aMat(k+1,j) -= p*y;
                                     aMat(k,j) -= p*x;
                                 }
                                 mmin = nn<k+3 ? nn : k+3;
                                 for (i=l;i<=mmin;i++) {
                                     p=x*aMat(i,k)+y*aMat(i,k+1);
                                     if (k != (nn-1)) {
                                         p += z*aMat(i,k+2);
                                         aMat(i,k+2) -= p*r;
                                     }
                                     aMat(i,k+1) -= p*q;
                                     aMat(i,k) -= p;
                                 }
                             }
                        }
                 }
             }
        } while (l < nn-1);
     }
	if (nn == 0) return 1;
		else return 0;
}


void TransShift::getRPhiBw (mytype *wr,  mytype *wi, mytype *radius,  mytype *phi ,mytype *bandwith)
{// 

  /* The following sorts the roots in wr and wi.  It is adapted from a matlab script that Reiner wrote */
  mytype	arc[MAX_NLPC], arc2[MAX_NLPC], wr2[MAX_NLPC], wi2[MAX_NLPC];
  mytype	wreal, wimag, warc, wmag, mag[MAX_NLPC], mag2[MAX_NLPC];
  int		numroots, i0, j0, nmark;
  
  /* calculate angles for all defined roots */

  numroots = 0;  

  for (i0=0; i0<p.nLPC; i0++)
	{          
    arc[i0] = atan2(wi[i0],wr[i0]);
    mag[i0] = sqrt(wi[i0]*wi[i0] + wr[i0]*wr[i0]);
    if  ( /*(arc[i0] > F1_min) && */ (wi[i0]>0) /*&& 
    	 (mag[i0] > 0.9) && (mag[i0] < 1.0) */ )
        /* only store positive arc root of conjugate pairs */

        {
            mag2[numroots]	 = mag[i0];
            arc2[numroots]   = arc[i0];  
            wr2[numroots]    = wr[i0];
            wi2[numroots++]  = wi[i0];
         }
  }

  /* sort according to arc using a stupid sort algorithm. */
 
  for (i0=0; i0<numroots; i0++)  /* look for minimal first */
    {
      nmark = i0;
      for (j0=i0+1; j0<numroots; j0++)  /* find smallest arc (frequency) */
            if (arc2[j0] < arc2[nmark]) nmark = j0;
      if (nmark != i0) /* switch places if smaller arc */
          {
            wreal = wr2[i0];
            wimag = wi2[i0];
            warc  = arc2[i0];
            wmag  = mag2[i0];
            wr2[i0] = wr2[nmark];
            wi2[i0] = wi2[nmark];
            arc2[i0] = arc2[nmark];
            mag2[i0] = mag2[nmark];
            wr2[nmark] = wreal;
            wi2[nmark] = wimag;
            arc2[nmark] = warc;
            mag2[nmark] = wmag;
          }
    }
  for (i0=0; i0<numroots; i0++)
  {
	  radius[i0]=mag2[i0];
	  bandwith[i0]=-log(mag2[i0])*mytype(p.sr)/M_PI;
	  phi[i0]=arc2[i0];
  }
}

void TransShift::trackPhi(mytype *r_ptr,mytype *phi_ptr,mytype time)
{// Dynamic programming based formant tracking (c.f., Xia and Espy-Wilson, 2000, ICSLP)
	mytype cum_Mat[MAX_FMT_TRACK_JUMP][MAX_NTRACKS];// cumulative cost matrix
	mytype cost_Mat[MAX_FMT_TRACK_JUMP][MAX_NTRACKS];// local cost Mat // just for debugging
	mytype fmts_min[MAX_NTRACKS]={0,350,1200,2000,3000}; // defines minimal value for each formant
	mytype fmts_max[MAX_NTRACKS]={1500,3500,4500,5000,7000};// defines maximal value for each formant
	mytype fn[MAX_NTRACKS]={500,1500,2500,3500,4500};// neutral formant values (start formants : here vowel [a])
	static mytype last_f[MAX_NTRACKS]={500,1500,2500,3500,4500};// last moving average estimated formants 
	int f_list[MAX_NTRACKS]={0,0,0,0,0};
	const int tri[MAX_NTRACKS]={0,1,2,3,4};
	mytype this_fmt,this_bw,this_cost, this_cum_cost, low_cost,min_cost,inf_cost=10000000;
	bool in_range=false;
	int k=0,i0=0,j0;
	int bound,indx=0,new_indx;

	mytype a_fact, g_fact, b_fact;

	k=0;
	i0=0;
	j0=0;

	fn[0]=p.fn1;
	fn[1]=p.fn2;

	a_fact=p.aFact;
	b_fact=p.bFact;
	g_fact=p.gFact;


	// loop builds the cumulative cost Matrix cum_Mat
	// each column represents the cumulative cost for each node which will be the cost entry value for the next column
	for (k=0;k<p.nTracks;k++)
	{
		low_cost=inf_cost;
			for(i0=k;i0<(p.nCands-p.nTracks+k+2);i0++)
			{
				//cum_Mat[i0-k][k]=inf_cost;
				this_cum_cost=inf_cost;
				this_fmt=phi_ptr[i0]*p.sr/(2*M_PI);
				this_bw=-log(r_ptr[i0])*p.sr/M_PI;
				if((this_fmt>fmts_min[k]) && (this_fmt<fmts_max[k]))// check if actual formant is in range
				{
					in_range=true;
					this_cost=a_fact*this_bw+b_fact*fabs(this_fmt-fn[k])+g_fact*fabs(this_fmt-last_f[k]);//calc local cost
					cost_Mat[i0-k][k]=this_cost;
					if (k==0)// build first column: cumulative cost = local cost
						this_cum_cost=this_cost;
					else// build all other columns: cumulative cost(this column) = cumulative cost(previous column)+local cost
						this_cum_cost=cum_Mat[i0-k][k-1]+this_cost;

					if (this_cum_cost<low_cost)
						low_cost=this_cum_cost;// low_cost is the lowest cumulative cost of all elements in this column until element [i0] (included)
											  // therefore, for each column :i0=0 low_cost=cumulative cost
											  //							:i0=n low_cost=min(cumulative_cost(from 0 to [n]))                           
				
				}
				if (k<p.nTracks-1)// for all columns except last
					cum_Mat[i0-k][k]=low_cost;//ATTENTION: represents the minimal cost that serves as entry cost for the next node (same row element, but next column)
				else// last column  
					cum_Mat[i0-k][k]=this_cum_cost;//shows the overall accumulated cost... from here will start the viterbi traceback
			}
	}
	
	bound=p.nCands-p.nTracks+2;// VERY IMPORTANT!!! because values of cum_Mat beyond this point are not referenced !!
	indx=0;
	// viterbi traceback updates index vector f_list
	// ATTENTION!!! f_list is not the definitive index vector.. has to be diagonalised
	for(k=0;k<p.nTracks;k++)
	{	
		min_cost=inf_cost;
		for(i0=0;i0<bound;i0++)
		{
			if(cum_Mat[i0][p.nTracks-k-1]<min_cost)
			{
				min_cost=cum_Mat[i0][p.nTracks-k-1];
				indx=i0;
			}
		}
		if(indx==0)
			break;
		else
		{
			bound=indx+1;
			f_list[p.nTracks-k-1]=indx;
		}
	}


	// update r, phi and last_f
	for(k=0;k<p.nTracks;k++)
	{
		new_indx=f_list[k]+k;// rediagonalize index vector
		r_ptr[k]=r_ptr[new_indx];
		phi_ptr[k]=phi_ptr[new_indx];
		last_f[k]=(1-p.trackFF)*phi_ptr[k]*p.sr/(2*M_PI)+p.trackFF*last_f[k];
	}

	
}

void TransShift::myFilt (mytype *xin_ptr, mytype* xout_ptr,mytype *oldPhi_ptr,mytype *newPhi_ptr,mytype *r_ptr,const int size)
	{// filter cascading two biquad IIR filters 

	// coefficients for the first filter (f1 shift) NOTE: b_filt1[0]=1 (see initilization)
	//SC The first filter is for correcting the first formant (?)
	b_filt1[1]=-2*r_ptr[0]*cos(oldPhi_ptr[0]);
	b_filt1[2]= r_ptr[0]*r_ptr[0]; 	
	a_filt1[0]=-2*r_ptr[0]*cos(newPhi_ptr[0]);  
	a_filt1[1]= r_ptr[0]*r_ptr[0];              
	
	// coefficients for the second filter (f2 shift) NOTE: b_filt2[0]=1 (see initilization)
	//SC The second filter is for correcting the second formant (?)
	b_filt2[1]=-2*r_ptr[1]*cos(oldPhi_ptr[1]);
	b_filt2[2]= r_ptr[1]*r_ptr[1]; 	
	a_filt2[0]=-2*r_ptr[1]*cos(newPhi_ptr[1]);  
	a_filt2[1]= r_ptr[1]*r_ptr[1];              

	DSPF_dp_biquad(xin_ptr, &b_filt1[0], &a_filt1[0], &filt_delay1[0], &filtbuf[0], size);//first formant
	DSPF_dp_biquad(&filtbuf[0],&b_filt2[0], &a_filt2[0], &filt_delay2[0], xout_ptr, size);//second formant
	
	}


// Calculate rms of buffer
mytype TransShift::calcRMS1(const mytype *xin_ptr, int size)
{
	//SC rmsFF: RMF forgetting factor, by default equals 0.9.
	ma_rms1=(1-p.rmsFF)*sqrt(DSPF_dp_vecsum_sq(xin_ptr,size)/(mytype)size)+p.rmsFF*ma_rms1;
	return ma_rms1 ;
}

mytype TransShift::calcRMS2(const mytype *xin_ptr, int size)
{
	ma_rms2=(1-p.rmsFF)*sqrt(DSPF_dp_vecsum_sq(xin_ptr,size)/(mytype)size)+p.rmsFF*ma_rms2;
	return ma_rms2 ;
}

int TransShift::sign(mytype x)
{
	if (x>0)
		return 1;
	else if (x<0)
		return -1;
	else
		return 0;
}


void TransShift::downSampSig(mytype *b,mytype *a, mytype *x, mytype *buffer, mytype *r,mytype  *d,const int nr, const int n_coeffs, const int downfact)
//SC Input parameters: 
//SC	b: IIR numerators;
//SC	a: IIR denominators; 
//SC	x: input frame;
//SC	buffer: filter buffer
//SC	r: final output
//SC	d: filter delay for downsampling
//SC	nr: frameLength after the downsampling
//SC	n_coeffs: number of coefficients
//SC	downfact: downsampling factor
{// filtering and decimation

	// filtering
	//SC implem: downSampSig(&srfilt_b[0], &srfilt_a[0],inFrame_ptr,&srfilt_buf[0], 
	//SC	&inFrameBuf[0],&srfilt_delay_down[0], p.frameLen , N_COEFFS_SRFILT , DOWNSAMP_FACT);
	iir_filt(b, a, x, buffer,d,nr*downfact, n_coeffs, 1);	//SC gain (g) is 1 here


	// decimation

	int i0;

	for(i0=0; i0 < nr; i0++)
	{
		r[i0] = buffer[downfact*i0]; 
	}
}

void TransShift::downSampSig_noFilt(mytype *b,mytype *a, mytype *x, mytype *buffer, mytype *r,mytype  *d,const int nr, const int n_coeffs, const int downfact)
//SC Input parameters: 
//SC	b: IIR numerators;
//SC	a: IIR denominators; 
//SC	x: input frame;
//SC	buffer: filter buffer
//SC	r: final output
//SC	d: filter delay for downsampling
//SC	nr: frameLength after the downsampling
//SC	n_coeffs: number of coefficients
//SC	downfact: downsampling factor
{// filtering and decimation

	// filtering
	//SC implem: downSampSig(&srfilt_b[0], &srfilt_a[0],inFrame_ptr,&srfilt_buf[0], 
	//SC	&inFrameBuf[0],&srfilt_delay_down[0], p.frameLen , N_COEFFS_SRFILT , DOWNSAMP_FACT);
	//iir_filt(b, a, x, buffer,d,nr*downfact, n_coeffs, 1);	//SC gain (g) is 1 here


	// decimation

	int i0;

	for(i0=0; i0 < nr; i0++)
	{
		r[i0] = x[downfact*i0]; 
	}
}


void TransShift::upSampSig (mytype  *b, mytype *a, mytype *x, mytype *buffer, mytype *r,mytype  *d,const int nr, const int n_coeffs, const int upfact,const mytype scalefact)
{//  interpolation and filtering

// interpolation
	int i0;

 for(i0=0; i0 < nr; i0++)
  {
	  buffer[i0]=0;
	if (i0 % upfact ==0)
	{
			buffer[i0]= x[int(i0/upfact)];
	}
  }

	// filtering
	iir_filt(b, a,buffer,r, d,nr, n_coeffs, upfact*scalefact);
}

void TransShift::iir_filt (mytype *b, mytype *a,  mytype *x, mytype *r,mytype  *d,const int nr, const int n_coeffs,  mytype g)
{// iir transposed II form
	// b : Numerator coeffs
	// a : Denominator coeffs
	// x : input frame
	// r : output frame
	// d : filter delays
	// nr: length of output frame
	// n_coeeffs : number of filter coeffs
	// !!! a and b must both!!! be of length (n_coeffs)
	// if you want a ~= b, fill with zeros
	// g : additional gain

	int m,k;

	for(m=0; m < nr; m++)
	{
		r[m]=g*b[0]*x[m] + d[0];
		for(k=0; k < n_coeffs-2; k++)
		{// start delay recursion
			d[k]=g*b[k+1]*x[m]+d[k+1]-a[k+1]*r[m];
		}
		d[n_coeffs-2] = g*b[n_coeffs-1]*x[m]-a[n_coeffs-1]*r[m]; 
	}

}

mytype TransShift::hz2mel(mytype hz){	// Convert frequency from Hz to mel
	mytype mel;
	mel=1127.01048*log(1+hz/700);
	return mel;
}

mytype TransShift::mel2hz(mytype mel){	// Convert frequency from mel to Hz
	mytype hz;
	hz=(exp(mel/1127.01048)-1)*700;
	return hz;
}

mytype TransShift::locateF2(mytype f2){	
//SC Locate the value of f2 in the pertF2 table, through a binary search.
//SC Usef for subsequent interpolation. 
	mytype loc;
	int k=1<<(PF_NBIT-1),n;

	for(n=0;n<PF_NBIT-1;n++){
		if (f2>=p.pertF2[k])
			k=k+(1<<(PF_NBIT-n-2));
		else
			k=k-(1<<(PF_NBIT-n-2));
	}	
	if (f2<p.pertF2[k])
		k--;

	loc=(mytype)k;

	loc+=(f2-p.pertF2[k])/(p.pertF2[k+1]-p.pertF2[k]);

	if (loc>=PF_NPOINTS-1){	//PF_NPOINTS=257, so locint not be greater than 255 (PF_NPOINTS-2)
		loc=PF_NPOINTS-1-0.000000000001;
	}
	if (loc<0){
		loc=0;
	}
	return loc;
}

void TransShift::DSPF_dp_cfftr2(int n, mytype * x, mytype * w, int n_min)	//SC Fast Fourier transform on complex numbers
{
	int n2, ie, ia, i, j, k, m;
	mytype rtemp, itemp, c, s;
	n2 = n;

	ie = 1;

	for(k = n; k > n_min; k >>= 1)
	{
		n2 >>= 1;
		ia = 0;
		for(j=0; j < ie; j++)
		{
			for(i=0; i < n2; i++)
			{
				c = w[2*i];
				s = w[2*i+1];
				m = ia + n2;
				rtemp = x[2*ia] - x[2*m];
				x[2*ia] = x[2*ia] + x[2*m];
				itemp = x[2*ia+1] - x[2*m+1];
				x[2*ia+1] = x[2*ia+1] + x[2*m+1];
				x[2*m] = c*rtemp - s*itemp;
				x[2*m+1] = c*itemp + s*rtemp;
				ia++;
			}
			ia += n2;
		}	
		ie <<= 1;
		w = w + k;
	}
}

void TransShift::DSPF_dp_icfftr2(int n, double * x, double * w, int n_min)	//SC Inverse Fast Fourier transform on complex numbers
{
	int n2, ie, ia, i, j, k, m;
	double rtemp, itemp, c, s;
	n2 = n;
	ie = 1;
	for(k = n; k > n_min; k >>= 1)
	{
		n2 >>= 1;
		ia = 0;
		for(j=0; j < ie; j++)
		{
			for(i=0; i < n2; i++)
			{
				c = w[2*i];
				s = w[2*i+1];
				m = ia + n2;
				rtemp = x[2*ia] - x[2*m];
				x[2*ia] = x[2*ia] + x[2*m];
				itemp = x[2*ia+1] - x[2*m+1];
				x[2*ia+1] = x[2*ia+1] + x[2*m+1];
				x[2*m] = c*rtemp + s*itemp;
				x[2*m+1] = c*itemp - s*rtemp;
				ia++;
			}
			ia += n2;
		}
		ie <<= 1;
		w = w + k;
	}
}

void TransShift::gen_w_r2(double* w, int n)		//SC An FFT subroutine
{
	int i, j=1;
	double pi = 4.0*atan(1.0);
	double e = pi*2.0/n;
	for(j=1; j < n; j <<= 1)
	{
	for(i=0; i < ( n>>1 ); i += j)
	{
		*w++ = cos(i*e);
		*w++ = -sin(i*e);
	}
}
}

void TransShift::bit_rev(double* x, int n)	//SC Bit reversal: an FFT subroutine
{
	int i, j, k;
	double rtemp, itemp;

	j = 0;
	for(i=1; i < (n-1); i++)
	{
		k = n >> 1;
		while(k <= j)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
		if(i < j)
		{
			rtemp = x[j*2];
			x[j*2] = x[i*2];
			x[i*2] = rtemp;
			itemp = x[j*2+1];
			x[j*2+1] = x[i*2+1];
			x[i*2+1] = itemp;
		}
	}
}

void TransShift::smbFft(mytype *fftBuffer, double fftFrame_Size, int sign)
/* 
	
 * FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
	Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
	time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
	and returns the cosine and sine parts in an interleaved manner, ie.
	fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
	must be a power of 2. It expects a complex input signal (see footnote 2),
	ie. when working with 'common' audio signals our input signal has to be
	passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
	of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
{
	mytype wr, wi, arg,  *p1, *p2, temp;
	mytype tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
	long i, bitm, j, k, le, le2, logN;
	
	logN = (long)(log(fftFrame_Size)/log(2.)+0.5);


	for (i = 2; i < 2*fftFrame_Size-2; i += 2) 
	{
		for (bitm = 2, j = 0; bitm < 2*fftFrame_Size; bitm <<= 1) 
		{
			if (i & bitm) j++;
			j <<= 1;
		}

		if (i < j) 
		{
			p1 = fftBuffer+i; p2 = fftBuffer+j;
			temp = *p1;
			*(p1++) = *p2;
			*(p2++) = temp; 
			temp = *p1;
			*p1 = *p2; 
			*p2 = temp;
		}
	}

	
	for (k = 0, le = 2; k < logN; k++) 
	{
		le <<= 1;
		le2 = le>>1;
		ur = 1.0;
		ui = 0.0;
		arg = M_PI /(le2>>1);
		wr = cos(arg);
		wi = sign*sin(arg);
		for (j = 0; j < le2; j += 2) 
		{
			p1r = fftBuffer+j; p1i = p1r+1;
			p2r = p1r+le2; p2i = p2r+1;
			for (i = j; i < 2*fftFrame_Size; i += le) 
			{
				tr = *p2r * ur - *p2i * ui;
				ti = *p2r * ui + *p2i * ur;
				*p2r = *p1r - tr; *p2i = *p1i - ti;
				*p1r += tr; *p1i += ti;
				p1r += le; p1i += le;
				p2r += le; p2i += le;
			}
			tr = ur*wr - ui*wi;
			ui = ur*wi + ui*wr;
			ur = tr;
		}


	}

}

