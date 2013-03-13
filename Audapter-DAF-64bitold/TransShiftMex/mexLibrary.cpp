#include "mex.h"
#include "TransShift.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "audioIO.h"

//#include "procWarpAtom.h"

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	int i0, j0, action;
	unsigned int k0;
	mytype  *signal_ptr,*data_ptr;
	const mytype  *algosignal_ptr, *algodata_ptr;
	mytype   retval = 0;
	mytype   value= 0;
	int size=1;
	int vecsize=1;
	int toPrompt=0;
	static TransShift algo;		//SC algo is a TransShift class object. Notice it is a static pointer.
	static audioIO audio_obj;	//SC Refer to the audioIO project
	static	DeviceParams devpar;//SC DeviceParams is a struct, fields: num, chans, fs, set
	static bool started = 0;	
	int frame_len; 
	int activeDeviceNum;

	action = (int)floor(*(double*)mxGetPr(prhs[0]));	//SC Figure out the input argument

	switch (action){
		
        case 0:		//SC enumerate all the audio devices
			printf("Information: \n");
			for(k0 = 0; k0 < audio_obj.devices.size(); k0++)
			{
				if (audio_obj.devices[k0].inputChannels > 0)	//SC input devices
					printf("[I]: %d : %s\n", k0, audio_obj.devices[k0].name.c_str());
				if (audio_obj.devices[k0].outputChannels > 0)	//SC input devices
					printf("[O]: %d : %s\n", k0, audio_obj.devices[k0].name.c_str());
			}

			frame_len = algo.getparams((void *)"framelen");
			audio_obj.setcallbackparams(frame_len * DOWNSAMP_FACT, (void *)&algoCallbackFunc, (void *)&algo);			
			break;

		case 1:			//SC set audio device parameters and start the action			
			frame_len = algo.getparams((void *)"framelen");
			/* printf("frame_len = %d\n", frame_len); */

			audio_obj.setcallbackparams(frame_len * DOWNSAMP_FACT, (void *)&algoCallbackFunc, (void *)&algo);
			devpar.chans = 2;
			devpar.fs = algo.getparams((void *)"srate")*DOWNSAMP_FACT;	//SC device sampling rate
			
			//for (k0 = 0; k0 < audio_obj.devices.size(); k0++)
			//{
			//	if (audio_obj.devices[k0].inputChannels > 0 && audio_obj.devices[k0].outputChannels > 0)	//SC input devices
			//		break;
			//}
			for(i0 = 0; i0 < audio_obj.devices.size(); i0++)			
				if ((audio_obj.devices[i0].inputChannels > 0) && (audio_obj.devices[i0].outputChannels > 0)) { // Find the first active output device
					activeDeviceNum = i0 + 1;
					break;
				}

			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,1);
			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,2);

			if (!started)
			{
				printf("Action Start:  %d\n", action);
				algo.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;

		case 2:			//SC Stop the audio device
			if (started){
				printf("Action End:  %d\n", action);
				if (audio_obj.started)
					audio_obj.stopdev();
			}
			else
				printf("Not started\n");
			started = 0;
			break;
		
		case 3:			//SC Set paramters of the algo object
			if (nrhs >= 3){				
				retval=algo.setparams((void *)mxArrayToString(prhs[1]), (void *)mxGetPr(prhs[2]));
			}
			toPrompt=1;
			if (nrhs >=4){
				toPrompt=(int)floor(*(double*)mxGetPr(prhs[3]));
			}
			if (toPrompt==1){
				printf("Setting param: %s = %4.0f \n", mxArrayToString(prhs[1]), *(mytype *)mxGetPr(prhs[2]));
			}
			if (retval)
				printf("!!! UNKNOWN PARAMETER !!! Could not set param: %s  \n" , mxArrayToString(prhs[1]));						

			break;
			
		case 4:			//SC Manually get three outputs: 1
			// signal
			algosignal_ptr = algo.getsignal(size);	//SC size is updated in this process: size = frame_counter*p.frameLen
			if (size>0)
			{

				plhs[0] = mxCreateDoubleMatrix(size,2, mxREAL);
				signal_ptr = mxGetPr(plhs[0]);
				for(i0=0;i0<size;i0++)
				{
					signal_ptr[i0] = algosignal_ptr[i0];
					signal_ptr[size+i0] = algosignal_ptr[MAX_REC_SIZE+i0];
				}

				// data 
				algodata_ptr = algo.getdata(size,vecsize);
				plhs[1] = mxCreateDoubleMatrix(size,vecsize, mxREAL);
				data_ptr = mxGetPr(plhs[1]);

				for(i0=0;i0<size;i0++)
				{
					for(j0=0;j0<vecsize;j0++)
						data_ptr[j0*size+i0] = algodata_ptr[j0*MAX_DATA_SIZE+i0];
				}				
			}
			else
			{
				plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
				plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
				plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
			}

			break;

		case 5:				//SC manually supply the buffer and run the process
			data_ptr = (double*)mxGetPr(prhs[1]);	//SC pointer to buffer
			algoCallbackFunc((char *)data_ptr, algo.getparams((void *)"framelen")*DOWNSAMP_FACT, (void *)&algo);
			break;

		case 6:				//SC manually reset algo object
			algo.reset();
			break;

		case 11:			//SC Sine wave generator
			audio_obj.setcallbackparams(algo.getparams((void *)"framelen")*DOWNSAMP_FACT, (void *)&algoCallbackFuncSineGen, (void *)&algo);
			devpar.chans = 2;
			
			devpar.fs = algo.getparams((void *)"srate")*DOWNSAMP_FACT;	//SC device sampling rate
			for(i0 = 0; i0 < audio_obj.devices.size(); i0++)			
				if ((audio_obj.devices[i0].inputChannels > 0) && (audio_obj.devices[i0].outputChannels > 0)) { // Find the first active output device
					activeDeviceNum = i0 + 1;
					break;
				}

			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,1);
			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,2);

			if (!started)
			{
				printf("Action Start:  %d\n", action);
				algo.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;

		case 12:			//SC wave playback
			audio_obj.setcallbackparams(algo.getparams((void *)"framelen")*DOWNSAMP_FACT, (void *)&algoCallbackFuncWavePB, (void *)&algo);
			devpar.chans = 2;

			devpar.fs = algo.getparams((void *)"srate")*DOWNSAMP_FACT;	//SC device sampling rate

			for(i0 = 0; i0 < audio_obj.devices.size(); i0++)			
				if ((audio_obj.devices[i0].inputChannels > 0) && (audio_obj.devices[i0].outputChannels > 0)) { // Find the first active output device
					activeDeviceNum = i0 + 1;
					break;
				}

			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,1);
			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,2);

			if (!started)
			{
				printf("Action Start:  %d\n", action);
				algo.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;


		case 21:	//SC(2012/02/29) Write data to binary file
			algosignal_ptr = algo.getsignal(size);	//SC size is updated in this process: size = frame_counter*p.frameLen
			if (size > 0){
				printf("Writing data to file data0.bin: size = %d (%d bytes)\n", size, size * (sizeof mytype));
				ofstream dataFileCl ("data0.bin");
				dataFileCl.close();

				fstream dataFile("data0.bin", ios::binary | ios::out);
				if (!dataFile){
					printf("WARNING: Cannot open file data0.bin\n");
				}

				dataFile.write((char *) algosignal_ptr, size * (sizeof mytype));
				dataFile.close();
			}
			break;

		case 22:			//SC Manually get the entire recording buffer
			// signal
			//algosignal_ptr = algo.getsignal(size);	//SC size is updated in this process: size = frame_counter*p.frameLen
			plhs[0] = mxCreateDoubleMatrix(MAX_REC_SIZE, 1, mxREAL);
			signal_ptr = mxGetPr(plhs[0]);
			algosignal_ptr = algo.getsignal(size);
			for(i0=0; i0 < MAX_REC_SIZE;i0++)
			{
				signal_ptr[i0] = algosignal_ptr[i0];
			}

			
			break;

		case 23:			// SC Manually get the entire playback buffer
			plhs[0] = mxCreateDoubleMatrix(MAX_REC_SIZE, 1, mxREAL);
			signal_ptr = mxGetPr(plhs[0]);
			algosignal_ptr = algo.getsignal(size);
			for(i0=0; i0 < MAX_REC_SIZE;i0++)
			{
				signal_ptr[i0] = algosignal_ptr[MAX_REC_SIZE+i0];
			}

			
			break;

		default:
			printf("Action Unknown:  %d\n", action);
			
	}

}