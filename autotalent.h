/* 
 
 autotalent.h
 
 autotalent
 
 An auto-tuning plugin, based on 
 
 Autotalent An auto-tuning LADSPA plugin.
 
 Free software by Thomas A. Baran.
 http://web.mit.edu/tbaran/www/autotalent.html
 VERSION 0.2
 March 20, 2010
 
 port to IPlug by oli larkin
 
 http://www.olilarkin.co.uk
 
 version 0.2
 
 This program is free software; you can redistribute it and/or modify        
 it under the terms of the GNU General Public License as published by        
 the Free Software Foundation; either version 2 of the License, or           
 (at your option) any later version.                                         
 
 This program is distributed in the hope that it will be useful,             
 but WITHOUT ANY WARRANTY; without even the implied warranty of              
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               
 GNU General Public License for more details.                                
 
 You should have received a copy of the GNU General Public License           
 along with this program; if not, write to the Free Software                 
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  
 

*/

#ifndef __AUTOTALENT__
#define __AUTOTALENT__

#include "IPlug_include_in_plug_hdr.h"
#define L2SC (float)3.32192809488736218171

#include "mayer_fft.c"
#include "fftsetup.h"

class autotalent : public IPlug
{
public:
	
	autotalent(IPlugInstanceInfo instanceInfo);
	~autotalent();
	
	// Implement these if your audio or GUI logic requires doing something, when params change or when audio processing stops/starts.
	void Reset();
	void OnParamChange(int paramIdx);

	void ProcessDoubleReplacing(double** inputs, double** outputs, int SampleCount);
	void init(unsigned long SampleRate);
	fft_vars* mfmembvars; // member variables for fft routine

private:
	
	
	unsigned long mfs; // Sample rate
	
	unsigned long mcbsize; // size of circular buffer
	unsigned long mcorrsize; // cbsize/2 + 1
	unsigned long mcbiwr;
	unsigned long mcbord;
	float* mcbi; // circular input buffer
	float* mcbf; // circular formant correction buffer
	float* mcbo; // circular output buffer
	
	float* mcbwindow; // hann of length N/2, zeros for the rest
	float* macwinv; // inverse of autocorrelation of window
	float* mhannwindow; // length-N hann
	int mnoverlap;
	
	float* mffttime;
	float* mfftfreqre;
	float* mfftfreqim;
	
	// VARIABLES FOR LOW-RATE SECTION
	float maref; // A tuning reference (Hz)
	float minpitch; // Input pitch (semitones)
	float mconf; // Confidence of pitch period estimate (between 0 and 1)
	float moutpitch; // Output pitch (semitones)
	float mvthresh; // Voiced speech threshold
	
	float mpmax; // Maximum allowable pitch period (seconds)
	float mpmin; // Minimum allowable pitch period (seconds)
	unsigned long mnmax; // Maximum period index for pitch prd est
	unsigned long mnmin; // Minimum period index for pitch prd est
	
	float mlrshift; // Shift prescribed by low-rate section
	int mptarget; // Pitch target, between 0 and 11
	float msptarget; // Smoothed pitch target
	
	float mlfophase;
	
	// VARIABLES FOR PITCH SHIFTER
	float mphprdd; // default (unvoiced) phase period
	double minphinc; // input phase increment
	double moutphinc; // input phase increment
	double mphincfact; // factor determining output phase increment
	double mphasein;
	double mphaseout;
	float* mfrag; // windowed fragment of speech
	unsigned long mfragsize; // size of fragment in samples
	
	// VARIABLES FOR FORMANT CORRECTOR
	int mford;
	float mfalph;
	float mflamb;
	float* mfk;
	float* mfb;
	float* mfc;
	float* mfrb;
	float* mfrc;
	float* mfsig;
	float* mfsmooth;
	float mfhp;
	float mflp;
	float mflpa;
	float** mfbuff;
	float* mftvec;
	float mfmute;
	float mfmutealph;
	
	
	float m_pfTune;
	float m_pfFixed;
	float m_pfPull;
	int m_pfA;
	int m_pfBb;
	int m_pfB;
	int m_pfC;
	int m_pfDb;
	int m_pfD;
	int m_pfEb;
	int m_pfE;
	int m_pfF;
	int m_pfGb;
	int m_pfG;
	int m_pfAb;
	float m_pfAmount;
	float m_pfSmooth;
	float m_pfShift;
	int m_pfScwarp;
	float m_pfLfoamp;
	float m_pfLforate;
	float m_pfLfoshape;
	float m_pfLfosymm;
	int m_pfLfoquant;
	int m_pfFcorr;
	float m_pfFwarp;
	float m_pfMix;
	float m_pfPitch;
	float m_pfConf;
	float m_pfInputBuffer1;
	float m_pfOutputBuffer1;
	float m_pfLatency;
};

////////////////////////////////////////

#endif
