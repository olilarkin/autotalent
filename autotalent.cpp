/* 
 
 autotalent.cpp
 
 autotalent vst/au
 
 An auto-tuning vst/au plugin, based on 
 
 autotalent An auto-tuning LADSPA plugin.
 
 Free software by Thomas A. Baran.
 http://web.mit.edu/tbaran/www/autotalent.html
 VERSION 0.2
 March 20, 2010
 
 port by oli larkin
 
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

#include "autotalent.h"
#include "IPlug_include_in_plug_src.h"
#include "IControl.h"
#include "resource.h"
#include <math.h>

#define AT_TUNE 0
#define AT_FIXED 1
#define AT_PULL 2
#define AT_A 3
#define AT_Bb 4
#define AT_B 5
#define AT_C 6
#define AT_Db 7
#define AT_D 8
#define AT_Eb 9
#define AT_E 10
#define AT_F 11
#define AT_Gb 12
#define AT_G 13
#define AT_Ab 14
#define AT_AMOUNT 15
#define AT_SMOOTH 16
#define AT_SHIFT 17
#define AT_SCWARP 18
#define AT_LFOAMP 19
#define AT_LFORATE 20
#define AT_LFOSHAPE 21
#define AT_LFOSYMM 22
#define AT_LFOQUANT 23
#define AT_FCORR 24
#define AT_FWARP 25
#define AT_MIX 26
//#define AT_PITCH 27
//#define AT_CONF 28
//#define AT_INPUT1  29
//#define AT_OUTPUT1 30
//#define AT_LATENCY 31

#define SQUARE 0
#define SINE 1
#define TRI 2

// The number of presets/programs
const int kNumPrograms = 1;




autotalent::autotalent(IPlugInstanceInfo instanceInfo)
:	IPLUG_CTOR(27, kNumPrograms, instanceInfo), mfs(44100)
{
	TRACE;

	// Define parameter ranges, display units, labels.
	//arguments are: name, defaultVal, minVal, maxVal, step, label
	//GetParam(kGain)->InitDouble("Gain", 0.0, -70.0, 12.0, 0.1, "dB");
	
	init(mfs);
	
	GetParam(AT_TUNE)->InitDouble("Concert A ref", 440.0, 430.0, 450.0, 0.1, "Hz");
	GetParam(AT_FIXED)->InitDouble("Fixed pitch", 0.0, -36, 36.0, 0.01, "semitones");
	GetParam(AT_PULL)->InitDouble("Pull to fixed pitch", 0.0, 0., 1.0, 0.01, "");
	
	GetParam(AT_A)->InitEnum("A", 1, 3);
	GetParam(AT_Bb)->InitEnum("Bb", 1, 3);
	GetParam(AT_B)->InitEnum("B", 1, 3);
	GetParam(AT_C)->InitEnum("C", 1, 3);
	GetParam(AT_Db)->InitEnum("Db", 1, 3);
	GetParam(AT_D)->InitEnum("D", 1, 3);
	GetParam(AT_Eb)->InitEnum("Eb", 1, 3);
	GetParam(AT_E)->InitEnum("E", 1, 3);
	GetParam(AT_F)->InitEnum("F", 1, 3);
	GetParam(AT_Gb)->InitEnum("Gb", 1, 3);
	GetParam(AT_G)->InitEnum("G", 1, 3);
	GetParam(AT_Ab)->InitEnum("Ab", 1, 3);
	
	for(int i=3; i<15;i++)
	{
		GetParam(i)->SetDisplayText(0, "out");
		GetParam(i)->SetDisplayText(1, "in");
		GetParam(i)->SetDisplayText(2, "snap");
	}
	
	GetParam(AT_AMOUNT)->InitDouble("Correction strength", 0.0, 0., 1.0, 0.01, "");
	GetParam(AT_SMOOTH)->InitDouble("Correction smoothness", 0.0, 0., 1.0, 0.01, "");
	GetParam(AT_SHIFT)->InitDouble("Pitch shift", 0.0, -12.0, 12.0, 0.01, "semitones");
	GetParam(AT_SCWARP)->InitInt("Output scale rotate", 0, -3, 3,  "steps");
	GetParam(AT_LFOAMP)->InitDouble("LFO depth", 0.0, 0., 1.0, 0.01, "");
	GetParam(AT_LFORATE)->InitDouble("LFO rate", 0.05, 0.01, 10.0, 0.0001, "Hz");
	GetParam(AT_LFOSHAPE)->InitEnum("LFO shape", SINE, 3);
	GetParam(AT_LFOSHAPE)->SetDisplayText(SQUARE, "square");
	GetParam(AT_LFOSHAPE)->SetDisplayText(SINE, "sine");
	GetParam(AT_LFOSHAPE)->SetDisplayText(TRI, "triangle");
	GetParam(AT_LFOSYMM)->InitDouble("LFO symmetry", 0.0, -1.0, 1.0, 0.01, "");
	GetParam(AT_LFOQUANT)->InitBool("LFO quantise", 0, "");
	GetParam(AT_FCORR)->InitBool("Formant correction", 1, "");
	GetParam(AT_FWARP)->InitDouble("Formant warp", 0.0, -1.0, 1.0, 0.01, "");
	GetParam(AT_MIX)->InitDouble("Mix", 1., 0., 1., 0.01, "");
	
	
	//MakePreset("preset 1", -5.0, 5.0, 17, kReversed);
	MakeDefaultPreset("-", kNumPrograms);
}

autotalent::~autotalent() 
{
	fft_des(mfmembvars);
	free(mcbi);
	free(mcbf);
	free(mcbo);
	free(mcbwindow);
	free(mhannwindow);
	free(macwinv);
	free(mfrag);
	free(mffttime);
	free(mfftfreqre);
	free(mfftfreqim);
	free(mfk);
	free(mfb);
	free(mfc);
	free(mfrb);
	free(mfrc);
	free(mfsmooth);
	free(mfsig);
	
	for (int ti=0; ti<mford; ti++) 
	{
		free(mfbuff[ti]);
	}
	
	free(mfbuff);
	free(mftvec);
}

void autotalent::ProcessDoubleReplacing(double** inputs, double** outputs, int SampleCount)
{
	double* pfInput;
	double* pfOutputL;
	double* pfOutputR;
	float fAmount;
	float fSmooth;
	int iNotes[12];
	int iPitch2Note[12];
	int iNote2Pitch[12];
	int numNotes;
	float fTune;
	float fFixed;
	float fPull;
	float fShift;
	int iScwarp;
	float fLfoamp;
	float fLforate;
	float fLfoshape;
	float fLfosymm;
	int iLfoquant;
	int iFcorr;
	float fFwarp;
	float fMix;

	unsigned long lSampleIndex;
	
	long int N;
	long int Nf;
	long int fs;
	float pmin;
	float pmax;
	unsigned long nmin;
	unsigned long nmax;
	
	long int ti;
	long int ti2;
	long int ti3;
	long int ti4;
	float tf;
	float tf2;
	
	// Variables for cubic spline interpolator
	float indd;
	int ind0;
	int ind1;
	int ind2;
	int ind3;
	float vald;
	float val0;
	float val1;
	float val2;
	float val3;
	
	int lowersnap;
	int uppersnap;
	float lfoval;
	
	float pperiod;
	float inpitch;
	float conf;
	float outpitch;
	float aref;
	float fa;
	float fb;
	float fc;
	float fk;
	float flamb;
	float frlamb;
	float falph;
	float foma;
	float f1resp;
	float f0resp;
	float flpa;
	int ford;
	
	pfInput = inputs[0];
	pfOutputL = outputs[0];
	pfOutputR = outputs[1];
	fAmount = (float) m_pfAmount;
	fSmooth = (float) m_pfSmooth * 0.8; // Scales max to a more reasonable value
	fTune = (float) m_pfTune;
	iNotes[0] = (int) m_pfA;
	iNotes[1] = (int) m_pfBb;
	iNotes[2] = (int) m_pfB;
	iNotes[3] = (int) m_pfC;
	iNotes[4] = (int) m_pfDb;
	iNotes[5] = (int) m_pfD;
	iNotes[6] = (int) m_pfEb;
	iNotes[7] = (int) m_pfE;
	iNotes[8] = (int) m_pfF;
	iNotes[9] = (int) m_pfGb;
	iNotes[10] = (int) m_pfG;
	iNotes[11] = (int) m_pfAb;
	fFixed = (float) m_pfFixed;
	fPull = (float) m_pfPull;
	fShift = (float) m_pfShift;
	iScwarp = (int) m_pfScwarp;
	fLfoamp = (float) m_pfLfoamp;
	fLforate = (float) m_pfLforate;
	fLfoshape = (float) m_pfLfoshape;
	fLfosymm = (float) m_pfLfosymm;
	iLfoquant = (int) m_pfLfoquant;
	iFcorr = (int) m_pfFcorr;
	fFwarp = (float) m_pfFwarp;
	fMix = (float) m_pfMix;
	
	// Some logic for the semitone->scale and scale->semitone conversion
	// If no notes are selected as being in the scale, instead snap to all notes
	ti2 = 0;
	for (ti=0; ti<12; ti++) {
		if (iNotes[ti]>=0) {
			iPitch2Note[ti] = ti2;
			iNote2Pitch[ti2] = ti;
			ti2 = ti2 + 1;
		}
		else {
			iPitch2Note[ti] = -1;
		}
	}
	numNotes = ti2;
	while (ti2<12) {
		iNote2Pitch[ti2] = -1;
		ti2 = ti2 + 1;
	}
	if (numNotes==0) {
		for (ti=0; ti<12; ti++) {
			iNotes[ti] = 1;
			iPitch2Note[ti] = ti;
			iNote2Pitch[ti] = ti;
		}
		numNotes = 12;
	}
	iScwarp = (iScwarp + numNotes*5)%numNotes;
	
	ford = mford;
	falph = mfalph;
	foma = (float)1 - falph;
	flpa = mflpa;
	flamb = mflamb;
	tf = pow((float)2,fFwarp/2)*(1+flamb)/(1-flamb);
	frlamb = (tf - 1)/(tf + 1);
	
	maref = (float)fTune;
	
	N = mcbsize;
	Nf = mcorrsize;
	fs = mfs;
	
	pmax = mpmax;
	pmin = mpmin;
	nmax = mnmax;
	nmin = mnmin;
	
	aref = maref;
	pperiod = mpmax;
	inpitch = minpitch;
	conf = mconf;
	outpitch = moutpitch;
	
	
	/*******************
	 *  MAIN DSP LOOP  *
	 *******************/
	for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++)  {
		
		// load data into circular buffer
		tf = (float) *(pfInput++);
		ti4 = mcbiwr;
		mcbi[ti4] = tf;
		
		if (iFcorr>=1) {
			// Somewhat experimental formant corrector
			//  formants are removed using an adaptive pre-filter and
			//  re-introduced after pitch manipulation using post-filter
			// tf is signal input
			fa = tf - mfhp; // highpass pre-emphasis filter
			mfhp = tf;
			fb = fa;
			for (ti=0; ti<ford; ti++) {
				mfsig[ti] = fa*fa*foma + mfsig[ti]*falph;
				fc = (fb-mfc[ti])*flamb + mfb[ti];
				mfc[ti] = fc;
				mfb[ti] = fb;
				fk = fa*fc*foma + mfk[ti]*falph;
				mfk[ti] = fk;
				tf = fk/(mfsig[ti] + 0.000001);
				tf = tf*foma + mfsmooth[ti]*falph;
				mfsmooth[ti] = tf;
				mfbuff[ti][ti4] = tf;
				fb = fc - tf*fa;
				fa = fa - tf*fc;
			}
			mcbf[ti4] = fa;
			// Now hopefully the formants are reduced
			// More formant correction code at the end of the DSP loop
		}
		else {
			mcbf[ti4] = tf;
		}
		
		
		// Input write pointer logic
		mcbiwr++;
		if (mcbiwr >= N) {
			mcbiwr = 0;
		}
		
		
		// ********************
		// * Low-rate section *
		// ********************
		
		// Every N/noverlap samples, run pitch estimation / manipulation code
		if ((mcbiwr)%(N/mnoverlap) == 0) {
			
			// ---- Obtain autocovariance ----
			
			// Window and fill FFT buffer
			ti2 = mcbiwr;
			for (ti=0; ti<N; ti++) {
				mffttime[ti] = (float)(mcbi[(ti2-ti+N)%N]*mcbwindow[ti]);
			}
			
			// Calculate FFT
			fft_forward(mfmembvars, mffttime, mfftfreqre, mfftfreqim);
			
			// Remove DC
			mfftfreqre[0] = 0;
			mfftfreqim[0] = 0;
			
			// Take magnitude squared
			for (ti=1; ti<Nf; ti++) {
				mfftfreqre[ti] = (mfftfreqre[ti])*(mfftfreqre[ti]) + (mfftfreqim[ti])*(mfftfreqim[ti]);
				mfftfreqim[ti] = 0;
			}
			
			// Calculate IFFT
			fft_inverse(mfmembvars, mfftfreqre, mfftfreqim, mffttime);
			
			// Normalize
			tf = (float)1/mffttime[0];
			for (ti=1; ti<N; ti++) {
				mffttime[ti] = mffttime[ti] * tf;
			}
			mffttime[0] = 1;
			
			//  ---- END Obtain autocovariance ----
			
			
			//  ---- Calculate pitch and confidence ----
			
			// Calculate pitch period
			//   Pitch period is determined by the location of the max (biased)
			//     peak within a given range
			//   Confidence is determined by the corresponding unbiased height
			tf2 = 0;
			pperiod = pmin;
			for (ti=nmin; ti<nmax; ti++) {
				ti2 = ti-1;
				ti3 = ti+1;
				if (ti2<0) {
					ti2 = 0;
				}
				if (ti3>Nf) {
					ti3 = Nf;
				}
				tf = mffttime[ti];
				
				if (tf>mffttime[ti2] && tf>=mffttime[ti3] && tf>tf2) {
					tf2 = tf;
					ti4 = ti;
				}
			}
			if (tf2>0) {
				conf = tf2*macwinv[ti4];
				if (ti4>0 && ti4<Nf) {
					// Find the center of mass in the vicinity of the detected peak
					tf = mffttime[ti4-1]*(ti4-1);
					tf = tf + mffttime[ti4]*(ti4);
					tf = tf + mffttime[ti4+1]*(ti4+1);
					tf = tf/(mffttime[ti4-1] + mffttime[ti4] + mffttime[ti4+1]);
					pperiod = tf/fs;
				}
				else {
					pperiod = (float)ti4/fs;
				}
			}
			
			// Convert to semitones
			tf = (float) -12*log10((float)aref*pperiod)*L2SC;
			if (conf>=mvthresh) {
				inpitch = tf;
				minpitch = tf; // update pitch only if voiced
			}
			mconf = conf;
			
			// think this is LADSPA display OL
			//      *(m_pfPitch) = (LADSPA_Data) inpitch;
			//     *(m_pfConf) = (LADSPA_Data) conf;
			
			//  ---- END Calculate pitch and confidence ----
			
			
			//  ---- Modify pitch in all kinds of ways! ----
			
			outpitch = inpitch;
			
			// Pull to fixed pitch
			outpitch = (1-fPull)*outpitch + fPull*fFixed;
			
			// -- Convert from semitones to scale notes --
			ti = (int)(outpitch/12 + 32) - 32; // octave
			tf = outpitch - ti*12; // semitone in octave
			ti2 = (int)tf;
			ti3 = ti2 + 1;
			// a little bit of pitch correction logic, since it's a convenient place for it
			if (iNotes[ti2%12]<0 || iNotes[ti3%12]<0) { // if between 2 notes that are more than a semitone apart
				lowersnap = 1;
				uppersnap = 1;
			}
			else {
				lowersnap = 0;
				uppersnap = 0;
				if (iNotes[ti2%12]==1) { // if specified by user
					lowersnap = 1;
				}
				if (iNotes[ti3%12]==1) { // if specified by user
					uppersnap = 1;
				}
			}
			// (back to the semitone->scale conversion)
			// finding next lower pitch in scale
			while (iNotes[(ti2+12)%12]<0) {
				ti2 = ti2 - 1;
			}
			// finding next higher pitch in scale
			while (iNotes[ti3%12]<0) {
				ti3 = ti3 + 1;
			}
			tf = (tf-ti2)/(ti3-ti2) + iPitch2Note[(ti2+12)%12];
			if (ti2<0) {
				tf = tf - numNotes;
			}
			outpitch = tf + numNotes*ti;
			// -- Done converting to scale notes --
			
			// The actual pitch correction
			ti = (int)(outpitch+128) - 128;
			tf = outpitch - ti - 0.5;
			ti2 = ti3-ti2;
			if (ti2>2) { // if more than 2 semitones apart, put a 2-semitone-like transition halfway between
				tf2 = (float)ti2/2;
			}
			else {
				tf2 = (float)1;
			}
			if (fSmooth<0.001) {
				tf2 = tf*tf2/0.001;
			}
			else {
				tf2 = tf*tf2/fSmooth;
			}
			if (tf2<-0.5) tf2 = -0.5;
			if (tf2>0.5) tf2 = 0.5;
			tf2 = 0.5*sin(PI*tf2) + 0.5; // jumping between notes using horizontally-scaled sine segment
			tf2 = tf2 + ti;
			if ( (tf<0.5 && lowersnap) || (tf>=0.5 && uppersnap) ) {
				outpitch = fAmount*tf2 + ((float)1-fAmount)*outpitch;
			}
			
			// Add in pitch shift
			outpitch = outpitch + fShift;
			
			// LFO logic
			tf = fLforate*N/(mnoverlap*fs);
			if (tf>1) tf=1;
			mlfophase = mlfophase + tf;
			if (mlfophase>1) mlfophase = mlfophase-1;
			lfoval = mlfophase;
			tf = (fLfosymm + 1)/2;
			if (tf<=0 || tf>=1) {
				if (tf<=0) lfoval = 1-lfoval;
			}
			else {
				if (lfoval<=tf) {
					lfoval = lfoval/tf;
				}
				else {
					lfoval = 1 - (lfoval-tf)/(1-tf);
				}
			}
			if (fLfoshape>=0) {
				// linear combination of cos and line
				lfoval = (0.5 - 0.5*cos(lfoval*PI))*fLfoshape + lfoval*(1-fLfoshape);
				lfoval = fLfoamp*(lfoval*2 - 1);
			}
			else {
				// smoosh the sine horizontally until it's squarish
				tf = 1 + fLfoshape;
				if (tf<0.001) {
					lfoval = (lfoval - 0.5)*2/0.001;
				}
				else {
					lfoval = (lfoval - 0.5)*2/tf;
				}
				if (lfoval>1) lfoval = 1;
				if (lfoval<-1) lfoval = -1;
				lfoval = fLfoamp*sin(lfoval*PI*0.5);
			}
			// add in quantized LFO
			if (iLfoquant>=1) {
				outpitch = outpitch + (int)(numNotes*lfoval + numNotes + 0.5) - numNotes;
			}
			
			
			// Convert back from scale notes to semitones
			outpitch = outpitch + iScwarp; // output scale rotate implemented here
			ti = (int)(outpitch/numNotes + 32) - 32;
			tf = outpitch - ti*numNotes;
			ti2 = (int)tf;
			ti3 = ti2 + 1;
			outpitch = iNote2Pitch[ti3%numNotes] - iNote2Pitch[ti2];
			if (ti3>=numNotes) {
				outpitch = outpitch + 12;
			}
			outpitch = outpitch*(tf - ti2) + iNote2Pitch[ti2];
			outpitch = outpitch + 12*ti;
			outpitch = outpitch - (iNote2Pitch[iScwarp] - iNote2Pitch[0]); //more scale rotation here
			
			// add in unquantized LFO
			if (iLfoquant<=0) {
				outpitch = outpitch + lfoval*2;
			}
			
			
			if (outpitch<-36) outpitch = -48;
			if (outpitch>24) outpitch = 24;
			
			moutpitch = outpitch;
			
			//  ---- END Modify pitch in all kinds of ways! ----
			
			// Compute variables for pitch shifter that depend on pitch
			minphinc = aref*pow(2,inpitch/12)/fs;
			moutphinc = aref*pow(2,outpitch/12)/fs;
			mphincfact = moutphinc/minphinc;
		}
		// ************************
		// * END Low-Rate Section *
		// ************************
		
		
		
		// *****************
		// * Pitch Shifter *
		// *****************
		
		// Pitch shifter (kind of like a pitch-synchronous version of Fairbanks' technique)
		//   Note: pitch estimate is naturally N/2 samples old
		mphasein = mphasein + minphinc;
		mphaseout = mphaseout + moutphinc;
		
		//   When input phase resets, take a snippet from N/2 samples in the past
		if (mphasein >= 1) {
			mphasein = mphasein - 1;
			ti2 = mcbiwr - N/2;
			for (ti=-N/2; ti<N/2; ti++) {
				mfrag[(ti+N)%N] = mcbf[(ti + ti2 + N)%N];
			}
		}
		
		//   When output phase resets, put a snippet N/2 samples in the future
		if (mphaseout >= 1) {
			mfragsize = mfragsize*2;
			if (mfragsize > N) {
				mfragsize = N;
			}
			mphaseout = mphaseout - 1;
			ti2 = mcbord + N/2;
			ti3 = (long int)(((float)mfragsize) / mphincfact);
			if (ti3>=N/2) {
				ti3 = N/2 - 1;
			}
			for (ti=-ti3/2; ti<(ti3/2); ti++) {
				tf = mhannwindow[(long int)N/2 + ti*(long int)N/ti3];
				// 3rd degree polynomial interpolator - based on eqns from Hal Chamberlin's book
				indd = mphincfact*ti;
				ind1 = (int)indd;
				ind2 = ind1+1;
				ind3 = ind1+2;
				ind0 = ind1-1;
				val0 = mfrag[(ind0+N)%N];
				val1 = mfrag[(ind1+N)%N];
				val2 = mfrag[(ind2+N)%N];
				val3 = mfrag[(ind3+N)%N];
				vald = 0;
				vald = vald - (float)0.166666666667 * val0 * (indd - ind1) * (indd - ind2) * (indd - ind3);
				vald = vald + (float)0.5 * val1 * (indd - ind0) * (indd - ind2) * (indd - ind3);
				vald = vald - (float)0.5 * val2 * (indd - ind0) * (indd - ind1) * (indd - ind3);
				vald = vald + (float)0.166666666667 * val3 * (indd - ind0) * (indd - ind1) * (indd - ind2);
				mcbo[(ti + ti2 + N)%N] = mcbo[(ti + ti2 + N)%N] + vald*tf;
			}
			mfragsize = 0;
		}
		mfragsize++;
		
		//   Get output signal from buffer
		tf = mcbo[mcbord]; // read buffer
		
		mcbo[mcbord] = 0; // erase for next cycle
		mcbord++; // increment read pointer
		if (mcbord >= N) {
			mcbord = 0;
		}
		
		// *********************
		// * END Pitch Shifter *
		// *********************
		
		ti4 = (mcbiwr + 2)%N;
		if (iFcorr>=1) {
			// The second part of the formant corrector
			// This is a post-filter that re-applies the formants, designed
			//   to result in the exact original signal when no pitch
			//   manipulation is performed.
			// tf is signal input
			// gotta run it 3 times because of a pesky delay free loop
			//  first time: compute 0-response
			tf2 = tf;
			fa = 0;
			fb = fa;
			for (ti=0; ti<ford; ti++) {
				fc = (fb-mfrc[ti])*frlamb + mfrb[ti];
				tf = mfbuff[ti][ti4];
				fb = fc - tf*fa;
				mftvec[ti] = tf*fc;
				fa = fa - mftvec[ti];
			}
			tf = -fa;
			for (ti=ford-1; ti>=0; ti--) {
				tf = tf + mftvec[ti];
			}
			f0resp = tf;
			//  second time: compute 1-response
			fa = 1;
			fb = fa;
			for (ti=0; ti<ford; ti++) {
				fc = (fb-mfrc[ti])*frlamb + mfrb[ti];
				tf = mfbuff[ti][ti4];
				fb = fc - tf*fa;
				mftvec[ti] = tf*fc;
				fa = fa - mftvec[ti];
			}
			tf = -fa;
			for (ti=ford-1; ti>=0; ti--) {
				tf = tf + mftvec[ti];
			}
			f1resp = tf;
			//  now solve equations for output, based on 0-response and 1-response
			tf = (float)2*tf2;
			tf2 = tf;
			tf = ((float)1 - f1resp + f0resp);
			if (tf!=0) {
				tf2 = (tf2 + f0resp) / tf;
			}
			else {
				tf2 = 0;
			}
			//  third time: update delay registers
			fa = tf2;
			fb = fa;
			for (ti=0; ti<ford; ti++) {
				fc = (fb-mfrc[ti])*frlamb + mfrb[ti];
				mfrc[ti] = fc;
				mfrb[ti] = fb;
				tf = mfbuff[ti][ti4];
				fb = fc - tf*fa;
				fa = fa - tf*fc;
			}
			tf = tf2;
			tf = tf + flpa*mflp;  // lowpass post-emphasis filter
			mflp = tf;
			// Bring up the gain slowly when formant correction goes from disabled
			// to enabled, while things stabilize.
			if (mfmute>0.5) {
				tf = tf*(mfmute - 0.5)*2;
			}
			else {
				tf = 0;
			}
			tf2 = mfmutealph;
			mfmute = (1-tf2) + tf2*mfmute;
			// now tf is signal output
			// ...and we're done messing with formants
		}
		else {
			mfmute = 0;
		}
		
		// Write audio to output of plugin
		// Mix (blend between original (delayed) =0 and processed =1)
		*(pfOutputL++) = *(pfOutputR++) = fMix*tf + (1-fMix)*mcbi[ti4];
		
	}
	
	// Tell the host the algorithm latency
	SetLatency(N-1);
	
}	

void autotalent::Reset()
{
	TRACE;
	IMutexLock lock(this);
	
	//mSampleRate = GetSampleRate();
	//mSamplePeriod = 1./mSampleRate;
	
	unsigned long sr = GetSampleRate();
	
	if( mfs != sr ) init(sr);
}


void autotalent::OnParamChange(int paramIdx)
{
	IMutexLock lock(this);
	
	switch (paramIdx)
	{
			
		case AT_TUNE:
			m_pfTune = GetParam(AT_TUNE)->Value();
			break;
		case AT_FIXED:
			m_pfFixed = GetParam(AT_FIXED)->Value();
			break;
		case AT_PULL:
			m_pfPull = GetParam(AT_PULL)->Value();
			break;
		case AT_A:
			m_pfA = GetParam(AT_A)->Value() - 1;
			break;
		case AT_Bb:
			m_pfBb = GetParam(AT_Bb)->Value() - 1;
			break;
		case AT_B:
			m_pfB = GetParam(AT_B)->Value() - 1;
			break;
		case AT_C:
			m_pfC = GetParam(AT_C)->Value() - 1;
			break;
		case AT_Db:
			m_pfDb = GetParam(AT_Db)->Value() - 1;
			break;
		case AT_D:
			m_pfD = GetParam(AT_D)->Value() - 1;
			break;
		case AT_Eb:
			m_pfEb = GetParam(AT_Eb)->Value() - 1;
			break;
		case AT_E:
			m_pfE = GetParam(AT_E)->Value() - 1;
			break;
		case AT_F:
			m_pfF = GetParam(AT_F)->Value() - 1;
			break;
		case AT_Gb:
			m_pfGb = GetParam(AT_Gb)->Value() - 1;
			break;
		case AT_G:
			m_pfG = GetParam(AT_G)->Value() - 1;
			break;
		case AT_Ab:
			m_pfAb = GetParam(AT_Ab)->Value() - 1;
			break;
		case AT_AMOUNT:
			m_pfAmount = GetParam(AT_AMOUNT)->Value();
			break;
		case AT_SMOOTH:
			m_pfSmooth = GetParam(AT_SMOOTH)->Value();
			break;
		case AT_SHIFT:
			m_pfShift = GetParam(AT_SHIFT)->Value();
			break;
		case AT_SCWARP:
			m_pfScwarp = GetParam(AT_SCWARP)->Value();
			break;
		case AT_LFOAMP:
			m_pfLfoamp = GetParam(AT_LFOAMP)->Value();
			break;
		case AT_LFORATE:
			m_pfLforate = GetParam(AT_LFORATE)->Value();
			break;
		case AT_LFOSHAPE:
			m_pfLfoshape = GetParam(AT_LFOSHAPE)->Value() - 1;
			break;
		case AT_LFOSYMM:
			m_pfLfosymm = GetParam(AT_LFOSYMM)->Value();
			break;
		case AT_LFOQUANT:
			m_pfLfoquant = GetParam(AT_LFOQUANT)->Value();
			break;
		case AT_FCORR:
			m_pfFcorr = GetParam(AT_FCORR)->Value();
			break;
		case AT_FWARP:
			m_pfFwarp = GetParam(AT_FWARP)->Value();
			break;
		case AT_MIX:
			m_pfMix = GetParam(AT_MIX)->Value();
			break;
			
		default:
			break;
	}
}

void autotalent::init(unsigned long SampleRate)
{
	unsigned long ti;
	
	//Autotalent* membvars = malloc(sizeof(Autotalent));
	
	maref = 440;
	
	mfs = SampleRate;
	
	if (SampleRate>=88200) {
		mcbsize = 4096;
	}
	else {
		mcbsize = 2048;
	}
	mcorrsize = mcbsize / 2 + 1;
	
	mpmax = 1/(float)70;  // max and min periods (ms)
	mpmin = 1/(float)700; // eventually may want to bring these out as sliders
	
	mnmax = (unsigned long)(SampleRate * mpmax);
	if (mnmax > mcorrsize) {
		mnmax = mcorrsize;
	}
	mnmin = (unsigned long)(SampleRate * mpmin);
	
	mcbi = (float*) calloc(mcbsize, sizeof(float));
	mcbf = (float*) calloc(mcbsize, sizeof(float));
	mcbo = (float*) calloc(mcbsize, sizeof(float));
	
	mcbiwr = 0;
	mcbord = 0;
	
	mlfophase = 0;
	
	// Initialize formant corrector
	mford = 7; // should be sufficient to capture formants
	mfalph = pow(0.001f, (float) 80 / (SampleRate));
	mflamb = -(0.8517*sqrt(atan(0.06583*SampleRate))-0.1916); // or about -0.88 @ 44.1kHz
	mfk = (float*) calloc(mford, sizeof(float));
	mfb = (float*) calloc(mford, sizeof(float));
	mfc = (float*) calloc(mford, sizeof(float));
	mfrb = (float*) calloc(mford, sizeof(float));
	mfrc = (float*) calloc(mford, sizeof(float));
	mfsig = (float*) calloc(mford, sizeof(float));
	mfsmooth = (float*) calloc(mford, sizeof(float));
	mfhp = 0;
	mflp = 0;
	mflpa = pow(0.001f, (float) 10 / (SampleRate));
	mfbuff = (float**) malloc((mford)*sizeof(float*));
	for (ti=0; ti<mford; ti++) 
	{
		mfbuff[ti] = (float*) calloc(mcbsize, sizeof(float));
	}
	mftvec = (float*) calloc(mford, sizeof(float));
	mfmute = 1;
	mfmutealph = pow(0.001f, (float)1 / (SampleRate));
	
	// Standard raised cosine window, max height at N/2
	mhannwindow = (float*) calloc(mcbsize, sizeof(float));
	for (ti=0; ti<mcbsize; ti++) {
		mhannwindow[ti] = -0.5*cos(2*PI*ti/mcbsize) + 0.5;
	}
	
	// Generate a window with a single raised cosine from N/4 to 3N/4
	mcbwindow = (float*) calloc(mcbsize, sizeof(float));
	for (ti=0; ti<(mcbsize / 2); ti++) {
		mcbwindow[ti+mcbsize/4] = -0.5*cos(4*PI*ti/(mcbsize - 1)) + 0.5;
	}
	
	mnoverlap = 4;
	
	mfmembvars = fft_con(mcbsize);
	
	mffttime = (float*) calloc(mcbsize, sizeof(float));
	mfftfreqre = (float*) calloc(mcorrsize, sizeof(float));
	mfftfreqim = (float*) calloc(mcorrsize, sizeof(float));
	
	
	// ---- Calculate autocorrelation of window ----
	macwinv = (float*) calloc(mcbsize, sizeof(float));
	for (ti=0; ti<mcbsize; ti++) {
		mffttime[ti] = mcbwindow[ti];
	}
	fft_forward(mfmembvars, mcbwindow, mfftfreqre, mfftfreqim);
	for (ti=0; ti<mcorrsize; ti++) {
		mfftfreqre[ti] = (mfftfreqre[ti])*(mfftfreqre[ti]) + (mfftfreqim[ti])*(mfftfreqim[ti]);
		mfftfreqim[ti] = 0;
	}
	fft_inverse(mfmembvars, mfftfreqre, mfftfreqim, mffttime);
	for (ti=1; ti<mcbsize; ti++) {
		macwinv[ti] = mffttime[ti]/mffttime[0];
		if (macwinv[ti] > 0.000001) {
			macwinv[ti] = (float)1/macwinv[ti];
		}
		else {
			macwinv[ti] = 0;
		}
	}
	macwinv[0] = 1;
	// ---- END Calculate autocorrelation of window ----
	
	
	mlrshift = 0;
	mptarget = 0;
	msptarget = 0;
	
	mvthresh = 0.7;  //  The voiced confidence (unbiased peak) threshold level
	
	// Pitch shifter initialization
	mphprdd = 0.01; // Default period
	minphinc = (float)1/(mphprdd * SampleRate);
	mphincfact = 1;
	mphasein = 0;
	mphaseout = 0;
	mfrag = (float*) calloc(mcbsize, sizeof(float));
	mfragsize = 0;
	
	SetLatency(mcbsize-1);
}


