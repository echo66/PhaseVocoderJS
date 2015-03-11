function PhaseVocoder(winSize, sampleRate) {

	var _sampleRate = sampleRate;

	var _RS = 0;
	var _RA = 0;
	var _omega;

	var _previousInputPhase;
	var _previousOutputPhase;
	var _framingWindow;
	var _squaredFramingWindow;
	
	var _winSize = winSize;

	var _overlapBuffers;

	var _owOverlapBuffers;

	var _first = true;


	this.init = function() {

		var _ = this;

		_omega = _.create_omega_array(winSize);

		_previousInputPhase = _.create_constant_array(winSize/2, 0);
		_previousOutputPhase = _.create_constant_array(winSize/2, 0);
		_framingWindow = _.create_sin_beta_window_array(winSize, 1);
		_squaredFramingWindow = _framingWindow.map(function(x,i){ return x*x; });

		_overlapBuffers = _.create_constant_array(winSize, 0);

		_owOverlapBuffers = _.create_constant_array(winSize, 0);

		_.set_alpha(1);
	}


	this.get_previous_input_phase = function() {
		return _previousInputPhase;
	}

	this.get_previous_output_phase = function() {
		return _previousOutputPhase;
	}


	/*
	 *  TODO
	 *  
	 *  Method ported from the Java Phase Vocoder (https://github.com/groakley/phase-vocoder-java).
	 *  
	 *  @param RS: integer describing the Resynthesis Hop Size.
	 *  @param frame: Float32Array array with 'windowSize' time samples.
	 *  @param overlapBuffer: Float32Array array with samples from previous iterations.
	 *  @param windowSize: integer describing the size of the window.
	 *
	 *  @returns Float32Array array with 'RS' size.
	 */
	this.overlap_and_slide = function(RS, frame, overlapBuffer, windowSize) {
		//var finishedBytes = new Float32Array(RS);
		var finishedBytes = new Array(RS);

		for (var i=0; i<RS; i++) {
			finishedBytes[i] = overlapBuffer.shift();
			while(overlapBuffer.length > windowSize - 1 && overlapBuffer.length >= 0)
				overlapBuffer.shift();
			overlapBuffer.push.apply(overlapBuffer, [0.0]);
		}

		var outBytes = [].concat(overlapBuffer);

		for (var i=0; i<outBytes.length; i++) 
			outBytes[i] = frame[i] + overlapBuffer[i];

		while ((overlapBuffer.length > windowSize - outBytes.length) && overlapBuffer.length >= 0) {
			overlapBuffer.shift();
	    }

	    overlapBuffer.push.apply(overlapBuffer, outBytes);
		
		return finishedBytes;
	}


	/*
	 *  TODO
	 *
	 *  [CHECKED IN MATLAB]
	 * 
	 *  @param size: number (integer) of frequency bins.
	 *
	 *  @returns Phase advances per sample for the frequency bins.
	 */
	this.create_omega_array = function(size) {
		return Array.apply(null, Array(size/2 + 1)).map(function (x, i) { 
			return 2 * Math.PI * i / size;
		});
	}


	/*
	 *  TODO
	 *
	 *  [CHECKED IN MATLAB]
	 *  
	 *  @param size: TODO
	 *  @param beta: TODO
	 *
	 *  @returns TODO
	 */
	this.create_sin_beta_window_array = function(size, beta) {
		return Array.apply(null, Array(size)).map(function(x,i){
			return Math.pow(Math.sin(Math.PI*i/size), beta);
		});
	}


	this.create_constant_array = function(size, constant) {
		return Array.apply(null, Array(size)).map(function () { 
			return constant; 
		});
	}


	/*
	 *  An index in spec is considered a peak if its value is the largest among its 4 nearest neighbours.
	 *
	 *  [CHECKED IN MATLAB]
	 *
	 *  @param magFrame: object array with the magnitude of a STFT frame.
	 *
	 *  @returns object with fields 'peaks', the bins' indexes of the peaks, 'inflRegionStart', 
	 * the bins' indexes where the influence regions start, and 'inflRegionEnd', the bins' indexes 
	 * where the influence regions end.
	 */
	this.find_peaks_v3 = function(magFrame) {
		var magSpecPad = [0,0].concat(magFrame).concat([0,0]);

		var peaks = magSpecPad.slice(2,magSpecPad.length-2).map(function(x,i){
			var I = i + 2;
			if(x > magSpecPad[I-2] && x > magSpecPad[I-1] && x > magSpecPad[I+1] && x > magSpecPad[I+2]) {
				return i;
			}
		}).filter(function(x){ 
			return x!=undefined && x!=null; 
		});

		// var inflRegStart = []; var inflRegEnd = [];
		var inflRegStart = new Array(peaks.length); 
		var inflRegEnd = new Array(peaks.length);

		inflRegStart[0] = 0;
		for (var i=0; i<peaks.length-1; i++) {
			inflRegStart[i+1] = Math.ceil((peaks[i] + peaks[i+1])/2); 
		}

		for (var i=1; i<inflRegStart.length; i++) {
			inflRegEnd[i-1] = inflRegStart[i]-1;
		}
		inflRegEnd[inflRegEnd.length] = inflRegEnd.length-1;

		return {peaks: peaks, inflRegionStart: inflRegStart, inflRegionEnd: inflRegEnd};
	}

	this.find_peaks_v4 = function(magFrame) {
		var magSpecPad = [0,0].concat(magFrame).concat([0,0]);
		var peaks = [];

		for (var i=2, I=0; i<=magSpecPad.length-2; i++, I++) {
			x = magSpecPad[i];
			if (x > magSpecPad[i-2] && x > magSpecPad[i-1] && x > magSpecPad[i+1] && x > magSpecPad[i+2]) {
				peaks = peaks.concat(I);
			}
		}

		var inflRegStart = new Array(peaks.length); 

		inflRegStart[0] = 0;
		for (var i=0; i<peaks.length-1; i++) {
			inflRegStart[i+1] = Math.ceil((peaks[i] + peaks[i+1])/2); 
		}

		var inflRegEnd = new Array(peaks.length);
		for (var i=1; i<inflRegStart.length; i++) {
			inflRegEnd[i-1] = inflRegStart[i]-1;
		}
		inflRegEnd[inflRegEnd.length] = inflRegEnd.length-1;

		// var influenceRegions = new Array(inflRegionStart.length);
		// for (var i=0; i<influenceRegions.length; i++) 
		// 	influenceRegions[i] = Math.max(0, inflRegionEnd[i] - inflRegionStart[i] + 1);

		return {
			peaks: peaks, 
			inflRegionStart: inflRegStart, 
			inflRegionEnd: inflRegEnd,
			// influenceRegions: influenceRegions
		};
	}

	/*
	 *  Returns the instantaneous phase advances per synthesis hopsize.
	 * 
	 *  [CHECKED IN MATLAB]
	 *
	 *  @param currentInputPhase, Float32Array vector with the phases of the current input frame.
	 *  @param omega, phase advances per sample for the frequencies k.
	 *  @param previousInputPhase, phases of the last input frame.
	 *  @param RA, analysis hopsize. Currently, the code doesn't allow 
	 * the usage of anchor points.
	 *  @param RS, synthesis hopsize.
	 *
	 *  @returns Float32Array array with the instantaneous phase advances.
	 */
	this.get_phase_advances_v4 = function(currentInputPhase, previousInputPhase, omega, RA, RS) {
		var twoPI = 2 * Math.PI;
		var instPhaseAdvHop = new Array(omega.length);

		for (var i=0; i<omega.length; i++) {
			var expectedPhaseAdv = omega[i] * RA;

			var auxheterodynedPhaseIncr = (currentInputPhase[i] - previousInputPhase[i]) - expectedPhaseAdv;
			var heterodynedPhaseIncr = auxheterodynedPhaseIncr - twoPI * Math.round(auxheterodynedPhaseIncr/twoPI);

			var instPhaseAdvPerSampleHop = omega[i] + heterodynedPhaseIncr / RA;

			instPhaseAdvHop[i] = instPhaseAdvPerSampleHop * RS;
		}

		return instPhaseAdvHop;
	}

	/*
	 *  TODO
	 *
	 *  [CHECKED IN MATLAB]
	 *
	 *  @param currentInputPhase: Float32Array with the phases of the current input frame.
	 *  @param previousOutputPhase: Float32Array with the phases of the last output frame.
	 *  @param instPhaseAdv: Float32Array with the instantaneous phases advance.
	 *  @param frequencyBins: integer array, each value representing a frequency 
	 * bin to be used in the phasor estimation.
	 *  @param influenceRegions: sorted integer array, each value representing 
	 * the ID of the frequency region.
	 *
	 *	@returns TODO
	 */
	this.get_phasor_theta_v2 = function(currentInputPhase, previousOutputPhase, instPhaseAdv, frequencyBins, influenceRegions) {
		// Get the peaks in the spectrum together with their regions of influence.

		//var theta = new Float32Array(currentInputPhase.length);
		var theta = new Array(currentInputPhase.length);

		var theta_idx = 0;
		for (var i=0; i<frequencyBins.length; i++) {
			var bin = frequencyBins[i];
			for (var j=0; j<influenceRegions[i]; j++, theta_idx++) {
				theta[theta_idx] = previousOutputPhase[bin] + instPhaseAdv[bin] - currentInputPhase[bin];
			}
		}

		var remaining_length = theta.length - theta_idx;
		for (var i=0; i<remaining_length; i++, theta_idx++)
			theta[theta_idx] = 0;

		return theta;
	}

	/*
	 *  Compute a phasor that rotates the phase angles of the current
	 * input frame by angles theta such that no phase discontinuities occur 
	 * when resynthesizing the resulting spectrogram with the synthesis 
	 * hopsize.
	 *
	 *  @param currentInputMagnitude: Float32Array holding the magnitude values 
	 * of the current input frame.
	 *  @param currentInputMagnitude: Float32Array holding the phase values of 
	 * the current input frame.
	 *  @param previousOutputPhase: Float32Array array with the phases of the last output frame.
	 *  @param instPhaseAdv: Float32Array array with the instantaneous phases advance.
	 *
	 *  @returns Float32Array array with the phasor angles.
	 */
	this.identity_phase_locking_v2 = function(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {
		var _ = this;

		var r = _.find_peaks_v4(currentInputMagnitude);

		var influenceRegions = new Array(r.inflRegionStart.length);

		for (var i=0; i<influenceRegions.length; i++) 
			influenceRegions[i] = Math.max(0, r.inflRegionEnd[i] - r.inflRegionStart[i] + 1);

		var phasor_theta = _.get_phasor_theta_v2(currentInputPhase, previousOutputPhase, instPhaseAdv, r.peaks, influenceRegions);

		return phasor_theta;
	}

	this.no_phase_locking = function(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {
		var theta = new Array(currentInputPhase.length);

		for (var i=0; i<currentInputPhase.length; i++) {
			theta[i] = previousOutputPhase[i] + instPhaseAdv[i] - currentInputPhase[i];
		}

		return theta;
	}

	/**
	 *  TODO
	 *
	 *  @param fftObject, an object with the following fields: 
	 * 		'real', a Float32Array vector with the real part of the FFT frame;
	 *		'imag', a Float32Array vector with the imaginary part of the FFT frame;
	 *		'magnitude', a Float32Array vector with the magnitude of the FFT frame;
	 *		'phase', a Float32Array vector with the phase/angle of the FFT frame.
	 *  @param previousInputPhase: TODO
	 *  @param previousOutputPhase: TODO
	 *  @param omega: TODO
	 *  @param RA: Analysis Hop Size
	 *  @param RS: Resynthesis Hop Size
	 *
	 *  @returns an object with the following fields:
	 *		'real': a Float32Array vector with the real part of the output frame.
	 *		'imag': a Float32Array vector with the imaginary part of the output frame.
	 *		'phase': a Float32Array vector with the phase/angle of the output frame.
	 */
	this.pv_step_v2 = function(fftObject, previousInputPhase, previousOutputPhase, omega, RA, RS) {

		var _ = this;

		var currentInputPhase = fftObject.phase;

		var instPhaseAdv = _.get_phase_advances_v4(currentInputPhase, previousInputPhase, omega, RA, RS);

		var currentInputMag = fftObject.magnitude;

		// var phasor_theta = _.no_phase_locking(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);
		var phasor_theta = _.identity_phase_locking_v2(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);

		var out_real = new Array((phasor_theta.length-1)*2);
		var out_imag = new Array((phasor_theta.length-1)*2);
		var out_phase = new Array((phasor_theta.length-1)*2);
		var out_magnitude = new Array((phasor_theta.length-1)*2);
		var doubleSize = (phasor_theta.length-1)*2;
		var sqrt	= Math.sqrt;
		var cos		= Math.cos;
		var sin		= Math.sin;
		var atan2	= Math.atan2;

		
		for (var i=0; i<phasor_theta.length; i++) {
			var theta = phasor_theta[i];

			var phasor_theta_real = cos(theta);
			var phasor_theta_imag = sin(theta);
			out_real[i] = phasor_theta_real * fftObject.real[i] - phasor_theta_imag * fftObject.imag[i];
			out_imag[i] = phasor_theta_real * fftObject.imag[i] + phasor_theta_imag * fftObject.real[i];
			out_phase[i] = atan2(out_imag[i], out_real[i]);
			out_magnitude[i] = sqrt(out_imag[i]*out_imag[i] + out_real[i]*out_real[i]);

			if (i>0) {
				out_real[doubleSize-i] = out_real[i];
				out_imag[doubleSize-i] = -out_imag[i];
				out_phase[doubleSize-i] = atan2(out_imag[doubleSize-i], out_real[doubleSize-i]);
				out_magnitude[doubleSize-i] = sqrt(out_imag[doubleSize-i]*out_imag[doubleSize-i] + out_real[doubleSize-i]*out_real[doubleSize-i]);
			}

		}
		
		return {real: out_real, imag: out_imag, phase: out_phase, magnitude: out_magnitude};
	}


	this.process = function(inputFrame) {

		var _ = this;

		var __RS = _RS;
		var __RA = _RA;

		var fftObject = (_first)? _.STFT(inputFrame, _framingWindow, _winSize) : _.STFT(inputFrame, _framingWindow, Math.round(_winSize/2)+1);

		var out = (_first)? null : _.pv_step_v2(fftObject, _previousInputPhase, _previousOutputPhase, _omega, __RA, __RS);
		_previousOutputPhase = (_first)? fftObject.phase : out.phase;
		_previousInputPhase = fftObject.phase;
		var processedFrame = (_first)? _.ISTFT(fftObject.real, fftObject.imag, _framingWindow) : _.ISTFT(out.real, out.imag, _framingWindow);
		_first = false;
		var outputFrame = _.overlap_and_slide(__RS, processedFrame, _overlapBuffers, _winSize);

		var owFrame = _.overlap_and_slide(__RS, _squaredFramingWindow, _owOverlapBuffers, _winSize);

		outputFrame = outputFrame.map(function(sample, i){
			return sample / ((owFrame[i]<10e-3)? 1 : owFrame[i]);
		});

		return outputFrame;

	}

	this.process_debug = function(inputFrame, params) {

		var _ = this;

		var oldFirst = _first;
		_first = false;

		var fftObject = {};

		if (params.stft.do_stft) {
			fftObject = (oldFirst)? _.STFT(inputFrame, _framingWindow, _winSize) : _.STFT(inputFrame, _framingWindow, Math.round(_winSize/2)+1);
		} else {
			fftObject.real = params.stft.real;
			fftObject.imag = params.stft.imag;
			fftObject.magnitude = params.stft.magnitude;
			fftObject.phase = params.stft.phase;
		}

		var out = (oldFirst)? null : _.pv_step_v2(fftObject, _previousInputPhase, _previousOutputPhase, _omega, _RA, _RS);
		_previousOutputPhase = (oldFirst)? fftObject.phase : out.phase;
		_previousInputPhase = fftObject.phase;
		

		if (!params.istft.do_istft) {
			return {
				correctedSpectrum: (oldFirst)? fftObject : out
			};
		}

		var processedFrame = (oldFirst)? _.ISTFT(fftObject.real, fftObject.imag, _framingWindow, true) : _.ISTFT(out.real, out.imag, _framingWindow, true);

		if (!params.overlap_and_slide.do_overlap_and_slide)
			return {
				istft: processedFrame,
				correctedSpectrum: (oldFirst)? fftObject : out
			};

		var outputFrame = _.overlap_and_slide(_RS, processedFrame, _overlapBuffers, _winSize);

		var owFrame = _.overlap_and_slide(_RS, _squaredFramingWindow, _owOverlapBuffers, _winSize);

		outputFrame = outputFrame.map(function(sample, i){
			return sample / ((owFrame[i]<10e-3)? 1 : owFrame[i]);
		});

		return {
			outputFrame: outputFrame,
			istft: processedFrame,
			correctedSpectrum: fftObject
		};
	}

	this.reset = function() {

		var _ = this;

		_previousInputPhase = _.create_constant_array(winSize/2, 0);
		_previousOutputPhase = _.create_constant_array(winSize/2, 0);

		_overlapBuffers = _.create_constant_array(winSize, 0);
		_owOverlapBuffers = _.create_constant_array(winSize, 0);

		_first = true;
	}

	this.reset2 = function() {

		var _ = this;

		_previousInputPhase = _.create_constant_array(winSize/2, 0);
		_previousOutputPhase = _.create_constant_array(winSize/2, 0);

		_first = true;
	}

	

	this.STFT = function(inputFrame, windowFrame, wantedSize) {
		var winSize = windowFrame.length;
		var _inputFrame = new Array(winSize);
		var fftFrame = new Array(2*winSize);

		for (var i=0; i<winSize; i++) {
			_inputFrame[i] = inputFrame[i] * windowFrame[i];
		}
		var fft = new FFT.complex(winSize, false);
		fft.simple(fftFrame, _inputFrame, 'real');

		var real = new Array(Math.min(winSize,wantedSize));
		var imag = new Array(Math.min(winSize,wantedSize));
		var magnitude = new Array(Math.min(winSize,wantedSize));
		var phase = new Array(Math.min(winSize,wantedSize));

		for (var p=0; p<winSize && p<wantedSize; p++) {
			real[p] = fftFrame[2*p];
			imag[p] = fftFrame[2*p+1];
			magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
			phase[p] = Math.atan2(imag[p], real[p]);
		}

		return {
			real: real,
			imag: imag,
			magnitude: magnitude,
			phase: phase
		}
	}

	this.ISTFT = function(real, imaginary, windowFrame, restoreEnergy) {
		var input = new Array(2 * real.length);
		var output1 = new Array(2 * real.length);
		var output2 = new Array(real.length);

		for (var i=0; i<real.length; i++) {
			input[2*i] = real[i];
			input[2*i+1] = imaginary[i];
		}

		// var ifft = new FFT.complex(real.length, true);

		// ifft.simple(output1, input);

		var ifft = new FFT.complex(real.length, true);

		ifft.simple(output1, input, 'complex');

		if (restoreEnergy) {
			var energy1 = 0;
			var energy2 = 0;
			var eps = 2.2204e-16;
			for (var i=0; i<windowFrame.length; i++) {
				energy1 += Math.abs(output1[2*i]);
				output2[i] = output1[2*i] / windowFrame.length;
				output2[i] *= windowFrame[i];
				energy2 += Math.abs(output1[2*i]);
				output2[i] *= energy1/(energy2+eps);
			}
		} else if (windowFrame) {
			for (var i=0; i<windowFrame.length; i++) {
				output2[i] = output1[2*i] / windowFrame.length;
				output2[i] *= windowFrame[i];
			}
		} else {
			for (var i=0; i<real.length; i++) {
				output2[i] = output1[2*i] / real.length;
			}
		}

		return output2;
	}

	this.get_analysis_hop = function() {
		return _RA;
	}

	this.get_synthesis_hop = function() {
		return _RS;
	}

	this.get_alpha = function() {
		return _RS / _RA;
	}

	this.get_framing_window = function() {
		return _framingWindow;
	}

	this.get_squared_framing_window = function() {
		return _squaredFramingWindow;
	}

	this.set_alpha = function(newAlpha) {
		_RA = _winSize/4;
		_RS = Math.round(newAlpha * _RA);
		// _RS = Math.round(_winSize/2);
		// _RA = Math.round(_RS / newAlpha);
	}

	this.get_alpha_step = function() {
		return 1/_RA;
	}

	this.set_hops = function(RA, RS) {
		_RA = RA;
		_RS = RS;
	}
}