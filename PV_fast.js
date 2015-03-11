function PhaseVocoder(winSize, sampleRate) {

	var _sampleRate = sampleRate; var _RS = 0; var _RA = 0; var _omega;

	var _previousInputPhase; var _previousOutputPhase; var _framingWindow;
	
	var _squaredFramingWindow; var _winSize = winSize;

	var _overlapBuffers; var _owOverlapBuffers;

	var _first = true;


	
	function overlap_and_slide(RS, frame, overlapBuffer, windowSize, finishedBytes) {

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
		
		return;
	}


	function find_peaks(magFrame, out) {
		var magSpecPad = [0,0].concat(magFrame).concat([0,0]);
		out.peaks = [];

		for (var i=2, I=0; i<=magSpecPad.length-2; i++, I++) {
			x = magSpecPad[i];
			if (x > magSpecPad[i-2] && x > magSpecPad[i-1] && x > magSpecPad[i+1] && x > magSpecPad[i+2]) {
				out.peaks = out.peaks.concat(I);
			}
		}

		out.inflRegStart = new Array(out.peaks.length); 

		out.inflRegStart[0] = 0;
		for (var i=0; i<out.peaks.length-1; i++) {
			out.inflRegStart[i+1] = Math.ceil((out.peaks[i] + out.peaks[i+1])/2); 
		}

		out.inflRegEnd = new Array(out.peaks.length);
		for (var i=1; i<out.inflRegStart.length; i++) {
			out.inflRegEnd[i-1] = out.inflRegStart[i]-1;
		}
		out.inflRegEnd[out.inflRegEnd.length] = out.inflRegEnd.length-1;

		out.influenceRegions = new Array(out.inflRegStart.length);
		for (var i=0; i<out.influenceRegions.length; i++) 
			out.influenceRegions[i] = Math.max(0, out.inflRegEnd[i] - out.inflRegStart[i] + 1);

		return;
	}

	
	function get_phase_advances(currInPh, prevInPh, omega, RA, RS, instPhaseAdvHop) {
		var twoPI = 2 * Math.PI;

		for (var i=0; i<omega.length; i++) {
			var expectedPhaseAdv = omega[i] * RA;

			var auxHeterodynedPhaseIncr = (currInPh[i] - prevInPh[i]) - expectedPhaseAdv;
			var heterodynedPhaseIncr = auxHeterodynedPhaseIncr - twoPI * Math.round(auxHeterodynedPhaseIncr/twoPI);

			var instPhaseAdvPerSampleHop = omega[i] + heterodynedPhaseIncr / RA;

			instPhaseAdvHop[i] = instPhaseAdvPerSampleHop * RS;
		}

		return;
	}

	
	function get_phasor_theta(currInPh, prevOutPh, instPhaseAdv, frequencyBins, inflRegs, theta) {
		// Get the peaks in the spectrum together with their regions of influence.

		var theta_idx = 0;
		for (var i=0; i<frequencyBins.length; i++) {
			var bin = frequencyBins[i];
			for (var j=0; j<inflRegs[i]; j++, theta_idx++) {
				theta[theta_idx] = prevOutPh[bin] + instPhaseAdv[bin] - currInPh[bin];
			}
		}

		var remainingLength = theta.length - theta_idx;
		for (var i=0; i<remainingLength; i++, theta_idx++)
			theta[theta_idx] = 0;

		return;
	}

	
	function identity_phase_locking(currInMag, currInPh, prevOutPh, instPhaseAdv, phTh) {
		var _ = this; var r = {};

		find_peaks(currInMag, r);

		get_phasor_theta(currInPh, prevOutPh, instPhaseAdv, r.peaks, r.influenceRegions, phTh);

		return;
	}


	function pv_step(fftObj, prevInPh, prevOutPh, omega, RA, RS, out) {

		var currInPh = fftObj.phase;

		var instPhaseAdv = new Float32Array(omega.length);
		get_phase_advances(currInPh, prevInPh, omega, RA, RS, instPhaseAdv);

		var currInMag = fftObj.magnitude;

		var PhTh = new Float32Array(currInPh.length);
		identity_phase_locking(currInMag, currInPh, prevOutPh, instPhaseAdv, PhTh);

		var dblSize = (PhTh.length-1)*2;
		// out.real = new Float32Array(dblSize);
		// out.imag = new Float32Array(dblSize);
		// out.phase = new Float32Array(dblSize);
		// out.magnitude = new Float32Array(dblSize);
		var sqrt = Math.sqrt; var cos = Math.cos;
		var sin = Math.sin; var atan2 = Math.atan2;

		
		for (var i=0; i<PhTh.length; i++) {
			var theta = PhTh[i];

			var PhThRe = cos(theta);
			var PhThIm = sin(theta);
			out.real[i] = PhThRe * fftObj.real[i] - PhThIm * fftObj.imag[i];
			out.imag[i] = PhThRe * fftObj.imag[i] + PhThIm * fftObj.real[i];
			out.phase[i] = atan2(out.imag[i], out.real[i]);
			out.magnitude[i] = sqrt(out.imag[i]*out.imag[i] + out.real[i]*out.real[i]);

			if (i>0) {
				var idx = dblSize - 1;
				out.real[idx] = out.real[i];
				out.imag[idx] = -out.imag[i];
				out.phase[idx] = atan2(out.imag[idx], out.real[idx]);
				out.magnitude[idx] = sqrt(out.imag[idx]*out.imag[idx] + out.real[idx]*out.real[idx]);
			}

		}
		
		return;
	}


	this.process = function(inputFrame) {

		var _ = this;

		var __RS = _RS;
		var __RA = _RA;

		// ----------------------------------
		// ----------ANALYSIS STEP-----------
		// ----------------------------------
		
		var processedFrame = [];

		if (_first) {
			var fftObj = {
				real: new Float32Array(_winSize), 
				imag: new Float32Array(_winSize), 
				magnitude: new Float32Array(_winSize), 
				phase: new Float32Array(_winSize)
			};
			_.STFT(inputFrame, _framingWindow, _winSize, fftObj);
			_previousOutputPhase = fftObj.phase;
			_previousInputPhase = fftObj.phase;
			processedFrame = new Float32Array(fftObj.real.length);
			_.ISTFT(fftObj.real, fftObj.imag, _framingWindow, false, processedFrame);
		} else {
			var hlfSize = Math.round(_winSize/2)+1;
			var fftObj = {
				real: new Float32Array(hlfSize), 
				imag: new Float32Array(hlfSize), 
				magnitude: new Float32Array(hlfSize), 
				phase: new Float32Array(hlfSize)
			};
			var pvOut = {
				real: new Float32Array(_winSize), 
				imag: new Float32Array(_winSize), 
				magnitude: new Float32Array(_winSize), 
				phase: new Float32Array(_winSize)
			};
			_.STFT(inputFrame, _framingWindow, hlfSize, fftObj);
			pv_step(fftObj, _previousInputPhase, _previousOutputPhase, _omega, __RA, __RS, pvOut);
			_previousOutputPhase = pvOut.phase;
			_previousInputPhase = fftObj.phase;
			processedFrame = new Float32Array(pvOut.real);
			_.ISTFT(pvOut.real, pvOut.imag, _framingWindow, false, processedFrame);
		}

		_first = false;


		// ----------------------------------
		// ------OVERLAP AND SLIDE STEP------
		// ----------------------------------
		var outputFrame = [];
		overlap_and_slide(__RS, processedFrame, _overlapBuffers, _winSize, outputFrame);
		var owFrame = [];
		overlap_and_slide(__RS, _squaredFramingWindow, _owOverlapBuffers, _winSize, owFrame);

		for (var i=0; i<outputFrame.length; i++)
			outputFrame[i] = outputFrame[i] / ((owFrame[i]<10e-3)? 1 : owFrame[i]);

		return outputFrame;

	}

	
	this.STFT = function(inputFrame, windowFrame, wantedSize, out) {
		var winSize = windowFrame.length;
		var _inputFrame = new Array(winSize);
		var fftFrame = new Array(2*winSize);

		for (var i=0; i<winSize; i++) {
			_inputFrame[i] = inputFrame[i] * windowFrame[i];
		}
		var fft = new FFT.complex(winSize, false);
		fft.simple(fftFrame, _inputFrame, 'real');

		// out.real = new Array(Math.min(winSize,wantedSize));
		// out.imag = new Array(Math.min(winSize,wantedSize));
		// out.magnitude = new Array(Math.min(winSize,wantedSize));
		// out.phase = new Array(Math.min(winSize,wantedSize));

		for (var p=0; p<winSize && p<wantedSize; p++) {
			var real = out.real; var imag = out.imag;
			var phase = out.phase; var magnitude = out.magnitude;
			real[p] = fftFrame[2*p];
			imag[p] = fftFrame[2*p+1];
			magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
			phase[p] = Math.atan2(imag[p], real[p]);
		}

		return;
	}

	this.ISTFT = function(real, imaginary, windowFrame, restoreEnergy, output2) {
		var input = new Array(2 * real.length);
		var output1 = new Array(2 * real.length);

		for (var i=0; i<real.length; i++) {
			input[2*i] = real[i];
			input[2*i+1] = imaginary[i];
		}

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

		return;
	}


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

	this.create_omega_array = function(size) {
		return Array.apply(null, Array(size/2 + 1)).map(function (x, i) { 
			return 2 * Math.PI * i / size;
		});
	}
	
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


	this.get_previous_input_phase = function() {
		return _previousInputPhase;
	}

	this.get_previous_output_phase = function() {
		return _previousOutputPhase;
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