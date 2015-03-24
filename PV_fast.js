function PhaseVocoder(winSize, sampleRate) {

	var _sampleRate = sampleRate; var _RS = 0; var _RA = 0; var _omega;

	var _previousInputPhase; var _previousOutputPhase; var _framingWindow;
	
	var _squaredFramingWindow; var _winSize = winSize;

	var _overlapBuffers; var _owOverlapBuffers;

	var _first = true;

	var _fftProcessor = new FFT.complex(winSize, false);

	var _ifftProcessor = new FFT.complex(winSize, true);

	var _overlapFactor = 4;

	var _lastInputAlpha = 1;


	
	function overlap_and_slide(RS, inF, squaredWinF, oBuf, owOBuf, windowSize, outF) {

		for (var i=0; i<RS; i++) {
			var owSample = owOBuf.shift();
			outF[i] = oBuf.shift() / ((owSample<10e-3)? 1 : owSample);
			oBuf[oBuf.length] = owOBuf[owOBuf.length] = 0;
		}

		for (var i=0; i<windowSize; i++) {
			oBuf[oBuf.length-1] = inF[i] + oBuf.shift();
			owOBuf[owOBuf.length-1] = squaredWinF[i] + owOBuf.shift();
		}
	}


	function find_peaks(magFrame, out) {
		var msp = [0,0].concat(magFrame).concat([0,0]);
		out.peaks = [];
		var aux = new Float32Array(magFrame.length);

		for (var i=2, I=0; i<=msp.length-2; i++, I++) {
			x = msp[i];
			if (x > msp[i-2] && x > msp[i-1] && x > msp[i+1] && x > msp[i+2]) {
				out.peaks = out.peaks.concat(I);
				aux[I] = I;
			}
		}

		out.inflRegStart = new Array(out.peaks.length); 
		out.inflRegEnd = new Array(out.peaks.length);
		out.inflRegs = new Array(out.peaks.length);

		out.inflRegStart[0] = 0;
		for (var i=1; i<out.peaks.length; i++) {
			out.inflRegStart[i] = Math.ceil((out.peaks[i-1] + out.peaks[i])/2); 
			out.inflRegEnd[i-1] = out.inflRegStart[i]-1;
			out.inflRegs[i-1] = Math.max(0, out.inflRegEnd[i-1] - out.inflRegStart[i-1] + 1);
		}
		out.inflRegEnd[out.inflRegEnd.length] = out.inflRegEnd.length-1;

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

		return;
	}

	
	function identity_phase_locking(currInMag, currInPh, prevOutPh, instPhaseAdv, phTh) {
		var r = {};

		find_peaks(currInMag, r);

		get_phasor_theta(currInPh, prevOutPh, instPhaseAdv, r.peaks, r.inflRegs, phTh);

		return;
	}


	function pv_step(fftObj, prevInPh, prevOutPh, omega, RA, RS, out) {

		var currInPh = fftObj.phase;
		var currInMag = fftObj.magnitude;
		var instPhaseAdv = new Float32Array(omega.length);
		var phTh = new Float32Array(currInPh.length);

		get_phase_advances(currInPh, prevInPh, omega, RA, RS, instPhaseAdv);

		identity_phase_locking(currInMag, currInPh, prevOutPh, instPhaseAdv, phTh);

		var dblSize = (phTh.length-1)*2;
		var sqrt = Math.sqrt; var cos = Math.cos;
		var sin = Math.sin; var atan2 = Math.atan2;

		
		for (var i=0; i<phTh.length; i++) {
			var theta = phTh[i];

			var phThRe = cos(theta);
			var phThIm = sin(theta);
			out.real[i] = phThRe * fftObj.real[i] - phThIm * fftObj.imag[i];
			out.imag[i] = phThRe * fftObj.imag[i] + phThIm * fftObj.real[i];
			out.phase[i] = atan2(out.imag[i], out.real[i]);

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
			// IF I USE Float32Array FOR THE fftObj, I GET "PHASEY" ARTIFACTS.
			var fftObj = {
				real: new Array(_winSize), 
				imag: new Array(_winSize), 
				magnitude: new Array(_winSize), 
				phase: new Array(_winSize)
			};
			_.STFT(inputFrame, _framingWindow, _winSize, fftObj);
			_previousOutputPhase = fftObj.phase;
			_previousInputPhase = fftObj.phase;
			processedFrame = new Array(fftObj.real.length);
			_.ISTFT(fftObj.real, fftObj.imag, _framingWindow, false, processedFrame);
		} else {
			var hlfSize = Math.round(_winSize/2)+1;
			// IF I USE Float32Array for the fftObj, I get "phasey" artifacts.
			var fftObj = {
				real: new Array(hlfSize), 
				imag: new Array(hlfSize), 
				magnitude: new Array(hlfSize), 
				phase: new Array(hlfSize)
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
			processedFrame = new Array(pvOut.real);
			_.ISTFT(pvOut.real, pvOut.imag, _framingWindow, false, processedFrame);
		}

		_first = false;


		// ----------------------------------
		// ------OVERLAP AND SLIDE STEP------
		// ----------------------------------
		var outputFrame = new Array(__RS);

		overlap_and_slide(__RS, processedFrame, _squaredFramingWindow, _overlapBuffers, _owOverlapBuffers, _winSize, outputFrame);

		return outputFrame;

	}


	this.processv2 = function(fftObj) {

		var _ = this;

		var __RS = _RS;
		var __RA = _RA;

		// ----------------------------------
		// ----------ANALYSIS STEP-----------
		// ----------------------------------
		
		var processedFrame = [];

		if (_first) {
			// IF I USE Float32Array FOR THE fftObj, I GET "PHASEY" ARTIFACTS.
			fftObj.magnitude = new Array(_winSize);
			fftObj.phase = new Array(_winSize);
			for (var i=0; i<_winSize; i++) {
				var real = fftObj.real; var imag = fftObj.imag;
				var phase = fftObj.phase; var magnitude = fftObj.magnitude;
				magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
				phase[p] = Math.atan2(imag[p], real[p]);
			}
			_previousOutputPhase = fftObj.phase;
			_previousInputPhase = fftObj.phase;
			processedFrame = new Array(fftObj.real.length);
			_.ISTFT(fftObj.real, fftObj.imag, _framingWindow, false, processedFrame);
		} else {
			var hlfSize = Math.round(_winSize/2)+1;
			// IF I USE Float32Array for the fftObj, I get "phasey" artifacts.
			fftObj.magnitude = new Array(hlfSize);
			fftObj.phase = new Array(hlfSize);
			for (var i=0; i<hlfSize; i++) {
				var real = fftObj.real; var imag = fftObj.imag;
				var phase = fftObj.phase; var magnitude = fftObj.magnitude;
				magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
				phase[p] = Math.atan2(imag[p], real[p]);
			}
			var pvOut = {
				real: new Float32Array(_winSize), 
				imag: new Float32Array(_winSize), 
				magnitude: new Float32Array(_winSize), 
				phase: new Float32Array(_winSize)
			};
			pv_step(fftObj, _previousInputPhase, _previousOutputPhase, _omega, __RA, __RS, pvOut);
			_previousOutputPhase = pvOut.phase;
			_previousInputPhase = fftObj.phase;
			processedFrame = new Array(pvOut.real);
			_.ISTFT(pvOut.real, pvOut.imag, _framingWindow, false, processedFrame);
		}

		_first = false;


		// ----------------------------------
		// ------OVERLAP AND SLIDE STEP------
		// ----------------------------------
		var outputFrame = new Array(__RS);

		overlap_and_slide(__RS, processedFrame, _squaredFramingWindow, _overlapBuffers, _owOverlapBuffers, _winSize, outputFrame);

		return outputFrame;

	}

	
	this.STFT = function(inputFrame, windowFrame, wantedSize, out) {
		var winSize = windowFrame.length;
		var _inputFrame = new Array(winSize);
		var fftFrame = new Array(2*winSize);

		for (var i=0; i<winSize; i++) {
			_inputFrame[i] = inputFrame[i] * windowFrame[i];
		}
		
		_fftProcessor.simple(fftFrame, _inputFrame, 'real');

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

	this.STFTv2 = function(inputFrame, windowFrame, wantedSize, out) {
		var winSize = windowFrame.length;
		var _inputFrame = new Array(winSize);
		var fftFrame = new Array(2*winSize);

		for (var i=0; i<winSize; i++) {
			_inputFrame[i] = inputFrame[i] * windowFrame[i];
		}
		
		_fftProcessor.simple(fftFrame, _inputFrame, 'real');

		for (var p=0; p<winSize && p<wantedSize; p++) {
			var real = out.real; var imag = out.imag;
			real[p] = fftFrame[2*p];
			imag[p] = fftFrame[2*p+1];
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

		_ifftProcessor.simple(output1, input, 'complex');

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

		_omega = create_omega_array(winSize);

		_previousInputPhase = create_constant_array(winSize/2, 0);
		_previousOutputPhase = create_constant_array(winSize/2, 0);

		_framingWindow = create_sin_beta_window_array(winSize, 1);

		_squaredFramingWindow = _framingWindow.map(function(x,i){ return x*x; });

		_overlapBuffers = create_constant_array(winSize, 0);
		_owOverlapBuffers = create_constant_array(winSize, 0);

		this.set_alpha(1);
	}

	function create_omega_array(size) {
		return Array.apply(null, Array(size/2 + 1)).map(function (x, i) { 
			return 2 * Math.PI * i / size;
		});
	}
	
	function create_sin_beta_window_array(size, beta) {
		return Array.apply(null, Array(size)).map(function(x,i){
			return Math.pow(Math.sin(Math.PI*i/size), beta);
		});
	}

	function create_constant_array(size, constant) {
		return Array.apply(null, Array(size)).map(function () { 
			return constant; 
		});
	}

	this.reset_phases_and_overlap_buffers = function() {

		_previousInputPhase = create_constant_array(winSize/2, 0);
		_previousOutputPhase = create_constant_array(winSize/2, 0);

		_overlapBuffers = create_constant_array(winSize, 0);
		_owOverlapBuffers = create_constant_array(winSize, 0);

		_first = true;
	}

	this.reset_phases = function() {

		_previousInputPhase = create_constant_array(winSize/2, 0);
		_previousOutputPhase = create_constant_array(winSize/2, 0);

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
		_lastInputAlpha = newAlpha;
		_RA = Math.round(_winSize/_overlapFactor);
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

	this.get_specified_alpha = function() {
		return _lastInputAlpha;
	}

	this.set_overlap_factor = function(overlapFactor) {
		_overlapFactor = overlapFactor;
		this.set_alpha(_lastInputAlpha);
	}
}