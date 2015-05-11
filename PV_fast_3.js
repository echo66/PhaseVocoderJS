function PhaseVocoder(winSize, sampleRate) {

	var _sampleRate = sampleRate; var _RS = 0; var _RA = 0; var _omega;

	var _previousInputPhase; var _previousOutputPhase; var _framingWindow;
	
	var _squaredFramingWindow; var _winSize = winSize;

	var _overlapBuffers; var _owOverlapBuffers;

	var _first = true;

	var _overlapFactor = 4;

	var _lastInputAlpha = 1;

	var stdlib = {
	    Math: Math,
	    Float32Array: Float32Array,
	    Float64Array: Float64Array
	};

	var dromHeap = new ArrayBuffer(32*_winSize);
	var dromFFT = fourier.custom["fft_f32_"+_winSize+"_asm"](stdlib, null, dromHeap);
	dromFFT.init();



	//--------------------------------------------------
	//------------PRE-ALLOCATED MEMORY------------------
	//--------------------------------------------------

	//find_peaks
	var _hlfSize = Math.round(_winSize/2)+1;

	var _find_peaks = {
		msp : create_constant_array(_hlfSize+4, 0, Float32Array), 
		aux : create_constant_array(_hlfSize, 0, Uint8Array), 
		peaks : create_constant_array(_hlfSize, 0.1, Float32Array), 
		inflRegStart : create_constant_array(_hlfSize, 1, Uint8Array), 
		inflRegEnd : create_constant_array(_hlfSize, 1, Uint8Array), 
		inflRegs : create_constant_array(_hlfSize, 1, Uint8Array)
	};
	

	// // process
	var _process = {
		fftObj : {
			real: new Float32Array(_hlfSize), 
			imag: new Float32Array(_hlfSize), 
			magnitude: new Float32Array(_hlfSize), 
			phase: new Float32Array(_hlfSize)
		}, 
		pvOut : {
			real: create_constant_array(_winSize, 0, Float32Array), 
			imag: create_constant_array(_winSize, 0, Float32Array), 
			magnitude: create_constant_array(_winSize, 0, Float32Array), 
			phase: create_constant_array(_winSize, 0, Float32Array)
		},
		processedFrame : new Float32Array(_winSize)
	};

	var _pv_step = {
		instPhaseAdv : new Float32Array(_hlfSize), 
		phTh : new Float32Array(_hlfSize)
	};

	var _STFT = {
		_inputFrame : new Float32Array(_winSize),
		_zeros: new Float32Array(_winSize)
	}
	//--------------------------------------------------
	//--------------------------------------------------
	//--------------------------------------------------


	
	function overlap_and_slide(RS, inF, squaredWinF, oBuf, owOBuf, windowSize, outF) {

		var i = 0;
		var owSample = 0;

		for (var i = 0; i < RS; i++) {
			owSample = owOBuf.shift();
			outF[i] = oBuf.shift() / ((owSample<10e-3)? 1 : owSample);
			oBuf[oBuf.length] = owOBuf[owOBuf.length] = 0;
		}

		for (var i = 0; i < windowSize; i++) {
			oBuf[oBuf.length-1] = inF[i] + oBuf.shift();
			owOBuf[owOBuf.length-1] = squaredWinF[i] + owOBuf.shift();
		}
	}

	
	function get_phase_advances(currInPh, prevInPh, omega, RA, RS, instPhaseAdvHop) {
		var twoPI = 2 * Math.PI;

		for (var i = 0; i < omega.length; i++) {
			var expectedPhaseAdv = omega[i] * RA;

			var auxHeterodynedPhaseIncr = (currInPh[i] - prevInPh[i]) - expectedPhaseAdv;
			var heterodynedPhaseIncr = auxHeterodynedPhaseIncr - twoPI * Math.round(auxHeterodynedPhaseIncr/twoPI);

			var instPhaseAdvPerSampleHop = omega[i] + heterodynedPhaseIncr / RA;

			instPhaseAdvHop[i] = instPhaseAdvPerSampleHop * RS;
		}

		return;
	}

	
	function identity_phase_locking(currInMag, currInPh, prevOutPh, instPhaseAdv, phTh) {

		var peaksLength = 0;

		var msp = _find_peaks.msp;
		msp.set(currInMag, 2);
		for (var i=2, I=0; i<=msp.length-2; i++, I++) {
			x = msp[i];
			if (x > msp[i-2] && x > msp[i-1] && x > msp[i+1] && x > msp[i+2]) {
				_find_peaks.peaks[peaksLength++] = I;
				_find_peaks.aux[I] = I;
			}
		}


		_find_peaks.inflRegStart[0] = 0;
		for (var i=1; i<peaksLength; i++) {
			_find_peaks.inflRegStart[i] = Math.ceil((_find_peaks.peaks[i-1] + _find_peaks.peaks[i])/2); 
			_find_peaks.inflRegEnd[i-1] = _find_peaks.inflRegStart[i]-1;
			_find_peaks.inflRegs[i-1] = Math.max(0, _find_peaks.inflRegEnd[i-1] - _find_peaks.inflRegStart[i-1] + 1);
		}


		var phTh_idx = 0;
		for (var i=0; i<peaksLength; i++) {
			var bin = _find_peaks.peaks[i];
			for (var j=0; j<_find_peaks.inflRegs[i]; j++, phTh_idx++) {
				phTh[phTh_idx] = prevOutPh[bin] + instPhaseAdv[bin] - currInPh[bin];
			}
		}


		return;
	}


	function pv_step(fftObj, prevInPh, prevOutPh, omega, RA, RS, out) {

		var currInPh = fftObj.phase;
		var currInMag = fftObj.magnitude;
		var instPhaseAdv = _pv_step.instPhaseAdv;
		var phTh = _pv_step.phTh;

		get_phase_advances(currInPh, prevInPh, omega, RA, RS, instPhaseAdv);

		identity_phase_locking(currInMag, currInPh, prevOutPh, instPhaseAdv, phTh);

		var dblSize = (phTh.length-1)*2;
		var sqrt = Math.sqrt; var cos = Math.cos;
		var sin = Math.sin; var atan2 = Math.atan2;

		var theta = phTh[1];
		var phThRe = cos(theta);
		var phThIm = sin(theta);
		out.real[1] = phThRe * fftObj.real[1] - phThIm * fftObj.imag[1];
		out.imag[1] = phThRe * fftObj.imag[1] + phThIm * fftObj.real[1];
		out.phase[1] = atan2(out.imag[1], out.real[1]);

		var aux = phTh.length-1;
		var theta = phTh[aux];
		var phThRe = cos(theta);
		var phThIm = sin(theta);
		out.real[aux] = phThRe * fftObj.real[aux] - phThIm * fftObj.imag[aux];
		out.imag[aux] = phThRe * fftObj.imag[aux] + phThIm * fftObj.real[aux];
		out.phase[aux] = atan2(out.imag[aux], out.real[aux]);

		
		for (var i=1; i<phTh.length-1; i++) {
			var theta = phTh[i];

			var phThRe = cos(theta);
			var phThIm = sin(theta);
			
			out.real[i] = phThRe * fftObj.real[i] - phThIm * fftObj.imag[i];
			out.imag[i] = phThRe * fftObj.imag[i] + phThIm * fftObj.real[i];
			out.phase[i] = atan2(out.imag[i], out.real[i]);

			// var aux = dblSize-i;
			// out.real[aux] = out.real[i];
			// out.imag[aux] = -out.imag[i];
			// out.phase[aux] = atan2(out.imag[aux], out.real[aux]);

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
		
		var processedFrame = _process.processedFrame;;

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
			_.ISTFT(fftObj.real, fftObj.imag, _framingWindow, false, processedFrame);
		} else {
			var fftObj = _process.fftObj;
			// FOR SOME REASON, IF I DON'T CREATE A NEW "phase" ARRAY, I GET ARTIFACTS.
			// fftObj.phase = new Float32Array(_hlfSize); 
			var pvOut = _process.pvOut;
			_.STFT(inputFrame, _framingWindow, _hlfSize, fftObj);
			pv_step(fftObj, _previousInputPhase, _previousOutputPhase, _omega, __RA, __RS, pvOut);
			_previousOutputPhase = pvOut.phase;
			// The "phase" issue mentioned above is related to this line. 
			// If I create a new Float array using the phase array, I get no issues.
			_previousInputPhase = new Float32Array(fftObj.phase); 
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
		this.STFT_drom(inputFrame, windowFrame, wantedSize, out);
	}

	this.STFT_drom = function(inputFrame, windowFrame, wantedSize, out) {
		var winSize = windowFrame.length;
		var _inputFrame = _STFT._inputFrame;

		for (var i=0; i<winSize; i++) {
			_inputFrame[i] = inputFrame[i] * windowFrame[i];
		}
		
		// fourier.js forward FFT 
		(new Float32Array(dromHeap)).set(_inputFrame);
		fourier.custom.array2heap(_STFT._zeros, new Float32Array(dromHeap), winSize, winSize);
		dromFFT.transform();
		fourier.custom.heap2array(new Float32Array(dromHeap), out.real, wantedSize, 0);
		fourier.custom.heap2array(new Float32Array(dromHeap), out.imag, wantedSize, winSize);

		for (var p=0; p<winSize && p<wantedSize; p++) {
			var R = out.real; var I = out.imag;
			var P = out.phase; var M = out.magnitude;
			M[p] = Math.sqrt(I[p]*I[p] + R[p]*R[p]);
			P[p] = Math.atan2(I[p], R[p]);
		}

		return;
	}



	this.ISTFT = function(real, imag, windowFrame, restoreEnergy, timeFrame) {
		this.ISTFT_drom(real, imag, windowFrame, restoreEnergy, timeFrame);
	}

	this.ISTFT_drom = function(real, imag, windowFrame, restoreEnergy, timeFrame) {

		var size = windowFrame.length;

		(new Float32Array(dromHeap, 0, size)).set(imag);
		fourier.custom.array2heap(real, new Float32Array(dromHeap), size, size);
		dromFFT.transform();
		timeFrame.set(new Float32Array(dromHeap, size, size));

		for (var i=0; i<size; i++) {
			timeFrame[i] = timeFrame[i] / windowFrame.length;
			timeFrame[i] *= windowFrame[i];
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

	function create_constant_array(size, constant, ArrayType) {
		var arr = new ((ArrayType)?ArrayType:Array)(size);
		for (var i=0; i<size; i++) 
			arr[i] = constant;
		return arr;
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
		if (newAlpha <= 0.8)
			_overlapFactor = 2;
		else if (newAlpha <= 1)
			_overlapFactor = 4;
		else
			_overlapFactor = 5;
		_RA = Math.round(_winSize/_overlapFactor);
		_RS = Math.round(newAlpha * _RA);
		// _RS = _RA;
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