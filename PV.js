function PhaseVocoder(winSize, sampleRate, alpha) {

	var _fftProcessor = new FFT(winSize, sampleRate);

	var _RS = winSize/4;
	var _RA = Math.round(_RS / alpha);	
	// var _RA = winSize/4;
	// var _RS = Math.round(alpha * _RA);
	var _omega = createOmegaArray(winSize);

	var _previousInputPhase = createConstantArray(winSize/2, 0);
	var _previousOutputPhase = createConstantArray(winSize/2, 0);
	var _framingWindow = createSinBetaWindowArray(winSize, 1);
	var _squaredFramingWindow = _framingWindow.map(function(x,i){ return x*x; });
	var _winSize = winSize;

	var _overlapBuffers = []; 
	_overlapBuffers[0] = []; // current frame
	_overlapBuffers[1] = createConstantArray(3*winSize/4, 0);; // -1 previous frame 
	_overlapBuffers[2] = createConstantArray(winSize/2, 0);; // -2 previous frame 
	_overlapBuffers[3] = createConstantArray(winSize/4, 0); // -3 previous frame 

	// RUN THIS FOR THE JAVA-PORT OVERLAPPING METHOD
	//_overlapBuffers = createConstantArray(winSize, 0);

	var _owOverlapBuffers = createConstantArray(winSize, 0);
	_overlapBuffers[0] = []; // current frame
	_overlapBuffers[1] = createConstantArray(3*winSize/4, 0);; // -1 previous frame 
	_overlapBuffers[2] = createConstantArray(winSize/2, 0);; // -2 previous frame 
	_overlapBuffers[3] = createConstantArray(winSize/4, 0); // -3 previous frame 

	var _first = true;



	function concat(a,b) {
		return a.push.apply(a,b);
	}


	/*
	 *  TODO
	 *
	 *  @param buffers: TODO
	 *
	 *  @returns TODO
	 */
	function processOverlapBuffers(RS, frame, buffers) {
		var output = createConstantArray(RS, 0);

		var output = output.map(function(x, i) {
			var vb1 = buffers[1].shift();
			var vb2 = buffers[2].shift();
			var vb3 = buffers[3].shift();
			buffers[1].push(frame[i]);
			buffers[2].push(vb1);
			buffers[3].push(vb2);
			return frame[i] + vb1 + vb2 + vb3;
		});

		return output;
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
	function overlapAndSlide(RS, frame, overlapBuffer, windowSize) {
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



	function fftShift(data) {
		if (data.length % 2 == 0) {
			fftShiftEven(data);
		} else {
		  	fftShiftOdd(data);
		}
	}	


	function fftShiftOdd(data) {
		var shiftAmt = data.length / 2;
		var remaining = data.length;
		var curr = 0;
		while (remaining >= 0) {
			var next = data[(curr + shiftAmt) % data.length];
			data[(curr+shiftAmt) % data.length] = save;
			save = next;
			curr = (curr + shiftAmt) % data.length;
			remaining--;
		}
	}


	function fftShiftEven(data) {
		for (var i = 0; i < data.length / 2; i++) {
	      var tmp = data[i];
	      data[i] = data[i + data.length / 2];
	      data[i + data.length / 2] = tmp;
	    }
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
	function createOmegaArray(size) {
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
	function createSinBetaWindowArray(size, beta) {
		return Array.apply(null, Array(size)).map(function(x,i){
			return Math.pow(Math.sin(Math.PI*i/size), beta);
		});
	}


	function createConstantArray(size, constant) {
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
	function findPeaks(magFrame) {

		var magSpecPad = [0,0].concat(Array.prototype.slice.call(magFrame)).concat([0,0]);

		var peaks = magSpecPad.slice(2,magSpecPad.length-2).map(function(x,i){
			var I = i + 2;
			if(x > magSpecPad[I-2] && x > magSpecPad[I-1] && x > magSpecPad[I+1] && x > magSpecPad[I+2]) {
				return i;
			}
		}).filter(function(x){ 
			return x!=undefined && x!=null; 
		});

		var inflRegStart = []; var inflRegEnd = [];

		inflRegStart = inflRegStart.concat(0).concat(peaks.slice(0,peaks.length-1).map(function(x, i){ 
			return Math.ceil((x + peaks[i+1])/2); 
		}));

		inflRegEnd = inflRegEnd.concat(inflRegStart.slice(1,inflRegStart.length).map(function(x,i){
			return x - 1;
		})).concat(inflRegEnd.length-1);

		return {peaks: peaks, inflRegionStart: inflRegStart, inflRegionEnd: inflRegEnd};
	}

	function findPeaksV2(magFrame, numberOfNeighbors) {
		var zeros = createConstantArray(Math.round(numberOfNeighbors/2), 0);
		var magSpecPad = zeros.concat(Array.prototype.slice.call(magFrame)).concat(zeros);

		var peaks = magSpecPad.slice(zeros.length,magSpecPad.length-zeros.length).map(function(x,i){
			var I = i + zeros.length;
			for (var j=-zeros.length; j<=zeros.length; j++) 
				if (x < magSpecPad[I+j])
					return null;
			return i;
		}).filter(function(x){ 
			return x!=undefined && x!=null; 
		});

		var inflRegStart = []; var inflRegEnd = [];

		inflRegStart = inflRegStart.concat(0).concat(peaks.slice(0,peaks.length-1).map(function(x, i){ 
			return Math.ceil((x + peaks[i+1])/2); 
		}));

		inflRegEnd = inflRegEnd.concat(inflRegStart.slice(1,inflRegStart.length).map(function(x,i){
			return x - 1;
		})).concat(inflRegEnd.length-1);

		return {peaks: peaks, inflRegionStart: inflRegStart, inflRegionEnd: inflRegEnd};
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
	function get_phase_advances(currentInputPhase, previousInputPhase, omega, RA, RS) {

		// Expected phase advances from the last to the current input frame.
		var dphi = omega.map(function(x,i){ return x * RA;  });

		// Heterodyned phase increments.
		var twoPI = 2 * Math.PI; // Just to avoid a lot of duplicate multiplications.
		//var hpi = dphi.map(function(x,i){ return (currentInputPhase[i] - previousInputPhase[i]) * x; });
		// Reduce to the range -pi:pi
		//var hpi = hpi.map(function(x,i){ return x - twoPI * Math.round(x/twoPI); });

		// Heterodyned phase increments and reduce hpi to the range -pi:pi.
		var hpi = dphi.map(function(_dphi,i){
			var _hpi = (currentInputPhase[i] - previousInputPhase[i]) - _dphi;
			return _hpi - twoPI * Math.round(_hpi/twoPI);
		});

		// Instantaneous phase advances per sample.
		var ipa_sample = hpi.map(function(x,i){ return omega[i] + x / RA; });

		// Instantaneous phase advances per synthesis hopsize.
		var ipa_hop = ipa_sample.map(function(x,i){ return x * RS; });

		return ipa_hop;
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
	function get_phasor_theta(currentInputPhase, previousOutputPhase, instPhaseAdv, frequencyBins, influenceRegions) {
		// Get the peaks in the spectrum together with their regions of influence.
		var theta = [];

		frequencyBins.map(function(bin,i){
			var new_theta = Array.apply(null, Array(influenceRegions[i])).map(function (x, i) { 
				return previousOutputPhase[bin] + instPhaseAdv[bin] - currentInputPhase[bin];
			});
			theta = theta.concat(new_theta);
		});

		theta = theta.concat(createConstantArray(currentInputPhase.length - theta.length, 0));

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
	function identity_phase_locking(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {
		//var r = findPeaks(currentInputMagnitude);
		var r = findPeaksV2(currentInputMagnitude,6);

		var influenceRegions = r.inflRegionStart.map(function(inflRegStart,i){ return Math.max(0, r.inflRegionEnd[i] - inflRegStart + 1); });

		var phasor_theta = get_phasor_theta(currentInputPhase, previousOutputPhase, instPhaseAdv, r.peaks, influenceRegions);

		return phasor_theta;
	}


	/*
	 *  TODO
	 *
	 *  @param currentInputMagnitude: Float32Array holding the magnitude values 
	 * of the current input frame.
	 *  @param currentInputMagnitude: Float32Array holding the phase values of 
	 * the current input frame.
	 *  @param previousOutputPhase: Float32Array array with the phases of the last output frame.
	 *  @param instPhaseAdv: Float32Array array with the instantaneous phases advance.
	 *
	 *  @returns phasor: float array with the phasor angles.
	 */
	function no_phase_locking(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {

		var phasor_theta = Array.prototype.slice.call(currentInputPhase).map(function(x,i){
			return previousOutputPhase[i] + instPhaseAdv[i] - x;
		});

		return phasor_theta;
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
	function pv_step_v2(fftObject, previousInputPhase, previousOutputPhase, omega, RA, RS) {

		var currentInputPhase = fftObject.angle;

		var instPhaseAdv = get_phase_advances(currentInputPhase, previousInputPhase, omega, RA, RS);

		var currentInputMag = fftObject.spectrum;

		var phasor_theta = identity_phase_locking(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);
		//var phasor_theta = no_phase_locking(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);

		// Multiplication of two vectors of complex numbers
		// var out_real = new Float32Array(phasor_theta.length);
		// var out_imag = new Float32Array(phasor_theta.length);
		var out_real = new Array(phasor_theta.length);
		var out_imag = new Array(phasor_theta.length);
		var out_phase = new Array(phasor_theta.length);
		var doubleSize = phasor_theta.length*2-1;
		for (var i=0; i<phasor_theta.length; i++) {
			var theta = phasor_theta[i];

			var phasor_theta_real = Math.cos(theta);
			var phasor_theta_imag = Math.sin(theta);
			out_real[i] = phasor_theta_real * fftObject.real[i] - phasor_theta_imag * fftObject.imag[i];
			out_imag[i] = phasor_theta_real * fftObject.imag[i] + phasor_theta_imag * fftObject.real[i];

			out_real[doubleSize-i] = out_real[i];
			out_imag[doubleSize-i] = -out_imag[i];

			out_phase[i] = Math.atan2(out_imag[i], out_real[i]);
		}
		
		return {real: out_real, imag: out_imag, phase: out_phase};
	}


	this.process = function(inputFrame) {

		var _inputFrame = inputFrame.map(function(sample,i){
			return sample * _framingWindow[i];
		});

		var outputFrame = [];

		_fftProcessor.forward(_inputFrame);
		_fftProcessor.calculateAngle();

		var processedFrame = [];

		if (_first) {

			_previousOutputPhase = _fftProcessor.angle;
			_previousInputPhase = _fftProcessor.angle;

			processedFrame = _fftProcessor.inverse();

			_first = !_first;

		} else {

			var out = pv_step_v2(_fftProcessor, _previousInputPhase, _previousOutputPhase, _omega, _RA, _RS);

			_previousOutputPhase = out.phase;
			_previousInputPhase = _fftProcessor.angle;

			processedFrame = _fftProcessor.inverse(out.real, out.imag);

		}

		processedFrame = Array.prototype.slice.call(processedFrame).map(function(sample,i){
			return sample * _framingWindow[i];
		});

		outputFrame = overlapAndSlide(_RS, processedFrame, _overlapBuffers, _winSize);

		owFrame = overlapAndSlide(_RS, _squaredFramingWindow, _owOverlapBuffers, _winSize);

		// outputFrame = processOverlapBuffers(_RS, processedFrame, _overlapBuffers);

		// owFrame = processOverlapBuffers(_RS, _squaredFramingWindow, _owOverlapBuffers);

		outputFrame.map(function(sample, i){
			return sample / (owFrame[i]<10e-3)?1:owFrame[i];
		});

		return outputFrame;

	}

	this.reset = function() {
		// TODO
	}
}