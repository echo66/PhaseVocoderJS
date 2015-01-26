function PhaseVocoder(winSize, sampleRate, alpha) {

	var _position = 0;

	var _fftProcessor = new FFT(winSize, sampleRate);

	var _RS = 0;
	var _RA = 0;
	var _omega = createOmegaArray(winSize);

	var _previousInputPhase = createConstantArray(winSize/2, 0);
	var _previousOutputPhase = createConstantArray(winSize/2, 0);
	var _framingWindow = createSinBetaWindowArray(winSize, 1);
	var _squaredFramingWindow = _framingWindow.map(function(x,i){ return x*x; });
	
	var _winSize = winSize;

	var _overlapBuffers = createConstantArray(winSize, 0);

	var _owOverlapBuffers = createConstantArray(winSize, 0);

	var _first = true;


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
	function findPeaksV3(magFrame) {
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
	function get_phase_advances_v4(currentInputPhase, previousInputPhase, omega, RA, RS) {
		var twoPI = 2 * Math.PI;
		var ipa_hop = new Array(omega.length);

		for (var i=0; i<omega.length; i++) {
			var dphi = omega[i] * RA;

			var auxHpi = (currentInputPhase[i] - previousInputPhase[i]) - dphi;
			var hpi = auxHpi - twoPI * Math.round(auxHpi/twoPI);

			var ipa_sample = omega[i] + hpi / RA;

			ipa_hop[i] = ipa_sample * RS;
		}

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
	function get_phasor_theta_v2(currentInputPhase, previousOutputPhase, instPhaseAdv, frequencyBins, influenceRegions) {
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
	function identity_phase_locking_v2(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {
		var r = findPeaks(currentInputMagnitude);

		var influenceRegions = new Array(r.inflRegionStart.length);

		for (var i=0; i<influenceRegions.length; i++) 
			influenceRegions[i] = Math.max(0, r.inflRegionEnd[i] - r.inflRegionStart[i] + 1);

		var phasor_theta = get_phasor_theta_v2(currentInputPhase, previousOutputPhase, instPhaseAdv, r.peaks, influenceRegions);

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

		var instPhaseAdv = get_phase_advances_v4(currentInputPhase, previousInputPhase, omega, RA, RS);

		var currentInputMag = fftObject.spectrum;

		var phasor_theta = identity_phase_locking_v2(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);

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

		var fftObject = this.STFT(inputFrame, _framingWindow);
		//_position += _RA;

		var outputFrame = [];

		var processedFrame = [];

		var out = (_first)? null : pv_step_v2(_fftProcessor, _previousInputPhase, _previousOutputPhase, _omega, _RA, _RS);
		_previousOutputPhase = (_first)? _fftProcessor.angle : out.phase;
		_previousInputPhase = _fftProcessor.angle;
		processedFrame = (_first)? _fftProcessor.inverse() : _fftProcessor.inverse(out.real, out.imag);
		_first = false;

		processedFrame = _framingWindow.process(processedFrame);

		//fftShift(processedFrame);

		outputFrame = overlapAndSlide(_RS, processedFrame, _overlapBuffers, _winSize);

		owFrame = overlapAndSlide(_RS, _squaredFramingWindow, _owOverlapBuffers, _winSize);

		outputFrame = outputFrame.map(function(sample, i){
			return sample / ((owFrame[i]<10e-3)? 10e-3 : owFrame[i]);
		});

		return outputFrame;

	}

	/*
	 * TODO: I NEED TO TEST THIS.
	 */
	this.processAll = function(data) {
		var Y = [];
		var OW = [];
		var win = _framingWindow.process(createConstantArray(winSize, 1));
		var numFrames = Math.floor((data.length - _winSize)/_RA + 1);

		for (var i=0, anaPos=0, synPos=0; i<numFrames; i++, anaPos+=_RA, synPos+=_RS) {
			console.log('Frame ' + (i+1) + '/' + numFrames);
			// Analysis Step
			var frame = _framingWindow.process(data.subarray(anaPos,anaPos+_winSize));
			_fftProcessor.forward(frame);
			_fftProcessor.calculateSpectrumAndPhase();
			// -------------

			// Phase Correction Step
			var out = (_first)? null : pv_step_v2(_fftProcessor, _previousInputPhase, _previousOutputPhase, _omega, _RA, _RS);
			_previousOutputPhase = (_first)? _fftProcessor.angle : out.phase;
			_previousInputPhase = _fftProcessor.angle;
			processedFrame = (_first)? _fftProcessor.inverse() : _fftProcessor.inverse(out.real, out.imag);
			processedFrame = _framingWindow.process(processedFrame);
			// ---------------------

			// Overlapping  Step
			processedFrame = Array.prototype.slice.call(processedFrame);
			if (_first) {
				Y = Y.concat(processedFrame);
				OW = OW.concat(win);
			} else {
				var overlappedFrames = new Array(processedFrame.length);
				var overlappedWindows = new Array(processedFrame.length);
				for (var j=0; j<processedFrame.length; j++) {
					overlappedWindows[j] = win[j] + ((OW[synPos+j])? OW[synPos+j] : 0);
					overlappedFrames[j] = (processedFrame[j] + ((Y[synPos+j])? Y[synPos+j] : 0)) / overlappedWindows[j];
				}
				Y.splice(synPos);
				Y = Y.concat(overlappedFrames);
				OW.splice(synPos);
				OW = OW.concat(overlappedWindows);
			}
			// -----------------

			_first = false;

		}

		return Y;
	}

	this.reset = function() {

		_previousInputPhase = createConstantArray(winSize/2, 0);
		_previousOutputPhase = createConstantArray(winSize/2, 0);

		_overlapBuffers = createConstantArray(winSize, 0);
		_owOverlapBuffers = createConstantArray(winSize, 0);

		_first = true;
	}

	this.setAlpha = function(newAlpha) {
		this.reset();
		_RA = _winSize/4;
		_RS = Math.round(newAlpha * _RA);
	}

	this.setPosition = function(newPosition) {
		this.reset();
		_position = newPosition;
	}

	this.STFT = function(inputFrame, windowFrame) {
		var winSize = windowFrame.length;
		var _inputFrame = new Array(winSize);
		var fftFrame = new Array(2*winSize);

		for (var i=0; i<winSize; i++) {
			myData1[i] = inputFrame[i] * windowFrame[i];
		}
		var fft = new FFT.complex(winSize, false);
		fft.simple(fftFrame, _inputFrame, 'real');

		var real = new Array(winSize);
		var imag = new Array(winSize);
		var magnitude = new Array(winSize);
		var phase = new Array(winSize);

		for (var i=0; i<winSize; i++) {
			real[p] = fftFrame[2*p];
			imag[p] = fftFrame[2*p+1];
			magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
			phase[p] = Math.atan2(imag[p], real[p]);
		}

		return {
			real: real,
			imag: imag,
			magnitude: magnitude,
			angle: phase
		}
	}

	this.ISTFT = function(real, imaginary, windowFrame) {
		var input = new Array(2 * windowFrame.length);
		var output1 = new Array(2 * windowFrame.length);
		var output2 = new Array(2 * windowFrame.length);

		for (var i=0; i<2048; i++) {
			input[2*i] = real[i];
			input[2*i+1] = imaginary[i];
		}

		var ifft = new FFT.complex(2048, true);

		ifft.simple(output1, input);

		for (var i=0; i<2048; i++) {
			output2[i] = output1[2*i] * windowFrame[i] / windowFrame.length;
		}

		return output2;
	}
}