	var eps = 2.2204e-16;

	function restoreEnergy(){
		// TODO
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
		// var inflRegStart = new Array(peaks.length); 
		// var inflRegEnd = new Array(peaks.length);

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

	function get_phase_advances_v2(currentInputPhase, previousInputPhase, omega, RA, RS) {

		var ipa_hop = new Array(omega.length);

		for (var i=0; i<omega.length; i++) {
			var twoPI = 2 * Math.PI;

			var dphi = omega[i] * RA;

			var auxHpi = (currentInputPhase[i] - previousInputPhase[i]) - dphi;
			var hpi = auxHpi - twoPI * Math.round(auxHpi/twoPI);

			var ipa_sample = omega[i] + hpi / RA;

			ipa_hop[i] = ipa_sample * RS;
		}

		return ipa_hop;
	}

	function get_phase_advances_v3(currentInputPhase, previousInputPhase, omega, RA, RS) {
		var twoPI = 2 * Math.PI;

		var ipa_hop = omega.map(function(omegaI, i){
			var dphi = omegaI * RA;

			var auxHpi = (currentInputPhase[i] - previousInputPhase[i]) - dphi;
			var hpi = auxHpi - twoPI * Math.round(auxHpi/twoPI);

			var ipa_sample = omegaI + hpi / RA;

			return ipa_sample * RS;

		});

		return ipa_hop;
	}

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
	function identity_phase_locking(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {
		// var r = findPeaks(currentInputMagnitude);
		// var r = findPeaksV2(currentInputMagnitude,10);
		var r = findPeaksV3(currentInputMagnitude);

		var influenceRegions = r.inflRegionStart.map(function(inflRegStart,i){ return Math.max(0, r.inflRegionEnd[i] - inflRegStart + 1); });

		var phasor_theta = get_phasor_theta_v2(currentInputPhase, previousOutputPhase, instPhaseAdv, r.peaks, influenceRegions);

		return phasor_theta;
	}

	function identity_phase_locking_v2(currentInputMagnitude, currentInputPhase, previousOutputPhase, instPhaseAdv) {
		var r = findPeaks(currentInputMagnitude);

		var influenceRegions = new Array(r.inflRegionStart.length);

		for (var i=0; i<influenceRegions.length; i++) 
			influenceRegions[i] = Math.max(0, r.inflRegionEnd[i] - r.inflRegionStart[i] + 1);

		var phasor_theta = get_phasor_theta_v2(currentInputPhase, previousOutputPhase, instPhaseAdv, r.peaks, influenceRegions);

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

		var instPhaseAdv = get_phase_advances_v4(currentInputPhase, previousInputPhase, omega, RA, RS);

		var currentInputMag = fftObject.spectrum;

		var phasor_theta = identity_phase_locking_v2(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);
		//var phasor_theta = no_phase_locking(currentInputMag, currentInputPhase, previousOutputPhase, instPhaseAdv);

		// return generateY(phasor_theta, fftObject, 1025);

		// Multiplication of two vectors of complex numbers
		// var out_real = new Float32Array((phasor_theta.length-1)*2);
		// var out_imag = new Float32Array((phasor_theta.length-1)*2);
		// var out_phase = new Float32Array((phasor_theta.length-1)*2);
		// var out_magnitude = new Float32Array((phasor_theta.length-1)*2);

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


			// out[2047] = out[1]
			// out[2046] = out[2]
			// ...
			// out[1026] = out[1024]
		}
		
		return {real: out_real, imag: out_imag, phase: out_phase, magnitude: out_magnitude};
	}

	// For some (divine) reason, the forloop inside this function is A LOT slower than the one
	// in pv_step_v2. Go figure...
	function generateY(theta, complexFrame, length) {
		// phasor = e^(j*theta) = cos(theta) + j*sin(theta)
		// X = x_real + j*x_imag
		// Y = phasor .* X =
		//   = [cos(theta) + j*sin(theta)] * [x_real + j*x_imag] =
		//	 = [cos(theta) * x_real - sin(theta) * x_imag] + j*[cos(theta) * x_imag + sin(theta) * x_real]
		//   = Y_real + j*Y_imag

		// var Y_real  = new Float32Array(theta.length);
		// var Y_imag  = new Float32Array(theta.length);
		// var Y_phase = new Float32Array(theta.length);
		// var Y_mag   = new Float32Array(theta.length);
		var Y_real  = new Array(theta.length);
		var Y_imag  = new Array(theta.length);
		var Y_phase = new Array(theta.length);
		var Y_mag   = new Array(theta.length);
		var sqrt	= Math.sqrt;
		var cos		= Math.cos;
		var sin		= Math.sin;
		var atan2	= Math.atan2;

		var dLen = (length-1)*2;

		for (var i=0; i<theta.length; i++) {
			var phasor_theta_real = cos(theta[i]);
			var phasor_theta_imag = sin(theta[i]);
			var X_real  = complexFrame.real[i];
			var X_imag  = complexFrame.imag[i];

			Y_real[i]  = phasor_theta_real * X_real - phasor_theta_imag * X_imag;
			Y_imag[i]  = phasor_theta_real * X_imag + phasor_theta_imag * X_real;
			Y_phase[i] = atan2(Y_real, Y_imag);
			Y_mag[i]   = sqrt(Y_real * Y_real + Y_imag * Y_imag);

			if (i>0) {
				Y_real[dLen-i] = Y_real[i];
				Y_imag[dLen-i] = -Y_imag[i];
				Y_phase[dLen-i] = atan2(Y_imag[dLen-i], Y_real[dLen-i]);
				Y_mag[dLen-i] = sqrt(Y_imag[dLen-i]*Y_imag[dLen-i] + Y_real[dLen-i]*Y_real[dLen-i]);
			}
		}

		return {real: Y_real, imag: Y_imag, magnitude: Y_mag, phase: Y_phase};

	}


	function STFT(frame, windowProcessor, returnSize, fftProcessor) {
		var windowedFrame = windowProcessor.process(frame);
		fftProcessor.forward(windowedFrame);
		//fftProcessor.calculateSpectrumAndPhase(returnSize);
		return {
			real: fftProcessor.real.subarray(0,returnSize),
			imag: fftProcessor.imag.subarray(0,returnSize),
			magnitude: fftProcessor.spectrum.subarray(0,returnSize),
			phase: fftProcessor.angle.subarray(0,returnSize)
		}
	}

	function STFTv2(frame, windowProcessor, returnSize, fftProcessor) {
		var windowedFrame = windowProcessor.process(frame);
		fftProcessor.imag = new Float32Array(frame.length);
		fftProcessor.forwardV2(windowedFrame);
		//fftProcessor.calculateSpectrumAndPhase(returnSize);
		return {
			real: fftProcessor.real.subarray(0,returnSize),
			imag: fftProcessor.imag.subarray(0,returnSize),
			magnitude: fftProcessor.spectrum.subarray(0,returnSize),
			phase: fftProcessor.angle.subarray(0,returnSize)
		}
	}