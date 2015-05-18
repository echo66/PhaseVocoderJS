function PhaseVocoder2(_fftFrameSize, _sampleRate) {

	var k = i = peaksLength = phTh_idx = 0;
	var twoPI = 2 * Math.PI;
	var j, mag, phase, tmp, qpd, index, signal;
	var prevPeak, reg, regStart, prevRegEnd, prevRegStart, 
		d, thRe, thIm;
	var expectedPhaseAdv, auxHeterodynedPhaseIncr, 
		heterodynedPhaseIncr, instPhaseAdvPerSampleHop, 
		instPhaseAdv, prevInstPhaseAdv;

	
	var fftSize = _fftFrameSize;
	var sampleRate = _sampleRate;
	var fftSizeHalf = fftSize/2 + 1;

	var win  = new Float32Array(fftSize);
	var win2 = new Float32Array(fftSize);
	for (k = 0; k < fftSize; k++) {
		// win[k]= WindowFunction.Hann(fftSize, k);
		win[k] = Math.pow(Math.sin(Math.PI * k / fftSize), 1);
		// win[k]  = 0.5 * (1 - Math.cos(twoPI * k / (fftSize - 1)));
		win2[k] = win[k] * win[k];
	}
	var omega = new Float32Array(fftSizeHalf);
	for (k = 0; k < fftSizeHalf; k++) {
		omega[k] = twoPI * k / fftSize;
	}

	var RA = RS = 0;

	var gRover = false;
	// This has to go.
	this.MAX_FRAME_LENGTH = 8192;

	var inputFIFO = new Float32Array(this.MAX_FRAME_LENGTH);
	var outputFIFO = new Float32Array(this.MAX_FRAME_LENGTH);
	var prevInPhase = new Float32Array(fftSizeHalf);
	var prevOutPhase = new Float32Array(fftSizeHalf);
	var phasorTheta = new Float32Array(fftSizeHalf);
	var outputAccum = new Array(2 * this.MAX_FRAME_LENGTH);
	var owOutputAccum = new Array(2 * this.MAX_FRAME_LENGTH);
	for (k = 0; k < outputAccum.length; k++) {
		outputAccum[k] = 0;
		owOutputAccum[k] = 0;
	}
	// Not two, 'cos we haven't to fill phases with 0's.
	var frame = new Float32Array(fftSize);

	// Real and imaginary parts of the resynthesized signal
	var RE = new Float32Array(fftSize);
	var IM = new Float32Array(fftSize);

	var osamp = 4;
	var overlapFactor = 4;
	var RA = RS = fftSize / overlapFactor;

	var stepSize = fftSize / this.osamp;
	var freqPerBin = sampleRate / fftSize;
	var expct = twoPI * RA / fftSize;
	var inFifoLatency = fftSize - RA;


	// FFT fields
	var stdlib = { Math: Math, Float32Array: Float32Array, Float64Array: Float64Array };
	var heap = new ArrayBuffer(32*fftSize);
	var fft = fourier.custom["fft_f32_"+fftSize+"_asm"](stdlib, null, heap);
	fft.init();
	var zeros = new Float32Array(fftSize);


	this.process = function (numSampsToProcess, inData, inDataOffset, outData, outDataOffset) {
		
		if (gRover === false) 
	        gRover = inFifoLatency;
	    
		/* main processing loop */
		for (j = 0; j < numSampsToProcess; j++){
			// /* As long as we have not yet collected enough data just read in */
			inputFIFO[gRover] = inData[j + inDataOffset];
			outData[j + outDataOffset] = outputFIFO[gRover - inFifoLatency];
			// console.log(gRover - inFifoLatency)
			gRover++;

			// /* now we have enough data for processing */
			if (gRover >= fftSize) {

				gRover = inFifoLatency;

				/* Do the windowing */
				for (k = 0 ; k < fftSize ; k++) {
					frame[k] = inputFIFO[k] * win[k];
				}

				/* Forward FFT */
				(new Float32Array(heap)).set(frame);
				fourier.custom.array2heap(zeros, new Float32Array(heap), fftSize, fftSize);
				fft.transform();
				fourier.custom.heap2array(new Float32Array(heap), RE, fftSizeHalf, 0);
				fourier.custom.heap2array(new Float32Array(heap), IM, fftSizeHalf, fftSize);


				/****************** ANALYSIS & PROCESSING & SYNTHESIS *******************/

				for (i = 0; i < fftSizeHalf; i++) {

					mag = 2 * Math.sqrt (RE[i] * RE[i] + IM[i] * IM[i]) * 1000;
					currInPhase = Math.atan2 (IM[i], RE[i]);

					expectedPhaseAdv = omega[i] * RA;
					auxHeterodynedPhaseIncr = (currInPhase - prevInPhase[i]) - expectedPhaseAdv;
					heterodynedPhaseIncr = auxHeterodynedPhaseIncr - twoPI * Math.round(auxHeterodynedPhaseIncr/twoPI);
					instPhaseAdvPerSampleHop = omega[i] + heterodynedPhaseIncr / RA;
					instPhaseAdv = instPhaseAdvPerSampleHop * RS;

					if (mag[i] > (mag[i-2]|0) && mag[i] > (mag[i-1]|0) && mag[i] > (mag[i+1]|0) && mag[i] > (mag[i+2]|0)) {
						regStart = Math.ceil((prevPeak + i)/2) | 0; 
						prevRegEnd = regStart-1;
						reg = Math.max(0, prevRegEnd - prevRegStart + 1);
						prevRegStart = regStart;
						for (d = 0; d < reg; d++, phTh_idx++) {
							phasorTheta[phTh_idx] = prevOutPhase[prevPeak] + prevInstPhaseAdv - currInPh[prevPeak];
						}
						prevPeak = i;
						prevInstPhaseAdv = instPhaseAdv;
					}

					prevInPhase[i] = currInPhase;


					thRe = Math.cos(phasorTheta[i]);
					thIm = Math.sin(phasorTheta[i]);

					RE[i] = thRe * RE[i] - thIm * IM[i];
					IM[i] = thRe * IM[i] + thIm * RE[i];
					prevOutPhase[i] = Math.atan2(IM[i], RE[i]);

				}

				/* Inverse FFT */
				(new Float32Array(heap, 0, fftSize)).set(IM);
				fourier.custom.array2heap(RE, new Float32Array(heap), fftSize, fftSize);
				fft.transform();
				frame.set(new Float32Array(heap, fftSize, fftSize));

				// Do inverse windowing and add to output accumulator

				for (k = 0; k < fftSize; k++) {
					outputAccum[k] += win[k] * frame[k] / fftSize;
				}

				for (k = 0; k < RS; k++) {
					outputFIFO[k] = outputAccum[k];
				}
				

				// Shift the output accumulator.
				// Rough memmove implementation.

				// var tempArray = outputAccum.slice (RS, RS + fftSize);
				for (k = 0; k < fftSize; k++) {
					// outputAccum[k] = tempArray[k];
					outputAccum[k] = outputAccum[k + RS + fftSize];
				}

				// Shift the input FIFO
				// These memory shifts have to be optimized.

				for (k = 0; k < inFifoLatency; k++) {
					inputFIFO[k] = inputFIFO[k + RS];
				}
			}
		}
	}

	this.set_alpha = function (alpha) {
		if (alpha <= 0.8)
			overlapFactor = 2;
		else if (alpha <= 1)
			overlapFactor = 4;
		else
			overlapFactor = 5;
		RA = Math.round(fftSize / overlapFactor);
		RS = Math.round(alpha * RA);
	}

	this.set_pitch = function (newPitchFactor) {

	}

  	this.get_alpha = function () {
  		return RS / RA;
  	}

  	this.get_pitch_factor = function() {

  	}

  	this.get_analysis_hop = function() {
  		return RA;
  	}

  	this.get_synthesis_hop = function() {
  		return RS;
  	}

	// Initialize analysis & synthesis hops.
	this.set_alpha(1);
}