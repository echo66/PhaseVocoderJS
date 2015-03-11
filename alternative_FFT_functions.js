	this.STFT = function(inputFrame, windowFrame, wantedSize, out) {
		var winSize = Math.min(windowFrame.length, wantedSize);
		var data = new complex_array.ComplexArray(winSize);

		data.map(function(value, i, n) {
			value.real = inputFrame[i] * windowFrame[i];
	    });

	    data.FFT();

	    out.real = data.real.subarray(0, winSize);
		out.imag = data.imag.subarray(0, winSize);
		out.magnitude = data.magnitude(winSize);
		out.phase = data.angle(winSize);
	 
	 	return;   
	}

	this.ISTFT = function(real, imag, windowFrame, restoreEnergy, output2) {
		var data = new complex_array.ComplexArray(real.length);

		data.map(function(freq, i, n) {
			freq.real = real[i];
			freq.imag = imag[i];
		});

		output2.set(data.InvFFT().real);

		if (restoreEnergy) {
			// TODO
		} else if (windowFrame) {
			for (var i=0; i<windowFrame.length; i++) 
				output2[i] = (output2[i] / windowFrame.length) * windowFrame[i];
		} else {
			for (var i=0; i<real.length; i++) 
				output2[i] /= real.length;
		}

		return;
	}




	this.STFT = function(real, windowFrame, wantedSize, out) {
		var winSize = Math.min(windowFrame.length, wantedSize);
		var imag = new Float32Array(real.length);

		transform(real, imag);

		out.real = real.subarray(0, winSize);
		out.imag = imag.subarray(0, winSize);
		out.phase = new Float32Array(winSize);
		out.magnitude = new Float32Array(winSize);

		for (var p=0; p<winSize; p++) {
			out.magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
			out.phase[p] = Math.atan2(imag[p], real[p]);
		}

		return;
	}

	this.ISTFT = function(real, imag, windowFrame, restoreEnergy, output2) {
		
		inverseTransform(real, imag);

		for (var i=0; i<windowFrame.length; i++) 
			output2[i] = (real[i] / windowFrame.length) * windowFrame[i];
	}




	this.STFT = function(inputFrame, windowFrame, wantedSize, out) {
		var fft = new FFT(windowFrame.length, 44100);

		var _window = new WindowFunction(DSP.HAMMING, 1);
		_window.process(inputFrame);

		// for (var i=0; i<wantedSize; i++)
		// 	inputFrame[i] *= windowFrame[i];

		fft.forward(inputFrame);
		fft.calculateSpectrum();
		// fft.calculateAngle();

		out.real = new Float32Array(wantedSize);
		out.imag = new Float32Array(wantedSize);
		out.magnitude = new Float32Array(wantedSize);
		out.phase = new Float32Array(wantedSize);

		out.real.set(fft.real.subarray(0,wantedSize));
		out.imag.set(fft.imag.subarray(0,wantedSize));
		out.magnitude.set(fft.spectrum.subarray(0,wantedSize));
		out.phase.set(fft.angle.subarray(0,wantedSize));

		return;
	}

	this.ISTFT = function(real, imag, windowFrame, restoreEnergy, output2) {
		var fft = new FFT(windowFrame.length, 44100);
		
		var buffer = fft.inverse(real,imag);

		var _window = new WindowFunction(DSP.HAMMING, 1);
		_window.process(buffer);

		for (var i=0; i<windowFrame.length; i++) {
			output2[i] = buffer[i];
		}

		// for (var i=0; i<windowFrame.length; i++) {
		// 	output2[i] = buffer[i] * windowFrame[i];
		// }

		return;
	}






	this.STFT = function(inputFrame, windowFrame, wantedSize, out) {
		var fftasm = new FftModule(2048, false);
		var re = new Float32Array(inputFrame.length);
		var im = new Float32Array(inputFrame.length);
		out.magnitude = new Float32Array(wantedSize);
		out.phase = new Float32Array(wantedSize);

		re.set(inputFrame);

		for(var i=0; i<windowFrame.length; i++)
			re[i] *= windowFrame[i];

		fftasm.fft(re, im, false);

		out.real = re.subarray(0,wantedSize);
		out.imag = im.subarray(0,wantedSize);

		for (var i=0; i<wantedSize; i++) {
			out.magnitude[i] = Math.sqrt(re*re + im*im);
			out.phase[i] = Math.atan2(im, re);
		}

		return
	}

	this.ISTFT = function(real, imag, windowFrame, restoreEnergy, output2) {
		var fftasm = new FftModule(2048, false);

		fftasm.fft(imag, real, true);

		for (var i=0; i<windowFrame.length; i++) {
			output2[i] = imag[i] * windowFrame[i];
		}
	}