	this.STFT = function(inputFrame, windowFrame, wantedSize, out) {
		var winSize = windowFrame.length;
		var _inputFrame = new Array(winSize);
		var fftFrame = new Array(2*winSize);

		for (var i=0; i<winSize; i++) {
			_inputFrame[i] = inputFrame[i] * windowFrame[i];
		}
		var fft = new FFT.complex(winSize, false);
		fft.simple(fftFrame, _inputFrame, 'real');

		out.real = new Array(Math.min(winSize,wantedSize));
		out.imag = new Array(Math.min(winSize,wantedSize));
		out.magnitude = new Array(Math.min(winSize,wantedSize));
		out.phase = new Array(Math.min(winSize,wantedSize));

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