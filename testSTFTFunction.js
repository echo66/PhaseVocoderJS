function doSTFTTest(data,i) {

	function do_diff(a,b,length) {
		var c = new Array(length);
		for (var i=0; i<length; i++) {
			c[i] = (a[i] - b[i])*(a[i] - b[i]);
			c[i] = (c[i]<10e-6 && c[i]>-10e-6)? 0 : c[i]; // to very marginal errors, due to numerical issues.
		}
		return Math.sqrt(c.reduce(function(a,b){ return a+b;}));
	}

	function mag(real,imag) {
		return Math.sqrt(real*real+imag*imag);
	}

	function ang(real,imag) {
		return Math.atan2(imag,real);
	}

	var frameNumber = data[i].frameNumber._ArrayData_;
	var inputFrame = data[i].inputFrame._ArrayData_;
	var STFTInputFrame = {};
	STFTInputFrame.real = data[i].STFTInputFrame.real._ArrayData_;
	STFTInputFrame.imag = data[i].STFTInputFrame.imag._ArrayData_;
	STFTInputFrame.magnitude = data[i].STFTInputFrame.magnitude._ArrayData_;
	STFTInputFrame.phase = data[i].STFTInputFrame.phase._ArrayData_;

	var fftProcessor = new FFT(2048, 44100);
	var winProcessor = new WindowFunction(DSP.SINBETA);
	var meuSTFTInputFrame = STFT(inputFrame, winProcessor, 1025, fftProcessor);

	var diff_real = do_diff(STFTInputFrame.real, meuSTFTInputFrame.real, 1025);

	var diff_imag = do_diff(STFTInputFrame.imag, meuSTFTInputFrame.imag, 1025);

	var diff_magnitude = do_diff(STFTInputFrame.magnitude, meuSTFTInputFrame.magnitude, 1025);

	var diff_phase = do_diff(STFTInputFrame.phase, meuSTFTInputFrame.phase, 1025);

	return {
		real: diff_real,
		imag: diff_imag,
		magnitude: diff_magnitude,
		phase: diff_phase
	};
}

function doSTFTTests(data) {
	var totalError = {};
	totalError.real = 0;
	totalError.imag = 0;
	totalError.magnitude = 0;
	totalError.phase = 0;
	for (var i=0; i<data.length; i++) {
		var diffs = doSTFTTest(data,i);
		console.log("Time Frame "+i+": ");
		console.log(diffs);
		console.log("-------------");
		totalError.real += diffs.real;
		totalError.imag += diffs.imag;
		totalError.magnitude += diffs.magnitude;
		totalError.phase += diffs.phase;
	}
	return totalError;
}