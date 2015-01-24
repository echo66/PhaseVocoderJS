function doTimeToSTFTTest(data,i) {

	var frameNumber = data[i].frameNumber._ArrayData_;
	var inputFrame = data[i].inputFrame._ArrayData_;
	var windowedInputFrame = data[i].windowedInputFrame._ArrayData_;
	var STFTInputFrame = {};
	STFTInputFrame.real = data[i].STFTInputFrame.real._ArrayData_;
	STFTInputFrame.imag = data[i].STFTInputFrame.imag._ArrayData_;
	STFTInputFrame.magnitude = data[i].STFTInputFrame.magnitude._ArrayData_;
	STFTInputFrame.phase = data[i].STFTInputFrame.phase._ArrayData_;

	// jsfft.js
	// var meuSTFTInputFrame = new complex_array.ComplexArray(2048);
	// meuSTFTInputFrame.map(function(value,i,n){
	// 	value.real = inputFrame[i];
	// }); 
	// meuSTFTInputFrame.FFT();


	// dsp.js
	var fftprocessor = new FFT(2048,44100);
	var meuSTFTInputFrame = {};
	meuSTFTInputFrame.buffer = new Float32Array(2048);
	meuSTFTInputFrame.real = new Float32Array(2048);
	meuSTFTInputFrame.imag = new Float32Array(2048);
	meuSTFTInputFrame.magnitude = new Float32Array(2048);
	meuSTFTInputFrame.angle = new Float32Array(2048);
	for (var p=0; p<2048; p++) {
		meuSTFTInputFrame.buffer[p] = inputFrame[p];
	}
	fftprocessor.forward(meuSTFTInputFrame.buffer);
	fftprocessor.calculateSpectrum();
	fftprocessor.calculateAngle();
	for (var p=0; p<2048; p++) {
		meuSTFTInputFrame.real[p] = fftprocessor.real[p];
		meuSTFTInputFrame.imag[p] = fftprocessor.imag[p];
		meuSTFTInputFrame.magnitude[p] = fftprocessor.spectrum[p];
		meuSTFTInputFrame.angle[p] = fftprocessor.angle[p];
	}

	console.log('Calculating diff_real');
	var diff_real = do_diff_v2(STFTInputFrame.real, Array.prototype.slice.call(meuSTFTInputFrame.real), 2048, i);
	console.log('Calculating diff_imag');
	var diff_imag = do_diff_v2(STFTInputFrame.imag, Array.prototype.slice.call(meuSTFTInputFrame.imag), 2048, i);
	console.log('Calculating diff_magnitude');
	var diff_magnitude = do_diff_v2(STFTInputFrame.magnitude, Array.prototype.slice.call(meuSTFTInputFrame.magnitude), 2048, i);
	console.log('Calculating diff_phase');
	var diff_phase = do_diff_v2(STFTInputFrame.phase, Array.prototype.slice.call(meuSTFTInputFrame.angle), 2048, i);

	return {
		real : diff_real.meanL2norm,
		imag : diff_imag.meanL2norm,
		magnitude : diff_magnitude.meanL2norm,
		phase : diff_phase.meanL2norm
	};

}

function doTimeToSTFTTests(data) {
	var totalError = {};
	totalError.real = 0;
	totalError.imag = 0;
	totalError.magnitude = 0;
	totalError.phase = 0;
	for (var i=0; i<data.length; i++) {
		console.log("Time Frame "+i+": ");
		var diffs = doTimeToSTFTTest(data,i);
		console.log(diffs);
		console.log("-------------");
		totalError.real += diffs.real;
		totalError.imag += diffs.imag;
		totalError.magnitude += diffs.magnitude;
		totalError.phase += diffs.phase;
	}
	return totalError;
}





function doSTFTToTimeTest(data,i) {

	var frameNumber = data[i].frameNumber._ArrayData_;
	var processedSTFTFrame = {};
	processedSTFTFrame.real = data[i].processedSTFTFrame.real._ArrayData_;
	processedSTFTFrame.imag = data[i].processedSTFTFrame.imag._ArrayData_;
	processedSTFTFrame.magnitude = data[i].processedSTFTFrame.magnitude._ArrayData_;
	processedSTFTFrame.phase = data[i].processedSTFTFrame.phase._ArrayData_;
	var inverseSTFTFrame = data[i].inverseSTFTFrame._ArrayData_;

	// jsfft.js
	// var meuSTFTInputFrame = new complex_array.ComplexArray(2048);
	// meuSTFTInputFrame.map(function(value,p,n){
	// 	value.real = processedSTFTFrame.real[p];
	// 	value.imag = -processedSTFTFrame.imag[p];
	// }); 
	// meuSTFTInputFrame.InvFFT();


	// dsp.js
	var fftprocessor = new FFT(2048,44100);
	for (var p=0; p<2048; p++) {
		fftprocessor.real[p] = processedSTFTFrame.real[p];
		fftprocessor.imag[p] = processedSTFTFrame.imag[p];
	}
	var meuSTFTInputFrame = {};
	meuSTFTInputFrame.real = fftprocessor.inverse();

	// fft.js
	// var meuSTFTInputFrame = {};
	// meuSTFTInputFrame.real = new Float32Array(2048);
	// meuSTFTInputFrame.imag = new Float32Array(2048);
	// for (var p=0; p<2048; p++) {
	// 	meuSTFTInputFrame.real[p] = processedSTFTFrame.real[p];
	// 	meuSTFTInputFrame.imag[p] = processedSTFTFrame.imag[p];
	// }
	// inverseTransform(meuSTFTInputFrame.real, meuSTFTInputFrame.imag);


	var diff_real = do_diff_v2(inverseSTFTFrame, Array.prototype.slice.call(meuSTFTInputFrame.real), 2048, i);

	return {
		real : diff_real
	};

}

function doSTFTToTimeTests(data) {
	var totalError = {};
	totalError.real = 0;
	totalError.imag = 0;
	totalError.magnitude = 0;
	totalError.phase = 0;
	for (var i=0; i<data.length; i++) {
		console.log("STFT Frame "+i+": ");
		var diffs = doSTFTToTimeTest(data,i);
		console.log(diffs);
		console.log("-------------");
		totalError.real += diffs.real;
	}
	return totalError;
}