/**
 *	opts.libraryName
 */
function do_simple_stft_to_time_tests(data, sampleRate, winSize, opts) {
	for (var frameNumber=0; frameNumber<data.length; frameNumber++) {
		console.log("STFT Frame "+frameNumber+": ");
		var diffs = do_simple_stft_to_time_test(data, frameNumber, sampleRate, winSize, opts);
		console.log(diffs);
		console.log("-------------");
	}
}

/**
 *	opts.libraryName
 */
function do_simple_stft_to_time_test(data, i, sampleRate, winSize, opts) {

	//var frameNumber = data[i].frameNumber._ArrayData_;
	var processedSTFTFrame = {};
	processedSTFTFrame.real = data[i].processedSTFTFrame.real._ArrayData_;
	processedSTFTFrame.imag = data[i].processedSTFTFrame.imag._ArrayData_;
	processedSTFTFrame.magnitude = data[i].processedSTFTFrame.magnitude._ArrayData_;
	processedSTFTFrame.phase = data[i].processedSTFTFrame.phase._ArrayData_;
	var inverseSTFTFrame = data[i].inverseSTFTFrame._ArrayData_;

	var dataToBeTested = [];

	if (opts.libraryName=='dsp.js') {
		dataToBeTested = process_with_dsp_js(processedSTFTFrame, sampleRate, winSize);
	} else if (opts.libraryName=='jsfft.js') {
		dataToBeTested = process_with_jsfft_js(processedSTFTFrame, sampleRate, winSize);
	} else if (opts.libraryName=='fft.js') {
		dataToBeTested = process_with_fft_js(processedSTFTFrame, sampleRate, winSize);
	}

	var diff_real = do_diff_v2(inverseSTFTFrame, Array.prototype.slice.call(dataToBeTested.real), winSize, i, {showWrongSign: opts.showWrongSign, valName: 'real'});

	return {
		diff_real : diff_real.meanL2norm
	};

}

function process_with_dsp_js(data, sampleRate, winSize) {
	var fftprocessor = new FFT(winSize, sampleRate);

	for (var p=0; p<winSize; p++) {
		fftprocessor.real[p] = data.real[p];
		fftprocessor.imag[p] = data.imag[p];
	}

	var myData = {};

	myData.real = fftprocessor.inverse();

	return myData;
}

function process_with_jsfft_js(data, sampleRate, winSize) {
	var myData = new complex_array.ComplexArray(winSize);
	myData.map(function(value,p,n){
		value.real = data.real[p];
		value.imag = data.imag[p];
	}); 
	myData.InvFFT();

	return {
		real: myData.real
	}
}

function process_with_fft_js(data, sampleRate, winSize) {
	var myData = {};
	meuSTFTInputFrame.real = new Float32Array(winSize);
	meuSTFTInputFrame.imag = new Float32Array(winSize);
	for (var p=0; p<winSize; p++) {
		meuSTFTInputFrame.real[p] = data.real[p];
		meuSTFTInputFrame.imag[p] = data.imag[p];
	}
	inverseTransform(myData.real, myData.imag);

	return {
		real: myData.real
	}
}