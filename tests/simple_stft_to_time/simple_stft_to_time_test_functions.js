/**
 *	opts.libraryName
 *	opts.useWindowedFrame
 *  opts.showWrongSign
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
 *	opts.useWindowedFrame
 *  opts.showWrongSign
 */
function do_simple_stft_to_time_test(data, i, sampleRate, winSize, opts) {

	// var PRECISION = 20;
	//var frameNumber = data[i].frameNumber._ArrayData_;
	var processedSTFTFrame = {};
	processedSTFTFrame.real = data[i].processed_stft_frame.real._ArrayData_;
	processedSTFTFrame.imag = data[i].processed_stft_frame.imag._ArrayData_;
	processedSTFTFrame.magnitude = data[i].processed_stft_frame.magnitude._ArrayData_;
	processedSTFTFrame.phase = data[i].processed_stft_frame.phase._ArrayData_;
	var inverseSTFTFrame = data[i].processed_time_frame._ArrayData_;
	if (opts.useWindowedFrame)
		inverseSTFTFrame = data[i].processed_windowed_time_frame._ArrayData_;

	var dataToBeTested = [];

	if (opts.libraryName=='dsp.js') {
		dataToBeTested = process_with_dsp_js(processedSTFTFrame, sampleRate, winSize, opts.useWindowedFrame);
	} else if (opts.libraryName=='jsfft.js') {
		dataToBeTested = process_with_jsfft_js(processedSTFTFrame, sampleRate, winSize, opts.useWindowedFrame);
	} else if (opts.libraryName=='fft.js') {
		dataToBeTested = process_with_fft_js(processedSTFTFrame, sampleRate, winSize, opts.useWindowedFrame);
	} else if (opts.libraryName=='complex.js') {
		dataToBeTested = process_with_complex_js(processedSTFTFrame, sampleRate, winSize, opts.useWindowedFrame);
	}

	var diff_real = do_diff_v2(inverseSTFTFrame, Array.prototype.slice.call(dataToBeTested.real), winSize, i, {showWrongSign: opts.showWrongSign, valName: 'real'});

	return {
		diff_real : diff_real.meanL2norm
	};

}

function process_with_dsp_js(data, sampleRate, winSize, useWindowedFrame) {
	var fftprocessor = new FFT(winSize, sampleRate);

	for (var p=0; p<winSize; p++) {
		fftprocessor.real[p] = data.real[p];
		fftprocessor.imag[p] = data.imag[p];
	}

	var myData = {};

	myData.real = fftprocessor.inverse();

	return myData;
}

function process_with_jsfft_js(data, sampleRate, winSize, useWindowedFrame) {
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

function process_with_fft_js(data, sampleRate, winSize, useWindowedFrame) {
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

function process_with_complex_js(data, sampleRate, winSize, useWindowedFrame) {
	var input = new Array(2 * winSize);
	var output1 = new Array(2 * winSize);
	var output2 = new Array(winSize);
	var windowFrame = createSinBetaWindowArray(winSize, 1);

	for (var i=0; i<winSize; i++) {
		input[2*i] = data.real[i];
		input[2*i+1] = data.imag[i];
	}

	var ifft = new FFT.complex(winSize, true);

	ifft.simple(output1, input);

	for (var i=0; i<winSize; i++) {
		output2[i] = output1[2*i] / winSize;
		if (useWindowedFrame)
			output2[i] *= windowFrame[i];
	}

	return {
		real: output2
	}
}