/**
 *	opts.libraryName
 *	opts.useWindowedFrame
 *  opts.showWrongSign
 */
function do_simple_time_to_stft_tests(data, sampleRate, winSize, opts) {
	for (var frameNumber=0; frameNumber<data.length; frameNumber++) {
		console.log("Time Frame "+frameNumber+": ");
		var diffs = do_simple_time_to_stft_test(data, frameNumber, sampleRate, winSize, opts);
		console.log(diffs);
		console.log("-------------");
	}
}


/**
 *	opts.libraryName
 *	opts.useWindowedFrame
 *  opts.showWrongSign
 */
function do_simple_time_to_stft_test(data, i, sampleRate, winSize, opts) {

	//var frameNumber = data[i].frameNumber._ArrayData_;
	var inputFrame = data[i].inputFrame._ArrayData_;
	var windowedInputFrame =  data[i].windowedInputFrame._ArrayData_;
	var STFTFrame = {};
	if (opts.useWindowedFrame) {
		STFTFrame.real = data[i].STFTWindowedInputFrame.real._ArrayData_;
		STFTFrame.imag = data[i].STFTWindowedInputFrame.imag._ArrayData_;
		STFTFrame.magnitude = data[i].STFTWindowedInputFrame.magnitude._ArrayData_;
		STFTFrame.phase = data[i].STFTWindowedInputFrame.phase._ArrayData_;
	} else {
		STFTFrame.real = data[i].STFTInputFrame.real._ArrayData_;
		STFTFrame.imag = data[i].STFTInputFrame.imag._ArrayData_;
		STFTFrame.magnitude = data[i].STFTInputFrame.magnitude._ArrayData_;
		STFTFrame.phase = data[i].STFTInputFrame.phase._ArrayData_;
	}
	

	var dataToBeTested = [];

	if (opts.libraryName=='dsp.js') {
		dataToBeTested = process_with_dsp_js(inputFrame, sampleRate, winSize, {useWindowedFrame: opts.useWindowedFrame});
	} else if (opts.libraryName=='jsfft.js') {
		dataToBeTested = process_with_jsfft_js(inputFrame, sampleRate, winSize, {useWindowedFrame: opts.useWindowedFrame});
	} else if (opts.libraryName=='fft.js') {
		//TODO
	}


	var diff_real = do_diff_v2(STFTFrame.real, Array.prototype.slice.call(dataToBeTested.real), winSize, i, {showWrongSign: opts.showWrongSign});

	var diff_imag = do_diff_v2(STFTFrame.imag, Array.prototype.slice.call(dataToBeTested.imag), winSize, i, {showWrongSign: opts.showWrongSign});

	var diff_magnitude = do_diff_v2(STFTFrame.magnitude, Array.prototype.slice.call(dataToBeTested.magnitude), winSize, i, {showWrongSign: opts.showWrongSign});
	
	var diff_phase = do_diff_v2(STFTFrame.phase, Array.prototype.slice.call(dataToBeTested.angle), winSize, i, {showWrongSign: opts.showWrongSign});

	return {
		diff_real : diff_real.meanL2norm,
		diff_imag : diff_imag.meanL2norm,
		diff_magnitude : diff_magnitude.meanL2norm,
		diff_phase : diff_phase.meanL2norm
	};
}

function process_with_dsp_js(data, sampleRate, winSize, useWindowedFrame) {
	var fftprocessor = new FFT(winSize,sampleRate);
	var myData = {};
	myData.buffer = new Float32Array(winSize);
	myData.real = new Float32Array(winSize);
	myData.imag = new Float32Array(winSize);
	myData.magnitude = new Float32Array(winSize);
	myData.angle = new Float32Array(winSize);
	for (var p=0; p<winSize; p++) {
		myData.buffer[p] = data[p];
	}
	if (useWindowedFrame) {
		var winprocessor = new WindowFunction(DSP.SINBETA, 1);
		myData.buffer = winprocessor.process(myData.buffer);
	}
	fftprocessor.forward(myData.buffer);
	fftprocessor.calculateSpectrum();
	fftprocessor.calculateAngle();
	for (var p=0; p<winSize; p++) {
		myData.real[p] = fftprocessor.real[p];
		myData.imag[p] = fftprocessor.imag[p];
		myData.magnitude[p] = fftprocessor.spectrum[p];
		myData.angle[p] = fftprocessor.angle[p];
	}

	return {
		real: myData.real,
		imag: myData.imag,
		magnitude: myData.magnitude,
		angle: myData.angle
	};
}

function process_with_jsfft_js(data, sampleRate, winSize, useWindowedFrame) {
	var myData1 = new complex_array.ComplexArray(winSize);
	myData1.buffer = new Float32Array(winSize);
	for (var p=0; p<winSize; p++) {
		myData1.buffer[p] = data[p];
	}
	if (useWindowedFrame) {
		var winprocessor = new WindowFunction(DSP.SINBETA, 1);
		myData1.buffer = winprocessor.process(myData1.buffer);
	}
	myData1.map(function(value,i,n){
		value.real = myData1.buffer[i];
	}); 
	myData1.FFT();
	var myData2 = {};
	myData2.real = new Float32Array(winSize);
	myData2.imag = new Float32Array(winSize);
	myData2.magnitude = new Float32Array(winSize);
	myData2.angle = new Float32Array(winSize);
	var auxMag = myData1.magnitude();
	var auxAng = myData1.angle();
	for (var p=0; p<winSize; p++) {
		myData2.real[p] = myData1.real[p];
		myData2.imag[p] = myData1.imag[p];
		myData2.magnitude[p] = auxMag[p];
		myData2.angle[p] = auxAng[p];
	}

	return {
		real: myData2.real,
		imag: myData2.imag,
		magnitude: myData2.magnitude,
		angle: myData2.angle
	};
}