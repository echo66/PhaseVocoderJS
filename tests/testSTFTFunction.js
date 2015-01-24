'use strict'

function doSTFTTest(data,i) {

	function do_diff(mine, test, length, frame) {
		var c = new Array(length);
		for (var i=0; i<length; i++) {
			c[i] = (mine[i] - test[i])*(mine[i] - test[i]);
			c[i] = (c[i]<10e-6 && c[i]>-10e-6)? 0 : c[i]; // to very marginal errors, due to numerical issues.
			// if (c[i]>0) {
			// 	console.log('**** Issue at index '+i+' at frame '+frame+' ****');
			// 	console.log('* my value:   '+mine[i]);
			// 	console.log('* test value: '+test[i]);
			// 	if (Math.sign(mine[i])!=Math.sign(test[i])) {
			// 		console.log('* WRONG SIGN *');
			// 	}
			// 	console.log('');
			// }
		}
		return Math.sqrt(c.reduce(function(a,b){ return a+b;}));
	}

	function do_diff_v2(mine, test, length, frame) {

		var amine = mine.slice(0,length);
		var atest = test.slice(0,length);
		var c = new Array(length);

		var avgMine = amine.reduce(function(a, b) { return a + b })/amine.length;
		var avgTest = atest.reduce(function(a, b) { return a + b })/atest.length;

		var maxMine = Math.max.apply(Math, amine);
		var minMine = Math.min.apply(Math, amine);

		var maxTest = Math.max.apply(Math, atest);
		var minTest = Math.min.apply(Math, atest);

		var normMine = amine.map(function(x,i){
			return (x-minMine)/(maxMine-minMine);
		});

		var normTest = mine.map(function(x,i){
			return (x-minTest)/(maxTest-minTest);
		});

		for(var j=0; j<c.length; j++) {
			c[j] = (normMine[i]-normTest[i])*(normMine[i]-normTest[i]);
		}
		
		return Math.sqrt(c.reduce(function(a,b){ return a+b;}));

	}

	function mag(real,imag) { return Math.sqrt(real*real+imag*imag);}

	function ang(real,imag) { return Math.atan2(imag,real); }

	var frameNumber = data[i].frameNumber._ArrayData_;
	var inputFrame = data[i].inputFrame._ArrayData_;
	var STFTInputFrame = {};
	STFTInputFrame.real = data[i].STFTInputFrame.real._ArrayData_;
	STFTInputFrame.imag = data[i].STFTInputFrame.imag._ArrayData_;
	STFTInputFrame.magnitude = data[i].STFTInputFrame.magnitude._ArrayData_;
	STFTInputFrame.phase = data[i].STFTInputFrame.phase._ArrayData_;

	var fftProcessor = new FFT(2048, 44100);
	var winProcessor = new WindowFunction(DSP.SINBETA);

	// DSP.js v1
	// var meuSTFTInputFrame = STFT(inputFrame, winProcessor, 1025, fftProcessor);

	// DSP.js v2
	// var meuSTFTInputFrame = STFTv2(inputFrame, winProcessor, 1025, fftProcessor);

	// fft.js
	var meuSTFTInputFrame = {};
	meuSTFTInputFrame.real = winProcessor.process(inputFrame);
	// If I use Float32Array instead of Array, I get numerical errors
	// meuSTFTInputFrame.imag = new Float32Array(2048); 
	// If I initialize 
	// meuSTFTInputFrame.imag = new Array(2048);
	// meuSTFTInputFrame.imag = createConstantArray(2048,0);
	meuSTFTInputFrame.imag = [];
	meuSTFTInputFrame.magnitude = new Float32Array(2048);
	meuSTFTInputFrame.phase = new Float32Array(2048); 
	transform(meuSTFTInputFrame.real, meuSTFTInputFrame.imag);
	for (var j=0; j<2048; j++) {
		meuSTFTInputFrame.magnitude[j] = mag(meuSTFTInputFrame.real[j], meuSTFTInputFrame.imag[j]);
		meuSTFTInputFrame.phase[j] = ang(meuSTFTInputFrame.real[j], meuSTFTInputFrame.imag[j]);
	}



	var diff_real = do_diff_v2(STFTInputFrame.real, meuSTFTInputFrame.real, 1025, i);

	var diff_imag = do_diff_v2(STFTInputFrame.imag, meuSTFTInputFrame.imag, 1025, i);

	var diff_magnitude = do_diff_v2(STFTInputFrame.magnitude, Array.prototype.slice.call(meuSTFTInputFrame.magnitude), 1025, i);

	var diff_phase = do_diff_v2(STFTInputFrame.phase, Array.prototype.slice.call(meuSTFTInputFrame.phase), 1025, i);

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
		console.log("Time Frame "+i+": ");
		var diffs = doSTFTTest(data,i);
		console.log(diffs);
		console.log("-------------");
		totalError.real += diffs.real;
		totalError.imag += diffs.imag;
		totalError.magnitude += diffs.magnitude;
		totalError.phase += diffs.phase;
	}
	return totalError;
}