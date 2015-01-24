function doISTFTTest(data, i) {
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



	var processedSTFTFrame = {};
	processedSTFTFrame.real = data[i].processedSTFTFrame.real._ArrayData_;
	processedSTFTFrame.imag = data[i].processedSTFTFrame.imag._ArrayData_;
	processedSTFTFrame.magnitude = data[i].processedSTFTFrame.magnitude._ArrayData_;
	processedSTFTFrame.phase = data[i].processedSTFTFrame.phase._ArrayData_;

	var inverseSTFTFrame = data[i].inverseSTFTFrame._ArrayData_;
	var inverseSTFTWindowedFrame = data[i].inverseSTFTWindowedFrame._ArrayData_;
	var inverseSTFTWindowedFrameEnergyR = data[i].inverseSTFTWindowedFrameEnergyR._ArrayData_;

	var fftProcessor = new FFT(2048, 44100);
	var winProcessor = new WindowFunction(DSP.SINBETA);

	// DSP.js v1
	// var meuInverseSTFTFrame = fftProcessor.inverse(processedSTFTFrame.real, processedSTFTFrame.imag);
	// var meuInverseSTFTWindowedFrame = winProcessor.process(meuInverseSTFTFrame);

	// fft.js
	// inverseTransform(processedSTFTFrame.real, processedSTFTFrame.imag);
	// var meuInverseSTFTFrame = processedSTFTFrame.real;
	// var meuInverseSTFTWindowedFrame = winProcessor.process(processedSTFTFrame.real);

	// jsfft.js
	var toProcess = new complex_array.ComplexArray(2048);
	toProcess.map(function(freq, i, n) {
		freq.real = processedSTFTFrame.real[i];
		freq.imag = processedSTFTFrame.imag[i];
	});
	var meuInverseSTFTFrame = Array.prototype.slice.call(toProcess.InvFFT().real);
	var meuInverseSTFTWindowedFrame = Array.prototype.slice.call(winProcessor.process(processedSTFTFrame.real));


	var diff_frame = do_diff_v2(meuInverseSTFTFrame, inverseSTFTFrame, 2048, i);

	var diff_windowed_frame = do_diff_v2(meuInverseSTFTWindowedFrame, inverseSTFTWindowedFrame, 2048, i);

	return {
		diff_frame: diff_frame,
		diff_windowed_frame: diff_windowed_frame
	}
}

function doISTFTTests(data) {
	var totalError = {};
	totalError.diff_frame = 0;
	totalError.diff_windowed_frame = 0;
	for (var i=0; i<data.length; i++) {
		console.log("Time Frame "+i+": ");
		var diffs = doISTFTTest(data,i);
		console.log(diffs);
		console.log("-------------");
		totalError.diff_frame += diffs.diff_frame;
		totalError.diff_windowed_frame += diffs.diff_windowed_frame;
	}
	return totalError;
}