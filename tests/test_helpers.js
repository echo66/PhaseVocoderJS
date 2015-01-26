function do_diff_v2(test, mine, length, frame, opt) {
	var amine = mine.slice(0,length);
	var atest = test.slice(0,length);
	var c = new Array(length);
	var wrongSign = [];

	var avgMine = amine.reduce(function(a, b) { return a + b })/amine.length;
	var avgTest = atest.reduce(function(a, b) { return a + b })/atest.length;

	var maxMine = Math.max.apply(Math, amine);
	var minMine = Math.min.apply(Math, amine);

	var maxTest = Math.max.apply(Math, atest);
	var minTest = Math.min.apply(Math, atest);

	var normalizedMine = amine.map(function(x,i){
		return (x-minMine)/(maxMine-minMine);
	});

	var normalizedTest = mine.map(function(x,i){
		return (x-minTest)/(maxTest-minTest);
	});

	for(var j=0; j<c.length; j++) {
		if(Math.sign(amine[j])!=Math.sign(atest[j]) && opt.showWrongSign) {
			wrongSign.push({
				valName: opt.valName, 
				frame: frame, 
				index: j,
				mineVal: amine[j],
				testVal: atest[j]
			});
		}
		var l2norm = Math.sqrt((normalizedMine[j]-normalizedTest[j])*(normalizedMine[j]-normalizedTest[j]));
		//c[j] = (l2norm>1e-10)? l2norm : 0;
		c[j] = l2norm;
	}
	
	return {
		meanL2norm: c.reduce(function(a,b){ return a+b;}) / length,
		L2normVector: c,
		wrongSignAttrs: wrongSign
	}
}

function mag(real,imag) { return Math.sqrt(real*real+imag*imag);}

function ang(real,imag) { return Math.atan2(imag,real); }

function forward_stft(time_frame, sampleRate, winSize, useWindowedFrame, libraryName) {
	if (libraryName=='dsp.js')
		return forward_stft_with_dsp_js(time_frame, sampleRate, winSize, useWindowedFrame);
	else if (libraryName=='jsfft.js')
		return forward_stft_with_jsfft_js(time_frame, sampleRate, winSize, useWindowedFrame);
	else if (libraryName=='complex.js')
		return forward_stft_with_complex_js(time_frame, sampleRate, winSize, useWindowedFrame);
	else if (libraryName=='fft.js')
		return null;

	function forward_stft_with_dsp_js(data, sampleRate, winSize, useWindowedFrame) {
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
			frame: myData.buffer,
			real: myData.real,
			imag: myData.imag,
			magnitude: myData.magnitude,
			angle: myData.angle
		};
	}

	function forward_stft_with_jsfft_js(data, sampleRate, winSize, useWindowedFrame) {
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
			frame: myData1.buffer,
			real: myData2.real,
			imag: myData2.imag,
			magnitude: myData2.magnitude,
			angle: myData2.angle
		};
	}

	function forward_stft_with_complex_js(data, sampleRate, winSize, useWindowedFrame) {
		var myData1 = new Float32Array(winSize);
		var myData2 = new Float32Array(2*winSize);
		if (useWindowedFrame) {
			var winprocessor = new WindowFunction(DSP.SINBETA, 1);
			myData1 = winprocessor.process(data);
		} else {
			for (var p=0; p<winSize; p++) {
				myData1[p] = data[p];
			}
		}
		var fft = new FFT.complex(winSize, false);
		fft.simple(myData2, myData1, 'real');

		var real = new Float32Array(winSize);
		var imag = new Float32Array(winSize);
		var magnitude = new Float32Array(winSize);
		var phase = new Float32Array(winSize);

		for (var p=0; p<winSize; p++) {
			real[p] = myData2[2*p];
			imag[p] = myData2[2*p+1];
			magnitude[p] = Math.sqrt(imag[p]*imag[p] + real[p]*real[p]);
			phase[p] = Math.atan2(imag[p], real[p]);
		}

		return {
			frame: myData1,
			real: real,
			imag: imag,
			magnitude: magnitude,
			angle: phase
		}
	}
}

function inverse_stft(fft_frame_real, fft_frame_imag, libraryName) {
	if (libraryName=='dsp.js')
		return null;
	else if (libraryName=='jsfft.js')
		return null;
	else if (libraryName=='complex.js')
		return null;
	else if (libraryName=='fft.js')
		return null;
}

