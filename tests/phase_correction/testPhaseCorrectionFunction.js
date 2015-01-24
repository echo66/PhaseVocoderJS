function doPhaseCorrectionTest(data, i) {

	function do_diff(mine, test, length, frame) {
		var c = new Array(length);
		for (var i=0; i<length; i++) {
			c[i] = (mine[i] - test[i])*(mine[i] - test[i]);
			c[i] = (c[i]<10e-6 && c[i]>-10e-6)? 0 : c[i]; // to very marginal errors, due to numerical issues.
			if (c[i]>0) {
				console.log('**** Issue at index '+i+' at frame '+frame+' ****');
				console.log('* my value:   '+mine[i]);
				console.log('* test value: '+test[i]);
				if (Math.sign(mine[i])!=Math.sign(test[i])) {
					console.log('* WRONG SIGN *');
				}
				console.log('');
			}
		}
		return Math.sqrt(c.reduce(function(a,b){ return a+b;}));
	}

	function do_diff_v2(mine, test, length, frame) {
		var c = new Array(length);

		var avgMine = mine.reduce(function(a, b) { return a + b })/mine.length;
		var avgTest = test.reduce(function(a, b) { return a + b })/test.length;

		var maxMine = Math.max.apply(Math, mine);
		var minMine = Math.min.apply(Math, mine);

		var maxTest = Math.max.apply(Math, test);
		var minTest = Math.min.apply(Math, test);

		var normMine = mine.map(function(x,i){
			return (x-minMine)/(maxMine-minMine);
		});

		var normTest = mine.map(function(x,i){
			return (x-minTest)/(maxTest-minTest);
		});

		for(var j=0; j<c.length; j++) {
			c[j] = normMine[i]*normMine[i] - normTest[i]*normTest[i];
		}
		
		return Math.sqrt(c.reduce(function(a,b){ return a+b;}));

	}

	var RA = data[i].RA._ArrayData_;
	var RS = data[i].RS._ArrayData_;
	var currentInputPhase = data[i].currentInputPhase._ArrayData_;
	var previousInputPhase = data[i].previousInputPhase._ArrayData_;
	var fftObject = {
		real: data[i].fftObject.real._ArrayData_, 
		imag: data[i].fftObject.imag._ArrayData_, 
		spectrum: data[i].fftObject.magnitude._ArrayData_, 
		angle: data[i].fftObject.phase._ArrayData_
	};
	var previousOutputPhase = data[i].previousOutputPhase._ArrayData_;
	var omega = data[i].omega._ArrayData_ ;

	var meuOutput = pv_step_v2(fftObject, previousInputPhase, previousOutputPhase, omega, RA, RS);

	var testOutput = {
		real: data[i].outputSTFTFrame.real._ArrayData_,
		imag: data[i].outputSTFTFrame.imag._ArrayData_
	};

	var	diff_real = do_diff_v2(meuOutput.real, testOutput.real, 1025, i);

	var	diff_imag = do_diff_v2(meuOutput.imag, testOutput.imag, 1025, i);

	return {
		real: diff_real,
		imag: diff_imag
	};

}

function doPhaseCorrectionTests(data) {
	var totalError = {};
	totalError.real = 0;
	totalError.imag = 0;
	for (var i=0; i<data.length; i++) {
		console.log("STFT Frame "+i+": ");
		var diffs = doPhaseCorrectionTest(data,i);
		console.log(diffs);
		console.log("-------------");
		totalError.real += diffs.real;
		totalError.imag += diffs.imag;
	}
	return totalError;
}