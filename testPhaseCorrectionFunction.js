function doPhaseCorrectionTest(data, i) {

	function do_diff(a,b,length) {
		var c = new Array(length);
		for (var i=0; i<length; i++) {
			c[i] = (a[i] - b[i])*(a[i] - b[i]);
			c[i] = (c[i]<10e-6 && c[i]>-10e-6)? 0 : c[i]; // to very marginal errors, due to numerical issues.
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

	var	diff_real = do_diff(meuOutput.real, testOutput.real, 1025);

	var	diff_imag = do_diff(meuOutput.imag, testOutput.imag, 1025);

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
		var diffs = doPhaseCorrectionTest(data,i);
		console.log("STFT Frame "+i+": ");
		console.log(diffs);
		console.log("-------------");
		totalError.real += diffs.real;
		totalError.imag += diffs.imag;
	}
	return totalError;
}