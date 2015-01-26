 /*
  *  opts.showWrongSign
  */
function do_phase_correction_tests(data, opts) {
	for (var i=0; i<data.length; i++) {
		console.log("STFT Frame "+i+": ");
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
		var testOutput = {
			real: data[i].outputSTFTFrame.real._ArrayData_,
			imag: data[i].outputSTFTFrame.imag._ArrayData_
		};
		var diffs = do_phase_correction_test(fftObject, previousInputPhase, previousOutputPhase, omega, RA, RS, testOutput, 1025, i, opts);
		console.log(diffs);
		console.log("-------------");
	}
}

 /*
  *  opts.showWrongSign
  */
function do_phase_correction_test(fftObject, previousInputPhase, previousOutputPhase, omega, RA, RS, testOutput, winSize, i, opts) {
	var myData = pv_step_v2(fftObject, previousInputPhase, previousOutputPhase, omega, RA, RS);

	var	diff_real = do_diff_v2(myData.real, testOutput.real, winSize, i, {showWrongSign: opts.showWrongSign});

	var	diff_imag = do_diff_v2(myData.imag, testOutput.imag, winSize, i, {showWrongSign: opts.showWrongSign});

	return {
		real: diff_real.meanL2norm,
		imag: diff_imag.meanL2norm
	};
}