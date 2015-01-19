function doTests(data, i) {
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

	var diff_real = meuOutput.real.map(function(x,i){
		return testOutput.real[i] - x;
	});

	var diff_imag = meuOutput.imag.map(function(x,i){
		return testOutput.imag[i] - x;
	});

	return {
		real: diff_real,
		imag: diff_imag
	};

}