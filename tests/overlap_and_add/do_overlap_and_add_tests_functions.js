function do_overlap_and_add_tests(data, sampleRate, winSize, opts) {
	// var x_buffer = createConstantArray(winSize, 0);
	// var ow_buffer = createConstantArray(winSize, 0);
	var x_buffer = data[0].overlapped_x_after._ArrayData_;
	var ow_buffer = data[0].overlapped_window_after._ArrayData_;
	for (var i=1; i<data.length; i++) {
		console.log("Time Frame "+i+": ");
		var ow_frame = data[i].overlapped_window_before._ArrayData_;
		var x_frame  = data[i].processed_windowed_time_frame_energy_restored._ArrayData_;

		inputData = {};
		inputData.x_buffer = x_buffer;
		inputData.ow_buffer = ow_buffer;
		inputData.x_frame = x_frame;
		inputData.ow_frame = createSinBetaWindowArray(winSize, 1).map(function(x,i){ return x*x; });

		testData = {};
		testData.overlapped_x_after  = data[i].overlapped_x_after._ArrayData_;
		testData.overlapped_ow_after = data[i].overlapped_window_after._ArrayData_;

		var diffs = do_overlap_and_add_test(inputData, testData, data[i].RS._ArrayData_[0], winSize, i);

		console.log(diffs);
	}
}

function do_overlap_and_add_test(inputData, testData, RS, winSize, frame) {
	var x_output = overlapAndSlide(RS, inputData.x_frame, inputData.x_buffer, winSize);

	var ow_output = overlapAndSlide(RS, inputData.ow_frame, inputData.ow_buffer, winSize);

	var diff_x  = do_diff_v2(testData.overlapped_x_after, inputData.x_buffer, winSize, frame, {});

	var diff_ow = do_diff_v2(testData.overlapped_ow_after, inputData.ow_buffer, winSize, frame, {});

	return {
		diff_x: diff_x.meanL2norm,
		diff_ow: diff_ow.meanL2norm
	}
}


function overlapAndSlide(RS, frame, overlapBuffer, windowSize) {
	//var finishedBytes = new Float32Array(RS);
	var finishedBytes = new Array(RS);

	for (var i=0; i<RS; i++) {
		finishedBytes[i] = overlapBuffer.shift();
		while(overlapBuffer.length > windowSize - 1 && overlapBuffer.length >= 0)
			overlapBuffer.shift();
		overlapBuffer.push.apply(overlapBuffer, [0.0]);
	}

	var outBytes = [].concat(overlapBuffer);

	for (var i=0; i<outBytes.length; i++) 
		outBytes[i] = frame[i] + overlapBuffer[i];

	while ((overlapBuffer.length > windowSize - outBytes.length) && overlapBuffer.length >= 0) {
		overlapBuffer.shift();
    }

    overlapBuffer.push.apply(overlapBuffer, outBytes);
	
	return finishedBytes;
}