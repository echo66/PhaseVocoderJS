function do_pipeline_tests(data, winSize, sampleRate) {
	var phase_vocoder = new PhaseVocoder(winSize, sampleRate);
	phase_vocoder.init();

	for (var i=0; i<data.length; i++) {
		console.log("Time Frame "+i+": ");

		var input_time_frame = data[i].input_time_frame;

		var processed_windowed_time_frame = data[i].processed_windowed_time_frame;

		var RA = data[i].RA;

		var RS = data[i].RS;

		var diffs = do_pipeline_test(input_time_frame, processed_windowed_time_frame, RA, RS, phase_vocoder, i);

		console.log(diffs);

		console.log("-------------");
	}
}

function do_pipeline_test(input_time_frame, processed_windowed_time_frame, RA, RS, phasevocoder, frameNumber) {
	phasevocoder.set_hops(RA, RS);
	
	var myOutputFrame = phasevocoder.process_without_overlap_and_add(input_time_frame);

	var diff = do_diff_v2(processed_windowed_time_frame, myOutputFrame, myOutputFrame.length, frameNumber, {});

	return diff.meanL2norm;
}

function build_data(data1, data2) {
	var data = new Array(Math.min(data1.length, data2.length));
	for (var i=0; i<data1.length && i<data2.length; i++) {
		data[i] = {};
		data[i].input_time_frame = data1[i].input_time_frame._ArrayData_;
		data[i].RA = data1[i].RA._ArrayData_[0];
		data[i].RS = data2[i].RS._ArrayData_[0];
		data[i].processed_windowed_time_frame = data2[i].processed_windowed_time_frame._ArrayData_;
	}
	return data;
}