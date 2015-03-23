importScripts("jensnockert-fft.js/lib/complex.js", "jensnockert-fft.js/lib/real.js", "PV_fast.js");   

var BUFFER_SIZE = 2048;
var div = 1;
var sample_rate = 44100;
var pvL = new PhaseVocoder(BUFFER_SIZE/div, sample_rate); pvL.init();
var pvR = new PhaseVocoder(BUFFER_SIZE/div, sample_rate); pvR.init();
var pv3 = new PhaseVocoder(BUFFER_SIZE/div, sample_rate); pv3.init();
var audioDataL = [];
var audioDataR = [];
var myInterval;
var outBufferL = [];
var outBufferR = [];
var position = 0;
var alphas = [];
var last_alpha = 1;

/*
 * Message format: [type, channel, op-data1, ..., op-dataN]
 */
onmessage = function(e) {
	var type = e.data.type;
	// console.log(e);

	if (type=="set-alpha") {

		alphas.push(e.data.alpha);
		console.log("added new alpha: " + e.data.alpha);

	} else if (type=="play") {

		myInterval = setInterval(function () {

			// console.log("processing");

			do {
				if (alphas.length!=0) {
					last_alpha = alphas.splice(0,1)[0];
					console.log("using new alpha: " + last_alpha);
				}
				pvL.set_alpha(last_alpha);
				pvR.set_alpha(last_alpha);

		        var bufL = new Float32Array(BUFFER_SIZE);
		        var bufR = new Float32Array(BUFFER_SIZE);
		        bufL = audioDataL.subarray(position, position+BUFFER_SIZE);
		        bufR = audioDataR.subarray(position, position+BUFFER_SIZE);

		        position += pvL.get_analysis_hop();

		        // Process left input channel
		        outBufferL = outBufferL.concat(pvL.process(bufL));

		        // Process right input channel
		        outBufferR = outBufferR.concat(pvR.process(bufR));

		    } while(outBufferL.length < BUFFER_SIZE);

		    postMessage({
		    	left: outBufferL.splice(0,BUFFER_SIZE),
		    	right: outBufferR.splice(0,BUFFER_SIZE)
		    })

		    // console.log("done");
		}, 1/(1024/sample_rate));

	} else if (type=="pause") {

		clearInterval(myInterval);

	} else if (type=="position") {

		outBufferL = [];
    	outBufferR = [];
    	position = 0;
    	pvL.reset_phases_and_overlap_buffers();
		pvR.reset_phases_and_overlap_buffers();

	} else {

		audioDataL = new Float32Array(e.data.left);
		audioDataR = new Float32Array(e.data.right);
		console.log("loaded into the worker");

	}
}