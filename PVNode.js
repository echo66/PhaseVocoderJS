function PVNode(context, bufferSize, div) {

	var BUFFER_SIZE = bufferSize;

	var ctx = context;

	var node = context.createScriptProcessor(BUFFER_SIZE, 2, 2);

	var outBufferL = []; 
	var outBufferR = [];
	var position = 0;

	var phasevocoderL = new PhaseVocoder(BUFFER_SIZE/div, 44100); 
	phasevocoderL.init();
	var phasevocoderR = new PhaseVocoder(BUFFER_SIZE/div, 44100); 
	phasevocoderR.init();

	this.setAudioData = function(buffer) {
		this.buffer = buffer;
	}
	
	this.onaudioprocess = function(e) {
		var il = buffer.getChannelData(0); 
		var ir = buffer.getChannelData(1);

	    var ol = e.outputBuffer.getChannelData(0); 
	    var or = e.outputBuffer.getChannelData(1);

	    // Fill output buffers (left & right) until the system has 
	    // enough processed samples to reproduce.
	    do {

	        var bufL = new Float32Array(BUFFER_SIZE);
	        var bufR = new Float32Array(BUFFER_SIZE);
	        bufL = il.subarray(position,position+BUFFER_SIZE);
	        bufR = ir.subarray(position,position+BUFFER_SIZE);

	        position += phasevocoderL.get_analysis_hop();

	        // Process left input channel
	        outBufferL = outBufferL.concat(phasevocoderL.process(bufL));

	        // Process right input channel
	        outBufferR = outBufferR.concat(phasevocoderR.process(bufR));

	    } while(outBufferL.length < BUFFER_SIZE);

	    ol.set(outBufferL.splice(0,BUFFER_SIZE));
	    or.set(outBufferR.splice(0,BUFFER_SIZE));
	}

	this.play = function() {
		node.connect(ctx.destination);
		node.onaudioprocess = this.onaudioprocess;
	}

	this.pause = function() {
		node.disconnect();
	}

	this.setAlpha = function(newAlpha) {
		phasevocoderL.set_alpha(newAlpha);
    	phasevocoderR.set_alpha(newAlpha);
	}

	this.setPosition = function(newPosition) {
		outBufferL = [];
    	outBufferR = [];
    	position = 0;
    	this.reset_all();
	}

	this.reset_phases = function() {
		phasevocoderL.reset_phases();
		phasevocoderR.reset_phases();
	}

	this.reset_all = function() {
		phasevocoderL.reset_phases_and_overlap_buffers();
		phasevocoderR.reset_phases_and_overlap_buffers();
	}
}