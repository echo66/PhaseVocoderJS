var BUFFER_SIZE = 2048;

var winLenHalf = Math.round(BUFFER_SIZE/2);

var context = new AudioContext();

var buffer = context.createBuffer(2, BUFFER_SIZE, context.sampleRate);

var node = context.createScriptProcessor(BUFFER_SIZE, 2, 2);

var alpha = 1;

var phasevocoderL = new PhaseVocoder(BUFFER_SIZE/2, 44100); phasevocoderL.setAlpha(1);
var phasevocoderR = new PhaseVocoder(BUFFER_SIZE/2, 44100); phasevocoderR.setAlpha(1);

loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            // buffer = decodedData;

            buffer = context.createBuffer(2, decodedData.length+BUFFER_SIZE*3+winLenHalf, context.sampleRate);
            var bufL = buffer.getChannelData(0);
            var bufR = buffer.getChannelData(1);
            bufL.set(decodedData.getChannelData(0),winLenHalf);
            bufR.set(decodedData.getChannelData(1),winLenHalf);

            buffer = context.createBuffer(2, decodedData.length+winLenHalf, context.sampleRate);
            var bufL = buffer.getChannelData(0);
            var bufR = buffer.getChannelData(1);
            bufL.set(decodedData.getChannelData(0),winLenHalf);
            bufR.set(decodedData.getChannelData(1),winLenHalf);
        });
    }

    console.log('reading url');
    request.send();
}

function createBuffer(arrayBuffer) {
    offset = 0;
    startTime = 0;
    var start = new Date();
    // NOTE the second parameter is required, or a TypeError is thrown
    buffer = context.createBuffer(2, arrayBuffer.byteLength, context.sampleRate);

    console.log('loaded audio in ' + (new Date() - start));
}

loadSample('../soundtouchjs/4.mp3');

var position = 0;

var buffers = [];

// node.onaudioprocess = function(e) {
//     var il = buffer.getChannelData(0);
//     var ir = buffer.getChannelData(1);
//     var i = new Array(2);
//     i[0] = il;
//     i[1] = ir;

//     var ol = e.outputBuffer.getChannelData(0);
//     var or = e.outputBuffer.getChannelData(1);
//     var o = new Array(2);
//     o[0] = ol;
//     o[1] = or;

//     process_samples(position, BUFFER_SIZE, i, 0, o, alpha);

//     position += BUFFER_SIZE*alpha;
// }

node.onaudioprocess = function (e) {

    // var outBufferL = new Float32Array(BUFFER_SIZE);
    // var outBufferR = new Float32Array(BUFFER_SIZE);

    var outBufferL = [];
    var outBufferR = [];

    var il = buffer.getChannelData(0);
    var ir = buffer.getChannelData(1);

    var ol = e.outputBuffer.getChannelData(0);
    var or = e.outputBuffer.getChannelData(1);

    // Fill output buffers (left & right) until the system has 
    // enough processed samples to reproduce.
    do {

        var bufL = new Float32Array(BUFFER_SIZE);
        var bufR = new Float32Array(BUFFER_SIZE);

        for (var i = 0; i < BUFFER_SIZE; i++) {
            // bufL = bufL.concat(il[i + position]);
            // bufR = bufR.concat(ir[i + position]);
            bufL[i] = il[i + position];
            bufR[i] = ir[i + position];
        }

        position += phasevocoderL.getAnalysisHop();
        // position += BUFFER_SIZE/4;
        // position += Math.round(BUFFER_SIZE/8);

        // Process left input channel
        outBufferL = outBufferL.concat(phasevocoderL.process(bufL));

        // Process right input channel
        outBufferR = outBufferR.concat(phasevocoderR.process(bufR));

    } while(outBufferL.length < BUFFER_SIZE);

    for (var i = 0, LEN = outBufferL.length; i < BUFFER_SIZE; i++) {
        // ol[i] = bufL.shift();
        // or[i] = bufR.shift();
        ol[i] = outBufferL.shift();
        or[i] = outBufferR.shift()
    }
    
};

function setAlpha(newAlpha) {
    phasevocoderL.setAlpha(newAlpha);
    phasevocoderR.setAlpha(newAlpha);
}

function reset() {
    phasevocoderL.reset2();
    phasevocoderR.reset2();
}

function play() {
    node.connect(context.destination);
}

function pause() {
    node.disconnect();
}

function process_samples(input_start, buffer_size, input_channels, output_start, output_channels, rate) {
    var beat, destination_offset, sample_l, sample_r, source_offset, source_offset_float;
    while (--buffer_size >= 0) {
        source_offset_float = input_start + (buffer_size * rate);
        source_offset = Math.round(source_offset_float);
        destination_offset = output_start + buffer_size;
        sample_l = input_channels[0][source_offset];
        sample_r = input_channels[1][source_offset];
        output_channels[0][destination_offset] = sample_l;
        output_channels[1][destination_offset] = sample_r;
    }
    return null;
};