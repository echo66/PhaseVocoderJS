var BUFFER_SIZE = 2048;

var context = new AudioContext();

var buffer = context.createBuffer(2, BUFFER_SIZE, context.sampleRate);

var node = context.createScriptProcessor(BUFFER_SIZE, 2, 2);

var alpha = 1;

var phasevocoderL = new PhaseVocoder(BUFFER_SIZE, 44100, alpha);
var phasevocoderR = new PhaseVocoder(BUFFER_SIZE, 44100, alpha);

loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            buffer = decodedData
        })
        // createBuffer(request.response);
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

//loadSample('badromance.mp3')
//loadSample('http://localhost/annotator/1.mp3')
loadSample('../soundtouchjs/2.mp3');

var position = 0;

var buffers = [];

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

        var bufL = new Array(BUFFER_SIZE);
        var bufR = new Array(BUFFER_SIZE);

        for (var i = 0; i < BUFFER_SIZE; i++) {
            // bufL = bufL.concat(il[i + position]);
            // bufR = bufR.concat(ir[i + position]);
            bufL[i] = il[i + position];
            bufR[i] = ir[i + position];
        }

        position += Math.round(BUFFER_SIZE/4 / alpha);
        // position += BUFFER_SIZE/4;
        // position += Math.round(BUFFER_SIZE/8);

        // Process left input channel
        outBufferL = outBufferL.concat(phasevocoderL.process(bufL));

        // Process right input channel
        //outBufferR = outBufferR.concat(phasevocoderR.process(bufR));

    } while(outBufferL.length < BUFFER_SIZE);

    for (var i = 0, LEN = outBufferL.length; i < BUFFER_SIZE; i++) {
        // ol[i] = bufL.shift();
        // or[i] = bufR.shift();
        ol[i] = outBufferL.shift();
        or[i] = ol[i];
    }
    
};

function play() {
    node.connect(context.destination);
}

function pause() {
    node.disconnect();
}