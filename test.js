var BUFFER_SIZE = 2048;

var context = new AudioContext();

var buffer = context.createBuffer(2, BUFFER_SIZE, context.sampleRate);

var node = context.createScriptProcessor(BUFFER_SIZE, 2, 2);

var alpha = 1;

var phasevocoderL = new PhaseVocoder(BUFFER_SIZE, 44100); phasevocoderL.init();
var phasevocoderR = new PhaseVocoder(BUFFER_SIZE, 44100); phasevocoderR.init();

loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            buffer = decodedData;
        });
    }

    console.log('reading url');
    request.send();
}

// loadSample('../soundtouchjs/4.mp3');

var position = 0;

var outBufferL = [];
var outBufferR = [];

// node.onaudioprocess = function (e) {
//     var il = buffer.getChannelData(0);
//     var ir = buffer.getChannelData(1);

//     var ol = e.outputBuffer.getChannelData(0);
//     var or = e.outputBuffer.getChannelData(1);

//     for (var i=0; i<BUFFER_SIZE; i++) {
//         ol[i] = il[position+i];
//         or[i] = ir[position+i];
//     }

//     position += BUFFER_SIZE;
// };

node.onaudioprocess = function (e) {

    var il = buffer.getChannelData(0);
    var ir = buffer.getChannelData(1);

    var ol = e.outputBuffer.getChannelData(0);
    var or = e.outputBuffer.getChannelData(1);

    // Fill output buffers (left & right) until the system has 
    // enough processed samples to reproduce.
    do {

        var bufL = new Float32Array(BUFFER_SIZE);
        var bufR = new Float32Array(BUFFER_SIZE);

        // for (var i = 0; i < BUFFER_SIZE; i++) {
        //     bufL[i] = il[i + position];
        //     bufR[i] = ir[i + position];
        // }
        bufL = il.subarray(position,position+BUFFER_SIZE);
        bufR = ir.subarray(position,position+BUFFER_SIZE);

        position += phasevocoderL.get_analysis_hop();

        // Process left input channel
        outBufferL = outBufferL.concat(phasevocoderL.process(bufL));

        // Process right input channel
        // outBufferR = outBufferR.concat(phasevocoderR.process(bufR));

    } while(outBufferL.length < BUFFER_SIZE);


    for (var i = 0, LEN = outBufferL.length; i < BUFFER_SIZE; i++) {
        ol[i] = outBufferL.shift();
        // or[i] = outBufferR.shift();
    }

    // ol = outBufferL.splice(0,BUFFER_SIZE);
    // or = outBufferR.splice(0,BUFFER_SIZE);
    
};

function setAlpha(newAlpha) {
    phasevocoderL.set_alpha(newAlpha);
    phasevocoderR.set_alpha(newAlpha);
}

function setPosition(v) {
    resetPVs2();
    outBufferL = [];
    outBufferR = [];
    position = Math.round(buffer.length * v);
}

function resetPVs() {
    phasevocoderL.reset();
    phasevocoderR.reset();
}

function resetPVs2() {
    phasevocoderL.reset2();
    phasevocoderR.reset2();
}

function play() {
    node.connect(context.destination);
}

function pause() {
    node.disconnect();
}

document.addEventListener('DOMContentLoaded', function () {
    var toggleActive = function (e, toggle) {
        e.stopPropagation();
        e.preventDefault();
        // toggle ? e.target.classList.add('wavesurfer-dragover') :
        //     e.target.classList.remove('wavesurfer-dragover');
    };

    var handlers = {
        // Drop event
        drop: function (e) {
            toggleActive(e, false);

            // Load the file into wavesurfer
            if (e.dataTransfer.files.length) {
                pause();
                position = 0;
                resetPVs();

                var my = this;
                // Create file reader
                var reader = new FileReader();
                reader.addEventListener('progress', function (e) {
                    console.log(e);
                });
                reader.addEventListener('load', function (e) {
                    document.getElementById('filename').innerHTML = "<b>" + filename + "</b> loaded";
                    context.decodeAudioData(e.target.result, function(decodedData) {
                        buffer = decodedData;
                    });
                });
                reader.addEventListener('error', function () {
                    console.error('Error reading file');
                });

                var filename = e.dataTransfer.files[0].name;
                reader.readAsArrayBuffer(e.dataTransfer.files[0].slice());

            } else {
                console.error('Not a file');
            }
        },

        // Drag-over event
        dragover: function (e) {
            toggleActive(e, true);
        },

        // Drag-leave event
        dragleave: function (e) {
            toggleActive(e, false);
        }
    };

    var dropTarget = document.querySelector('#drop');
    Object.keys(handlers).forEach(function (event) {
        dropTarget.addEventListener(event, handlers[event]);
    });
});