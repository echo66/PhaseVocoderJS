var audioContext = new AudioContext();

var app = {};

var BUFFER_SIZE = 2048;

function fetch(url) {
    return new Promise(function(resolve, reject) {
        var xhr = new XMLHttpRequest();
        if (url.indexOf(".mp3") !== -1) {
            xhr.responseType = "arraybuffer";
        }
        xhr.open("GET", url);
        xhr.onload = function() {
            resolve({
                text: function() {
                    return xhr.response;
                },
                arrayBuffer: function() {
                    return xhr.response;
                }
            });
        };
        xhr.onerror = reject;
        xhr.send();
    });
}

// audio from YouTube Audio Library
fetch("/soundtouchjs/2.mp3").then(function(res) {
    return audioContext.decodeAudioData(res.arrayBuffer());
}).then(function(audioBuffer) {
    app.audioBuffer = audioBuffer;
});



var bufferSource, phaseVocoderNode, intermediaryNode, finalNode;
var position = 0;

function start() {
    // bufferSource = audioContext.createBufferSource();

    // bufferSource.buffer = app.audioBuffer;
    // bufferSource.start(audioContext.currentTime);
    // bufferSource.onended = function() {
    //   app.stop();
    // };

    intermediaryNode = audioContext.createScriptProcessor(BUFFER_SIZE, 2, 2);

    finalNode = audioContext.createScriptProcessor(BUFFER_SIZE, 2, 2);

    phaseVocoderNode = audioContext.createAudioWorker("pv_worker3.js", 2, 2);
    
    phaseVocoderNode.addParameter("alpha", 1);


    intermediaryNode.onaudioprocess = function (e) {
        var il = app.audioBuffer.getChannelData(0);
        var ir = app.audioBuffer.getChannelData(1);

        var ol = e.outputBuffer.getChannelData(0);
        var or = e.outputBuffer.getChannelData(1);

        for (var i = 0; i < BUFFER_SIZE; i++) {
            ol[i] = il[position+i];
            or[i] = ir[position+i];
        }

        position += BUFFER_SIZE/4;

        // console.log("next");
    };

    finalNode.onaudioprocess = function (e) {
        var il = e.inputBuffer.getChannelData(0);
        var ir = e.inputBuffer.getChannelData(1);

        var ol = e.outputBuffer.getChannelData(0);
        var or = e.outputBuffer.getChannelData(1);

        for (var i = 0; i < BUFFER_SIZE; i++) {
            ol[i] = il[i];
            or[i] = ir[i];
        }

    };

    phaseVocoderNode.alpha.value = 1;

    intermediaryNode.connect(phaseVocoderNode);

    phaseVocoderNode.connect(finalNode);

    finalNode.connect(audioContext.destination);
}

function stop() {
    bufferSource.stop(audioContext.currentTime);
    bufferSource.disconnect();
    phaseVocoderNode.disconnect();
}

function updateParameters(alpha) {
    phaseVocoderNode.alpha.value = alpha;
}