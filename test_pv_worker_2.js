var BUFFER_SIZE = 2048;

var context = new AudioContext();

var buffer = context.createBuffer(2, BUFFER_SIZE, context.sampleRate);

var node = context.createScriptProcessor(BUFFER_SIZE, 2, 2);

var myWorker = new Worker('pv_worker_2.js');

var alpha = 1;

myWorker.onmessage = function(e) {
    bufL = bufL.concat(e.data.left);
    bufR = bufR.concat(e.data.right);
    node.connect(context.destination);
}

var bufL = [];
var bufR = [];


loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            var o = {
                left: decodedData.getChannelData(0),
                right: decodedData.getChannelData(1)
            }
            myWorker.postMessage(o, [o.left.buffer, o.right.buffer]);
        });
    }

    console.log('reading url');
    request.send();
}

loadSample('../soundtouchjs/2.mp3');

node.onaudioprocess = function (e) {

    // console.log("MAIN THREAD PROCESSING");

    if (bufL.length<BUFFER_SIZE) {
        // pause();
        myWorker.postMessage({
            type: "process", 
            alpha: alpha
        });
        return;
    }

    var ol = e.outputBuffer.getChannelData(0);
    var or = e.outputBuffer.getChannelData(1);

    ol.set(bufL.splice(0, BUFFER_SIZE));
    or.set(bufR.splice(0, BUFFER_SIZE));
    
};

function play() {
    node.connect(context.destination);
}

function pause() {
    node.disconnect();
}

function setAlpha(newAlpha) {
    alpha = newAlpha;
}