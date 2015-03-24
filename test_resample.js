var BUFFER_SIZE = 2048;

var context = new AudioContext();

var buffer = context.createBuffer(2, BUFFER_SIZE, context.sampleRate);

var node = context.createScriptProcessor(BUFFER_SIZE, 2, 2);

var alpha = 1;

var meuOut;

loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            buffer = decodedData
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

loadSample('../soundtouchjs/2.mp3');

var position = 0;

var bufL = [];
var bufR = [];

node.onaudioprocess = function (e) {

    var il = buffer.getChannelData(0);
    var ir = buffer.getChannelData(1);
    var inputC = [il,ir];

    var ol = e.outputBuffer.getChannelData(0);
    var or = e.outputBuffer.getChannelData(1);
    var outputC = [ol, or];

    if (alpha > 1) {
        while (bufL.length<BUFFER_SIZE) {
            bufL = bufL.concat(resample(il.subarray(position, position+BUFFER_SIZE), alpha));
            bufR = bufR.concat(resample(ir.subarray(position, position+BUFFER_SIZE), alpha));
            position += BUFFER_SIZE;
        }
    } else {
        bufL = [];
        bufR = [];
        // process_samples(position, BUFFER_SIZE, inputC, 0, outputC, alpha);
        bufL = bufL.concat(resample(il.subarray(position, position+BUFFER_SIZE), alpha));
        bufR = bufR.concat(resample(ir.subarray(position, position+BUFFER_SIZE), alpha));
        position += BUFFER_SIZE*alpha;
    }

    // while (bufL.length<BUFFER_SIZE) {
    //     bufL = bufL.concat(resample(il.subarray(position, position+BUFFER_SIZE), alpha));
    //     bufR = bufR.concat(resample(ir.subarray(position, position+BUFFER_SIZE), alpha));

    //     position += BUFFER_SIZE;
    // }

    for (var j=0; j<BUFFER_SIZE; j++) {
        ol[j] = bufL.shift();
        or[j] = bufR.shift();
    }

    var p = 0;
    
};

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

function resample(buffer, fromRate /* or speed */, fromFrequency /* or toRate */, toRate, toFrequency) {
    var argc        = arguments.length,
        speed       = (argc === 2 ? fromRate : (argc === 3 ? fromRate / fromFrequency : toRate / fromRate * toFrequency / fromFrequency)),
        l       = buffer.length,
        length      = Math.ceil(l / speed),
        newBuffer   = new Array(length),
        i, n;

    for (i=0, n=0; i<l; i += speed) {
        newBuffer[n++] = linear_interpolation(buffer, i);
    }

    return newBuffer;
};

function nearest_interpolation(arr, pos) {
    return pos >= arr.length - 0.5 ? arr[0] : arr[Math.round(pos)];
};

function linear_interpolation(arr, pos) {
    var first   = Math.floor(pos),
        second  = first + 1,
        frac    = pos - first;
        second  = second < arr.length ? second : 0;

    return arr[first] * (1 - frac) + arr[second] * frac;
};