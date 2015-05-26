var BUFFER_SIZE = 2048;

var div = 1;

var context = new AudioContext();

var masterNode = context.createScriptProcessor(BUFFER_SIZE*4, 2, 2);

var nodes = [];

var eqs = [];

loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            console.log('decoded data');
            var I = nodes.length;
            var N = nodes[I] = context.createScriptProcessor(BUFFER_SIZE, 2, 2);
            masterNode.connect(nodes[I]);
            eqs[I] = new Equalizer(nodes[I], context);
            N.buffer = decodedData;
            N.pvL = new PhaseVocoder(BUFFER_SIZE/div, 44100); N.pvL.init();
            N.pvR = new PhaseVocoder(BUFFER_SIZE/div, 44100); N.pvR.init();
            N.outBufferL = [];
            N.outBufferR = [];
            N.position = 0;
            N.pitch = 1;
            N.onaudioprocess = function (e) {

                var il = this.buffer.getChannelData(0);
                var ir = this.buffer.getChannelData(1);

                var ol = e.outputBuffer.getChannelData(0);
                var or = e.outputBuffer.getChannelData(1);

                // Fill output buffers (left & right) until the system has 
                // enough processed samples to reproduce.
                do {

                    // var bufL = new Float64Array(BUFFER_SIZE/div);
                    // var bufR = new Float64Array(BUFFER_SIZE/div);
                    var bufL = il.subarray(this.position, this.position+BUFFER_SIZE/div);
                    var bufR = ir.subarray(this.position, this.position+BUFFER_SIZE/div);

                    this.position += this.pvL.get_analysis_hop();

                    // Process left input channel
                    this.outBufferL = this.outBufferL.concat(this.pvL.process(bufL));

                    // Process right input channel
                    this.outBufferR = this.outBufferR.concat(this.pvR.process(bufR));

                } while(this.outBufferL.length < BUFFER_SIZE);

                ol.set(this.outBufferL.splice(0, BUFFER_SIZE));
                or.set(this.outBufferR.splice(0, BUFFER_SIZE));
            };
        });
    }

    console.log('reading url');
    request.send();
}

// loadSample('../soundtouchjs/4.mp3');
// loadSample('../soundtouchjs/2.mp3');
// loadSample('../soundtouchjs/3.mp3');


function set_pitch(newPitch) {
    pitch = phasevocoderL1.get_synthesis_hop()*newPitch / phasevocoderL1.get_analysis_hop();
    phasevocoderL1.set_overlap_factor(pitch);
    phasevocoderR1.set_overlap_factor(pitch);
}

function set_alpha(ids, newFactor) {
    for (var i=0; i<ids.length; i++) {
        nodes[i].pvL.set_alpha(newFactor);
        nodes[i].pvR.set_alpha(newFactor);
    }
}

function set_position(ids, newPosition) {
    for (var i=0; i<ids.length; i++) {
        nodes[i].position = newPosition;
    }
}

function play(ids) {
    for (var i=0; i<ids.length; i++) 
        eqs[ids[i]].connect();
}

function pause(ids) {
    for (var i=0; i<ids.length; i++) 
        eqs[ids[i]].disconnect();
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







function linearInterpolation (a, b, t) {
    return a + (b - a) * t;
};
function hannWindow (length) {

    var window = new Float32Array(length);
    for (var i = 0; i < length; i++) {
        window[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (length - 1)));
    }
    return window;
};









grainSize = 512;
overlapRatio = 0.70;
pitchShifterProcessor = context.createScriptProcessor(grainSize, 1, 1);
pitchShifterProcessor.buffer = new Float32Array(grainSize * 2);
pitchShifterProcessor.grainWindow = hannWindow(grainSize);
pitchRatio = 1;

pitchShifterProcessor.onaudioprocess = function (event) {

    var inputData = event.inputBuffer.getChannelData(0);
    var outputData = event.outputBuffer.getChannelData(0);

    for (i = 0; i < inputData.length; i++) {

        // Apply the window to the input buffer
        inputData[i] *= this.grainWindow[i];

        // Shift half of the buffer
        this.buffer[i] = this.buffer[i + grainSize];

        // Empty the buffer tail
        this.buffer[i + grainSize] = 0.0;
    }

    // Calculate the pitch shifted grain re-sampling and looping the input
    var grainData = new Float32Array(grainSize * 2);
    for (var i = 0, j = 0.0;
         i < grainSize;
         i++, j += pitchRatio) {

        var index = Math.floor(j) % grainSize;
        var a = inputData[index];
        var b = inputData[(index + 1) % grainSize];
        grainData[i] += linearInterpolation(a, b, j % 1.0) * this.grainWindow[i];
    }

    // Copy the grain multiple times overlapping it
    for (i = 0; i < grainSize; i += Math.round(grainSize * (1 - overlapRatio))) {
        for (j = 0; j <= grainSize; j++) {
            this.buffer[i + j] += grainData[j];
        }
    }

    // Output the first half of the buffer
    for (i = 0; i < grainSize; i++) {
        outputData[i] = this.buffer[i];
    }
};