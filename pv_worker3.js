// importScript("PV.js");
// importScript("jensnockert-fft.js/lib/complex.js");
// importScript("jensnockert-fft.js/lib/real.js");

var BUFFER_SIZE = 2048;

var phasevocoderL = new PhaseVocoder(BUFFER_SIZE, 44100); 
phasevocoderL.init();
var phasevocoderR = new PhaseVocoder(BUFFER_SIZE, 44100); 
phasevocoderR.init();

var outBufferL = [];
var outBufferR = [];

console.log("configurado");

// Parameters: alpha
onaudioprocess = function (e) {
  var il = e.inputBuffers[0];
  var ir = e.inputBuffers[0];

  var ol = e.outputBuffers[0];
  var or = e.outputBuffers[1];

  var alpha = e.parameters.alpha[0];

  phasevocoderL.set_alpha(alpha);
  phasevocoderR.set_alpha(alpha);

  do {

        outBufferL = outBufferL.concat(phasevocoderL.process(il));

        outBufferR = outBufferR.concat(phasevocoderR.process(ir));

  } while (outBufferL.length < BUFFER_SIZE);

  // console.log("done");

  for (var i = 0; i < BUFFER_SIZE; i++) {
      ol[i] = outBufferL.shift();
      or[i] = outBufferR.shift();
      // console.log([i, ol[i], or[i]]);
  }

};