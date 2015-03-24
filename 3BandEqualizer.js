function 3BandEqualizer(audioSource, context) {

  var gain_node = context.createGain();
  gain_node.gain.value = 1;

  var lowpass = createFilter("lowshelf", 100, 1, 0);

  var highpass = createFilter("highshelf", 4000, 0, 0);

  var mid = createFilter("peaking", 2000, 2, 0);

  lowpass.connect(mid);
  mid.connect(highpass);
  highpass.connect(gain_node);

  function createFilter(type, freq, Q, gain) {
    var filter = context.createBiquadFilter();
    mid.type = type;
    mid.frequency.value = freq;
    mid.Q.value = Q;
    mid.gain.value = gain;
    return filter;
  }

  /*
   *  Range: [0; 50]
   */
  function linear_to_db(linear) {
    if (linear < 50) {
      return -(db_to_linear(50.0 - linear));
    } else if (linear === 50) {
      return 0;
    } else {
      return db_to_linear(linear - 50.0);
    }
  }

  function db_to_linear(db) {
    return Math.pow(10, db / 30.0) - 1.0;
  }

  /*
   *  Range: [0, 2]
   *  Volume=1 produces sound without any gain.
   *  Volume=2 doubles the natural volume of the input sound.
   */ 
  this.set_volume = function(newVolume) {
    if (newVolume>=0 && newVolume<=2)
      gain_node.gain.value = newVolume;
  }

  /*
   *  Range: [0; 2]
   *  Gain=0 is the maximum attenuation.
   *  Gain=2 doubles the natural gain of the input sound.
   */
  this.set_low_gain = function(newGain) {
    if (newGain>=0 && newGain<=2)
      lowpass.gain.value = linear_to_db(50*newGain);
  }

  /*
   *  Range: [0; 2]
   *  Gain=0 is the maximum attenuation.
   *  Gain=2 doubles the natural gain of the input sound.
   */
  this.set_mid_gain = function(newGain) {
    if (newGain>=0 && newGain<=2)
      mid.gain.value = linear_to_db(50*newGain);
  }

  /*
   *  Range: [0; 2]
   *  Gain=0 is the maximum attenuation.
   *  Gain=2 doubles the natural gain of the input sound.
   */
  this.set_hi_gain = function(newGain) {
    if (newGain>=0 && newGain<=2)
      highpass.gain.value = linear_to_db(50*newGain);
  }



}