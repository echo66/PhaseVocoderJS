/*
params : {
	id: String, 
	audioContext (opt): AudioContext, 
	bufferSize (default 4096): Number, 
	sampleRate (default 44100): Number, 
	div (default 2): Number
}
 */
function AudioPlayer(params) {

	var _id = params.id;
	var _audioContext = (params.audioContext)? params.audioContext : new AudioContext();
	var _samplePosition = 0;
	var _synSamplePosition = 0;
	var _bufferSize = params.bufferSize || 4096;
	var _sampleRate = params.sampleRate || 44100;
	var _div = params.div || 2;
	var _pvL = new PhaseVocoder(_bufferSize/_div, _sampleRate); _pvL.init();
	var _pvR = new PhaseVocoder(_bufferSize/_div, _sampleRate); _pvR.init();
	var _outBufferL = [];
	var _outBufferR = [];
	var _node = context.createScriptProcessor(_bufferSize, 2, 2);
	var _stretchFactor = 1;
	var _pitchFactor = 1;
	var _playing = false;
	var _loaded = false;
	var _grid;
	var _buffer;
	var _connected = false;

	_node.onaudioprocess = function (e) {

		if (!_playing || !_loaded) return;

		var il = _buffer.getChannelData(0);
		var ir = _buffer.getChannelData(1);

		var ol = e.outputBuffer.getChannelData(0);
		var or = e.outputBuffer.getChannelData(1);

		// Fill output buffers (left & right) until the system has 
		// enough processed samples to reproduce.
		do {

			if (!_playing || !_loaded) return;

			if (_samplePosition/_sampleRate >= _buffer.duration || !_loaded) break;

			var bufL = il.subarray(_samplePosition, _samplePosition+_bufferSize/_div);
			var bufR = ir.subarray(_samplePosition, _samplePosition+_bufferSize/_div);

			_samplePosition += _pvL.get_analysis_hop();
			_synSamplePosition += _pvL.get_synthesis_hop();

			// Process left input channel
			_outBufferL = _outBufferL.concat(_pvL.process(bufL));

			// Process right input channel
			_outBufferR = _outBufferR.concat(_pvR.process(bufR));

		} while(_outBufferL.length < _bufferSize);

		ol.set(_outBufferL.splice(0, _bufferSize));
		or.set(_outBufferR.splice(0, _bufferSize));

	};


	Object.defineProperties(this, {
		'id' : {
			get: function() {
				return _id;
			}
		},
		'time' : {
			get: function() {
				return _synSamplePosition/_sampleRate;
			},
			set: function(position) {
				var __samplePosition = position*_sampleRate;
				if (!_loaded)
					throw "Load an audio buffer in AudioPlayer before changing the current time";
				if (__samplePosition >= _buffer.duration)
					_samplePosition = _buffer.length;
				else 
					_samplePosition = __samplePosition;
			}
		}, 
		'duration' : {
			get: function() {
				if (_loaded) 
					return _buffer.duration;
			}
		}, 
		'stretch' : {
			get: function() {
				return _stretchFactor;
			},
			set: function(factor) {
				_stretchFactor = factor;
				_pvL.set_alpha(_stretchFactor);
				_pvR.set_alpha(_stretchFactor);
			}
		},
		'pitch' : {
			get: function() {
				return _pitchFactor;
			},
			set: function(factor) {
				_pitchFactor = factor;
			}
		},
		'node' : {
			value: _node
		},
		'audio' : {
			get: function() {
				return _buffer;
			}
		},
		'isPlaying' : {
			get: function() { return _playing; }
		},
		'isLoaded' : {
			get: function() { return _loaded; }
		},
		'isConnected' : {
			get: function() { return _connected; }
		}
	});

	this.play   = function () {
		_playing = true;
	}

	this.stop   = function () {
		_playing = false;
	}

	// params: {buffer, grid}
	this.load   = function (params) {
		if (!_loaded && !_playing) {
			_buffer = params.buffer;
			_grid = params.grid;
			_loaded = true;
		} if (_playing)
			throw "Stop AudioPlayer before loading new audio";
	}

	// params: {start, end}
	this.loop = function(params) {
		// TODO: Set loop if params.start and params.end exist. Otherwise, unset the loop.
	}

	this.unload = function () {
		if (_loaded && !_playing) {
			_loaded = false;
			_buffer = undefined;
			_grid = undefined;
			_samplePosition = 0;
			_pitchFactor = 1;
			_stretchFactor = 1;
		} if (_playing)
			throw "Stop AudioPlayer before unloading the current audio";
	}

	this.connect = function (dest) {
		_node.connect(dest);
		_connected = true;
	}

	this.disconnect = function() {
		_node.disconnect();
		_connected = false;
	}

	

	var _emit = function(evenType, data) {
		for (var ci in _callbacks[evenType]) 
		_callbacks[evenType][ci](data);
	}

	this.on = function(observerID, eventType, callback) {

		if (!eventType || _callbacks[eventType]==undefined) 
			throw "Unsupported event type";

		if (observerID!=undefined && _callbacks[eventType][observerID]!=undefined) 
			throw "Illegal modification of callback";

		var __id = (observerID==undefined)? _id + "-associate-" + (_idCounter++) : observerID;
		_callbacks[eventType][__id] = callback;

		return __id;
	}

	this.off = function(observerID, eventType) {
		if (!eventType || _callbacks[eventType]==undefined) 
			throw "Unsupported event type";

		delete _callbacks[eventType][observerID];
	}

}

/*
AudioPlayer

Methods:

	AudioPlayer.play: () -> undefined

	AudioPlayer.stop: () -> undefined

	AudioPlayer.load: (data) -> undefined

	AudioPlayer.unload: () -> undefined

	AudioPlayer.time: (number) -> number

	AudioPlayer.pitch: (number) -> number

	AudioPlayer.stretch: (number) -> number

Events:

	"loaded"
	"unloaded"
	"progress"
	"changed-time"
	"changed-pitch"
	"changed-stretch"

*/