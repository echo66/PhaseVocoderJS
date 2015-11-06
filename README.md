# PhaseVocoder.js

A JavaScript implementation of the Phase Vocoder algorithm, with Identity Phase Locking, to perform time stretching. This implementation is independent of the Web Audio API but each time stretcher instance can be integrated in a ScriptProcessor or AudioWorker node.

demo: http://echo66.github.io/demos/PhaseVocoder.js/

# Constructor

*PhaseVocoder(Number frameSize, Number sampleRate)*: frameSize and sampleRate are integers.

# API

*process(Array inputFrame, CBuffer outputFrame)*: given a (mono) frame, performs a time stretching iteration and pushes H s samples in the output CBuffer.

*set_stft_fn(Function stftCallback)*: stftCallback(Array inputFrame, Array windowFrame, Number wantedSize, Object out), *inputFrame* is a sequence of samples; *windowFrame* is the discretization of the window function; *wantedSize* is the desired size for the real and imaginary arrays of the output 14 ; *out* is a JSON object with four arrays: real, imaginary, magnitude and phase, all of them describing the result of the forward Short Time Fourier Transform (STFT). This function will be invoked for each time frame being processed by PhaseVocoder.js

*set_istft_fn(Function istftCallback)*: istftCallback(Array real, Array imag, Array windowFrame, Array timeFrame), *real* and *imag* are the real and imaginary arrays describing a frequency frame; timeFrame is the result of the inverse STFT. This function will be invoked when synthesizing a frequency frame, after the phase adaptation.

*init_fft_fn(Function initCallback)*: this function is invoked after being added to a PhaseVocoder.js instance.

*clear_buffers()*: clears all internal buffers, like the overlapping buffer. This can be useful for audio players that need to create a noticeable stop in the transition to the next file in a playlist, in order to avoid using the phase of the previous song to adjust the phase of the next song.

*set_alpha(Number alpha, Number overlap, Number beta)*: given the new stretching factor, it computes the new values for Hs , Ha (both integers) and invokes the function pointed by *overlap_fn*.

*overlap_fn(Number alpha)*: public field pointing to a function that, given a stretching factor α, will return a new overlapping factor.

*alpha_step(Number alpha)*: TODO

*get_alpha()*: returns the last specified stretching factor.

*get_ha()*: returns the current analysis hop size. This function calculates the increment to the “read head” of the input signal, when playing an audio file.

*get_hs()*: returns the current synthesis hop size. This function calculates the increment to the output signal position which an be used to guide the cursor in the UI of an audio player using OLA-TS.js as time stretcher.

*get_overlap_factor()*: returns the current overlapping factor.



# Helper Classes

*BufferedTS*: it manages the intermediary frame buffering for two PhaseVocoder.js instances. This class offers two public writable fields, *alpha* and *position* to manipulate the stretching factor and the 'read head' of the input audio buffer, as well as two public methods:
* *process(AudioBuffer outputAudioBuffer)*: writes the next output frame in the provided output audio buffer.
* *set_audio_buffer(AudioBuffer newBuffer)*: defines the input audio buffer.

*WAAPlayer*: integrates a *BufferedTS* instance in a ScriptProcessor node, providing the usual functions *connect(AudioNode)* and *disconnect(AudioNode)*, two extra public methods, *play* and *stop*, as well as four public writable fields: *position*, *speed*, *audioContext* and *audioBuffer*.

# Notes

* The audio quality and computational costs are heavily dependent on two things (for both percussive and non-percussive sounds): (1) the overlapping factor of each frame and (2) the number of bins of each Fourier Transform. If we use a lower overlapping factor, we decrease the computational costs but but, for higher stretching factors, there will be some stuttering in the output. If we use a lower number of bins for the FFT, we decrease the computational costs of the algorithm but decrease the frequency resolution, meaning that we will, probably, miss several spectral peaks and, as such, the identity phase locking will not work so well. On the other hand, if we increase the number of bins, we will introduced additional transient smearing, resulting in undesired distortion in percussive sounds.

* Due to the computational costs of the phase vocoder algorithm (i.e.: the number of Fourier Transforms is equal to 2 x overlappingFactor), there will be issues like audio dropouts when using this implementation to play several songs, simultaneously.

# Roadmap

* Implement pitch shiting using a frequency-based method.

* Integrate new frequency-based effects like robotization, harmonization and whisperization.