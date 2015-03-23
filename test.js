var BUFFER_SIZE = 2048;

var context = new AudioContext();

var node1 = new PVNode(context, BUFFER_SIZE, 2);
var node2 = new PVNode(context, BUFFER_SIZE, 2);

loadSample = function(url) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    request.responseType = 'arraybuffer';

    request.onload = function() {
        console.log('url loaded');
        context.decodeAudioData(request.response, function(decodedData) {
            buffer = decodedData;
            node1.setAudioData(buffer);
            node2.setAudioData(buffer);
        });
    }

    console.log('reading url');
    request.send();
}

loadSample('../soundtouchjs/4.mp3');




function setAlpha(newAlpha) {
    node1.setAlpha(newAlpha);
    node2.setAlpha(newAlpha);
}

function setPosition(v) {
    node1.setPosition(v);
    node2.setPosition(v);
}

function resetPVs() {
    //TODO
}

function resetPVs2() {
    //TODO
}

function play() {
    node1.play();
    // node2.play();
}

function pause() {
    node1.pause();
    node2.pause();
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