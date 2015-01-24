function do_diff_v2(test, mine, length, frame, opt) {

	var amine = mine.slice(0,length);
	var atest = test.slice(0,length);
	var c = new Array(length);
	var wrongSign = [];

	var avgMine = amine.reduce(function(a, b) { return a + b })/amine.length;
	var avgTest = atest.reduce(function(a, b) { return a + b })/atest.length;

	var maxMine = Math.max.apply(Math, amine);
	var minMine = Math.min.apply(Math, amine);

	var maxTest = Math.max.apply(Math, atest);
	var minTest = Math.min.apply(Math, atest);

	var normalizedMine = amine.map(function(x,i){
		return (x-minMine)/(maxMine-minMine);
	});

	var normalizedTest = mine.map(function(x,i){
		return (x-minTest)/(maxTest-minTest);
	});

	for(var j=0; j<c.length; j++) {
		if(Math.sign(amine[j])!=Math.sign(atest[j]) && opt.showWrongSign) {
			wrongSign.push({
				valName: opt.valName, 
				frame: frame, 
				index: j,
				mineVal: amine[j],
				testVal: atest[j]
			});
		}
		var l2norm = Math.sqrt((normalizedMine[j]-normalizedTest[j])*(normalizedMine[j]-normalizedTest[j]));
		c[j] = (l2norm>1e-10)? l2norm : 0;
	}
	
	return {
		meanL2norm: c.reduce(function(a,b){ return a+b;}) / length,
		L2normVector: c,
		wrongSignAttrs: wrongSign
	}

}

function mag(real,imag) { return Math.sqrt(real*real+imag*imag);}

function ang(real,imag) { return Math.atan2(imag,real); }