
/* File: PlanetMaker.js */
/**
 * PlanetMaker JavaScript Library
 * http://planetmaker.wthr.us
 * 
 * Copyright 2013 Kevin M. Gill <kmsmgill@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */


if (!window.KMG) { window.KMG = {}; };
KMG.Util = {};


self.console = self.console || {

	info: function () {},
	log: function () {},
	debug: function () {},
	warn: function () {},
	error: function () {}

};

KMG.RAD_360 = 360 * (Math.PI / 180);
KMG.RAD_270 = 270 * (Math.PI / 180);
KMG.RAD_180 = 180 * (Math.PI / 180);
KMG.RAD_90 = 90 * (Math.PI / 180);
KMG.RAD_45 = 45 * (Math.PI / 180);
KMG.AU_TO_KM = 149597870.700;
KMG.PI_BY_180 = (Math.PI / 180);
KMG._180_BY_PI = (180 / Math.PI);


/* File: Math.js */
KMG.Math = {};


KMG.Math.sinh = function(a) {
	return (Math.exp(a) - Math.exp(-a)) / 2;
}

KMG.Math.cosh = function(a) {
	return (Math.pow(Math.E, a) + Math.pow(Math.E, -a)) / 2;
}

KMG.Math.sign = function(a) {
	return (a >= 0.0) ? 1 : -1;
}

KMG.Math.radians = function(d) {
	return d * (Math.PI / 180);
}

KMG.Math.degrees = function(r) {
	return r * (180 / Math.PI);
}

KMG.Math.sqr = function(v) {
	return v * v;
};

KMG.Math.clamp = function(v, within) {
	if (!within) {
		within = 360;
	}
	return v - within * Math.floor(v / within);
};




KMG.Math.dsin = function(v) {
	return Math.sin(v * Math.PI / 180);
};

KMG.Math.dcos = function(v) {
	return Math.cos(v * Math.PI / 180);
};

KMG.Math.dtan = function(v) {
	return Math.tan(v * Math.PI / 180);
};


KMG.Math.dasin = function(v) {
	return Math.asin(v) * 180 / Math.PI;
};

KMG.Math.dacos = function(v) {
	return Math.acos(v) * 180 / Math.PI;
};

KMG.Math.datan2 = function(y, x) {
	return Math.atan2(y, x) * 180 / Math.PI;
};

KMG.Math.datan = function(v) {
	return KMG.Math.dasin(v / Math.sqrt(Math.pow(v, 2) + 1));
};

KMG.Math.degToRad = function(v) {
	return v * KMG.PI_BY_180;
};

KMG.Math.radToDeg = function(v) {
	return v * KMG._180_BY_PI;
};

KMG.Math.round = function(v, factor) {
	if (!factor) {
		factor = 1000;
	}
	
	return Math.floor(v * factor) / factor;
};

KMG.Math.trimTo360Radians = function(x) {	
	if( x > 0.0 ) {
		while( x > KMG.RAD_360 )
			x = x-KMG.RAD_360;
	} else {
		while( x< 0.0 )
			x =x+ KMG.RAD_360;
	}
	return x;
}

KMG.Math.trimTo360 = function(x) {	
	if( x > 0.0 ) {
		while( x > 360.0 )
			x = x-360.0;
	} else {
		while( x< 0.0 )
			x =x+ 360.0;
	}
	return x;
}


KMG.Math.fixThetaDegrees = function(degrees){
	var limited;
    degrees /= 360.0;
    limited = 360.0 * (degrees - Math.floor(degrees));
    if (limited < 0)
		limited += 360.0;
	return limited;
}
        
KMG.Math.fixPhiDegrees = function(degrees) {
	degrees += 90.0;
	var limited;
	degrees /= 180.0;
	limited = 180.0 * (degrees - Math.floor(degrees));
	if (limited < 0)
			limited += 180.0;
	return limited - 90.0;
}




KMG.Math.getPoint3D = function(theta, // Longitude, in degrees
                               phi, // Latitude, in degrees
                               radius) {

	//theta += 90.0;
	theta = KMG.Math.fixThetaDegrees(theta) * KMG.PI_BY_180;
    phi = KMG.Math.fixPhiDegrees(phi) * KMG.PI_BY_180;

                
	var _y = Math.sqrt(KMG.Math.sqr(radius) - KMG.Math.sqr(radius * Math.cos(phi)));
    var r0 = Math.sqrt(KMG.Math.sqr(radius) - KMG.Math.sqr(_y));

    var _b = r0 * Math.cos(theta );
    var _z = Math.sqrt(KMG.Math.sqr(r0) - KMG.Math.sqr(_b));
    var _x = Math.sqrt(KMG.Math.sqr(r0) - KMG.Math.sqr(_z));

                
    if (theta <= KMG.RAD_90) {
		_z *= -1.0;
	} else if (theta  <= KMG.RAD_180) {
		_x *= -1.0;
		_z *= -1.0;
	} else if (theta  <= KMG.RAD_270) {
		_x *= -1.0;
	}

	if (phi >= 0) { 
		_y = Math.abs(_y);
	} else {
		_y = Math.abs(_y) * -1;
	}
	return new THREE.Vector3(_x, _y, _z);
}


// http://www.mathworks.com/help/aeroblks/radiusatgeocentriclatitude.html
KMG.Math.radiusAtGeocentricLatitude = function(equatorialRadius, latitude, flattening) {
	
	var R = equatorialRadius;
	var l = latitude;
	var f = flattening;
	
	var r = Math.sqrt(KMG.Math.sqr(R, 2) / (1 + (1 / KMG.Math.sqr(1 - f) - 1) * KMG.Math.sqr(Math.sin(l))));
	return r;
};





/* File: Util.js */

KMG.Util = {};

KMG.Util.isUserMobile = function()
{
	return /Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent);
}

KMG.Util.cardinalDirectionByValue = function(value, ifPos, ifNeg) {
	return (value >= 0) ? ifPos : ifNeg;
}

KMG.Util.formatDegrees = function(value, ifPos, ifNeg) {
	
	value = KMG.Math.round(value, 1000);
	
	var fmt = Math.abs(value) + "&deg;";
	if (ifPos && ifNeg) {
		fmt += KMG.Util.cardinalDirectionByValue(value, ifPos, ifNeg);
	}
	return fmt;
}

KMG.Util.formatDegreesToHours = function(value, ifPos, ifNeg) {
	
	var h = Math.floor(Math.abs(value));
	var m = Math.floor((Math.abs(value) - h) * 60);
	var s = KMG.Math.round(((Math.abs(value) - h) * 60 - m) * 60, 100);
	
	if (h < 10) {
		h = "0" + h;
	}
	
	if (m < 10) {
		m = "0" + m;
	}
	
	var fmt = h + "h " + m + "m " + s + "s";
	if (ifPos && ifNeg) {
		fmt += KMG.Util.cardinalDirectionByValue(value, ifPos, ifNeg);
	}
	return fmt;
}

KMG.Util.formatDegreesToMinutes = function(value, ifPos, ifNeg, skipSeconds) {
	
	var d = Math.floor(Math.abs(value));
	var m = Math.floor((Math.abs(value) - d) * 60);
	var s = KMG.Math.round(((Math.abs(value) - d) * 60 - m) * 60, 100);
	
	var sign = (value < 0) ? "-" : " ";
	
	var fmt = sign + d + "&deg; " + m + "\' " + ((!skipSeconds) ? (s + "\"") : "");
	if (ifPos && ifNeg) {
		fmt += KMG.Util.cardinalDirectionByValue(value, ifPos, ifNeg);
	}
	return fmt;
}


KMG.Util.intensityToWhiteColor = function(intensity)
{
	intensity = parseInt(intensity);
	var rgb = "rgb("+intensity+","+intensity+","+intensity+")";
	return new THREE.Color(rgb);
};

KMG.Util.arrayToColor = function(array)
{
	var r = parseInt(array[0]);
	var g = parseInt(array[1]);
	var b = parseInt(array[2]);
	var a = (array.length >= 4) ? parseInt(array[3]) : 255.0;
	var rgb = "rgb("+r+","+g+","+b+")";
	return new THREE.Color(rgb);
};

KMG.Util.rgbToArray = function(rgb)
{
	var c = new THREE.Color(rgb);
	return new THREE.Vector3(c.r, c.g, c.b);
}


KMG.Util.eyePosition = function(context)
{
	return context.camera.position;
};
	
KMG.Util.eyeDistanceToCenter = function(context)
{
	return context.primaryScene.position.distanceTo(KMG.Util.eyePosition(context));
};
	
KMG.Util.surfaceDistance = function(context, radius)
{
	return KMG.Util.eyeDistanceToCenter(context) - radius;
};


//farClipDistance
KMG.Util.horizonDistance = function(context, radius)
{
	var r = radius;
	var e = KMG.Util.surfaceDistance(context, radius);
	var f = Math.sqrt(e * (2 * r + e));
	return f;
};

// TODO: Strengthen this...
KMG.Util.clone = function(object) 
{
	if (typeof object !== "object") {
		return object;
	}

	var cloned;
	
	if (object instanceof Array) {
		cloned = new Array();
		for (var i = 0; i < object.length; i++) {
			cloned.push(KMG.Util.clone(object[i]));
		}
	} else {
		cloned = {};
		for(var key in object) {
			cloned[key] = KMG.Util.clone(object[key]);
		};
	}
	
	return cloned;
};

// TODO: Strengthen this...
KMG.Util.extend = function(target, source) {

	//var extended = KMG.Util.clone(target);
	var extended = target;
	for(var key in source) {
		if (extended[key] === undefined) {
			extended[key] = KMG.Util.clone(source[key]);
		}
		
		if (extended[key] && typeof extended[key] === "object") {
			extended[key] = KMG.Util.extend(extended[key], source[key]);
		}
	
	}
	
	return extended;
};

//http://stackoverflow.com/questions/1714786/querystring-encoding-of-a-javascript-object
KMG.Util.serialize = function(obj, prefix) {
	var str = [];
	for(var p in obj) {
		var k = prefix ? prefix + "[" + p + "]" : p, v = obj[p];
		str.push(typeof v == "object" ? KMG.Util.serialize(v, k) : encodeURIComponent(k) + "=" + encodeURIComponent(v));
	}
	return str.join("&");
};

/*
KMG.Util.getTimezoneOffset = function() {
	return (new Date()).getTimezoneOffset() * 1000;
};

KMG.Util.getTimezoneOffsetJulians = function() {
	return (new Date()).getTimezoneOffset() / 60 / 24;
};
*/


KMG.Util.replaceWithGreekLetters = function(str) {
	
	str = str.replace(/alpha/g, 'α');
	str = str.replace(/beta/g, 'β');
	
	str = str.replace(/gamma/g, 'γ');
	str = str.replace(/delta/g, 'δ');
	str = str.replace(/epsilon/g, 'ε');
	str = str.replace(/zeta/g, 'ζ');
	str = str.replace(/eta/g, 'η');
	str = str.replace(/theta/g, 'θ');
	str = str.replace(/iota/g, 'ι');
	str = str.replace(/kappa/g, 'κ');
	str = str.replace(/lambda/g, 'λ');
	str = str.replace(/mu/g, 'μ');
	str = str.replace(/nu/g, 'ν');
	str = str.replace(/xi/g, 'ξ');
	str = str.replace(/omicron/g, 'ο');
	str = str.replace(/pi/g, 'π');
	str = str.replace(/rho/g, 'ρ');
	str = str.replace(/sigma/g, 'σ');
	str = str.replace(/tau/g, 'τ');
	str = str.replace(/upsilon/g, 'υ');
	str = str.replace(/phi/g, 'φ');
	str = str.replace(/chi/g, 'χ');
	str = str.replace(/psi/g, 'ψ');
	str = str.replace(/omega/g, 'ω');
	
	
	
	return str;
};

KMG.Util.replaceWithGreekLettersAbbreviated = function(str) {
	
	str = str.replace(/Alp/g, 'α');
	str = str.replace(/Bet/g, 'β');
	
	str = str.replace(/Gam/g, 'γ');
	str = str.replace(/Del/g, 'δ');
	str = str.replace(/Eps/g, 'ε');
	str = str.replace(/Zet/g, 'ζ');
	str = str.replace(/Eta/g, 'η');
	str = str.replace(/The/g, 'θ');
	str = str.replace(/Iot/g, 'ι');
	str = str.replace(/Kap/g, 'κ');
	str = str.replace(/Lam/g, 'λ');
	str = str.replace(/Mu/g, 'μ');
	str = str.replace(/Nu/g, 'ν');
	str = str.replace(/Xi/g, 'ξ');
	str = str.replace(/Omi/g, 'ο');
	str = str.replace(/Pi/g, 'π');
	str = str.replace(/Rho/g, 'ρ');
	str = str.replace(/Sig/g, 'σ');
	str = str.replace(/Tau/g, 'τ');
	str = str.replace(/Ups/g, 'υ');
	str = str.replace(/Phi/g, 'φ');
	str = str.replace(/Chi/g, 'χ');
	str = str.replace(/Psi/g, 'ψ');
	str = str.replace(/Ome/g, 'ω');
	
	
	
	return str;
};

/** A best-effort attempt to convert a string to an intended data type when the intended type may not be known.
 * Returns the supplied parameter as-is if the type cannot be determined as anything other than a string.
 * 
 */ 
KMG.Util.stringToDataType = function(v) {
	if (/^true$/.test(v))
		return true;
	else if (/^false$/.test(v))
		return false;
	else if (/^-?\d+\.?\d*$/.test(v)) 
		return parseFloat(v);
	else if (/\,/.test(v)) {
		var a = v.split(",");
		var list = [];
		for (var i = 0; i < a.length; i++) {
			list.push(KMG.Util.stringToDataType(a[i]));
		}
		return list;
	} else {
		v = v.replace(/%20/g, " ");
		return v;
	}
}

/* File: AstroDate.js */

KMG.J2000 = 2451545.0;
KMG.J1900 = 2415020.0;

KMG.AstroDate = function(millis, jd, epoch) {
	
	if (!epoch) {
		epoch = KMG.J2000;
	}
	
	if (!millis && !jd) {
		jd = KMG.Util.julianNow();
	}
	
	if (millis != undefined && millis != null) {
		jd = KMG.Util.millisToJulian(millis);
	}
	
	// Julian day overrides millis
	if (jd != undefined && jd != null) {
		millis = KMG.Util.julianToDate(jd).getTime();
	}

	function getMillis() {
		return millis;
	}
	
	function getJulianDay() {
		return jd;
	}
	
	function getJulianCentury(_epoch) {
		if (!_epoch) {
			_epoch = epoch;
		}
		return (jd - _epoch) / 36525.0;
	}
	
	function getDayNumberNow() {
		return KMG.Astro.getDayNumberNow();
	}
	
	function getDayNumber(_jd) {
		if (!_jd) {
			_jd = jd;
		}
		return KMG.Astro.getDayNumber(_jd);
	}
	
	function getGMST(clampTo) {
		return KMG.Astro.getGMST(jd, clampTo);
	}
	
	function getLMST(lon) {
		return KMG.Astro.getLMST(jd, lon);
	}
	
	//http://www.satellite-calculations.com/TLETracker/scripts/tletracker.online.sat.calc
	function getGMST2(lon, _jd) {
		if (!_jd) {
			_jd = jd;
		}
		return KMG.Astro.getGMST2(lon, _jd);
	}
	
	function getDate(_jd) {
		if (!_jd) {
			_jd = jd;
		}
		return KMG.Util.julianToDate(_jd);
	}
	
	function getDayOfYear() {
		return KMG.Astro.getDayOfYear(jd);
	}
	
	function toString(format) {
		var d = getDate();
		return d.toString();
	}
		
	function toJSON() {
		return {
			millis : millis,
			jd : jd,
			date : toString() + " UTC"
		};
	}
	return {
		getMillis : getMillis,
		getJulianDay : getJulianDay,
		getGMST : getGMST,
		getLMST : getLMST,
		getGMST2 : getGMST2,
		getJulianCentury : getJulianCentury,
		getDayOfYear : getDayOfYear,
		getDayNumberNow : getDayNumberNow,
		getDayNumber : getDayNumber,
		getDate : getDate,
		toString : toString,
		toJSON : toJSON
	
	};
	
};
/* File: Astronomy.js */
/**
 * Astronomy Algorithms
 * http://www.apoapsys.com
 * 
 * Copyright 2014 Kevin M. Gill <kmsmgill@gmail.com>
 *
 * Uses algorithms from:
 *
 * 
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */

KMG.Astro = {};
 
KMG.Astro.J2000 = 2451545.0;
KMG.Astro.J1900 = 2415020.0;
 
// http://www.csgnetwork.com/juliangregcalconv.html
KMG.Astro.julianToDate = function(jd) {
	var _jd = jd;
		
	jd += 0.5;
	var z = Math.floor(jd);
	var f = jd - z;
	if (z < 2299161) {
		var A = z;
	} else {
		var omega = Math.floor((z-1867216.25)/36524.25);
		var A = z + 1 + omega - Math.floor(omega/4);
	}
	var B = A + 1524;
	var C = Math.floor((B-122.1)/365.25);
	var D = Math.floor(365.25*C);
	var Epsilon = Math.floor((B-D)/30.6001);
	var dayGreg = B - D - Math.floor(30.6001*Epsilon) + f;
	var monthGreg, yearGreg;
	if (Epsilon < 14) {
		monthGreg = Epsilon - 1;
	} else {
		monthGreg = Epsilon - 13;
	}
	if (monthGreg > 2) {
		yearGreg = C - 4716;
	} else {
		yearGreg = C - 4715;
	}
	
	var year = yearGreg;
	var month = monthGreg;
	var day = Math.floor(dayGreg);
	
	var dayMinutes = ((dayGreg - day) * 1440.0);
	var hour = Math.floor(dayMinutes / 60.0);
	var minute = Math.floor(dayMinutes - (hour * 60.0));
	var second = Math.floor(60.0 * (dayMinutes - (hour * 60.0) -minute));
	var millisecond =  (1000.0 * (60.0 * (dayMinutes - (hour * 60.0) -minute)- second) );
	
	return new Date(year, month - 1, day, hour, minute, second, millisecond);
};
KMG.Util.julianToDate = KMG.Astro.julianToDate;

// http://www.csgnetwork.com/juliandatetime.html
KMG.Astro.dateToJulian = function(year, month, day, hour, minute, second, millisecond, tz) {

	/*
	var extra = 100.0*year + month - 190002.5;
	var rjd = 367.0*year;
	rjd -= Math.floor(7.0*(year+Math.floor((month+9.0)/12.0))/4.0);
	rjd += Math.floor(275.0*month/9.0) ;
	rjd += day;
	rjd += (hour + (minute + second/60.0)/60.)/24.0;
	rjd += 1721013.5;
	rjd -= 0.5*extra/Math.abs(extra);
	rjd += 0.5;
	
	return rjd;*/
	if (!tz) {
		tz = 0;
	}
	if (!millisecond) {
		millisecond = 0;
	}
	var day_decimal, julian_day, a;

	day_decimal = day + (hour - tz + (minute + second / 60.0 + millisecond / 1000 / 60) / 60.0) / 24.0;

	if (month < 3) {
		month += 12;
		year--;
	}

	julian_day = Math.floor(365.25 * (year + 4716.0)) + Math.floor(30.6001 * (month + 1)) + day_decimal - 1524.5;
	if (julian_day > 2299160.0) {
		a = Math.floor(year / 100);
		julian_day += (2 - a + Math.floor(a / 4));
	}
	
	return julian_day;
};
KMG.Util.dateToJulian = KMG.Astro.dateToJulian;

KMG.Astro.millisToJulian = function(millis) {
	var d = new Date(millis);
	var jd =  KMG.Util.dateToJulian(d.getFullYear(),
							d.getMonth() + 1,
							d.getDate(),
							d.getHours(),
							d.getMinutes(),
							d.getSeconds(),
							d.getMilliseconds());
	return jd;
};
KMG.Util.millisToJulian = KMG.Astro.millisToJulian;

KMG.Astro.julianNow = function() {
	var d = new Date();
	return KMG.Util.dateToJulian(d.getUTCFullYear(),
						d.getUTCMonth() + 1,
						d.getUTCDate(),
						d.getUTCHours(),
						d.getUTCMinutes(),
						d.getUTCSeconds(),
						d.getUTCMilliseconds());
	
}
KMG.Util.julianNow = KMG.Astro.julianNow;

KMG.Astro.formatJulianDay = function(jd, isUtc, format) {
	if (!format) {
		format = "LLL";
	}
	var dt = KMG.Util.julianToDate(jd);
	if (isUtc) 
		return moment(dt).format(format);
	else 
		return moment(dt).utc().format(format);
};
KMG.Util.formatJulianDay = KMG.Astro.formatJulianDay;




KMG.Astro.convertCartesianEquatorialToEcliptic = function(equatorial) {
	var e = 23.4 * KMG.PI_BY_180;
	
	var m = new THREE.Matrix3();
	m.set(1, 0, 0, 
		  0, Math.cos(e), Math.sin(e), 
		  0, -Math.sin(e), Math.cos(e));
	
	var ecliptic = equatorial.clone().applyMatrix3(m);
	return ecliptic;
};
KMG.Math.convertCartesianEquatorialToEcliptic = KMG.Astro.convertCartesianEquatorialToEcliptic;



// Duffet-Smith, Peter: Practical Astronomy with Your Calculator, page 42
KMG.Astro.convertEquatorialToEcliptic = function(ra, dec) {
	var e = 23.4;
	
	var Y = KMG.Math.dsin(ra) * KMG.Math.dcos(e) + KMG.Math.dtan(dec) * KMG.Math.dsin(e);
	var X = KMG.Math.dcos(ra);
	
	var l = KMG.Math.datan2(Y, X);
	var b = KMG.Math.dasin(KMG.Math.dsin(dec) * KMG.Math.dcos(e) - KMG.Math.dcos(dec) * KMG.Math.dsin(e) * KMG.Math.dsin(ra));
	
	return {
		l : l,
		b : b
	};
};
KMG.Math.convertEquatorialToEcliptic = KMG.Astro.convertEquatorialToEcliptic;

// Adapted from Celestia customorbits.cpp:
// static double Obliquity(double t)
KMG.Astro.obliquity = function(t)
{
    // Parameter t represents the Julian centuries elapsed since 1900.
    // In other words, t = (jd - 2415020.0) / 36525.0

    return (2.345229444E1 - ((((-1.81E-3*t)+5.9E-3)*t+4.6845E1)*t)/3600.0) * KMG.PI_BY_180;
}
KMG.Math.obliquity = KMG.Astro.obliquity;

// Adapted from Celestia customorbits.cpp:
// static void Nutation(double t, double &deps, double& dpsi)
KMG.Astro.nutation = function(t)
{
    // Parameter t represents the Julian centuries elapsed since 1900.
    // In other words, t = (jd - 2415020.0) / 36525.0

    var ls, ld;	// sun's mean longitude, moon's mean longitude
    var ms, md;	// sun's mean anomaly, moon's mean anomaly
    var nm;	    // longitude of moon's ascending node
    var t2;
    var tls, tnm, tld;	// twice above
    var a, b;

    t2 = t*t;

    a = 100.0021358*t;
    b = 360.*(a-Math.floor(a));
    ls = 279.697+.000303*t2+b;

    a = 1336.855231*t;
    b = 360.*(a-Math.floor(a));
    ld = 270.434-.001133*t2+b;

    a = 99.99736056000026*t;
    b = 360.*(a-Math.floor(a));
    ms = 358.476-.00015*t2+b;

    a = 13255523.59*t;
    b = 360.*(a-Math.floor(a));
    md = 296.105+.009192*t2+b;

    a = 5.372616667*t;
    b = 360.*(a-Math.floor(a));
    nm = 259.183+.002078*t2-b;

    //convert to radian forms for use with trig functions.
    tls = 2*KMG.Math.degToRad(ls);
    nm = KMG.Math.degToRad(nm);
    tnm = 2*KMG.Math.degToRad(nm);
    ms = KMG.Math.degToRad(ms);
    tld = 2*KMG.Math.degToRad(ld);
    md = KMG.Math.degToRad(md);

    // find delta psi and eps, in arcseconds.
    var dpsi = (-17.2327-.01737*t)*Math.sin(nm)+(-1.2729-.00013*t)*Math.sin(tls)
        +.2088*Math.sin(tnm)-.2037*Math.sin(tld)+(.1261-.00031*t)*Math.sin(ms)
        +.0675*Math.sin(md)-(.0497-.00012*t)*Math.sin(tls+ms)
        -.0342*Math.sin(tld-nm)-.0261*Math.sin(tld+md)+.0214*Math.sin(tls-ms)
        -.0149*Math.sin(tls-tld+md)+.0124*Math.sin(tls-nm)+.0114*Math.sin(tld-md);
    var deps = (9.21+.00091*t)*Math.cos(nm)+(.5522-.00029*t)*Math.cos(tls)
        -.0904*Math.cos(tnm)+.0884*Math.cos(tld)+.0216*Math.cos(tls+ms)
        +.0183*Math.cos(tld-nm)+.0113*Math.cos(tld+md)-.0093*Math.cos(tls-ms)
        -.0066*Math.cos(tls-nm);

    // convert to radians.
    dpsi = (dpsi/3600) * KMG.PI_BY_180;
    deps = (deps/3600) * KMG.PI_BY_180;
    
    //double &deps, double& dpsi
    return {
		deps : deps,
		dpsi : dpsi
	};
}
KMG.Math.nutation = KMG.Astro.nutation;


// Adapted from Celestia customorbits.cpp:
// static void EclipticToEquatorial(double t, double fEclLat, double fEclLon,
//                                 double& RA, double& dec) 
KMG.Astro.convertEclipticToEquatorial = function(jd, fEclLat, fEclLon)
{
    // Parameter t represents the Julian centuries elapsed since 1900.
    // In other words, t = (jd - 2415020.0) / 36525.0

    var seps, ceps;	// sin and cos of mean obliquity
    var sx, cx, sy, cy, ty;
    var eps;
	var dec, ra;
	
	
	var t = (jd - 2415020.0) / 36525.0;
   // t = (2451545.0 - 2415020.0) / 36525.0;
   // t = 0;
    eps = KMG.Math.obliquity(t);		// mean obliquity for date
    var nut = KMG.Math.nutation(t);
    var deps = nut.deps;
    var dpsi = nut.dpsi;
    
    eps += deps;
    seps = Math.sin(eps);
    ceps = Math.cos(eps);

    sy = Math.sin(fEclLat);
    cy = Math.cos(fEclLat);				// always non-negative
    if (Math.abs(cy)<1e-20)
        cy = 1e-20;		// insure > 0
    ty = sy/cy;
    cx = Math.cos(fEclLon);
    sx = Math.sin(fEclLon);
    dec = Math.asin((sy*ceps)+(cy*seps*sx));
    
  //  ra = Math.atan(((sx*ceps)-(ty*seps))/cx);
	ra = Math.atan2(((sx*ceps)-(ty*seps)), cx);
    //if (cx<0)
   //     ra += Math.PI;		// account for atan quad ambiguity
	ra = KMG.Math.clamp(ra, 2 * Math.PI);
    
    
    return {
		ra : ra,
		dec : dec
	};
};
KMG.Math.convertEclipticToEquatorial = KMG.Astro.convertEclipticToEquatorial;

// Convert equatorial coordinates from one epoch to another.  Method is from
// Chapter 21 of Meeus's _Astronomical Algorithms_
// Actuall adapted from Celestia customorbits.cpp:
//void EpochConvert(double jdFrom, double jdTo,
//                  double a0, double d0,
//                  double& a, double& d)
KMG.Astro.epochConvert = function(jdFrom, jdTo, a0, d0)
{
	var a, d;
	
    var T = (jdFrom - 2451545.0) / 36525.0;
    var t = (jdTo - jdFrom) / 36525.0;

    var zeta = (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t +
        (0.30188 - 0.000344 * T) * t * t + 0.017998 * t * t * t;
    var z = (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t +
        (1.09468 + 0.000066 * T) * t * t + 0.018203 * t * t * t;
    var theta = (2004.3109 - 0.85330 * T - 0.000217 * T * T) * t -
        (0.42665 + 0.000217 * T) * t * t - 0.041833 * t * t * t;
    zeta  = KMG.Math.degToRad(zeta / 3600.0);
    z     = KMG.Math.degToRad(z / 3600.0);
    theta = KMG.Math.degToRad(theta / 3600.0);

    var A = Math.cos(d0) * Math.sin(a0 + zeta);
    var B = Math.cos(theta) * Math.cos(d0) * Math.cos(a0 + zeta) -
        Math.sin(theta) * Math.sin(d0);
    var C = Math.sin(theta) * Math.cos(d0) * Math.cos(a0 + zeta) +
        Math.cos(theta) * Math.sin(d0);

    a = Math.atan2(A, B) + z;
    d = Math.asin(C);
    
    return {
		ra : a,
		dec : d
	};
};
KMG.Math.epochConvert = KMG.Astro.epochConvert;



//http://www.satellite-calculations.com/TLETracker/scripts/tletracker.online.sat.calc
KMG.Astro.getGMST2 = function(lon, _jd) {
	if (!lon) {
		lon = 0.0; // Handle null or undefined
	}
	if (!_jd) {
		_jd = KMG.Astro.julianNow();
	}
	
	var dt = KMG.Astro.julianToDate(_jd);
	var day = dt.getDate();
	var month = dt.getMonth() + 1;
	var year = dt.getFullYear();
	var hour = dt.getHours();
	var minute  = dt.getMinutes();
	var second = dt.getSeconds();
	var ms = dt.getMilliseconds();
	if( month == 1 || month == 2 )
	{
	year = year - 1;
	month = month + 12;
	}

	var a = Math.floor( year/100 );
	var b = 2 - a + Math.floor( a/4 );

	var c = Math.floor(365.25 * year);
	var d = Math.floor(30.6001 * (month + 1));

	// days since J2000.0   
	var jd = b + c + d - 730550.5 + day + (hour + minute/60.0 + second/3600.0)/24.0;
	
	var jt   = (jd)/36525.0;                   // julian centuries since J2000.0         
	var GMST = 280.46061837 + 360.98564736629*jd + 0.000387933*jt*jt - jt*jt*jt/38710000 + lon;           
	if( GMST > 0.0 )
	{
		while( GMST > 360.0 )
			GMST -= 360.0;
	}
	else
	{
		while( GMST < 0.0 )
			GMST += 360.0;
	}
		
	return GMST;
};


KMG.Astro.div = function(a, b) {
	return ((a-a%b)/b);
};


KMG.Astro.getDayNumber = function(_jd) {
	var dt = KMG.Astro.julianToDate(_jd);
	var dd = dt.getDate();
	var mm = dt.getMonth() + 1;
	var yyyy = dt.getFullYear();
	var hh = dt.getHours();
	var min  = dt.getMinutes();
	var sec = dt.getSeconds();
	var ms = dt.getMilliseconds();

	var d=367.0*yyyy - KMG.Astro.div(  (7.0*(yyyy+(KMG.Astro.div((mm+9.0),12.0)))),4.0 ) + KMG.Astro.div((275.0*mm),9.0) + dd - 730530.0 ;
	d=d+ hh/24.0 + min/(60.0*24.0) + sec/(24.0*60.0*60.0);

	return d;
};

KMG.Astro.getDayNumberNow = function() {
	return KMG.Astro.getDayNumber(KMG.Astro.julianNow());
};


KMG.Astro.getGMST = function(jd, clampTo) {
	if (!clampTo) {
		clampTo = 24;
	}
	var jd0 = Math.floor(jd + 0.5) - .5;
	var H = (jd - jd0) * 24;
	var D = jd - 2451545.0;
	var D0 = jd0 - 2451545.0;
	var T = D / 36525.0;
	var gmst = (6.697374558 + 0.06570982441908 * D0 + 1.00273790935 * H + 0.000026 * Math.pow(T, 2));
	gmst = KMG.Math.clamp(gmst, clampTo);
	return gmst;
};

KMG.Astro.getLMST = function(jd, lon) {
	var gst = KMG.Astro.getGMST(jd);
	return gst + (lon / 15);
};

KMG.Astro.getDayOfYear = function(jd) {
	var date = KMG.Astro.julianToDate(jd);
	var start = new Date(date.getFullYear(), 0, 0);
	var diff = date - start;
	var oneDay = 1000 * 60 * 60 * 24;
	var day = Math.floor(diff / oneDay);
	return day;
};


KMG.Astro.vectorToHeliocentricLatitudeLongitude = function(vec, skipRotation) {
	/* Opposite of
	var x = Math.cos(l) * Math.sin(b) * r;
	var y = Math.cos(b) * r;
	var z = -Math.sin(l) * Math.sin(b) * r;
	*/
	
	var x, y, z;

	x = vec.x;
	y = vec.y;
	z = vec.z;
	
	var r = Math.sqrt(Math.pow(vec.x, 2) + Math.pow(vec.y, 2) + Math.pow(vec.z, 2));
	b = Math.abs(Math.acos(y / r)) * ((z < 0) ? -1 : 1);
	l = 2 * Math.PI - Math.acos((vec.x * (1 / Math.sin(b))) / r); 

	var rotB = (skipRotation) ? 0 : Math.PI / 2;
	var rotL = (skipRotation) ? 0 : -Math.PI;
	
	return {
		 r : r,
		 b : (b + rotB) * KMG._180_BY_PI,
		 l : (l + rotL) * KMG._180_BY_PI,
	 };
	
};



KMG.Astro.calculatePositionVector = function(date, orbit) {
	var position = orbit.positionAtTime(date.getJulianDay(), false);
	var E = position.E;
	var M = position.M;
	var trueAnomaly = position.trueAnomaly;
	
	return {
		x : position.x,
		y : position.y,
		z : position.z,
		E : position.E,
		M : position.M,
		trueAnomaly : position.trueAnomaly
	};
};

/**
 *
 * @param date
 * @param 
 */
//KMG.Astro.calculateSatellitePosition = function(date, p, E, M) {
KMG.Astro.calculateSatellitePosition = function(date, orbit) {

	var position = orbit.positionAtTime(date.getJulianDay(), false);
	var E = position.E;
	var M = position.M;
	var trueAnomaly = position.trueAnomaly;

	
	
	var semiMajorAxis = orbit.orbitProperties.semiMajorAxis;
	var arg_per = orbit.orbitProperties.argOfPeriapsis;
	var RAAN = orbit.orbitProperties.rightAscension;
	var i = orbit.orbitProperties.inclination * KMG.PI_BY_180;
	var e = orbit.orbitProperties.eccentricity;
	
	//var Epoch_now = date.getJulianDay() - apoapsys.J2000;
	var Epoch_now = date.getDayNumber();
	var Epoch_start = orbit.orbitProperties.epochStart;
	var Earth_equatorial_radius = 6378.135;
	

	var first_derative_mean_motion = orbit.orbitProperties.derivativeOfMeanMotion;
	var Satellite_rev_sidereal_day = orbit.orbitProperties.meanMotion;
	
	var TCdecimal=(1440/((1*Satellite_rev_sidereal_day)+(first_derative_mean_motion*(Epoch_now-Epoch_start )))) /60;   // Period in hours


	// bug here ?
	var RangeA=Math.pow(    (6028.9* (TCdecimal*60)), (2/3)   );
	
	
	var apogee =RangeA*(1+e*1);   // apogee
	var perigee=RangeA*(1-e*1);  //perigee
	var semimajoraxsis=(1*apogee+1*perigee)/2;  // semimajoraxsis
	

	
	perigee_perturbation=(Epoch_now-Epoch_start)*4.97*Math.pow((Earth_equatorial_radius/(1*semimajoraxsis)) , 3.5   )* (  5*Math.cos(i)*Math.cos(i) -1)/((1-e*e)*(1-e*e));

	// perturbation of ascending node

	ascending_node_perturbation=(Epoch_now-Epoch_start)*9.95*Math.pow((Earth_equatorial_radius/(1*semimajoraxsis)) , 3.5   )*   Math.cos(i)/((1-e*e)*(1-e*e));


	// perbutation of perigee
	arg_per=arg_per +perigee_perturbation;
	RAAN=RAAN-ascending_node_perturbation;

	arg_per *= KMG.PI_BY_180;
	RAAN *= KMG.PI_BY_180;

	var X0=1.0*semiMajorAxis*(Math.cos(E)-e);  //  = r*Cos(trueanomaly)
	var Y0=1.0*semiMajorAxis*Math.sqrt(1-e*e)*Math.sin(E);  // = r*sin (trueanomaly)
	var r=Math.sqrt(X0*X0+Y0*Y0); // distance
	
	X0 *= KMG.AU_TO_KM;
	Y0 *= KMG.AU_TO_KM;
	r *= KMG.AU_TO_KM;
	
	var Px = Math.cos(arg_per)*Math.cos(RAAN) - Math.sin(arg_per)*Math.sin(RAAN)*Math.cos(i);
	var Py = Math.cos(arg_per)*Math.sin(RAAN) + Math.sin(arg_per)*Math.cos(RAAN)*Math.cos(i);
	var Pz = Math.sin(arg_per)*Math.sin(i);

	var Qx=-Math.sin(arg_per)*Math.cos(RAAN)-Math.cos(arg_per)*Math.sin(RAAN)*Math.cos(i);
	var Qy=-Math.sin(arg_per)*Math.sin(RAAN)+Math.cos(arg_per)*Math.cos(RAAN)*Math.cos(i);
	var Qz=Math.cos(arg_per)*Math.sin(i);
	
	x=Px*X0+Qx*Y0;
	y=Py*X0+Qy*Y0;
	z=Pz*X0+Qz*Y0;

	
	var dec = Math.atan2(  z,Math.sqrt(x*x+y*y)  );
	var ra = Math.atan2(  y,x );

	
	ra = ra%(2 * Math.PI);

	var gmst = date.getGMST2(0);
	
	var lon = Math.atan2( y,x ) - (gmst * KMG.PI_BY_180);
	lon = KMG.Math.clamp(lon, 2 * Math.PI);
	
	if (lon > Math.PI) {
		lon = -1 * (2 * Math.PI - lon);
	}
	
	var lat = Math.atan2(  z,Math.sqrt(x*x+y*y)  );

	var radius = KMG.Math.radiusAtGeocentricLatitude(Earth_equatorial_radius, lat, 0.0033528);
	
	var altitude = Math.sqrt(x*x+y*y+z*z) - radius;
	var rh = (radius + altitude);
	var theta = Math.acos(radius / rh);
	var t = (6.284*rh)*(Math.sqrt(rh/398600))/60;
	
	var angularRange = radius * Math.tan(theta);
	var surfaceRange = 220 * (theta * KMG._180_BY_PI);
	var vis = ((2*theta*KMG._180_BY_PI)/360)*t;
	
	var maxAngularRange = .5 * Math.PI * 6378.135;
	if (angularRange > maxAngularRange) {
		angularRange = maxAngularRange;
	}
	
	//pos.ra_hms = apoapsys.Util.formatDegreesToHours(pos.ra / 15);
	//	pos.dec_dms = apoapsys.Util.formatDegreesToMinutes(pos.dec);
	return {
		ra : ra * KMG._180_BY_PI,
		ra_hms : KMG.Util.formatDegreesToHours(ra * KMG._180_BY_PI / 15),
		dec : dec * KMG._180_BY_PI,
		dec_dms : KMG.Util.formatDegreesToMinutes(dec * KMG._180_BY_PI),
		lat : lat * KMG._180_BY_PI,
		lon : lon * KMG._180_BY_PI,
		altitude : altitude,
		angularRange : angularRange,
		surfaceRange : surfaceRange,
		visibilityTimeMinutes : vis,
		eccentricAnomaly : E,
		meanAnomaly : M,
		trueAnomaly : trueAnomaly,
		
		xyz : {
			x : x, 
			y : y,
			z : z
		}
	};
};



KMG.Astro.getPositionAzimuthalSatellite = function(date, orbit, lat, lon) {
	//var equatorial = getPositionEquatorial(date, fromObject);
	//lat = 42.76537;
	//lon = -71.46757;
	
	var equatorial = KMG.Astro.calculateSatellitePosition(date, orbit);
	
	var dlon = lon;
	var dlat = lat;
	lat = lat * KMG.PI_BY_180;
	lon = lon * KMG.PI_BY_180;
	var ra = equatorial.ra * KMG.PI_BY_180;
	var dec = equatorial.dec * KMG.PI_BY_180;


	var gmst = date.getGMST2(0);
	var lst = gmst + (dlon / 15);
	
	var ha = lst - (((ra * KMG._180_BY_PI) / 15));
	var f = 0.0033528;
	var e2 = 2 * f - f * f;
	var C=1/Math.sqrt(1+0.0033528*(0.0033528-2)*Math.sin(lat)*Math.sin(lat)             );
	var omega=(1*gmst+1*dlon) * KMG.PI_BY_180;

	var Re = 6378.135;
	var C = 1 / Math.sqrt(1 - e2 * Math.pow(Math.sin(lat), 2));
	var S=(1-0.0033528)*(1-0.0033528)*C;
	var R =  Re * Math.cos(lat);
	var a = Re;
	
	/*
	console.info("omega: " + (omega * KMG._180_BY_PI));
	console.info("C: " + C);
	console.info("S: " + S);
	console.info("R: " + R);
	console.info("lat: " + (lat * KMG._180_BY_PI));
	console.info("lon: " + (lon * KMG._180_BY_PI));
	console.info("GMST: " + gmst);
	*/
	
	var x_ = a * C * Math.cos(lat) * Math.cos(omega);
	var y_ = a * C * Math.cos(lat) * Math.sin(omega);
	var z_=6378.135*S*Math.sin(lat);
	//console.info([x_, y_, z_]);

	var xs = equatorial.xyz.x;
	var ys = equatorial.xyz.y;
	var zs = equatorial.xyz.z;
	
	
	
	var xo = x_;
	var yo = y_;
	var zo = z_;
	
	//console.info(["sat", xs, ys, zs]);
	//console.info(["Obs", xo, yo, zo]);
	
	var rx=xs-xo;
	var ry=ys-yo;
	var rz=zs-zo;
	
	//console.info(["rxyz", rx, ry, rz]); 
	
	fi=(1*gmst+1*dlon)*KMG.PI_BY_180;
	var rS=Math.sin(lat)*Math.cos(fi)*rx+Math.sin(lat)*Math.sin(fi)*ry-Math.cos(lat)*rz;
	var rE=-Math.sin(fi)*rx+Math.cos(fi)*ry;
	var rZ=Math.cos(lat)*Math.cos(fi)*rx+Math.cos(lat)*Math.sin(fi)*ry+Math.sin(lat)*rz;
	//console.info("fi: " + (fi * KMG._180_BY_PI));
	//console.info(["rSEZ", rS, rE, rZ]);

	var range=Math.sqrt(rS*rS+rE*rE+rZ*rZ);
	//console.info("range: " + range);
	var Elevation=Math.asin(rZ/range);
	var Azimuth=Math.atan(-rE/rS);

	if (rS>0) Azimuth=Azimuth+Math.PI;
	if (Azimuth<0) Azimuth=Azimuth+ 2*Math.PI;

	var alt = Elevation;
	var az = Azimuth;

	if (az < 0) {
		az += 2 * Math.PI;
	}
	
	//console.info("Azimuth: " + (az * KMG._180_BY_PI));
	//console.info("Elevation: " + (alt * KMG._180_BY_PI));
	
	return {
		alt : alt * KMG._180_BY_PI,
		az : KMG.Math.clamp(az * KMG._180_BY_PI, 360),
		ha : ha,
		lst : lst,
		lst_hms : KMG.Util.formatDegreesToHours(lst),
		range : range,
		gmst : gmst
	};

};


KMG.Astro.getHeliocentricPosition = function(date, orbit) {

	var position = orbit.positionAtTime(date.getJulianDay(), false);
	var parentPosition, heliocentricPosition;

	if (orbit.parentOrbit) {
		parentPosition = orbit.parentOrbit.positionAtTime(date.getJulianDay(), false);
		heliocentricPosition = parentPosition.clone().add(position);
	} else {
		heliocentricPosition = position;
	}
	
	return heliocentricPosition;
}

KMG.Astro.getHeliocentricPositionEquatorial = function(date, orbit, fromObjectOrbit) {
	
	var heliocentricPosition = KMG.Astro.getHeliocentricPosition(date, orbit);
		
	// fromObject is assumed to be Earth
	if (fromObjectOrbit) {
		var fromHeliocentricPosition = KMG.Astro.getHeliocentricPosition(date, fromObjectOrbit);
		heliocentricPosition.sub(fromHeliocentricPosition);
	}
	
	var pos;
	
	var eclipticPos = KMG.Astro.vectorToHeliocentricLatitudeLongitude(heliocentricPosition);
	pos = KMG.Math.convertEclipticToEquatorial(date.getJulianDay(), eclipticPos.b *KMG.PI_BY_180, eclipticPos.l *KMG.PI_BY_180);
	pos = KMG.Math.epochConvert(date.getJulianDay(), 2451545.0, pos.ra, pos.dec);
	pos.ra *= KMG._180_BY_PI;
	pos.dec *= KMG._180_BY_PI;
	if (pos.dec < 0) {
		pos.ra += 180;
	}

	if (!pos.ra_hms) {
		pos.ra_hms = KMG.Util.formatDegreesToHours(pos.ra / 15);
	}
	if (!pos.dec_dms) {
		pos.dec_dms = KMG.Util.formatDegreesToMinutes(pos.dec);
	}

	return pos;
};
	
	
KMG.Astro.getPositionAzimuthalHeliocentricBody = function(date, orbit, fromObjectOrbit, lat, lon) {

	var equatorial = KMG.Astro.getHeliocentricPositionEquatorial(date, orbit, fromObjectOrbit);

	var dlon = lon;
	var dlat = lat;
	lat = lat * KMG.PI_BY_180;
	lon = lon * KMG.PI_BY_180;
	var ra = equatorial.ra * KMG.PI_BY_180;
	var dec = equatorial.dec * KMG.PI_BY_180;

	var lst = date.getLMST(dlon);

	var ha = lst - (((ra * KMG._180_BY_PI) / 15));

	var H = ha * 15 * KMG.PI_BY_180;
	var alt = Math.asin(Math.sin(dec) * Math.sin(lat) + Math.cos(dec) * Math.cos(lat) * Math.cos(H));
	var az = Math.asin(-Math.sin(H) * Math.cos(dec) / Math.cos(alt));
	//alt += (h * KMG.PI_BY_180);
	
	var Y = -Math.sin(H);
	var X = Math.cos(lat) * Math.tan(dec) - Math.sin(lat) * Math.cos(H);
	var az = Math.atan2(Y, X);

	if (az < 0) {
		az += 2 * Math.PI;
	}

	return {
		alt : alt * KMG._180_BY_PI,
		az : KMG.Math.clamp(az * KMG._180_BY_PI, 360),
		ha : ha,
		lst : lst,
		lst_hms : KMG.Util.formatDegreesToHours(lst)
	};

};

KMG.Astro.geocentricToECI = function(date, lat, lon, alt) {

	var F = 0.0033528;
	var mfactor = 2 * Math.PI * (1.00273790934 / 86400.0); 

	theta = KMG.Astro.getGMST2(0.0, date.getJulianDay());
	theta = KMG.Math.clamp(theta + lon, 360) * KMG.PI_BY_180;
	
	
	lat *= KMG.PI_BY_180;
	lon *= KMG.PI_BY_180;
	
	var c = 1.0 / Math.sqrt(1.0 + (1.0 / 298.26) * (F - 2.0) * KMG.Math.sqr(Math.sin(lat)));   
	var s = KMG.Math.sqr(1.0 - F) * c;   
	var achcp = (KMG.AU_TO_KM * c + alt) * Math.cos(lat);   
	
	var x = achcp * Math.cos(theta);
	var y = achcp * Math.sin(theta);
	var z = (KMG.AU_TO_KM * s + alt) * Math.sin(lat);
	var w = Math.sqrt(x*x + y*y + z*z);
	
	return {
		x : x,
		y : y,
		z : z,
		w : w
	};
};


/* File: UserLocationController.js */


KMG.Location = {};



	
KMG.Location.isPreciseLocationSupported = function() {
	if (navigator.geolocation)
		return true;
	else
		return false;
};
	
KMG.Location.getPreciseLocation = function(onLocation) {
	if (!KMG.Location.isPreciseLocationSupported()) {
		return false;
	}
	
	navigator.geolocation.getCurrentPosition(onLocation);
	return true;
};
	
	

/* File: StarFlares.js */



KMG.starFlares = [
	{
		name : "Yellow Star",
		texture : "/img/lensflare0.png"
	}
];
/* File: Backgrounds.js */


KMG.backgrounds = [
	{
		name : "Starfield",
		
		// I'm not sure the exact origin of this image. If you know, or you made it, 
		// please let me know at kevin@wthr.us so I can provide proper credit.
		texture : "/img/starfield2_1900x1250.jpg",
		enabled : true
	}
];

/* File: Textures.js */

KMG.textures = [
	{
		name : "Earth - Blue Marble / Oceans & Ice",
		texture : "/img/tx_composite_adjusted_4096x2048_1.jpg",
		bumpMap : "/img/tx_composite_adjusted_4096x2048_1.jpg",
		normalMap : "/img/tx_composite_adjusted_4096x2048_1.jpg",
		specularMap : "/img/tx_composite_adjusted_4096x2048_1.jpg",
		enabled : true
	},
	{
		name : "Mars",
		texture : "/img/tx_composite.adjusted_#resolution#.jpg",
		bumpMap : "/img/mars_mola_bumpmap_#resolution#.jpg",
		normalMap : "/img/mars_mola_normalmap_#resolution#.jpg",
		specularMap :"/img/flat_black_texture.jpg",
		enabled : true
	}
];


/* File: DefaultConfig.js */

KMG.DefaultConfig = {
	version : 2.0,
	initWithShadows : true,
	shadows : false,
	shadowDarkness : 0.5,
	radius : 200,
	textureResolution : "4096x4096",
	enableFps : false,
	postprocessingEnabled : true,
	
	useScript : true,
	
	// Light
	lightingType : "Directional", // or "Point"
	sunlightDirection : 60.0,
	realtimeSunlight : false,
	sunlightDate : (new Date()).getTime(),
	localStarDistance : 1.0,
	displayLocalStar : true,
	localStarTexture : KMG.starFlares[0].name,
	localStarColor : [ 255, 255, 255 ],
	starColorAffectsPlanetLighting : true,
	lensFlareEnabled : false,
	
	// Surface
	texture : KMG.textures[0].name,
	surfaceDetail : 0.0,
	elevationScale : 0,
	shininess : 0,
	diffuseIntensity : 170,
	specularIntensity : 4,
	ambientIntensity : 255,
	emissiveIntensity : 180,
	flattening : 0.0033528,
	axialTilt : 0.0,
	surfaceColorMode : "Normal",
	surfaceHue : 0.5,
	surfaceSaturation : 0.0,
	surfaceLightness : 0.75,
	surfaceWrapRGB : 0.031,
	surfaceRotation : 0.0,
	scaleSurface : 1.0, 

	
	
	// Background
	backgroundType : 'stars',
	backgroundImage : 'Starfield',
	backgroundImageType : 'flat',
	backgroundImageFitType : 'stretch',
	starQuantity : 6.5, // 0 - 10
	
	
	camera : {
		positionZ : 700,
		fieldOfView : 45,
		near : 0.01,
		far : 10000000,
		
		useSecondaryParameters : false,
		fieldOfViewSecondary : 45,
		nearSecondary : 0.01,
		farSecondary : 10000000
	},
	controls : {
		rotateSpeed : 0.5
	}
};

/* File: TextureMap.js */

KMG.TextureMap = {

	map : {},
	textureResolution : "2048x1024",
	texturesLoading : 0,
	sceneReadyCallback : null,
	resourceLoadingStart : null,
	resourceLoadingFinish : null,
	renderCallback : null,
	
	onResourceLoaded : function()
	{
		KMG.TextureMap.texturesLoading--;
		if (KMG.TextureMap.sceneReadyCallback && KMG.TextureMap.texturesLoading === 0) {
			KMG.TextureMap.sceneReadyCallback();
		}
		
		if (KMG.TextureMap.resourceLoadingFinish) {
			KMG.TextureMap.resourceLoadingFinish(true, KMG.TextureMap.texturesLoading);
		}
	},
	
	setupEncodedTexture : function(dat) {
		var img = new Image();
		var t = new THREE.Texture(img);
		t.wrapS = THREE.RepeatWrapping;

		img.onload = function() {
			t.needsUpdate = true;
			KMG.TextureMap.onResourceLoaded();
			if (KMG.TextureMap.renderCallback !== null) {
				KMG.TextureMap.renderCallback();
			}
		};
		img.src = dat;
		return t;
	},

	
	loadTexture : function(url, onload, noCache)
	{

		if (!url || url.length === 0) {
			return null;
		}
		
		if (!onload) {
			onload = KMG.TextureMap.onResourceLoaded;
		} else {
			var origOnload = onload;
			onload = function() {
				origOnload();
				KMG.TextureMap.onResourceLoaded();
			};
		}
		
		url = url.replace("#resolution#", KMG.TextureMap.textureResolution);
		
		if (KMG.TextureMap.map[url] !== undefined && !noCache) {
			return KMG.TextureMap.map[url];
		}
		
		if (KMG.TextureMap.resourceLoadingStart) {
			KMG.TextureMap.resourceLoadingStart(url);
		}
		
		KMG.TextureMap.texturesLoading++;
		
		var tex = null;
		if (/^data:/i.test(url)) {
			tex = KMG.TextureMap.setupEncodedTexture(url);
			onload();
		} else {
			tex = THREE.ImageUtils.loadTexture( url, {}, onload, onload );
		}
		
		tex.repeat.set( 0.998, 0.998 );
		tex.offset.set( 0.001, 0.001 )
		tex.wrapS = tex.wrapT = THREE.RepeatWrapping;
		tex.format = THREE.RGBFormat;
		//KMG.TextureMap.map[url].anisotropy = 4;
		
		if (!noCache) {
			KMG.TextureMap.map[url] = tex;
		} 
		
		return tex;
	},
	
	getDefinitionByName : function(list, name) 
	{
		for (var i = 0; i < list.length; i++) {
			if (list[i].name == name) {
				return list[i];
			}
		}
		return null;
	},
	
	getTextureDefinitionByName : function( name )
	{
		return KMG.TextureMap.getDefinitionByName(KMG.textures, name);
	},
	
	getBackgroundDefinitionByName : function( name )
	{
		return KMG.TextureMap.getDefinitionByName(KMG.backgrounds, name);
	},
	
	getFlareDefinitionByName : function( name )
	{
		return KMG.TextureMap.getDefinitionByName(KMG.starFlares, name);
	}

};

/* File: Procession.js */


KMG.Procession = {};

/*
0 Pc;
1 Qc;
2 Ps;
3 Qs;
4 period;
*/
KMG.Procession.EclipticPrecessionTerms =
[
    [   486.230527, 2559.065245, -2578.462809,   485.116645, 2308.98 ],
    [  -963.825784,  247.582718,  -237.405076,  -971.375498, 1831.25 ],
    [ -1868.737098, -957.399054,  1007.593090, -1930.464338,  687.52 ],
    [ -1589.172175,  493.021354,  -423.035168, -1634.905683,  729.97 ],
    [   429.442489, -328.301413,   337.266785,   429.594383,  492.21 ],
    [ -2244.742029, -339.969833,   221.240093, -2131.745072,  708.13 ]
];

/*
0 pc;
1 epsc;
2 ps;
3 epss;
4 period;
*/
KMG.Procession.PrecessionTerms =
[
    [ -6180.062400,   807.904635, -2434.845716, -2056.455197,  409.90 ],
    [ -2721.869299,  -177.959383,   538.034071,  -912.727303,  396.15 ],
    [  1460.746498,   371.942696, -1245.689351,   447.710000,  536.91 ],
    [ -1838.488899,  -176.029134,   529.220775,  -611.297411,  402.90 ],
    [   949.518077,   -89.154030,   277.195375,   315.900626,  417.15 ],
    [    32.701460,  -336.048179,   945.979710,    12.390157,  288.92 ],
    [   598.054819,   -17.415730,  -955.163661,   -15.922155, 4042.97 ],
    [  -293.145284,   -28.084479,    93.894079,  -102.870153,  304.90 ],
    [    66.354942,    21.456146,     0.671968,    24.123484,  281.46 ],
    [    18.894136,    30.917011,  -184.663935,     2.512708,  204.38 ]
];


KMG.Procession.EclipticPrecession_P03LP = function(T)
{

    var T2 = T * T;
    var T3 = T2 * T;

	var PA = (5750.804069
               +  0.1948311 * T
               -  0.00016739 * T2
               -  4.8e-8 * T3);
    var QA = (-1673.999018
               +   0.3474459 * T
               +   0.00011243 * T2
               -   6.4e-8 * T3);


    var nTerms = KMG.Procession.EclipticPrecessionTerms.length / KMG.Procession.EclipticPrecessionTerms[0].length;
		
    for (var i = 0; i < nTerms; i++) {
        var p = KMG.Procession.EclipticPrecessionTerms[i];
        var theta = 2.0 * Math.PI * T / p[4];
        var s = Math.sin(theta);
        var c = Math.cos(theta);
        PA += p[0] * c + p[2] * s;
        QA += p[1] * c + p[3] * s;
    }

    return {
		PA : PA,
		QA : QA
	};
};




KMG.Procession.PrecObliquity_P03LP = function(T)
{

    var T2 = T * T;
    var T3 = T2 * T;

    var pA   = (  7907.295950
                   + 5044.374034 * T
                   -    0.00713473 * T2
                   +    6e-9 * T3);
    var epsA = (  83973.876448
                   -     0.0425899 * T
                   -     0.00000113 * T2);

    var nTerms = KMG.Procession.PrecessionTerms.length / KMG.Procession.PrecessionTerms[0].length;
    

    for (var i = 0; i < nTerms; i++) {
		var p = KMG.Procession.PrecessionTerms[i];
        
        var theta = 2.0 * Math.PI * T / p[4];
        var s = Math.sin(theta);
        var c = Math.cos(theta);
        pA   += p[0] * c   + p[2] * s;
        epsA += p[1] * c + p[3] * s;
    }

    return {
		pA : pA,
		epsA : epsA
	};
};


KMG.Procession.EquatorialPrecessionAngles_P03 = function(T)
{
    var T2 = T * T;
    var T3 = T2 * T;
    var T4 = T3 * T;
    var T5 = T4 * T;
    
    var zetaA =  (     2.650545
                   + 2306.083227 * T
                   +    0.2988499 * T2
                   +    0.01801828 * T3
                   -    0.000005971 * T4
                   -    0.0000003173 * T5);
    var zA =     ( -    2.650545 
                    + 2306.077181 * T
                    +    1.0927348 * T2
                    +    0.01826837 * T3
                    -    0.000028596 * T4
                    -    0.0000002904 * T5);
    var thetaA = (   2004.191903 * T
                   -     0.4294934 * T2
                   -     0.04182264 * T3
                   -     0.000007089 * T4
                   -     0.0000001274 * T5);
    
    return {
		zetaA : zetaA,
		zA : zA,
		thetaA : thetaA
	};
};

KMG.Procession.EclipticPrecession_P03 = function(T)
{
    var T2 = T * T;
    var T3 = T2 * T;
    var T4 = T3 * T;
    var T5 = T4 * T;
    
    var PA = (  4.199094 * T
              + 0.1939873 * T2
              - 0.00022466 * T3
              - 0.000000912 * T4
              + 0.0000000120 * T5);
    var QA = (-46.811015 * T
              + 0.0510283 * T2
              + 0.00052413 * T3
              - 0.00000646 * T4
              - 0.0000000172 * T5);
    
    return {
		PA : PA,
		QA : QA
	};
};


KMG.Procession.EclipticPrecessionAngles_P03 = function(T)
{
    var T2 = T * T;
    var T3 = T2 * T;
    var T4 = T3 * T;
    var T5 = T4 * T;
    
    var piA = ( 46.998973 * T
               - 0.0334926 * T2
               - 0.00012559 * T3
               + 0.000000113 * T4
               - 0.0000000022 * T5);
    var PiA = (629546.7936
                - 867.95758 * T
                +   0.157992 * T2
                -   0.0005371 * T3
                -   0.00004797 * T4
                +   0.000000072 * T5);
    
    return {
		piA : piA,
		PiA : PiA
	};
};

KMG.Procession.eps0 = 84381.40889;

KMG.Procession.PrecObliquity_P03 = function(T)
{
    var T2 = T * T;
    var T3 = T2 * T;
    var T4 = T3 * T;
    var T5 = T4 * T;
    
    var epsA = (KMG.Procession.eps0
                - 46.836769 * T
                -  0.0001831 * T2
                +  0.00200340 * T3
                -  0.000000576 * T4
                -  0.0000000434 * T5);
    var pA   = ( 5028.796195 * T
                +   1.1054348 * T2
                +   0.00007964 * T3
                -   0.000023857 * T4
                -   0.0000000383 * T5);

    return {
		epsA : epsA,
		pA : pA
	};
};

/* File: Orbit.js */



KMG.SolveKeplerFunc1 = function(ecc, M) {
	this.solve = function(x) {
		return M + ecc * Math.sin(x);
	};
};

KMG.SolveKeplerFunc2 = function(ecc, M) {
	this.solve = function(x) {
		return x + (M + ecc * Math.sin(x) - x) / (1 - ecc * Math.cos(x));
	};
};

KMG.SolveKeplerLaguerreConway = function(ecc, M) {
	this.solve = function(x) {
		var s = ecc * Math.sin(x);
		var c = ecc * Math.cos(x);
		var f = x - s - M;
		var f1 = 1 - c;
		var f2 = s;

		x += -5 * f / (f1 + KMG.Math.sign(f1) * Math.sqrt(Math.abs(16 * f1 * f1 - 20 * f * f2)));
		return x;
	};
};

KMG.SolveKeplerLaguerreConwayHyp = function(ecc, M) {
	this.solve = function(x) {
		var s = ecc * KMG.Math.sinh(x);
		var c = ecc * KMG.Math.cosh(x);
		var f = x - s - M;
		var f1 = c - 1;
		var f2 = s;

		x += -5 * f / (f1 + KMG.Math.sign(f1) * Math.sqrt(Math.abs(16 * f1 * f1 - 20 * f * f2)));
		return x;
	};
};

KMG.Orbit = function() 
{



};

// Values valid from 1800 AD through 2050 AD
// See http://iau-comm4.jpl.nasa.gov/keplerformulae/kepform.pdf
KMG.OrbitDefinitions = {
	
	template : {
		semiMajorAxis : 0,
		longitudeOfPerihelion : 0,
		eccentricity : 0,
		inclination : 0,
		ascendingNode : 0, 
		argOfPeriapsis : 0,
		meanAnomalyAtEpoch : 0,
		period : 0
	},
	earth : {
		semiMajorAxis : 1.00000261,
		longitudeOfPerihelion : 102.947,
		eccentricity : 0.0167,
		inclination : 0.0001,
		ascendingNode : 348.73936,
		argOfPeriapsis : 114.20783,
		meanAnomalyAtEpoch : 357.51716,
		period : 1.000017421 * 365.25
	},

	mars : {
		semiMajorAxis : 1.52371034,
		longitudeOfPerihelion : -23.94362959,
		eccentricity : 0.09339410,
		inclination : 1.84969142,
		ascendingNode : 49.55953891, 
		argOfPeriapsis : 286.537,
		meanAnomalyAtEpoch : 19.3564,
		period : 1.8808 * 365.25
	}
};





/**
 * Adapted from http://sourceforge.net/p/celestia/code/5229/tree/trunk/celestia/src/celephem/orbit.cpp
 */
KMG.EllipticalOrbit = function(orbitProperties)
{
	KMG.Orbit.call( this );
	
	this.semiMajorAxis = orbitProperties.semiMajorAxis;
	this.eccentricity = orbitProperties.eccentricity;
	this.inclination = orbitProperties.inclination * (Math.PI / 180.0);
	this.ascendingNode = orbitProperties.ascendingNode * (Math.PI / 180.0);
	this.argOfPeriapsis = orbitProperties.argOfPeriapsis * (Math.PI / 180.0);
	this.meanAnomalyAtEpoch = orbitProperties.meanAnomalyAtEpoch;
	this.period = orbitProperties.period ;
	
	this.orbitProperties = orbitProperties;
	
	this.epoch = (orbitProperties.epoch) ? orbitProperties.epoch : 2451545;
	this.derivativeOfMeanMotion = (orbitProperties.derivativeOfMeanMotion) ? orbitProperties.derivativeOfMeanMotion : 0;
	
	//this.pericenterDistance = (orbitProperties.pericenterDistance) ? orbitProperties.pericenterDistance : this.semiMajorAxis * (1 - this.eccentricity);
	this.pericenterDistance = this.semiMajorAxis * (1 - this.eccentricity);
	//this.meanMotion = (orbitProperties.meanMotion) ? (orbitProperties.meanMotion) : (2.0 * Math.PI) / this.period;
	//this.meanMotion = 1 / this.period;
	this.meanMotion = (orbitProperties.meanMotion) ? (orbitProperties.meanMotion) : 1 / this.period;
	
	var ascendingNodeRotation = new THREE.Matrix4();
	ascendingNodeRotation.makeRotationZ(this.ascendingNode);
	
	var inclinationRotation = new THREE.Matrix4();
	inclinationRotation.makeRotationX(this.inclination);
	
	var argOfPeriapsisRotation = new THREE.Matrix4();
	argOfPeriapsisRotation.makeRotationZ(this.argOfPeriapsis );
	
	this.orbitPlaneRotation = new THREE.Matrix4();
	this.orbitPlaneRotation.identity();
	
	this.orbitPlaneRotation.multiplyMatrices( ascendingNodeRotation, inclinationRotation );
	this.orbitPlaneRotation.multiply( argOfPeriapsisRotation );
	
	
	var scope = this;
	
	function solveIterationFixed(f, x0, maxIter) {
		
		var x = 0;
		var x2 = x0;
		
		for (var i = 0; i < maxIter; i++) {
			x = x2;
			x2 = f.solve(x);
		}
		
		return [x2, x2 - x];
	}
	
	

	
	function eccentricAnomaly(M) {
		if (scope.eccentricity == 0.0) {
			return M;
		} else if (scope.eccentricity < 0.2) {
		
			var sol = solveIterationFixed(new KMG.SolveKeplerFunc1(scope.eccentricity, M), M, 5);

			return sol[0];
		
		} else if (scope.eccentricity < 0.9) {
		
			var sol = solveIterationFixed(new KMG.SolveKeplerFunc2(scope.eccentricity, M), M, 6);
			return sol[0];
		
		} else if (scope.eccentricity < 1.0) {
			var E = M + 0.85 * scope.eccentricity * ((Math.sin(M) >= 0.0) ? 1 : -1);
			
			var sol = solveIterationFixed(new KMG.SolveKeplerLaguerreConway(scope.eccentricity, M), E, 8);
			return sol[0];
			
		} else if (scope.eccentricity == 1.0) {
			return M;
		} else {
			var E = Math.log(2 * M / scope.eccentricity + 1.85);
			
			var sol = solveIterationFixed(new KMG.SolveKeplerLaguerreConwayHyp(scope.eccentricity, M), E, 30);
			return sol[0];
		}
	}
	
	this.positionAtE = function(E) {
		var x, y;
		
		if (this.eccentricity < 1.0) {
			var a = this.pericenterDistance / (1.0 - this.eccentricity);
			x = a * (Math.cos(E) - this.eccentricity);
			y = a * Math.sqrt(1 - this.eccentricity * this.eccentricity) * Math.sin(E);
		} else if (this.eccentricity > 1.0) {
			var a = this.pericenterDistance / (1.0 - this.eccentricity);
			x = -a * (this.eccentricity - KMG.Math.cosh(E));
			y = -a * Math.sqrt(this.eccentricity * this.eccentricity - 1) * KMG.Math.sinh(E);
		} else {
			x = 0.0;
			y = 0.0;
		}
		
		var pos = new THREE.Vector3(x, y, 0);
		pos.applyMatrix4(this.orbitPlaneRotation);
		
		return pos;
	};
	
	
	this.velocityAtE = function(E) {
		var x, y;

		if (this.eccentricity < 1.0) {
			var a = this.pericenterDistance / (1.0 - this.eccentricity);
			var sinE = Math.sin(E);
			var cosE = Math.cos(E);
        
			x = -a * sinE;
			y =  a * Math.sqrt(1 - KMG.Math.sqr(this.eccentricity)) * cosE;
		
			var meanMotion = 2.0 * Math.PI / this.period;
			var edot = meanMotion / (1 - this.eccentricity * cosE);
			x *= edot;
			y *= edot;
		} else if (this.eccentricity > 1.0) {
			var a = this.pericenterDistance / (1.0 - this.eccentricity);
			x = -a * (this.eccentricity - KMG.Math.cosh(E));
			y = -a * Math.sqrt(KMG.Math.sqr(this.eccentricity) - 1) * KMG.Math.sinh(E);
		} else {
			// TODO: Handle parabolic orbits
			x = 0.0;
			y = 0.0;
		}
		
		var v = new THREE.Vector3(x, y, 0);
		v.applyMatrix4(this.orbitPlaneRotation);
		
		return new THREE.Vector3(v.x, v.z, -v.y);
	};
	
	function trimTo360(x) {	
		if( x > 0.0 ) {
			while( x > 360.0 )
				x = x-360.0;
		} else {
			while( x< 0.0 )
				x =x+ 360.0;
		}
		return(x)
	}
	
	this.meanAnomalyAtTime = function(t) {
		var timeSinceEpoch = (t - this.epoch);
		var meanAnomaly = this.meanAnomalyAtEpoch*1+(360*(this.meanMotion*(timeSinceEpoch)+0.5*this.derivativeOfMeanMotion*(timeSinceEpoch)*(timeSinceEpoch))) ; 
		meanAnomaly = trimTo360(meanAnomaly) * (Math.PI / 180);
		return meanAnomaly; // in radians
	};
	
	this.trueAnomalyAtTime = function(t, meanAnomaly, E) {
		if (!meanAnomaly) {
			meanAnomaly = this.meanAnomalyAtTime(t);
		}
		if (!E) {
			E = eccentricAnomaly(meanAnomaly);
		}
		
		var true_anomaly=Math.acos((  Math.cos(E)-this.eccentricity)/(1-this.eccentricity*Math.cos(E))  ) ;
		return true_anomaly; // In radians
	};
	
	this.eccentricAnomalyAtTime = function(t) {
		var meanAnomaly = this.meanAnomalyAtTime(t);
		var E = eccentricAnomaly(meanAnomaly);
		return E;
	};
	
	this.velocityAtTime = function(t) {
		var E = this.eccentricAnomalyAtTime(t);
		var v = this.velocityAtE(E);
		v.multiplyScalar(1.0 / 86400.0);
		v.magnitude = Math.sqrt(KMG.Math.sqr(v.x) + KMG.Math.sqr(v.y) + KMG.Math.sqr(v.z));
		return v;
	};
	
	this.positionAtTime = function(t) {
		var meanAnomaly = this.meanAnomalyAtTime(t);
		var E = eccentricAnomaly(meanAnomaly);
		var pos = this.positionAtE(E);
		var v = new THREE.Vector3(pos.x, pos.z, -pos.y);
		v.E = E;
		v.M = meanAnomaly;
		v.trueAnomaly = this.trueAnomalyAtTime(t, meanAnomaly, E);
		return v;
		//return new THREE.Vector3(pos.x, pos.y, pos.z);
	};
	
	this.distanceAtTime = function(t) {
		var trueAnomaly = this.trueAnomalyAtTime(t);
		var p = this.semiMajorAxis * (1 - KMG.Math.sqr(this.eccentricity));
		var r = p / (1 + this.eccentricity * Math.cos(trueAnomaly));
		return r;
	
	};
	
	this.propertiesAtTime = function(t) {
		var timeSinceEpoch = (t - this.epoch);
		
		var meanAnomalyAtTime = this.meanAnomalyAtTime(t);
		meanAnomalyAtTime = KMG.Math.trimTo360Radians(meanAnomalyAtTime * (Math.PI / 180));
		
		var distanceAtTime = this.distanceAtTime(t);
		var velocityAtTime = this.velocityAtTime(t);
		var trueAnomalyAtTime = this.trueAnomalyAtTime(t);
		
		var eccentricAnomalyAtTime = this.eccentricAnomalyAtTime(meanAnomalyAtTime);
		var positionAtTime = this.positionAtE(eccentricAnomalyAtTime);
		var coordinatesAtTime = new THREE.Vector3(positionAtTime.x, positionAtTime.z, -positionAtTime.y);
		
		
		
		var props = {
			time : t,
			epoch : this.epoch,
			timeSinceEpoch : timeSinceEpoch,
			meanAnomalyAtEpoch : this.meanAnomalyAtEpoch,
			meanAnomalyAtTime : meanAnomalyAtTime,
			distanceAtTime : distanceAtTime, // To center of orbit
			velocityAtTime : velocityAtTime, // km per second
			velocityMagnitudeAtTime : velocityAtTime.magnitude, // km per second
			trueAnomalyAtTime : trueAnomalyAtTime,
			eccentricAnomalyAtTime : eccentricAnomalyAtTime,
			positionAtTime : positionAtTime,
			coordinatesAtTime : coordinatesAtTime,
			meanMotion : this.meanMotion,
			ephemeris : this.orbitProperties
		};
		return props;
	};
	
};
KMG.EllipticalOrbit.prototype = Object.create( KMG.Orbit.prototype );







/* File: VSOP87Orbits.js */


KMG.VSOPTerms = {};

KMG.VSOPTerms.earth_L0 = [
    [ 1.75347045673, 0, 0 ],
    [ 0.03341656453, 4.66925680415, 6283.07584999 ],
    [ 0.00034894275, 4.62610242189, 12566.1517 ],
    [ 3.417572e-05, 2.82886579754, 3.523118349 ],
    [ 3.497056e-05, 2.74411783405, 5753.3848849 ],
    [ 3.135899e-05, 3.62767041756, 77713.7714681 ],
    [ 2.676218e-05, 4.41808345438, 7860.41939244 ],
    [ 2.342691e-05, 6.13516214446, 3930.20969622 ],
    [ 1.273165e-05, 2.03709657878, 529.690965095 ],
    [ 1.324294e-05, 0.74246341673, 11506.7697698 ],
    [ 9.01854e-06, 2.04505446477, 26.2983197998 ],
    [ 1.199167e-05, 1.10962946234, 1577.34354245 ],
    [ 8.57223e-06, 3.50849152283, 398.149003408 ],
    [ 7.79786e-06, 1.17882681962, 5223.6939198 ],
    [ 9.9025e-06, 5.23268072088, 5884.92684658 ],
    [ 7.53141e-06, 2.53339052847, 5507.55323867 ],
    [ 5.05267e-06, 4.58292599973, 18849.22755 ],
    [ 4.92392e-06, 4.20505711826, 775.522611324 ],
    [ 3.56672e-06, 2.91954114478, 0.0673103028 ],
    [ 2.84125e-06, 1.89869240932, 796.298006816 ],
    [ 2.42879e-06, 0.34481445893, 5486.77784318 ],
    [ 3.17087e-06, 5.84901948512, 11790.6290887 ],
    [ 2.71112e-06, 0.31486255375, 10977.0788047 ],
    [ 2.06217e-06, 4.80646631478, 2544.31441988 ],
    [ 2.05478e-06, 1.86953770281, 5573.14280143 ],
    [ 2.02318e-06, 2.45767790232, 6069.77675455 ],
    [ 1.26225e-06, 1.08295459501, 20.7753954924 ],
    [ 1.55516e-06, 0.83306084617, 213.299095438 ],
    [ 1.15132e-06, 0.64544911683, 0.9803210682 ],
    [ 1.02851e-06, 0.63599845579, 4694.00295471 ],
    [ 1.01724e-06, 4.2667980198, 7.1135470008 ],
    [ 9.9206e-07, 6.20992926918, 2146.16541648 ],
    [ 1.32212e-06, 3.41118292683, 2942.46342329 ],
    [ 9.7607e-07, 0.68101342359, 155.420399434 ],
    [ 8.5128e-07, 1.29870764804, 6275.96230299 ],
    [ 7.4651e-07, 1.755089133, 5088.62883977 ],
    [ 1.01895e-06, 0.97569280312, 15720.8387849 ],
    [ 8.4711e-07, 3.67080093031, 71430.6956181 ],
    [ 7.3547e-07, 4.67926633877, 801.820931124 ],
    [ 7.3874e-07, 3.50319414955, 3154.6870849 ],
    [ 7.8757e-07, 3.03697458703, 12036.4607349 ],
    [ 7.9637e-07, 1.80791287082, 17260.1546547 ],
    [ 8.5803e-07, 5.9832263126, 161000.685738 ],
    [ 5.6963e-07, 2.78430458592, 6286.59896834 ],
    [ 6.1148e-07, 1.81839892984, 7084.89678112 ],
    [ 6.9627e-07, 0.83297621398, 9437.76293489 ],
    [ 5.6116e-07, 4.38694865354, 14143.4952424 ],
    [ 6.2449e-07, 3.97763912806, 8827.39026987 ],
    [ 5.1145e-07, 0.28306832879, 5856.47765912 ],
    [ 5.5577e-07, 3.47006059924, 6279.55273164 ],
    [ 4.1036e-07, 5.36817592855, 8429.24126647 ],
    [ 5.1605e-07, 1.33282739866, 1748.01641307 ],
    [ 5.1992e-07, 0.18914947184, 12139.5535091 ],
    [ 4.9e-07, 0.48735014197, 1194.44701022 ],
    [ 3.92e-07, 6.16833020996, 10447.3878396 ],
    [ 3.557e-07, 1.775968892, 6812.76681509 ],
    [ 3.677e-07, 6.04133863162, 10213.2855462 ],
    [ 3.6596e-07, 2.56957481827, 1059.38193019 ],
    [ 3.3296e-07, 0.59310278598, 17789.8456198 ],
    [ 3.5954e-07, 1.70875808777, 2352.86615377 ],
    [ 4.0938e-07, 2.39850938714, 19651.0484811 ]
    // 61 terms retained
];

KMG.VSOPTerms.earth_L1 = [
    [ 6283.07584999, 0, 0 ],
    [ 0.00206058863, 2.67823455808, 6283.07584999 ],
    [ 4.303419e-05, 2.63512233481, 12566.1517 ],
    [ 4.25264e-06, 1.59046982018, 3.523118349 ],
    [ 1.09017e-06, 2.96631010675, 1577.34354245 ],
    [ 9.3479e-07, 2.59211109542, 18849.22755 ],
    [ 1.19305e-06, 5.79555765566, 26.2983197998 ],
    [ 7.2121e-07, 1.13840581212, 529.690965095 ],
    [ 6.7784e-07, 1.87453300345, 398.149003408 ],
    [ 6.735e-07, 4.40932832004, 5507.55323867 ],
    [ 5.9045e-07, 2.88815790631, 5223.6939198 ],
    [ 5.5976e-07, 2.17471740035, 155.420399434 ],
    [ 4.5411e-07, 0.39799502896, 796.298006816 ],
    [ 3.6298e-07, 0.46875437227, 775.522611324 ],
    [ 2.8962e-07, 2.64732254645, 7.1135470008 ],
    [ 1.9097e-07, 1.84628376049, 5486.77784318 ],
    [ 2.0844e-07, 5.34138275149, 0.9803210682 ],
    [ 1.8508e-07, 4.96855179468, 213.299095438 ],
    [ 1.6233e-07, 0.03216587315, 2544.31441988 ],
    [ 1.7293e-07, 2.9911676063, 6275.96230299 ],
    [ 1.5832e-07, 1.43049301283, 2146.16541648 ],
    [ 1.4608e-07, 1.2046979369, 10977.0788047 ],
    [ 1.1877e-07, 3.25805082007, 5088.62883977 ],
    [ 1.1514e-07, 2.07502080082, 4694.00295471 ],
    [ 9.721e-08, 4.2392586526, 1349.86740966 ],
    [ 9.969e-08, 1.30263423409, 6286.59896834 ],
    [ 9.452e-08, 2.69956827011, 242.728603974 ],
    [ 1.2461e-07, 2.83432282119, 1748.01641307 ],
    [ 1.1808e-07, 5.27379760438, 1194.44701022 ],
    [ 8.577e-08, 5.6447608598, 951.718406251 ],
    [ 1.0641e-07, 0.76614722966, 553.569402842 ],
    [ 7.576e-08, 5.30056172859, 2352.86615377 ],
    [ 5.764e-08, 1.77228445837, 1059.38193019 ],
    [ 6.385e-08, 2.65034514038, 9437.76293489 ],
    [ 5.223e-08, 5.66135782131, 71430.6956181 ],
    [ 5.315e-08, 0.91110018969, 3154.6870849 ],
    [ 6.101e-08, 4.66633726278, 4690.47983636 ],
    [ 4.335e-08, 0.23934560382, 6812.76681509 ],
    [ 5.041e-08, 1.42489704722, 6438.49624943 ],
    [ 4.259e-08, 0.77355543889, 10447.3878396 ],
    [ 5.2e-08, 1.85528830215, 801.820931124 ]
    // 41 terms retained
];

KMG.VSOPTerms.earth_L2 = [
    [ 8.721859e-05, 1.07253635559, 6283.07584999 ],
    [ 9.9099e-06, 3.14159265359, 0 ],
    [ 2.94833e-06, 0.43717350256, 12566.1517 ],
    [ 2.7338e-07, 0.05295636147, 3.523118349 ],
    [ 1.6333e-07, 5.18820215724, 26.2983197998 ],
    [ 1.5745e-07, 3.68504712183, 155.420399434 ],
    [ 9.425e-08, 0.29667114694, 18849.22755 ],
    [ 8.938e-08, 2.05706319592, 77713.7714681 ],
    [ 6.94e-08, 0.82691541038, 775.522611324 ],
    [ 5.061e-08, 4.6624323168, 1577.34354245 ],
    [ 4.06e-08, 1.03067032318, 7.1135470008 ],
    [ 3.464e-08, 5.14021224609, 796.298006816 ],
    [ 3.172e-08, 6.05479318507, 5507.55323867 ],
    [ 3.02e-08, 1.19240008524, 242.728603974 ],
    [ 2.885e-08, 6.11705865396, 529.690965095 ],
    [ 3.809e-08, 3.44043369494, 5573.14280143 ],
    [ 2.719e-08, 0.30363248164, 398.149003408 ],
    [ 2.365e-08, 4.37666117992, 5223.6939198 ],
    [ 2.538e-08, 2.27966434314, 553.569402842 ],
    [ 2.078e-08, 3.75435095487, 0.9803210682 ],
    [ 1.675e-08, 0.90149951436, 951.718406251 ],
    [ 1.534e-08, 5.75895831192, 1349.86740966 ],
    [ 1.224e-08, 2.97285792195, 2146.16541648 ],
    [ 1.449e-08, 4.36401639552, 1748.01641307 ],
    [ 1.341e-08, 3.72019379666, 1194.44701022 ],
    [ 1.253e-08, 2.9488872631, 6438.49624943 ],
    [ 9.99e-09, 5.98665341008, 6286.59896834 ]
    // 27 terms retained
];

KMG.VSOPTerms.earth_L3 = [
    [ 2.89058e-06, 5.84173149732, 6283.07584999 ],
    [ 2.0712e-07, 6.0498393902, 12566.1517 ],
    [ 2.962e-08, 5.1956057957, 155.420399434 ],
    [ 2.527e-08, 3.14159265359, 0 ],
    [ 1.288e-08, 4.7219761197, 3.523118349 ]
    // 5 terms retained
];

KMG.VSOPTerms.earth_L4 = [
    [ 7.714e-08, 4.14117321449, 6283.07584999 ]
    // 1 terms retained
];

KMG.VSOPTerms.earth_L5 = [
    [ 0, 0, 0 ]
    // 0 terms retained
];

KMG.VSOPTerms.earth_B0 = [
    [ 2.7962e-06, 3.19870156017, 84334.6615813 ]
    // 1 terms retained
];

KMG.VSOPTerms.earth_B1 = [
    [ 0.00227777722, 3.4137662053, 6283.07584999 ],
    [ 3.805678e-05, 3.37063423795, 12566.1517 ],
    [ 3.619589e-05, 0, 0 ],
    [ 7.1542e-07, 3.32777549735, 18849.22755 ]
    // 4 terms retained
];

KMG.VSOPTerms.earth_B2 = [
    [ 9.721424e-05, 5.1519280992, 6283.07584999 ],
    [ 2.33002e-06, 3.14159265359, 0 ],
    [ 1.34188e-06, 0.64406212977, 12566.1517 ],
    [ 6.504e-08, 1.07333397797, 18849.22755 ]
    // 4 terms retained
];


KMG.VSOPTerms.earth_R0 = [
    [ 1.00013988784, 0, 0 ],
    [ 0.01670699632, 3.09846350258, 6283.07584999 ],
    [ 0.00013956024, 3.05524609456, 12566.1517 ],
    [ 3.08372e-05, 5.19846674381, 77713.7714681 ],
    [ 1.628463e-05, 1.17387558054, 5753.3848849 ],
    [ 1.575572e-05, 2.84685214877, 7860.41939244 ],
    [ 9.24799e-06, 5.45292236722, 11506.7697698 ],
    [ 5.42439e-06, 4.56409151453, 3930.20969622 ],
    [ 4.7211e-06, 3.66100022149, 5884.92684658 ],
    [ 3.2878e-06, 5.89983686142, 5223.6939198 ],
    [ 3.45969e-06, 0.96368627272, 5507.55323867 ],
    [ 3.06784e-06, 0.29867139512, 5573.14280143 ],
    [ 1.74844e-06, 3.01193636733, 18849.22755 ],
    [ 2.43181e-06, 4.2734953079, 11790.6290887 ],
    [ 2.11836e-06, 5.84714461348, 1577.34354245 ],
    [ 1.8574e-06, 5.02199710705, 10977.0788047 ],
    [ 1.09835e-06, 5.0551063586, 5486.77784318 ],
    [ 9.8316e-07, 0.88681311278, 6069.77675455 ],
    [ 8.65e-07, 5.68956418946, 15720.8387849 ],
    [ 8.5831e-07, 1.27079125277, 161000.685738 ],
    [ 6.2917e-07, 0.92177053978, 529.690965095 ],
    [ 5.7056e-07, 2.01374292245, 83996.8473181 ],
    [ 6.4908e-07, 0.27251341435, 17260.1546547 ],
    [ 4.9384e-07, 3.24501240359, 2544.31441988 ],
    [ 5.5736e-07, 5.2415979917, 71430.6956181 ],
    [ 4.252e-07, 6.01110257982, 6275.96230299 ],
    [ 4.6966e-07, 2.57799853213, 775.522611324 ],
    [ 3.8963e-07, 5.36063832897, 4694.00295471 ],
    [ 4.4666e-07, 5.53715663816, 9437.76293489 ],
    [ 3.5661e-07, 1.67447135798, 12036.4607349 ],
    [ 3.1922e-07, 0.18368299942, 5088.62883977 ],
    [ 3.1846e-07, 1.77775642078, 398.149003408 ],
    [ 3.3193e-07, 0.24370221704, 7084.89678112 ],
    [ 3.8245e-07, 2.39255343973, 8827.39026987 ],
    [ 2.8468e-07, 1.21344887533, 6286.59896834 ],
    [ 3.7486e-07, 0.82961281844, 19651.0484811 ],
    [ 3.6957e-07, 4.90107587287, 12139.5535091 ],
    [ 3.4537e-07, 1.84270693281, 2942.46342329 ],
    [ 2.6275e-07, 4.58896863104, 10447.3878396 ],
    [ 2.4596e-07, 3.78660838036, 8429.24126647 ],
    [ 2.3587e-07, 0.26866098169, 796.298006816 ],
    [ 2.7795e-07, 1.89934427832, 6279.55273164 ],
    [ 2.3927e-07, 4.99598548145, 5856.47765912 ],
    [ 2.0345e-07, 4.65282190725, 2146.16541648 ],
    [ 2.3287e-07, 2.80783632869, 14143.4952424 ],
    [ 2.2099e-07, 1.95002636847, 3154.6870849 ],
    [ 1.9509e-07, 5.38233922479, 2352.86615377 ],
    [ 1.7958e-07, 0.1987136996, 6812.76681509 ],
    [ 1.7178e-07, 4.43322156854, 10213.2855462 ],
    [ 1.619e-07, 5.23159323213, 17789.8456198 ],
    [ 1.7315e-07, 6.15224075188, 16730.4636896 ],
    [ 1.3814e-07, 5.18962074032, 8031.09226306 ],
    [ 1.8834e-07, 0.67280058021, 149854.400135 ],
    [ 1.833e-07, 2.25348717053, 23581.2581773 ],
    [ 1.3639e-07, 3.68511810757, 4705.73230754 ],
    [ 1.3142e-07, 0.65267698994, 13367.9726311 ],
    [ 1.0414e-07, 4.33285688501, 11769.8536932 ],
    [ 9.978e-08, 4.20126336356, 6309.37416979 ],
    [ 1.017e-07, 1.59366684542, 4690.47983636 ],
    [ 7.564e-08, 2.62560597391, 6256.77753019 ],
    [ 9.654e-08, 3.67583728703, 27511.4678735 ],
    [ 6.743e-08, 0.56269927047, 3340.6124267 ],
    [ 8.743e-08, 6.06359123461, 1748.01641307 ],
    [ 7.786e-08, 3.67371235367, 12168.0026966 ],
    [ 6.633e-08, 5.66149277789, 11371.7046898 ],
    [ 7.712e-08, 0.31242577788, 7632.94325965 ],
    [ 6.586e-08, 3.13580054586, 801.820931124 ],
    [ 7.46e-08, 5.6475806666, 11926.2544137 ],
    [ 6.933e-08, 2.92384586372, 6681.2248534 ],
    [ 6.805e-08, 1.42327153767, 23013.5395396 ],
    [ 6.118e-08, 5.13395999022, 1194.44701022 ],
    [ 6.477e-08, 2.64986648493, 19804.8272916 ]
    // 72 terms retained
];

KMG.VSOPTerms.earth_R1 = [
    [ 0.00103018607, 1.10748968172, 6283.07584999 ],
    [ 1.721238e-05, 1.06442300386, 12566.1517 ],
    [ 7.02217e-06, 3.14159265359, 0 ],
    [ 3.2345e-07, 1.02168583254, 18849.22755 ],
    [ 3.0801e-07, 2.84358443952, 5507.55323867 ],
    [ 2.4978e-07, 1.31906570344, 5223.6939198 ],
    [ 1.8487e-07, 1.42428709076, 1577.34354245 ],
    [ 1.0077e-07, 5.91385248388, 10977.0788047 ],
    [ 8.635e-08, 0.27158192945, 5486.77784318 ],
    [ 8.654e-08, 1.42046854427, 6275.96230299 ]
    // 10 terms retained
];

KMG.VSOPTerms.earth_R2 = [
    [ 4.359385e-05, 5.78455133808, 6283.07584999 ],
    [ 1.23633e-06, 5.57935427994, 12566.1517 ],
    [ 1.2342e-07, 3.14159265359, 0 ],
    [ 8.792e-08, 3.62777893099, 77713.7714681 ],
    [ 5.689e-08, 1.86958905084, 5573.14280143 ],
    [ 3.302e-08, 5.47034879713, 18849.22755 ]
    // 6 terms retained
];

KMG.VSOPTerms.earth_R3 = [
    [ 1.44595e-06, 4.27319433901, 6283.07584999 ],
    [ 6.729e-08, 3.91706261708, 12566.1517 ]
    // 2 terms retained
];

KMG.VSOPTerms.earth_R4 = [
    [ 3.858e-08, 2.56389016346, 6283.07584999 ]
    // 1 terms retained
];

KMG.VSOPTerms.earth_R5 = [
    [ 0, 0, 0 ]
    // 0 terms retained
];


KMG.VSOPTerms.mars_L0 = [
    [ 6.20347711581, 0, 0 ],
    [ 0.18656368093, 5.0503710027, 3340.6124267 ],
    [ 0.01108216816, 5.40099836344, 6681.2248534 ],
    [ 0.00091798406, 5.75478744667, 10021.8372801 ],
    [ 0.00027744987, 5.97049513147, 3.523118349 ],
    [ 0.00010610235, 2.93958560338, 2281.23049651 ],
    [ 0.00012315897, 0.84956094002, 2810.92146161 ],
    [ 8.926784e-05, 4.15697846427, 0.0172536522 ],
    [ 8.715691e-05, 6.11005153139, 13362.4497068 ],
    [ 6.797556e-05, 0.36462229657, 398.149003408 ],
    [ 7.774872e-05, 3.33968761376, 5621.84292321 ],
    [ 3.575078e-05, 1.6618650571, 2544.31441988 ],
    [ 4.161108e-05, 0.22814971327, 2942.46342329 ],
    [ 3.075252e-05, 0.85696614132, 191.448266112 ],
    [ 2.628117e-05, 0.64806124465, 3337.08930835 ],
    [ 2.937546e-05, 6.07893711402, 0.0673103028 ],
    [ 2.389414e-05, 5.03896442664, 796.298006816 ],
    [ 2.579844e-05, 0.02996736156, 3344.13554505 ],
    [ 1.528141e-05, 1.14979301996, 6151.5338883 ],
    [ 1.798806e-05, 0.65634057445, 529.690965095 ],
    [ 1.264357e-05, 3.62275122593, 5092.15195812 ],
    [ 1.286228e-05, 3.06796065034, 2146.16541648 ],
    [ 1.546404e-05, 2.91579701718, 1751.53953142 ],
    [ 1.024902e-05, 3.69334099279, 8962.45534991 ],
    [ 8.91566e-06, 0.18293837498, 16703.0621335 ],
    [ 8.58759e-06, 2.4009381194, 2914.01423582 ],
    [ 8.32715e-06, 2.46418619474, 3340.59517305 ],
    [ 8.3272e-06, 4.49495782139, 3340.62968035 ],
    [ 7.12902e-06, 3.66335473479, 1059.38193019 ],
    [ 7.48723e-06, 3.82248614017, 155.420399434 ],
    [ 7.23861e-06, 0.67497311481, 3738.76143011 ],
    [ 6.35548e-06, 2.92182225127, 8432.76438482 ],
    [ 6.55162e-06, 0.48864064125, 3127.31333126 ],
    [ 5.50474e-06, 3.81001042328, 0.9803210682 ],
    [ 5.5275e-06, 4.47479317037, 1748.01641307 ],
    [ 4.25966e-06, 0.55364317304, 6283.07584999 ],
    [ 4.15131e-06, 0.49662285038, 213.299095438 ],
    [ 4.72167e-06, 3.62547124025, 1194.44701022 ],
    [ 3.06551e-06, 0.38052848348, 6684.74797175 ],
    [ 3.12141e-06, 0.99853944405, 6677.70173505 ],
    [ 2.93198e-06, 4.22131299634, 20.7753954924 ],
    [ 3.02375e-06, 4.48618007156, 3532.06069281 ],
    [ 2.74027e-06, 0.54222167059, 3340.5451164 ],
    [ 2.81079e-06, 5.88163521788, 1349.86740966 ],
    [ 2.31183e-06, 1.28242156993, 3870.30339179 ],
    [ 2.83602e-06, 5.7688543494, 3149.16416059 ],
    [ 2.36117e-06, 5.75503217933, 3333.4988797 ],
    [ 2.74033e-06, 0.13372524985, 3340.679737 ],
    [ 2.99395e-06, 2.78323740866, 6254.62666252 ],
    [ 2.04162e-06, 2.82133445874, 1221.84856632 ],
    [ 2.38866e-06, 5.37153646326, 4136.91043352 ],
    [ 1.88648e-06, 1.4910406604, 9492.146315 ],
    [ 2.21228e-06, 3.50466812198, 382.896532223 ],
    [ 1.79196e-06, 1.00561962003, 951.718406251 ],
    [ 1.72117e-06, 0.43943649536, 5486.77784318 ],
    [ 1.93118e-06, 3.35716641911, 3.5904286518 ],
    [ 1.44304e-06, 1.41874112114, 135.065080035 ],
    [ 1.60016e-06, 3.94857092451, 4562.46099302 ],
    [ 1.74072e-06, 2.41361337725, 553.569402842 ],
    [ 1.30989e-06, 4.04491134956, 12303.0677766 ],
    [ 1.38243e-06, 4.30145122848, 7.1135470008 ],
    [ 1.28062e-06, 1.8066581622, 5088.62883977 ],
    [ 1.39898e-06, 3.32595559208, 2700.71514039 ],
    [ 1.28105e-06, 2.20807538189, 1592.59601363 ],
    [ 1.16944e-06, 3.12806863456, 7903.07341972 ],
    [ 1.10378e-06, 1.05194545948, 242.728603974 ],
    [ 1.13481e-06, 3.70070432339, 1589.07289528 ],
    [ 1.00099e-06, 3.24340223714, 11773.3768115 ],
    [ 9.5594e-07, 0.53950648295, 20043.6745602 ],
    [ 9.8947e-07, 4.84558326403, 6681.24210705 ],
    [ 1.04542e-06, 0.78532737699, 8827.39026987 ],
    [ 8.4186e-07, 3.98971116025, 4399.99435689 ],
    [ 8.6928e-07, 2.20183965407, 11243.6858464 ],
    [ 7.1438e-07, 2.80307223477, 3185.19202727 ],
    [ 7.2095e-07, 5.84669532401, 5884.92684658 ],
    [ 7.3482e-07, 2.18421190324, 8429.24126647 ],
    [ 9.8946e-07, 2.81481171439, 6681.20759975 ],
    [ 6.8413e-07, 2.73834597183, 2288.34404351 ],
    [ 8.6747e-07, 1.02091867465, 7079.37385681 ],
    [ 6.5316e-07, 2.68114882713, 28.4491874678 ],
    [ 8.3745e-07, 3.20254912006, 4690.47983636 ],
    [ 7.5031e-07, 0.76647765061, 6467.92575796 ],
    [ 6.8983e-07, 3.76403440528, 6041.32756709 ],
    [ 6.6706e-07, 0.73630288873, 3723.50895892 ],
    [ 6.3313e-07, 4.5277185022, 426.598190876 ],
    [ 6.1684e-07, 6.16831461502, 2274.11694951 ],
    [ 5.226e-07, 0.89938935091, 9623.68827669 ],
    [ 5.5485e-07, 4.60622447136, 4292.33083295 ],
    [ 5.1331e-07, 4.14823934301, 3341.59274777 ],
    [ 5.6633e-07, 5.06250402329, 15.252471185 ],
    [ 6.3376e-07, 0.91293637746, 3553.91152214 ],
    [ 4.5822e-07, 0.78790300125, 1990.74501704 ],
    [ 4.8553e-07, 3.95677994023, 4535.05943692 ],
    [ 4.1223e-07, 6.02013764154, 3894.18182954 ],
    [ 4.1941e-07, 3.58309124437, 8031.09226306 ],
    [ 5.6395e-07, 1.68727941626, 6872.67311951 ],
    [ 5.5907e-07, 3.46261441099, 263.083923373 ],
    [ 5.1677e-07, 2.81307639242, 3339.63210563 ],
    [ 4.0669e-07, 3.13838566327, 9595.23908922 ],
    [ 3.8111e-07, 0.73396370751, 10025.3603984 ],
    [ 3.9498e-07, 5.6322574136, 3097.88382273 ],
    [ 4.4175e-07, 3.19530118759, 5628.95647021 ],
    [ 3.6718e-07, 2.63750919104, 692.157601227 ],
    [ 4.5905e-07, 0.28717581576, 5614.72937621 ],
    [ 3.8351e-07, 5.82880639987, 3191.04922957 ],
    [ 3.8198e-07, 2.34832438823, 162.466636132 ],
    [ 3.2561e-07, 0.48401318272, 6681.2921637 ],
    [ 3.7135e-07, 0.68510839331, 2818.03500861 ],
    [ 3.1169e-07, 3.98160436995, 20.3553193988 ],
    [ 3.2561e-07, 0.89250965753, 6681.1575431 ],
    [ 3.7749e-07, 4.15481250779, 2803.8079146 ],
    [ 3.3626e-07, 6.11997987693, 6489.77658729 ],
    [ 2.9007e-07, 2.42707198395, 3319.83703121 ],
    [ 3.8794e-07, 1.35194224244, 10018.3141618 ],
    [ 3.3149e-07, 1.140241952, 5.5229243074 ],
    [ 2.7583e-07, 1.59721760699, 7210.91581849 ],
    [ 2.8699e-07, 5.7204755094, 7477.52286022 ],
    [ 3.4039e-07, 2.59525636978, 11769.8536932 ],
    [ 2.538e-07, 0.52092092633, 10.6366653498 ],
    [ 2.6355e-07, 1.34519007001, 3496.03282613 ],
    [ 2.4555e-07, 4.00321315879, 11371.7046898 ],
    [ 2.5637e-07, 0.24963503109, 522.577418094 ],
    [ 2.7275e-07, 4.55649766071, 3361.38782219 ],
    [ 2.3766e-07, 1.84063759173, 12832.7587417 ],
    [ 2.2814e-07, 3.52628452806, 1648.4467572 ],
    [ 2.2272e-07, 0.72111173236, 266.607041722 ]
    // 126 terms retained
];

KMG.VSOPTerms.mars_L1 = [
    [ 3340.61242701, 0, 0 ],
    [ 0.01457554523, 3.60433733236, 3340.6124267 ],
    [ 0.00168414711, 3.92318567804, 6681.2248534 ],
    [ 0.00020622975, 4.26108844583, 10021.8372801 ],
    [ 3.452392e-05, 4.7321039319, 3.523118349 ],
    [ 2.586332e-05, 4.60670058555, 13362.4497068 ],
    [ 8.41535e-06, 4.45864030426, 2281.23049651 ],
    [ 5.37567e-06, 5.01581256923, 398.149003408 ],
    [ 5.20948e-06, 4.99428054039, 3344.13554505 ],
    [ 4.32635e-06, 2.56070853083, 191.448266112 ],
    [ 4.29655e-06, 5.31645299471, 155.420399434 ],
    [ 3.81751e-06, 3.53878166043, 796.298006816 ],
    [ 3.2853e-06, 4.95632685192, 16703.0621335 ],
    [ 2.82795e-06, 3.15966768785, 2544.31441988 ],
    [ 2.05657e-06, 4.56889279932, 2146.16541648 ],
    [ 1.68866e-06, 1.3293655906, 3337.08930835 ],
    [ 1.57593e-06, 4.18519540728, 1751.53953142 ],
    [ 1.33686e-06, 2.23327245555, 0.9803210682 ],
    [ 1.16965e-06, 2.21414273762, 1059.38193019 ],
    [ 1.17503e-06, 6.02411290806, 6151.5338883 ],
    [ 1.13718e-06, 5.42753341019, 3738.76143011 ],
    [ 1.33565e-06, 5.97420357518, 1748.01641307 ],
    [ 9.1099e-07, 1.09626613064, 1349.86740966 ],
    [ 8.4256e-07, 5.29330740437, 6684.74797175 ],
    [ 1.13886e-06, 2.12863726524, 1194.44701022 ],
    [ 8.0823e-07, 4.42818326716, 529.690965095 ],
    [ 7.9847e-07, 2.24822372859, 8962.45534991 ],
    [ 7.2505e-07, 5.84203374239, 242.728603974 ],
    [ 7.2945e-07, 2.50193599662, 951.718406251 ],
    [ 7.149e-07, 3.85645759558, 2914.01423582 ],
    [ 8.534e-07, 3.90856932983, 553.569402842 ],
    [ 6.758e-07, 5.0233489507, 382.896532223 ],
    [ 6.506e-07, 1.01810963274, 3340.59517305 ],
    [ 6.5061e-07, 3.04888114328, 3340.62968035 ],
    [ 6.1478e-07, 4.15185188249, 3149.16416059 ],
    [ 4.8482e-07, 4.87339233007, 213.299095438 ],
    [ 4.6581e-07, 1.31461442691, 3185.19202727 ],
    [ 5.6642e-07, 3.88772102421, 4136.91043352 ],
    [ 4.7615e-07, 1.18228660215, 3333.4988797 ],
    [ 4.2052e-07, 5.30826745759, 20043.6745602 ],
    [ 4.133e-07, 0.71392238704, 1592.59601363 ],
    [ 4.028e-07, 2.72571311592, 7.1135470008 ],
    [ 3.304e-07, 5.40823104809, 6283.07584999 ],
    [ 2.8676e-07, 0.04305323493, 9492.146315 ],
    [ 2.2322e-07, 5.86718681699, 3870.30339179 ],
    [ 2.2432e-07, 5.46596961275, 20.3553193988 ],
    [ 2.2606e-07, 0.83782540818, 3097.88382273 ],
    [ 2.1416e-07, 5.37936489667, 3340.5451164 ],
    [ 2.3347e-07, 6.167744339, 3532.06069281 ],
    [ 2.6573e-07, 3.8900063113, 1221.84856632 ],
    [ 2.28e-07, 1.54501542908, 2274.11694951 ],
    [ 2.0474e-07, 2.3623686167, 1589.07289528 ],
    [ 2.0179e-07, 3.36390759347, 5088.62883977 ],
    [ 2.0013e-07, 2.57546546037, 12303.0677766 ],
    [ 1.992e-07, 0.44761063096, 6677.70173505 ],
    [ 2.655e-07, 5.11303525089, 2700.71514039 ],
    [ 2.1104e-07, 3.52541056271, 15.252471185 ],
    [ 2.1424e-07, 4.97083417225, 3340.679737 ],
    [ 1.8502e-07, 5.57854926842, 1990.74501704 ],
    [ 1.7805e-07, 6.12513609945, 4292.33083295 ],
    [ 1.6463e-07, 2.60307709195, 3341.59274777 ],
    [ 1.6592e-07, 1.25515357212, 3894.18182954 ],
    [ 1.9864e-07, 2.52765519587, 4399.99435689 ],
    [ 1.5002e-07, 1.03518790208, 2288.34404351 ],
    [ 2.0011e-07, 4.73112374598, 4690.47983636 ],
    [ 1.5431e-07, 2.46932776517, 4535.05943692 ],
    [ 2.0193e-07, 5.78561467842, 7079.37385681 ],
    [ 1.5298e-07, 2.26504738206, 3723.50895892 ],
    [ 1.5019e-07, 3.36690751539, 6681.24210705 ],
    [ 1.3219e-07, 5.61412860968, 10025.3603984 ],
    [ 1.3517e-07, 2.12392880454, 5486.77784318 ],
    [ 1.5019e-07, 1.33613594479, 6681.20759975 ],
    [ 1.2676e-07, 2.95036175206, 3496.03282613 ],
    [ 1.3644e-07, 1.97710249337, 5614.72937621 ],
    [ 1.3011e-07, 1.51458564766, 5628.95647021 ],
    [ 1.1353e-07, 6.23411904718, 135.065080035 ],
    [ 1.3508e-07, 3.42721826602, 5621.84292321 ],
    [ 1.0866e-07, 5.28165480979, 2818.03500861 ],
    [ 1.188e-07, 3.12847055823, 426.598190876 ],
    [ 1.0467e-07, 2.7359860705, 2787.04302386 ],
    [ 1.1131e-07, 5.84122566289, 2803.8079146 ],
    [ 1.177e-07, 2.58277425311, 8432.76438482 ],
    [ 1.1861e-07, 5.47552055459, 3553.91152214 ],
    [ 8.54e-08, 1.91739325491, 11773.3768115 ],
    [ 9.819e-08, 4.52958330672, 6489.77658729 ],
    [ 8.552e-08, 3.16147568714, 162.466636132 ],
    [ 1.0957e-07, 4.15775327007, 2388.89402045 ],
    [ 8.948e-08, 4.23164385777, 7477.52286022 ],
    [ 8.131e-08, 1.61308074119, 2957.71589448 ],
    [ 8.352e-08, 2.18475645206, 23.8784377478 ],
    [ 8.03e-08, 5.69889507906, 6041.32756709 ],
    [ 7.878e-08, 5.71359767892, 9623.68827669 ],
    [ 8.713e-08, 4.43300582398, 5092.15195812 ],
    [ 8.421e-08, 3.1635506725, 3347.7259737 ],
    [ 6.67e-08, 5.07423317095, 8031.09226306 ],
    [ 8.656e-08, 4.33239148117, 3339.63210563 ],
    [ 7.354e-08, 6.17934256606, 3583.34103067 ],
    [ 5.749e-08, 3.67719823582, 8429.24126647 ],
    [ 6.235e-08, 3.54003325209, 692.157601227 ],
    [ 5.458e-08, 1.05139431657, 4933.20844033 ],
    [ 6.132e-08, 1.66182646558, 6525.80445397 ],
    [ 5.197e-08, 1.14841109166, 28.4491874678 ],
    [ 4.95e-08, 5.28919125231, 6681.2921637 ],
    [ 5.516e-08, 6.12492946392, 2487.41604495 ],
    [ 4.89e-08, 3.10255139433, 5.5229243074 ],
    [ 5.354e-08, 0.37154896863, 12832.7587417 ],
    [ 4.751e-08, 0.2337468155, 36.0278666774 ],
    [ 6.362e-08, 2.11339432269, 5884.92684658 ],
    [ 4.996e-08, 2.44835744792, 5099.26550512 ],
    [ 4.952e-08, 5.69770765577, 6681.1575431 ],
    [ 4.678e-08, 0.27799012787, 10018.3141618 ],
    [ 4.746e-08, 0.00950199989, 7210.91581849 ],
    [ 4.862e-08, 5.60331599025, 6467.92575796 ],
    [ 5.544e-08, 2.00929051393, 522.577418094 ],
    [ 4.998e-08, 1.51094959078, 1744.42598442 ],
    [ 5.397e-08, 0.1884215497, 2942.46342329 ],
    [ 4.098e-08, 3.95776844736, 3.881335358 ],
    [ 5.414e-08, 5.66147396313, 23384.2869869 ],
    [ 5.467e-08, 0.19258681316, 7632.94325965 ],
    [ 4.305e-08, 2.8945229483, 2810.92146161 ],
    [ 4.118e-08, 1.59475420886, 7234.79425624 ],
    [ 4.489e-08, 4.16951490492, 2906.90068882 ],
    [ 5.277e-08, 2.22681020305, 3127.31333126 ],
    [ 3.882e-08, 2.26433789475, 2699.73481932 ],
    [ 3.544e-08, 1.76658498504, 1758.65307842 ],
    [ 3.408e-08, 2.65743533541, 4929.68532198 ],
    [ 4.336e-08, 4.43081904792, 640.877607382 ],
    [ 3.804e-08, 2.91373968131, 15643.6802033 ],
    [ 3.176e-08, 1.75893480952, 9595.23908922 ],
    [ 3.309e-08, 6.13831291678, 10419.9862835 ],
    [ 3.077e-08, 2.56194751001, 7064.12138562 ],
    [ 3.236e-08, 2.32387839882, 5085.03841111 ],
    [ 3.284e-08, 2.8621647971, 7740.60678359 ],
    [ 2.958e-08, 1.27767445188, 574.344798335 ],
    [ 2.805e-08, 0.43144651568, 5828.02847165 ],
    [ 2.851e-08, 0.98625869565, 3191.04922957 ],
    [ 3.324e-08, 2.5901098785, 2118.76386038 ],
    [ 3.039e-08, 1.86739127757, 7.046236698 ],
    [ 2.738e-08, 1.76460911547, 639.897286314 ],
    [ 2.757e-08, 3.70511041849, 10021.8545338 ],
    [ 3.376e-08, 1.53123149565, 6674.1113064 ],
    [ 2.757e-08, 1.67433972403, 10021.8200264 ],
    [ 2.67e-08, 3.11556212899, 6836.64525283 ],
    [ 2.583e-08, 3.77302627584, 2921.12778282 ],
    [ 2.51e-08, 0.30461555756, 3475.67750674 ],
    [ 2.288e-08, 2.81266012379, 7875.67186362 ],
    [ 2.411e-08, 0.97123911611, 3319.83703121 ],
    [ 2.41e-08, 2.95969382172, 6682.20517447 ],
    [ 2.211e-08, 0.61268074323, 10973.5556864 ],
    [ 2.246e-08, 4.12573972297, 59.3738619136 ],
    [ 2.183e-08, 2.17530786579, 15113.9892382 ],
    [ 2.445e-08, 5.91435376684, 5331.35744374 ]
    // 152 terms retained
];

KMG.VSOPTerms.mars_L2 = [
    [ 0.00058152577, 2.04961712429, 3340.6124267 ],
    [ 0.00013459579, 2.45738706163, 6681.2248534 ],
    [ 2.432575e-05, 2.79737979284, 10021.8372801 ],
    [ 4.01065e-06, 3.13581149963, 13362.4497068 ],
    [ 4.51384e-06, 0, 0 ],
    [ 2.22025e-06, 3.19437046607, 3.523118349 ],
    [ 1.20954e-06, 0.54327128607, 155.420399434 ],
    [ 6.2971e-07, 3.47765178989, 16703.0621335 ],
    [ 5.3644e-07, 3.54171478781, 3344.13554505 ],
    [ 3.4273e-07, 6.00208464365, 2281.23049651 ],
    [ 3.1659e-07, 4.14001980084, 191.448266112 ],
    [ 2.9839e-07, 1.9983873938, 796.298006816 ],
    [ 2.3172e-07, 4.33401932281, 242.728603974 ],
    [ 2.1663e-07, 3.44500841809, 398.149003408 ],
    [ 1.605e-07, 6.11000263211, 2146.16541648 ],
    [ 2.0369e-07, 5.42202383442, 553.569402842 ],
    [ 1.4924e-07, 6.09549588012, 3185.19202727 ],
    [ 1.6229e-07, 0.65685105422, 0.9803210682 ],
    [ 1.4317e-07, 2.61898820749, 1349.86740966 ],
    [ 1.4411e-07, 4.01941740099, 951.718406251 ],
    [ 1.1944e-07, 3.86196758615, 6684.74797175 ],
    [ 1.5655e-07, 1.22093822826, 1748.01641307 ],
    [ 1.1261e-07, 4.71857181601, 2544.31441988 ],
    [ 1.336e-07, 0.60151621438, 1194.44701022 ],
    [ 1.0395e-07, 0.25075540193, 382.896532223 ],
    [ 9.415e-08, 0.68050215057, 1059.38193019 ],
    [ 9.58e-08, 3.82256319681, 20043.6745602 ],
    [ 8.996e-08, 3.88272784458, 3738.76143011 ],
    [ 7.498e-08, 5.46428174266, 1751.53953142 ],
    [ 6.499e-08, 5.47802397833, 1592.59601363 ],
    [ 6.307e-08, 2.34134269478, 3097.88382273 ],
    [ 6.864e-08, 2.57523762859, 3149.16416059 ],
    [ 5.871e-08, 1.1486285578, 7.1135470008 ],
    [ 6.675e-08, 2.37862627319, 4136.91043352 ],
    [ 4.655e-08, 4.4310225149, 6151.5338883 ],
    [ 4.201e-08, 3.68638044545, 5614.72937621 ],
    [ 4.796e-08, 2.89378142432, 3333.4988797 ],
    [ 4.074e-08, 6.12610105396, 5628.95647021 ],
    [ 3.66e-08, 4.06581319964, 1990.74501704 ],
    [ 3.284e-08, 2.79214099721, 3894.18182954 ],
    [ 3.615e-08, 2.46526861067, 529.690965095 ],
    [ 3.214e-08, 0.68469193035, 8962.45534991 ],
    [ 3.087e-08, 4.56932030502, 3496.03282613 ],
    [ 2.918e-08, 5.41494777349, 2914.01423582 ],
    [ 2.925e-08, 1.23098223044, 2787.04302386 ],
    [ 2.808e-08, 1.38431632694, 4292.33083295 ],
    [ 2.652e-08, 1.05282528913, 3341.59274777 ],
    [ 2.92e-08, 3.41297158184, 3337.08930835 ],
    [ 2.423e-08, 0.9648433024, 4535.05943692 ],
    [ 2.311e-08, 4.84742235872, 9492.146315 ],
    [ 2.597e-08, 5.74792546254, 3340.59517305 ],
    [ 2.19e-08, 3.26596280325, 213.299095438 ],
    [ 2.598e-08, 1.49506860128, 3340.62968035 ],
    [ 2.365e-08, 4.1830384242, 10025.3603984 ],
    [ 2.63e-08, 4.67732434694, 3583.34103067 ],
    [ 2.606e-08, 2.64976204169, 2388.89402045 ],
    [ 1.822e-08, 0.97105743952, 1589.07289528 ],
    [ 2.397e-08, 1.04493547179, 4399.99435689 ],
    [ 2.203e-08, 0.16281603659, 6525.80445397 ],
    [ 2.373e-08, 4.26885534124, 7079.37385681 ],
    [ 2.366e-08, 0.00564620006, 4690.47983636 ],
    [ 1.623e-08, 4.95374152644, 5088.62883977 ],
    [ 2.143e-08, 0.47993241836, 2700.71514039 ],
    [ 1.646e-08, 4.94105214632, 1221.84856632 ],
    [ 1.588e-08, 1.11405721408, 12303.0677766 ],
    [ 1.518e-08, 0.11076145171, 2957.71589448 ],
    [ 1.774e-08, 3.80344931471, 3723.50895892 ],
    [ 1.364e-08, 3.86744855408, 6283.07584999 ],
    [ 1.764e-08, 2.51992889432, 2810.92146161 ],
    [ 1.394e-08, 2.7360876673, 7477.52286022 ],
    [ 1.28e-08, 5.47285286548, 6677.70173505 ],
    [ 1.447e-08, 2.97506973239, 6489.77658729 ],
    [ 1.248e-08, 3.77100223369, 2699.73481932 ],
    [ 1.527e-08, 2.92629955117, 640.877607382 ],
    [ 1.197e-08, 1.89205359446, 6681.24210705 ],
    [ 1.418e-08, 1.54599865534, 3347.7259737 ],
    [ 1.423e-08, 4.17063094406, 23384.2869869 ],
    [ 1.042e-08, 5.83283345776, 4933.20844033 ],
    [ 1.196e-08, 6.14479114175, 6681.20759975 ],
    [ 1.153e-08, 1.50265359557, 426.598190876 ],
    [ 1.099e-08, 3.80358943061, 3870.30339179 ],
    [ 9.09e-09, 3.81838122072, 5092.15195812 ],
    [ 1.071e-08, 5.04949161471, 5621.84292321 ],
    [ 8.46e-09, 3.82219998207, 3340.5451164 ],
    [ 1.075e-08, 3.81844135104, 3553.91152214 ],
    [ 8.56e-09, 3.42045045625, 3340.679737 ],
    [ 9.16e-09, 1.91472787569, 3532.06069281 ],
    [ 7.14e-09, 4.26169501052, 9623.68827669 ],
    [ 9.07e-09, 4.12943952579, 162.466636132 ],
    [ 6.53e-09, 3.10816357251, 7234.79425624 ],
    [ 7.92e-09, 5.20659969594, 87.3082045398 ],
    [ 6.54e-09, 1.57331630734, 2487.41604495 ],
    [ 6.49e-09, 2.78892909992, 574.344798335 ],
    [ 6.48e-09, 5.181110771, 12832.7587417 ],
    [ 7.07e-09, 5.8319586932, 3339.63210563 ],
    [ 5.2e-09, 4.64660657418, 6836.64525283 ],
    [ 6.6e-09, 0.24998045706, 8969.56889691 ],
    [ 6.4e-09, 1.70935421799, 7632.94325965 ],
    [ 5.28e-09, 0.3110540935, 8031.09226306 ],
    [ 5.1e-09, 4.63676288319, 10419.9862835 ],
    [ 6.04e-09, 3.85002715377, 5486.77784318 ],
    [ 5.14e-09, 1.38796992796, 7740.60678359 ]
    // 102 terms retained
];

KMG.VSOPTerms.mars_L3 = [
    [ 1.467867e-05, 0.4442983946, 3340.6124267 ],
    [ 6.92668e-06, 0.88679887123, 6681.2248534 ],
    [ 1.89478e-06, 1.28336839921, 10021.8372801 ],
    [ 4.1615e-07, 1.64210455567, 13362.4497068 ],
    [ 2.266e-07, 2.05278956965, 155.420399434 ],
    [ 8.126e-08, 1.99049724299, 16703.0621335 ],
    [ 1.0455e-07, 1.57992093693, 3.523118349 ],
    [ 4.902e-08, 2.8251687501, 242.728603974 ],
    [ 5.379e-08, 3.14159265359, 0 ],
    [ 3.782e-08, 2.01848153986, 3344.13554505 ],
    [ 3.181e-08, 4.59108786647, 3185.19202727 ],
    [ 3.133e-08, 0.65141319517, 553.569402842 ],
    [ 1.698e-08, 5.53803382831, 951.718406251 ],
    [ 1.525e-08, 5.71698515888, 191.448266112 ],
    [ 1.451e-08, 0.4606849022, 796.298006816 ],
    [ 1.473e-08, 2.33727441522, 20043.6745602 ],
    [ 1.314e-08, 5.36403056955, 0.9803210682 ],
    [ 1.178e-08, 4.14644990348, 1349.86740966 ],
    [ 1.138e-08, 2.37914351932, 6684.74797175 ],
    [ 1.046e-08, 1.76915268602, 382.896532223 ],
    [ 9.02e-09, 5.35475854699, 1194.44701022 ],
    [ 8.13e-09, 2.74852234414, 1748.01641307 ],
    [ 6.29e-09, 6.08292992203, 3496.03282613 ],
    [ 5.64e-09, 1.87914711325, 398.149003408 ],
    [ 5.66e-09, 5.8543921654, 7.1135470008 ],
    [ 6.46e-09, 3.17980126471, 3583.34103067 ]
    // 26 terms retained
];

KMG.VSOPTerms.mars_L4 = [
    [ 2.7242e-07, 5.6399774232, 6681.2248534 ],
    [ 2.5511e-07, 5.13956279086, 3340.6124267 ],
    [ 1.1147e-07, 6.03556608878, 10021.8372801 ],
    [ 3.19e-08, 3.56206901204, 155.420399434 ],
    [ 3.251e-08, 0.1291561646, 13362.4497068 ]
    // 5 terms retained
];

KMG.VSOPTerms.mars_L5 = [
    [ 7.62e-09, 4.03556368806, 6681.2248534 ],
    [ 5.11e-09, 4.4877039364, 10021.8372801 ],
    [ 3.6e-09, 5.07296615717, 155.420399434 ]
    // 3 terms retained
];

KMG.VSOPTerms.mars_B0 = [
    [ 0.03197134986, 3.76832042431, 3340.6124267 ],
    [ 0.00298033234, 4.10616996305, 6681.2248534 ],
    [ 0.00289104742, 0, 0 ],
    [ 0.00031365539, 4.4465105309, 10021.8372801 ],
    [ 3.4841e-05, 4.7881254926, 13362.4497068 ],
    [ 4.42999e-06, 5.65233014206, 3337.08930835 ],
    [ 4.43401e-06, 5.02642622964, 3344.13554505 ],
    [ 3.99109e-06, 5.13056816928, 16703.0621335 ],
    [ 2.92506e-06, 3.79290674178, 2281.23049651 ],
    [ 1.81982e-06, 6.13648041445, 6151.5338883 ],
    [ 1.63159e-06, 4.26399640691, 529.690965095 ],
    [ 1.59678e-06, 2.23194572851, 1059.38193019 ],
    [ 1.39323e-06, 2.41796458896, 8962.45534991 ],
    [ 1.49297e-06, 2.16501221175, 5621.84292321 ],
    [ 1.42686e-06, 1.18215016908, 3340.59517305 ],
    [ 1.42685e-06, 3.21292181638, 3340.62968035 ],
    [ 8.2544e-07, 5.36667920373, 6684.74797175 ],
    [ 7.3639e-07, 5.0918769577, 398.149003408 ],
    [ 7.266e-07, 5.53775735826, 6283.07584999 ],
    [ 8.6377e-07, 5.74429749104, 3738.76143011 ],
    [ 8.3276e-07, 5.98866355811, 6677.70173505 ],
    [ 6.0116e-07, 3.67960801961, 796.298006816 ],
    [ 6.3111e-07, 0.73049101791, 5884.92684658 ],
    [ 6.2338e-07, 4.8507212869, 2942.46342329 ]
    // 24 terms retained
];

KMG.VSOPTerms.mars_B1 = [
    [ 0.00217310991, 6.04472194776, 3340.6124267 ],
    [ 0.00020976948, 3.14159265359, 0 ],
    [ 0.00012834709, 1.60810667915, 6681.2248534 ],
    [ 3.320981e-05, 2.62947004077, 10021.8372801 ],
    [ 6.272e-06, 3.11898601248, 13362.4497068 ],
    [ 1.0199e-06, 3.52113557592, 16703.0621335 ],
    [ 7.5107e-07, 0.95983758515, 3337.08930835 ],
    [ 2.9264e-07, 3.4030768271, 3344.13554505 ],
    [ 2.3251e-07, 3.69342549027, 5621.84292321 ],
    [ 2.219e-07, 2.21703408598, 2281.23049651 ],
    [ 1.5454e-07, 3.89610159362, 20043.6745602 ],
    [ 1.1867e-07, 3.83861019788, 6684.74797175 ],
    [ 1.2038e-07, 2.13866775328, 6151.5338883 ],
    [ 9.697e-08, 5.48941186798, 3340.62968035 ],
    [ 9.697e-08, 3.45863925102, 3340.59517305 ],
    [ 1.1537e-07, 1.90395033905, 3532.06069281 ],
    [ 9.276e-08, 0.71941312462, 2942.46342329 ],
    [ 9.24e-08, 2.51747952408, 5884.92684658 ],
    [ 9.876e-08, 6.13507416822, 1059.38193019 ],
    [ 9.265e-08, 4.55759125226, 8962.45534991 ],
    [ 7.822e-08, 6.10932267009, 2810.92146161 ],
    [ 1.0369e-07, 0.60195347181, 529.690965095 ],
    [ 8.522e-08, 4.40106741096, 3496.03282613 ],
    [ 7.683e-08, 1.21169696624, 6677.70173505 ],
    [ 7.134e-08, 1.93610705535, 2544.31441988 ],
    [ 6.512e-08, 3.11636422105, 3738.76143011 ],
    [ 6.278e-08, 6.23176923902, 3185.19202727 ],
    [ 5.833e-08, 0.74324094343, 398.149003408 ],
    [ 5.033e-08, 2.28727456802, 3149.16416059 ],
    [ 4.958e-08, 1.54200127913, 6283.07584999 ]
    // 30 terms retained
];

KMG.VSOPTerms.mars_B2 = [
    [ 8.888446e-05, 1.06196052751, 3340.6124267 ],
    [ 2.595393e-05, 3.14159265359, 0 ],
    [ 9.18914e-06, 0.1153843119, 6681.2248534 ],
    [ 2.67883e-06, 0.78837893063, 10021.8372801 ],
    [ 6.6911e-07, 1.39435595847, 13362.4497068 ],
    [ 1.4267e-07, 1.87268116087, 16703.0621335 ],
    [ 7.948e-08, 2.58819177832, 3337.08930835 ],
    [ 2.709e-08, 2.29241371893, 20043.6745602 ],
    [ 2.911e-08, 1.36634316448, 3344.13554505 ],
    [ 2.528e-08, 6.00423798411, 3496.03282613 ],
    [ 1.617e-08, 5.72212771018, 5621.84292321 ],
    [ 1.625e-08, 4.63140305669, 3185.19202727 ]
    // 12 terms retained
];

KMG.VSOPTerms.mars_B3 = [
    [ 3.30418e-06, 2.04215300484, 3340.6124267 ],
    [ 9.3057e-07, 0, 0 ],
    [ 1.4546e-07, 5.38525967237, 10021.8372801 ],
    [ 8.731e-08, 4.90252313032, 6681.2248534 ],
    [ 5.215e-08, 5.97441462813, 13362.4497068 ],
    [ 1.422e-08, 0.21283650226, 16703.0621335 ]
    // 6 terms retained
];

KMG.VSOPTerms.mars_B4 = [
    [ 6.007e-08, 3.37637101191, 3340.6124267 ],
    [ 6.625e-08, 0, 0 ]
    // 2 terms retained
];

KMG.VSOPTerms.mars_B5 = [
    [ 0, 0, 0 ]
    // 0 terms retained
];

KMG.VSOPTerms.mars_R0 = [
    [ 1.53033488271, 0, 0 ],
    [ 0.1418495316, 3.47971283528, 3340.6124267 ],
    [ 0.00660776362, 3.81783443019, 6681.2248534 ],
    [ 0.00046179117, 4.15595316782, 10021.8372801 ],
    [ 8.109733e-05, 5.55958416318, 2810.92146161 ],
    [ 7.485318e-05, 1.77239078402, 5621.84292321 ],
    [ 5.523191e-05, 1.3643630377, 2281.23049651 ],
    [ 3.82516e-05, 4.49407183687, 13362.4497068 ],
    [ 2.306537e-05, 0.09081579001, 2544.31441988 ],
    [ 1.999396e-05, 5.36059617709, 3337.08930835 ],
    [ 2.484394e-05, 4.9254563992, 2942.46342329 ],
    [ 1.960195e-05, 4.74249437639, 3344.13554505 ],
    [ 1.167119e-05, 2.11260868341, 5092.15195812 ],
    [ 1.102816e-05, 5.00908403998, 398.149003408 ],
    [ 8.99066e-06, 4.40791133207, 529.690965095 ],
    [ 9.92252e-06, 5.83861961952, 6151.5338883 ],
    [ 8.07354e-06, 2.10217065501, 1059.38193019 ],
    [ 7.97915e-06, 3.44839203899, 796.298006816 ],
    [ 7.40975e-06, 1.49906336885, 2146.16541648 ],
    [ 6.92339e-06, 2.13378874689, 8962.45534991 ],
    [ 6.33144e-06, 0.89353283242, 3340.59517305 ],
    [ 7.25583e-06, 1.24516810723, 8432.76438482 ],
    [ 6.3314e-06, 2.92430446399, 3340.62968035 ],
    [ 5.74355e-06, 0.82896244455, 2914.01423582 ],
    [ 5.26166e-06, 5.38292991236, 3738.76143011 ],
    [ 6.29978e-06, 1.28737486495, 1751.53953142 ],
    [ 4.72775e-06, 5.19850522346, 3127.31333126 ],
    [ 3.48095e-06, 4.83219199976, 16703.0621335 ],
    [ 2.83713e-06, 2.90692064724, 3532.06069281 ],
    [ 2.79543e-06, 5.2574968538, 6283.07584999 ],
    [ 2.33857e-06, 5.10545987572, 5486.77784318 ],
    [ 2.19427e-06, 5.58340231744, 191.448266112 ],
    [ 2.69896e-06, 3.76393625127, 5884.92684658 ],
    [ 2.08335e-06, 5.25476078693, 3340.5451164 ],
    [ 2.75217e-06, 2.90817482492, 1748.01641307 ],
    [ 2.75506e-06, 1.21767950614, 6254.62666252 ],
    [ 2.39119e-06, 2.03669934656, 1194.44701022 ],
    [ 2.23189e-06, 4.19861535147, 3149.16416059 ],
    [ 1.82689e-06, 5.08062725665, 6684.74797175 ],
    [ 1.86207e-06, 5.6987157241, 6677.70173505 ],
    [ 1.76e-06, 5.95341919657, 3870.30339179 ],
    [ 1.78617e-06, 4.18423004741, 3333.4988797 ],
    [ 2.0833e-06, 4.84626439637, 3340.679737 ],
    [ 2.28126e-06, 3.25526555588, 6872.67311951 ],
    [ 1.44312e-06, 0.2130621946, 5088.62883977 ],
    [ 1.63527e-06, 3.79888811958, 4136.91043352 ],
    [ 1.33126e-06, 1.53906679361, 7903.07341972 ],
    [ 1.41755e-06, 2.47792380112, 4562.46099302 ],
    [ 1.14927e-06, 4.31748869065, 1349.86740966 ],
    [ 1.18789e-06, 2.12168482244, 1589.07289528 ],
    [ 1.02094e-06, 6.18145185708, 9492.146315 ],
    [ 1.2857e-06, 5.49884728795, 8827.39026987 ],
    [ 1.11546e-06, 0.55346108403, 11243.6858464 ],
    [ 8.2498e-07, 1.62220096558, 11773.3768115 ],
    [ 8.3204e-07, 0.61551135046, 8429.24126647 ],
    [ 8.4463e-07, 0.62274409931, 1592.59601363 ],
    [ 8.6666e-07, 1.74984525176, 2700.71514039 ],
    [ 7.1813e-07, 2.4749406548, 12303.0677766 ],
    [ 8.5321e-07, 1.61634750496, 4690.47983636 ],
    [ 6.3641e-07, 2.67334163937, 426.598190876 ],
    [ 6.8601e-07, 2.40188234283, 4399.99435689 ],
    [ 5.8559e-07, 4.7205283999, 213.299095438 ],
    [ 6.2009e-07, 1.10068565926, 1221.84856632 ],
    [ 6.6499e-07, 2.21296335919, 6041.32756709 ],
    [ 5.581e-07, 1.2328806632, 3185.19202727 ],
    [ 5.4969e-07, 5.72695354791, 951.718406251 ],
    [ 5.243e-07, 3.0236809553, 4292.33083295 ],
    [ 5.5688e-07, 5.44688671707, 3723.50895892 ],
    [ 5.8959e-07, 3.26242460622, 6681.24210705 ],
    [ 4.4638e-07, 2.01459444131, 8031.09226306 ],
    [ 5.8959e-07, 1.2316529679, 6681.20759975 ],
    [ 4.2439e-07, 2.26554261514, 155.420399434 ],
    [ 3.8955e-07, 2.57760417339, 3341.59274777 ],
    [ 5.155e-07, 5.72324451485, 7079.37385681 ],
    [ 4.894e-07, 5.61613493545, 3553.91152214 ],
    [ 4.5406e-07, 5.43303278149, 6467.92575796 ],
    [ 3.6438e-07, 4.43922435395, 3894.18182954 ],
    [ 3.598e-07, 1.15972378713, 2288.34404351 ],
    [ 3.5268e-07, 5.49032233898, 1990.74501704 ],
    [ 4.2192e-07, 1.63254827838, 5628.95647021 ],
    [ 4.4292e-07, 5.00344221303, 5614.72937621 ],
    [ 3.3616e-07, 5.17029030468, 20043.6745602 ],
    [ 4.3256e-07, 1.03722397198, 11769.8536932 ],
    [ 3.9237e-07, 1.24237030858, 3339.63210563 ],
    [ 3.1949e-07, 4.59259676953, 2274.11694951 ],
    [ 3.0352e-07, 2.44163963455, 11371.7046898 ],
    [ 3.2269e-07, 2.38222363233, 4535.05943692 ],
    [ 3.1855e-07, 4.37536980289, 3.523118349 ],
    [ 2.9342e-07, 4.06035002188, 3097.88382273 ],
    [ 3.1967e-07, 1.93969979134, 382.896532223 ],
    [ 2.6164e-07, 5.58463559826, 9623.68827669 ],
    [ 2.7903e-07, 4.25809486053, 3191.04922957 ],
    [ 3.3044e-07, 0.85475620169, 553.569402842 ],
    [ 2.7544e-07, 1.5766864517, 9595.23908922 ],
    [ 2.5163e-07, 0.81337734264, 10713.9948813 ],
    [ 2.2045e-07, 0.85711201558, 3319.83703121 ],
    [ 2.4759e-07, 5.38993953923, 2818.03500861 ],
    [ 2.3352e-07, 6.0145897459, 3496.03282613 ],
    [ 2.4723e-07, 2.58025225634, 2803.8079146 ],
    [ 1.9361e-07, 5.18528881954, 6681.2921637 ],
    [ 1.9118e-07, 5.419693554, 10025.3603984 ],
    [ 1.9361e-07, 5.59378511334, 6681.1575431 ],
    [ 1.8331e-07, 5.7956572331, 7064.12138562 ],
    [ 1.8188e-07, 5.61299105522, 7.1135470008 ],
    [ 2.0393e-07, 4.53615443964, 6489.77658729 ],
    [ 2.1258e-07, 6.19174428363, 14054.607308 ],
    [ 1.7094e-07, 1.54988538094, 2957.71589448 ],
    [ 2.2794e-07, 3.41719468533, 7632.94325965 ],
    [ 2.0561e-07, 2.98654120324, 3361.38782219 ],
    [ 1.705e-07, 6.15529583629, 10404.7338123 ],
    [ 1.8007e-07, 2.81505100996, 4032.77002793 ],
    [ 1.6487e-07, 3.84534133372, 10973.5556864 ],
    [ 1.6056e-07, 0.92819026247, 14584.2982731 ],
    [ 2.1008e-07, 2.38506850221, 4989.0591839 ],
    [ 1.6291e-07, 1.92190075688, 7373.38245463 ],
    [ 1.6286e-07, 6.28252184173, 7210.91581849 ],
    [ 1.8575e-07, 4.07319565284, 2388.89402045 ],
    [ 1.5976e-07, 4.58379703739, 3264.34635542 ],
    [ 1.9909e-07, 2.73523951203, 5099.26550512 ],
    [ 1.9667e-07, 1.86294734899, 3443.70520092 ],
    [ 1.65e-07, 4.1406165717, 7477.52286022 ],
    [ 1.9492e-07, 6.03778625701, 10018.3141618 ],
    [ 1.5097e-07, 2.65433832872, 2787.04302386 ],
    [ 1.9099e-07, 0.22623513076, 13745.346239 ],
    [ 1.7164e-07, 3.1882629935, 3347.7259737 ],
    [ 1.3407e-07, 2.12775612449, 3344.20285535 ],
    [ 1.5407e-07, 2.20766468871, 2118.76386038 ],
    [ 1.7246e-07, 3.67064642858, 3205.54734666 ],
    [ 1.3091e-07, 4.27475419816, 14314.168113 ],
    [ 1.6437e-07, 2.86612474805, 14712.3171165 ],
    [ 1.6648e-07, 4.521351492, 6674.1113064 ],
    [ 1.3718e-07, 1.68586111426, 3337.02199805 ],
    [ 1.1824e-07, 0.19675650045, 3475.67750674 ],
    [ 1.1757e-07, 3.23020638064, 5828.02847165 ],
    [ 1.1884e-07, 4.82075035433, 7234.79425624 ],
    [ 1.0608e-07, 1.73995972784, 639.897286314 ],
    [ 1.1143e-07, 0.23833349966, 12832.7587417 ],
    [ 1.1028e-07, 0.4455568729, 10213.2855462 ],
    [ 1.0238e-07, 5.74731032428, 242.728603974 ],
    [ 1.0052e-07, 2.45096419672, 4929.68532198 ],
    [ 1.0061e-07, 0.78904152333, 9381.93999379 ],
    [ 1.0065e-07, 5.37509927353, 5085.03841111 ],
    [ 1.1897e-07, 0.79890074455, 3265.83082813 ],
    [ 8.983e-08, 0.96474320941, 4933.20844033 ],
    [ 8.976e-08, 4.18310051894, 9225.53927328 ],
    [ 8.982e-08, 1.98499607259, 15113.9892382 ],
    [ 8.325e-08, 1.93706224943, 1648.4467572 ],
    [ 7.832e-08, 2.04997038646, 1758.65307842 ],
    [ 7.964e-08, 3.92258783522, 2921.12778282 ],
    [ 1.0223e-07, 2.66509814753, 2487.41604495 ],
    [ 8.277e-08, 0.94860765545, 2906.90068882 ],
    [ 7.371e-08, 0.84436508721, 692.157601227 ],
    [ 7.529e-08, 5.68043313811, 13916.0191096 ],
    [ 7.907e-08, 2.81314645975, 15643.6802033 ],
    [ 6.956e-08, 3.32212696002, 3230.40610548 ],
    [ 7.426e-08, 6.09654676653, 3583.34103067 ],
    [ 6.402e-08, 4.19806999276, 5202.35827934 ],
    [ 6.523e-08, 6.11927838278, 135.065080035 ],
    [ 6.127e-08, 0.00122595969, 6836.64525283 ],
    [ 6.223e-08, 6.1065313699, 17256.6315363 ],
    [ 8.161e-08, 5.24822786208, 10575.4066829 ],
    [ 6.163e-08, 3.60026818309, 10021.8545338 ],
    [ 6.163e-08, 1.56949585888, 10021.8200264 ],
    [ 5.673e-08, 0.13638905291, 13524.9163429 ],
    [ 6.257e-08, 4.50450316951, 8425.65083781 ],
    [ 5.249e-08, 2.70116504868, 4459.3682188 ],
    [ 6.47e-08, 2.74232480124, 7740.60678359 ],
    [ 5.523e-08, 6.06378363783, 10419.9862835 ],
    [ 5.548e-08, 5.75002125481, 12168.0026966 ],
    [ 6.827e-08, 4.69340338938, 17654.7805397 ],
    [ 4.993e-08, 4.68464837021, 522.577418094 ],
    [ 6.32e-08, 3.3193809127, 3767.21061758 ],
    [ 4.735e-08, 0.00770324607, 3325.35995551 ],
    [ 5.025e-08, 2.33675441772, 1052.26838319 ],
    [ 4.656e-08, 5.15033151106, 1066.49547719 ],
    [ 4.728e-08, 5.77993082374, 9808.53818466 ],
    [ 5.128e-08, 1.57178942294, 6525.80445397 ],
    [ 4.523e-08, 1.44233177206, 3369.06161417 ],
    [ 6.205e-08, 4.48163731718, 22747.2907149 ],
    [ 6.169e-08, 4.59085555242, 6531.66165626 ],
    [ 5.329e-08, 4.55141789349, 1744.42598442 ],
    [ 4.514e-08, 5.94508421612, 6894.52394884 ],
    [ 4.33e-08, 3.10899106071, 4569.57454002 ],
    [ 5.367e-08, 5.08071026709, 2707.82868739 ],
    [ 5.138e-08, 1.28584065229, 8439.87793182 ],
    [ 4.12e-08, 5.48544036931, 2699.73481932 ],
    [ 5.398e-08, 5.21710209952, 5305.45105355 ],
    [ 4.45e-08, 5.56771154217, 16865.5287696 ],
    [ 3.898e-08, 1.48753002285, 9168.64089835 ],
    [ 3.858e-08, 1.23056079731, 16858.4825329 ],
    [ 3.764e-08, 0.27080818668, 17395.2197347 ],
    [ 4.687e-08, 3.0570907584, 5518.75014899 ],
    [ 4.264e-08, 2.79046663043, 3503.07906283 ],
    [ 3.864e-08, 0.37957786186, 10177.2576795 ],
    [ 3.992e-08, 1.84425142473, 3134.42687826 ],
    [ 3.658e-08, 2.95544843123, 6144.4203413 ],
    [ 3.65e-08, 1.58041651396, 6680.24453233 ],
    [ 3.945e-08, 1.98631850445, 8969.56889691 ]
    // 198 terms retained
];

KMG.VSOPTerms.mars_R1 = [
    [ 0.01107433345, 2.03250524857, 3340.6124267 ],
    [ 0.00103175887, 2.37071847807, 6681.2248534 ],
    [ 0.000128772, 0, 0 ],
    [ 0.0001081588, 2.70888095665, 10021.8372801 ],
    [ 1.19455e-05, 3.04702256206, 13362.4497068 ],
    [ 4.38582e-06, 2.88835054603, 2281.23049651 ],
    [ 3.957e-06, 3.42323670971, 3344.13554505 ],
    [ 1.82576e-06, 1.58427562964, 2544.31441988 ],
    [ 1.35851e-06, 3.38507063082, 16703.0621335 ],
    [ 1.28199e-06, 0.62991771813, 1059.38193019 ],
    [ 1.27059e-06, 1.95391155885, 796.298006816 ],
    [ 1.18443e-06, 2.99762091382, 2146.16541648 ],
    [ 1.28362e-06, 6.04343227063, 3337.08930835 ],
    [ 8.7534e-07, 3.42053385867, 398.149003408 ],
    [ 8.3021e-07, 3.85575072018, 3738.76143011 ],
    [ 7.5604e-07, 4.45097659377, 6151.5338883 ],
    [ 7.2002e-07, 2.76443992447, 529.690965095 ],
    [ 6.6545e-07, 2.5487838147, 1751.53953142 ],
    [ 5.4305e-07, 0.67754203387, 8962.45534991 ],
    [ 5.1043e-07, 3.72584855417, 6684.74797175 ],
    [ 6.6413e-07, 4.40596377334, 1748.01641307 ],
    [ 4.786e-07, 2.28524521788, 2914.01423582 ],
    [ 4.942e-07, 5.72961379219, 3340.59517305 ],
    [ 4.942e-07, 1.47720011103, 3340.62968035 ],
    [ 5.7519e-07, 0.5435613312, 1194.44701022 ],
    [ 4.832e-07, 2.58061402348, 3149.16416059 ],
    [ 3.6383e-07, 6.02729341698, 3185.19202727 ],
    [ 3.7161e-07, 5.81436290851, 1349.86740966 ],
    [ 3.6035e-07, 5.89515829011, 3333.4988797 ],
    [ 3.1111e-07, 0.97820401887, 191.448266112 ],
    [ 3.8957e-07, 2.31902442004, 4136.91043352 ],
    [ 2.7256e-07, 5.41369838171, 1592.59601363 ],
    [ 2.4302e-07, 3.75838444077, 155.420399434 ],
    [ 2.2808e-07, 1.74818178182, 5088.62883977 ],
    [ 2.2322e-07, 0.93941901193, 951.718406251 ],
    [ 2.1712e-07, 3.83569490817, 6283.07584999 ],
    [ 2.1302e-07, 0.78030571909, 1589.07289528 ],
    [ 2.1631e-07, 4.56903942095, 3532.06069281 ],
    [ 1.7957e-07, 4.21923537063, 3870.30339179 ],
    [ 1.8241e-07, 0.41334220202, 5486.77784318 ],
    [ 1.625e-07, 3.80772429678, 3340.5451164 ],
    [ 1.6803e-07, 5.54855432911, 3097.88382273 ],
    [ 1.6852e-07, 4.53696884484, 4292.33083295 ],
    [ 1.5749e-07, 4.75766175289, 9492.146315 ],
    [ 1.5747e-07, 3.72356261757, 20043.6745602 ],
    [ 2.0429e-07, 3.13541604634, 4690.47983636 ],
    [ 1.4699e-07, 5.95340513928, 3894.18182954 ],
    [ 1.6251e-07, 3.39910570757, 3340.679737 ],
    [ 1.4256e-07, 3.99914527335, 1990.74501704 ],
    [ 1.6529e-07, 0.96740368703, 4399.99435689 ],
    [ 1.3011e-07, 5.14215010082, 6677.70173505 ],
    [ 1.2482e-07, 1.03238555854, 3341.59274777 ],
    [ 1.6454e-07, 3.53827765951, 2700.71514039 ],
    [ 1.6167e-07, 2.3489111087, 553.569402842 ],
    [ 1.3169e-07, 0.41462220221, 5614.72937621 ],
    [ 1.127e-07, 1.02387117266, 12303.0677766 ],
    [ 1.241e-07, 6.23139144626, 5628.95647021 ],
    [ 1.2747e-07, 0.69046237163, 3723.50895892 ],
    [ 1.1828e-07, 6.25270937134, 2274.11694951 ],
    [ 1.0382e-07, 1.23229650709, 426.598190876 ],
    [ 1.1207e-07, 1.31732435116, 3496.03282613 ],
    [ 1.0345e-07, 0.90062869301, 4535.05943692 ],
    [ 1.2214e-07, 4.22347837212, 7079.37385681 ],
    [ 9.764e-08, 3.45310129694, 382.896532223 ],
    [ 8.583e-08, 1.1647889051, 2787.04302386 ],
    [ 7.879e-08, 5.73808303461, 2288.34404351 ],
    [ 9.192e-08, 1.81719352796, 6681.24210705 ],
    [ 7.752e-08, 4.15038634174, 6041.32756709 ],
    [ 9.192e-08, 6.06960723129, 6681.20759975 ],
    [ 9.008e-08, 2.58179552398, 2388.89402045 ],
    [ 6.77e-08, 0.240118807, 11773.3768115 ],
    [ 7.088e-08, 3.51428380287, 8031.09226306 ],
    [ 9.159e-08, 3.90203365989, 3553.91152214 ],
    [ 7.233e-08, 3.70260535699, 2818.03500861 ],
    [ 6.701e-08, 4.25537421062, 242.728603974 ],
    [ 6.534e-08, 0.04317593308, 2957.71589448 ],
    [ 8.783e-08, 2.19764346848, 1221.84856632 ],
    [ 6.54e-08, 2.11811131682, 8429.24126647 ],
    [ 6.835e-08, 4.04527289029, 10025.3603984 ],
    [ 7.279e-08, 4.26969778292, 2803.8079146 ],
    [ 7.679e-08, 1.00816125095, 8432.76438482 ],
    [ 5.736e-08, 3.13988802339, 213.299095438 ],
    [ 5.343e-08, 3.7818416468, 5092.15195812 ],
    [ 5.985e-08, 2.96429619989, 6489.77658729 ],
    [ 5.132e-08, 3.98288020531, 7.1135470008 ],
    [ 6.264e-08, 1.90345600186, 5621.84292321 ],
    [ 5.238e-08, 2.67050910776, 7477.52286022 ],
    [ 6.264e-08, 1.60046198142, 3347.7259737 ],
    [ 6.527e-08, 2.76220386403, 3339.63210563 ],
    [ 4.594e-08, 1.82031785094, 2810.92146161 ],
    [ 5.46e-08, 4.60869963415, 3583.34103067 ],
    [ 4.73e-08, 0.90611934427, 5099.26550512 ],
    [ 5.484e-08, 4.91405753832, 7632.94325965 ],
    [ 4.002e-08, 4.1410000521, 9623.68827669 ],
    [ 3.836e-08, 0.03411499404, 7234.79425624 ],
    [ 3.618e-08, 5.76553319747, 4933.20844033 ],
    [ 3.747e-08, 0.08776717073, 6525.80445397 ],
    [ 3.016e-08, 3.73804058695, 6681.2921637 ],
    [ 3.975e-08, 4.91286826343, 2942.46342329 ],
    [ 3.911e-08, 0.67457174687, 3127.31333126 ],
    [ 3.923e-08, 3.07698893109, 3.523118349 ],
    [ 3.943e-08, 0.53936955267, 5884.92684658 ],
    [ 2.902e-08, 4.66228680082, 7210.91581849 ],
    [ 2.803e-08, 1.00505054832, 7064.12138562 ],
    [ 3.152e-08, 4.54611126545, 2487.41604495 ],
    [ 2.797e-08, 0.05226680768, 639.897286314 ],
    [ 2.758e-08, 5.17057629399, 5828.02847165 ],
    [ 3.02e-08, 4.14658810846, 6681.1575431 ],
    [ 3e-08, 0.82762095066, 5085.03841111 ],
    [ 3.022e-08, 2.59437829291, 2906.90068882 ],
    [ 2.673e-08, 0.69433657973, 2699.73481932 ],
    [ 2.593e-08, 1.08691889359, 4929.68532198 ],
    [ 3.127e-08, 0.99947199034, 2118.76386038 ],
    [ 2.597e-08, 5.01157388627, 10018.3141618 ],
    [ 2.606e-08, 5.34395258978, 10973.5556864 ],
    [ 2.779e-08, 3.98360727591, 6467.92575796 ],
    [ 2.457e-08, 1.52659064342, 6836.64525283 ],
    [ 2.381e-08, 3.93615187831, 11371.7046898 ],
    [ 2.584e-08, 5.08907827632, 12832.7587417 ]
    // 119 terms retained
];

KMG.VSOPTerms.mars_R2 = [
    [ 0.00044242249, 0.47930604954, 3340.6124267 ],
    [ 8.138042e-05, 0.86998389204, 6681.2248534 ],
    [ 1.274915e-05, 1.22593985222, 10021.8372801 ],
    [ 1.87388e-06, 1.57298976045, 13362.4497068 ],
    [ 4.0745e-07, 1.97082077028, 3344.13554505 ],
    [ 5.2395e-07, 3.14159265359, 0 ],
    [ 2.6617e-07, 1.91665337822, 16703.0621335 ],
    [ 1.7828e-07, 4.43491476419, 2281.23049651 ],
    [ 1.1713e-07, 4.52509926559, 3185.19202727 ],
    [ 1.021e-07, 5.3914732206, 1059.38193019 ],
    [ 9.95e-08, 0.41865678448, 796.298006816 ],
    [ 9.236e-08, 4.53559625376, 2146.16541648 ],
    [ 7.299e-08, 3.1421451312, 2544.31441988 ],
    [ 7.214e-08, 2.29302335628, 6684.74797175 ],
    [ 6.81e-08, 5.26707245601, 155.420399434 ],
    [ 6.526e-08, 2.307724561, 3738.76143011 ],
    [ 7.783e-08, 5.93373461009, 1748.01641307 ],
    [ 5.84e-08, 1.0519182029, 1349.86740966 ],
    [ 6.75e-08, 5.30191763402, 1194.44701022 ],
    [ 4.695e-08, 0.76881032874, 3097.88382273 ],
    [ 5.39e-08, 1.0020006836, 3149.16416059 ],
    [ 4.406e-08, 2.45557331437, 951.718406251 ],
    [ 4.286e-08, 3.89642578846, 1592.59601363 ],
    [ 3.516e-08, 1.84991934524, 398.149003408 ],
    [ 3.699e-08, 2.26016989021, 20043.6745602 ],
    [ 3.378e-08, 3.81703201748, 1751.53953142 ],
    [ 4.585e-08, 0.80785643853, 4136.91043352 ],
    [ 3.201e-08, 2.11661594157, 5614.72937621 ],
    [ 3.62e-08, 1.32428600053, 3333.4988797 ],
    [ 2.915e-08, 1.19342490174, 529.690965095 ],
    [ 2.979e-08, 2.86468474914, 6151.5338883 ],
    [ 3.057e-08, 4.55288594507, 5628.95647021 ],
    [ 2.906e-08, 1.20300479533, 3894.18182954 ],
    [ 3.848e-08, 3.86071515455, 553.569402842 ],
    [ 2.819e-08, 2.48714583532, 1990.74501704 ],
    [ 2.657e-08, 6.07409846258, 4292.33083295 ],
    [ 2.698e-08, 2.92100135189, 3496.03282613 ],
    [ 2.396e-08, 5.94193484091, 2787.04302386 ],
    [ 2.263e-08, 2.56188049651, 191.448266112 ],
    [ 2.169e-08, 5.36834559071, 8962.45534991 ],
    [ 2.149e-08, 2.74919289456, 242.728603974 ],
    [ 2.218e-08, 1.85260509629, 3337.08930835 ],
    [ 1.998e-08, 5.76396921426, 3341.59274777 ],
    [ 1.999e-08, 3.82347205028, 2914.01423582 ],
    [ 1.835e-08, 5.68648448195, 1589.07289528 ],
    [ 1.81e-08, 3.32122811143, 5088.62883977 ],
    [ 1.968e-08, 4.17404480033, 3340.59517305 ],
    [ 2.411e-08, 4.68376177281, 4690.47983636 ],
    [ 1.967e-08, 6.2057036343, 3340.62968035 ],
    [ 1.626e-08, 5.67648778513, 4535.05943692 ],
    [ 2.161e-08, 1.07446445419, 2388.89402045 ],
    [ 1.965e-08, 3.10811453974, 3583.34103067 ],
    [ 1.985e-08, 5.75867975763, 4399.99435689 ],
    [ 1.504e-08, 4.95929390466, 382.896532223 ],
    [ 1.276e-08, 4.82147500391, 2957.71589448 ],
    [ 1.475e-08, 2.22614544794, 3723.50895892 ],
    [ 1.196e-08, 3.26743061042, 9492.146315 ],
    [ 1.349e-08, 4.87558985925, 6525.80445397 ],
    [ 1.436e-08, 2.6975402327, 7079.37385681 ],
    [ 1.223e-08, 2.61880227353, 10025.3603984 ],
    [ 1.402e-08, 5.19177439326, 2700.71514039 ],
    [ 1.202e-08, 0.93436294282, 2810.92146161 ],
    [ 8.7e-09, 5.81258009514, 12303.0677766 ],
    [ 8.67e-09, 2.20048756217, 2699.73481932 ],
    [ 8.31e-09, 2.01782919511, 5092.15195812 ],
    [ 8.56e-09, 5.96129932558, 426.598190876 ],
    [ 8.47e-09, 2.26415579047, 6283.07584999 ],
    [ 9.17e-09, 1.4025908126, 6489.77658729 ],
    [ 8.33e-09, 1.17376008822, 7477.52286022 ],
    [ 1.041e-08, 6.27097603149, 3347.7259737 ],
    [ 9.65e-09, 3.40293030184, 5621.84292321 ],
    [ 7.23e-09, 4.26276570887, 4933.20844033 ],
    [ 7.7e-09, 2.06490049164, 5486.77784318 ],
    [ 7.06e-09, 2.34080074294, 7.1135470008 ],
    [ 9.54e-09, 2.11093711712, 3870.30339179 ],
    [ 8.44e-09, 2.2379157639, 3553.91152214 ],
    [ 6.47e-09, 2.24565892529, 3340.5451164 ],
    [ 6.53e-09, 3.98464883505, 6677.70173505 ],
    [ 7.17e-09, 0.29523050971, 6681.24210705 ],
    [ 8.28e-09, 0.22887694811, 3532.06069281 ],
    [ 6.12e-09, 1.56040446304, 7234.79425624 ],
    [ 7.17e-09, 4.54583138124, 6681.20759975 ],
    [ 5.85e-09, 3.29614213819, 1221.84856632 ],
    [ 6.46e-09, 1.8361516834, 3340.679737 ],
    [ 5.6e-09, 5.05995427063, 8031.09226306 ],
    [ 6.51e-09, 0.16211451692, 7632.94325965 ]
    // 86 terms retained
];

KMG.VSOPTerms.mars_R3 = [
    [ 1.113108e-05, 5.14987305093, 3340.6124267 ],
    [ 4.24447e-06, 5.61343952053, 6681.2248534 ],
    [ 1.00044e-06, 5.99727457548, 10021.8372801 ],
    [ 1.9606e-07, 0.07631453783, 13362.4497068 ],
    [ 3.478e-08, 0.42912010211, 16703.0621335 ],
    [ 4.693e-08, 3.14159265359, 0 ],
    [ 2.87e-08, 0.44692002393, 3344.13554505 ],
    [ 2.428e-08, 3.02114808809, 3185.19202727 ]
    // 8 terms retained
];

KMG.VSOPTerms.mars_R4 = [
    [ 1.9551e-07, 3.58210746512, 3340.6124267 ],
    [ 1.6322e-07, 4.05115851142, 6681.2248534 ],
    [ 5.848e-08, 4.4638164658, 10021.8372801 ],
    [ 1.533e-08, 4.84332951095, 13362.4497068 ],
    [ 3.75e-09, 1.50951652931, 3185.19202727 ],
    [ 3.4e-09, 5.20519444932, 16703.0621335 ]
    // 6 terms retained
];

KMG.VSOPTerms.mars_R5 = [
    [ 4.75e-09, 2.47621038205, 6681.2248534 ],
    [ 2.7e-09, 2.90961348988, 10021.8372801 ]
    // 2 terms retained
];



KMG.VSOPSeries = {};

KMG.VSOPSeries.earth_L = [
	KMG.VSOPTerms.earth_L0,
	KMG.VSOPTerms.earth_L1,
	KMG.VSOPTerms.earth_L2,
	KMG.VSOPTerms.earth_L3,
	KMG.VSOPTerms.earth_L4,
	KMG.VSOPTerms.earth_L5
];


KMG.VSOPSeries.earth_B = [
	KMG.VSOPTerms.earth_B0,
	KMG.VSOPTerms.earth_B1,
	KMG.VSOPTerms.earth_B2
];

KMG.VSOPSeries.earth_R = [
	KMG.VSOPTerms.earth_R0,
	KMG.VSOPTerms.earth_R1,
	KMG.VSOPTerms.earth_R2,
	KMG.VSOPTerms.earth_R3,
	KMG.VSOPTerms.earth_R4,
	KMG.VSOPTerms.earth_R5
];


KMG.VSOPSeries.mars_L = [
	KMG.VSOPTerms.mars_L0,
	KMG.VSOPTerms.mars_L1,
	KMG.VSOPTerms.mars_L2,
	KMG.VSOPTerms.mars_L3,
	KMG.VSOPTerms.mars_L4,
	KMG.VSOPTerms.mars_L5
];


KMG.VSOPSeries.mars_B = [
	KMG.VSOPTerms.mars_B0,
	KMG.VSOPTerms.mars_B1,
	KMG.VSOPTerms.mars_B2,
	KMG.VSOPTerms.mars_B3,
	KMG.VSOPTerms.mars_B4,
	KMG.VSOPTerms.mars_B5
];

KMG.VSOPSeries.mars_R = [
	KMG.VSOPTerms.mars_R0,
	KMG.VSOPTerms.mars_R1,
	KMG.VSOPTerms.mars_R2,
	KMG.VSOPTerms.mars_R3,
	KMG.VSOPTerms.mars_R4,
	KMG.VSOPTerms.mars_R5
];


KMG.VSOP87Orbit = function(L, B, R, period) {

	function positionAtTime(jd, preserveDirection) {
		var t = (jd - 2451545) / 365250;

		var l = 0;
		var b = 0;
		var r = 0;
		
		var T = 1;
		for (var i = 0; i < L.length; i++) {
			var s = 0;
			for (var j = 0; j < L[i].length; j++) {
				s += L[i][j][0] * Math.cos(L[i][j][1] + L[i][j][2] * t);
			}
			l += s * T;
            T = t * T;
		}
		
		var T = 1;
		for (var i = 0; i < B.length; i++) {
			var s = 0;
			for (var j = 0; j < B[i].length; j++) {
				s += B[i][j][0] * Math.cos(B[i][j][1] + B[i][j][2] * t);
			}
			b += s * T;
            T = t * T;
		}
		
		var T = 1;
		for (var i = 0; i < R.length; i++) {
			var s = 0;
			for (var j = 0; j < R[i].length; j++) {
				s += R[i][j][0] * Math.cos(R[i][j][1] + R[i][j][2] * t);
			}
			r += s * T;
            T = t * T;
		}
		

		l = KMG.Math.clamp(l, 2 * Math.PI);

		//b -= Math.PI / 2;
       // l += Math.PI;
		
		var modB = (preserveDirection) ? 0 : -(Math.PI / 2);
		var modL = (preserveDirection) ? 0 : Math.PI;
		b += modB;
		l += modL;
		
		var x = Math.cos(l) * Math.sin(b) * r;
		var y = Math.cos(b) * r;
		var z = -Math.sin(l) * Math.sin(b) * r;
		
		
		
		var position = new THREE.Vector3(x, y, z);
		position.l = l;
		position.b = b;
		position.r = r;

		return position;
	}
	
	return {
		positionAtTime : positionAtTime,
		period : period,
		epoch : KMG.Util.julianNow()
	};
	
};
KMG.VSOP87Orbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.VSOP87OrbitRect = function(X, Y, Z, period) {

	
	function evaluate(terms, t) {
		var v = 0;
		
		var T = 1;
		for (var i = 0; i < terms.length; i++) {
			var s = 0;
			for (var j = 0; j < terms[i].length; j++) {
				s += terms[i][j][0] * Math.cos(terms[i][j][1] + terms[i][j][2] * t);
			}
			v += s * T;
            T = t * T;
		}
		
		return v;
	}

	function positionAtTime(jd) {
		var t = (jd - 2451545) / 365250;

		var x = 0;
		var y = 0; 
		var z = 0;
		
		var T = 1;
		for (var i = 0; i < X.length; i++) {
			var s = 0;
			for (var j = 0; j < X[i].length; j++) {
				s += X[i][j][0] * Math.cos(X[i][j][1] + X[i][j][2] * t);
			}
			x += s * T;
            T = t * T;
		}
		
		
		var T = 1;
		for (var i = 0; i < Y.length; i++) {
			var s = 0;
			for (var j = 0; j < Y[i].length; j++) {
				s += Y[i][j][0] * Math.cos(Y[i][j][1] + Y[i][j][2] * t);
			}
			y += s * T;
            T = t * T;
		}
		
		var T = 1;
		for (var i = 0; i < Z.length; i++) {
			var s = 0;
			for (var j = 0; j < Z[i].length; j++) {
				s += Z[i][j][0] * Math.cos(Z[i][j][1] + Z[i][j][2] * t);
			}
			z += s * T;
            T = t * T;
		}
		
		
		
		var pos = new THREE.Vector3(x, z, -y);

		return pos;
	};
	
	return {
		positionAtTime : positionAtTime,
		period : period,
		epoch : KMG.Util.julianNow()
	};
	
};
KMG.VSOP87OrbitRect.prototype = Object.create( KMG.Orbit.prototype );	

/* File: CustomOrbits.js */
/** Custom and Body-Specific Orbital Algorithms
 * http://www.apoapsys.com
 * 
 * Copyright 2014 Kevin M. Gill <kmsmgill@gmail.com>
 *
 * Uses algorithms from the VSOP87 theory for the orbits
 * of major planets.
 * ftp://ftp.bdl.fr/pub/ephem/planets/vsop87/
 * 
 * Code is also adapted from Celestia source
 * vsop87.cpp
 * customorbit.cpp
 * 
 * Uses algorithms from:
 * Meeus, Jean: Astronomical Algorithms.
 * Richmond, Virg.: Willmann-Bell, 2009.
 * ISBN: 978-0943396613
 * http://amzn.com/0943396611
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */


KMG.OrbitUtil = {};

KMG.OrbitUtil.anomaly = function(meanAnomaly, eccentricity)
{
    var e, delta, err;
    var tol = 0.00000001745;
    var iterations = 20;	// limit while() to maximum of 20 iterations.

    e = meanAnomaly - 2 * Math.PI *  Math.floor(meanAnomaly / (2*Math.PI));
    err = 1;
    while(Math.abs(err) > tol && iterations > 0)
    {
        err = e - eccentricity*Math.sin(e) - meanAnomaly;
        delta = err / (1 - eccentricity * Math.cos(e));
        e -= delta;
        iterations--;
    }

    var trueAnomaly = 2*Math.atan(Math.sqrt((1+eccentricity)/(1-eccentricity))*Math.tan(e/2));
    var eccentricAnomaly = e;
    
    
    return {
		trueAnomaly : trueAnomaly,
		eccentricAnomaly : eccentricAnomaly
	};
};


KMG.CustomIoOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	}
	
	function degToRad(v) {
		return v * KMG.PI_BY_180;
	}
	
	function positionAtTime(jd) {
		//var t = (jd - 2451545) / 36525;
		var t = jd - 2443000.5;
		//var t = (jd - 2443000.5) / 36525;
		var e = KMG.Jupiter.computeElements(t);
		
		var LPEJ = e.Π;
		
		// Calculate periodic terms for longitude
		var Σ1 = 0.47259*KMG.Math.dsin(2*(e.l1 - e.l2)) - 0.03478*KMG.Math.dsin(e.p3 - e.p4)
				+ 0.01081*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p3) + 7.38e-3*KMG.Math.dsin(e.Φ)
				+ 7.13e-3*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p2) - 6.74e-3*KMG.Math.dsin(e.p1 + e.p3 - 2*LPEJ - 2*e.G)
				+ 6.66e-3*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p4) + 4.45e-3*KMG.Math.dsin(e.l1 - e.p3)
				- 3.54e-3*KMG.Math.dsin(e.l1 - e.l2) - 3.17e-3*KMG.Math.dsin(2*(e.Ψ - LPEJ))
				+ 2.65e-3*KMG.Math.dsin(e.l1 - e.p4) - 1.86e-3*KMG.Math.dsin(e.G)
				+ 1.62e-3*KMG.Math.dsin(e.p2 - e.p3) + 1.58e-3*KMG.Math.dsin(4*(e.l1 - e.l2))
				- 1.55e-3*KMG.Math.dsin(e.l1 - e.l3) - 1.38e-3*KMG.Math.dsin(e.Ψ + e.w3 - 2*LPEJ - 2*e.G)
				- 1.15e-3*KMG.Math.dsin(2*(e.l1 - 2*e.l2 + e.w2)) + 8.9e-4*KMG.Math.dsin(e.p2 - e.p4)
				+ 8.5e-4*KMG.Math.dsin(e.l1 + e.p3 - 2*LPEJ - 2*e.G) + 8.3e-4*KMG.Math.dsin(e.w2 - e.w3)
				+ 5.3e-4*KMG.Math.dsin(e.Ψ - e.w2);
		Σ1 = KMG.Math.clamp(Σ1, 360.0);
		//Σ1 = degToRad(Σ1);
		var L = e.l1 + Σ1;

		// Calculate periodic terms for the tangent of the latitude
		var B = 6.393e-4*KMG.Math.dsin(L - e.w1) + 1.825e-4*KMG.Math.dsin(L - e.w2)
			+ 3.29e-5*KMG.Math.dsin(L - e.w3) - 3.11e-5*KMG.Math.dsin(L - e.Ψ)
			+ 9.3e-6*KMG.Math.dsin(L - e.w4) + 7.5e-6*KMG.Math.dsin(3*L - 4*e.l2 - 1.9927*Σ1 + e.w2)
			+ 4.6e-6*KMG.Math.dsin(L + e.Ψ - 2*LPEJ - 2*e.G);
		B = KMG.Math.datan(B);

		// Calculate the periodic terms for distance
		var R = -4.1339e-3*KMG.Math.dcos(2*(e.l1 - e.l2)) - 3.87e-5*KMG.Math.dcos(e.l1 - e.p3)
		  - 2.14e-5*KMG.Math.dcos(e.l1 - e.p4) + 1.7e-5*KMG.Math.dcos(e.l1 - e.l2)
		  - 1.31e-5*KMG.Math.dcos(4*(e.l1 - e.l2)) + 1.06e-5*KMG.Math.dcos(e.l1 - e.l3)
		  - 6.6e-6*KMG.Math.dcos(e.l1 + e.p3 - 2*LPEJ - 2*e.G);
		R = 5.90569 * KMG.Jupiter.radius * (1 + R) / KMG.AU_TO_KM;

		var T = (jd - 2433282.423) / 36525.0;
		var P = 1.3966626*T + 3.088e-4*T*T;
		L += P;

		//L += 22.203;
	
		
		L = L * KMG.PI_BY_180;
		B = B * KMG.PI_BY_180;
		
		B -= Math.PI / 2;
        L += Math.PI;
		
		var x = Math.cos(L) * Math.sin(B) * R;
		var y = Math.cos(B) * R;
		var z = -Math.sin(L) * Math.sin(B) * R;

		var position = new THREE.Vector3(x, y, z);
		position.l = L;
		position.b = B;
		position.r = R;
		return position;
	}
	

	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 1.769138,
		epoch : KMG.Util.julianNow()
	};
	
};
KMG.CustomIoOrbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.CustomEuropaOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	}
	
	function degToRad(v) {
		return v * KMG.PI_BY_180;
	}
	
	function positionAtTime(jd) {
		//var t = (jd - 2451545) / 36525;
		var t = (jd - 2443000.5);// / 36525;
		var e = KMG.Jupiter.computeElements(t);
		
		var LPEJ = e.Π;
		
		
		// Calculate periodic terms for lone.Gitude
		var Σ1 = 1.06476*KMG.Math.dsin(2*(e.l2 - e.l3)) + 0.04256*KMG.Math.dsin(e.l1 - 2*e.l2 + e.p3)
			  + 0.03581*KMG.Math.dsin(e.l2 - e.p3) + 0.02395*KMG.Math.dsin(e.l1 - 2*e.l2 + e.p4)
			  + 0.01984*KMG.Math.dsin(e.l2 - e.p4) - 0.01778*KMG.Math.dsin(e.Φ)
			  + 0.01654*KMG.Math.dsin(e.l2 - e.p2) + 0.01334*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p2)
			  + 0.01294*KMG.Math.dsin(e.p3 - e.p4) - 0.01142*KMG.Math.dsin(e.l2 - e.l3)
			  - 0.01057*KMG.Math.dsin(e.G) - 7.75e-3*KMG.Math.dsin(2*(e.Ψ - LPEJ))
			  + 5.24e-3*KMG.Math.dsin(2*(e.l1 - e.l2)) - 4.6e-3*KMG.Math.dsin(e.l1 - e.l3)
			  + 3.16e-3*KMG.Math.dsin(e.Ψ - 2*e.G + e.w3 - 2*LPEJ) - 2.03e-3*KMG.Math.dsin(e.p1 + e.p3 - 2*LPEJ - 2*e.G)
			  + 1.46e-3*KMG.Math.dsin(e.Ψ - e.w3) - 1.45e-3*KMG.Math.dsin(2*e.G)
			  + 1.25e-3*KMG.Math.dsin(e.Ψ - e.w4) - 1.15e-3*KMG.Math.dsin(e.l1 - 2*e.l3 + e.p3)
			  - 9.4e-4*KMG.Math.dsin(2*(e.l2 - e.w2)) + 8.6e-4*KMG.Math.dsin(2*(e.l1 - 2*e.l2 + e.w2))
			  - 8.6e-4*KMG.Math.dsin(5*e.G_ - 2*e.G + 0.9115) - 7.8e-4*KMG.Math.dsin(e.l2 - e.l4)
			  - 6.4e-4*KMG.Math.dsin(3*e.l3 - 7*e.l4 + 4*e.p4) + 6.4e-4*KMG.Math.dsin(e.p1 - e.p4)
			  - 6.3e-4*KMG.Math.dsin(e.l1 - 2*e.l3 + e.p4) + 5.8e-4*KMG.Math.dsin(e.w3 - e.w4)
			  + 5.6e-4*KMG.Math.dsin(2*(e.Ψ - LPEJ - e.G)) + 5.6e-4*KMG.Math.dsin(2*(e.l2 - e.l4))
			  + 5.5e-4*KMG.Math.dsin(2*(e.l1 - e.l3)) + 5.2e-4*KMG.Math.dsin(3*e.l3 - 7*e.l4 + e.p3 +3*e.p4)
			  - 4.3e-4*KMG.Math.dsin(e.l1 - e.p3) + 4.1e-4*KMG.Math.dsin(5*(e.l2 - e.l3))
			  + 4.1e-4*KMG.Math.dsin(e.p4 - LPEJ) + 3.2e-4*KMG.Math.dsin(e.w2 - e.w3)
			  + 3.2e-4*KMG.Math.dsin(2*(e.l3 - e.G - LPEJ));
		Σ1 = KMG.Math.clamp(Σ1, 360.0);
		//Σ1 = dee.GToRad(Σ1);
		var L = e.l2 + Σ1;

		// Calculate periodic terms for the tane.Gent of the latitude
		var B = 8.1004e-3*KMG.Math.dsin(L - e.w2) + 4.512e-4*KMG.Math.dsin(L - e.w3)
		  - 3.284e-4*KMG.Math.dsin(L - e.Ψ) + 1.160e-4*KMG.Math.dsin(L - e.w4)
		  + 2.72e-5*KMG.Math.dsin(e.l1 - 2*e.l3 + 1.0146*Σ1 + e.w2) - 1.44e-5*KMG.Math.dsin(L - e.w1)
		  + 1.43e-5*KMG.Math.dsin(L + e.Ψ - 2*LPEJ - 2*e.G) + 3.5e-6*KMG.Math.dsin(L - e.Ψ + e.G)
		  - 2.8e-6*KMG.Math.dsin(e.l1 - 2*e.l3 + 1.0146*Σ1 + e.w3);
		B = KMG.Math.datan(B);

		// Calculate the periodic terms for distance
		var R = 9.3848e-3*KMG.Math.dcos(e.l1 - e.l2) - 3.116e-4*KMG.Math.dcos(e.l2 - e.p3)
		  - 1.744e-4*KMG.Math.dcos(e.l2 - e.p4) - 1.442e-4*KMG.Math.dcos(e.l2 - e.p2)
		  + 5.53e-5*KMG.Math.dcos(e.l2 - e.l3) + 5.23e-5*KMG.Math.dcos(e.l1 - e.l3)
		  - 2.9e-5*KMG.Math.dcos(2*(e.l1 - e.l2)) + 1.64e-5*KMG.Math.dcos(2*(e.l2 - e.w2))
		  + 1.07e-5*KMG.Math.dcos(e.l1 - 2*e.l3 + e.p3) - 1.02e-5*KMG.Math.dcos(e.l2 - e.p1)
		  - 9.1e-6*KMG.Math.dcos(2*(e.l1 - e.l3));
		R = 9.39657 * KMG.Jupiter.radius * (1 + R) / KMG.AU_TO_KM;

		var T = (jd - 2433282.423) / 36525.0;
		var P = 1.3966626*T + 3.088e-4*T*T;
		L += P;
		//L += dee.GToRad(P);
		//L += 22.203;
		
		
		//console.info([L, B, R]);
		
		L = L * KMG.PI_BY_180;
		B = B * KMG.PI_BY_180;
		
		B -= Math.PI / 2;
        L += Math.PI;
                   
		var x = Math.cos(L) * Math.sin(B) * R;
		var y = Math.cos(B) * R;
		var z = -Math.sin(L) * Math.sin(B) * R;

		var position = new THREE.Vector3(x, y, z);
		position.l = L;
		position.b = B;
		position.r = R;
		return position;
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 3.5511810791,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomEuropaOrbit.prototype = Object.create( KMG.Orbit.prototype );





KMG.CustomGanymedeOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	}
	
	function degToRad(v) {
		return v * KMG.PI_BY_180;
	}
	
	function positionAtTime(jd) {
		//var t = (jd - 2451545) / 36525;
		var t = (jd - 2443000.5);// / 36525;
		var e = KMG.Jupiter.computeElements(t);
		
		var LPEJ = e.Π;
		var psi = e.Ψ;
		var phi = e.Φ;

		  
		  
		//Calculate periodic terms for lone.Gitude
		var Σ1 = 0.1649*KMG.Math.dsin(e.l3 - e.p3) + 0.09081*KMG.Math.dsin(e.l3 - e.p4)
			  - 0.06907*KMG.Math.dsin(e.l2 - e.l3) + 0.03784*KMG.Math.dsin(e.p3 - e.p4)
			  + 0.01846*KMG.Math.dsin(2*(e.l3 - e.l4)) - 0.01340*KMG.Math.dsin(e.G)
			  - 0.01014*KMG.Math.dsin(2*(psi - LPEJ)) + 7.04e-3*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p3)
			  - 6.2e-3*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p2) - 5.41e-3*KMG.Math.dsin(e.l3 - e.l4)
			  + 3.81e-3*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p4) + 2.35e-3*KMG.Math.dsin(psi - e.w3)
			  + 1.98e-3*KMG.Math.dsin(psi - e.w4) + 1.76e-3*KMG.Math.dsin(phi)
			  + 1.3e-3*KMG.Math.dsin(3*(e.l3 - e.l4)) + 1.25e-3*KMG.Math.dsin(e.l1 - e.l3)
			  - 1.19e-3*KMG.Math.dsin(5*e.G_ - 2*e.G + 0.9115) + 1.09e-3*KMG.Math.dsin(e.l1 - e.l2)
			  - 1.0e-3*KMG.Math.dsin(3*e.l3 - 7*e.l4 + 4*e.p4) + 9.1e-4*KMG.Math.dsin(e.w3 - e.w4)
			  + 8.0e-4*KMG.Math.dsin(3*e.l3 - 7*e.l4 + e.p3 + 3*e.p4) - 7.5e-4*KMG.Math.dsin(2*e.l2 - 3*e.l3 + e.p3)
			  + 7.2e-4*KMG.Math.dsin(e.p1 + e.p3 - 2*LPEJ - 2*e.G) + 6.9e-4*KMG.Math.dsin(e.p4 - LPEJ)
			  - 5.8e-4*KMG.Math.dsin(2*e.l3 - 3*e.l4 + e.p4) - 5.7e-4*KMG.Math.dsin(e.l3 - 2*e.l4 + e.p4)
			  + 5.6e-4*KMG.Math.dsin(e.l3 + e.p3 - 2*LPEJ - 2*e.G) - 5.2e-4*KMG.Math.dsin(e.l2 - 2*e.l3 + e.p1)
			  - 5.0e-4*KMG.Math.dsin(e.p2 - e.p3) + 4.8e-4*KMG.Math.dsin(e.l3 - 2*e.l4 + e.p3)
			  - 4.5e-4*KMG.Math.dsin(2*e.l2 - 3*e.l3 + e.p4) - 4.1e-4*KMG.Math.dsin(e.p2 - e.p4)
			  - 3.8e-4*KMG.Math.dsin(2*e.G) - 3.7e-4*KMG.Math.dsin(e.p3 - e.p4 + e.w3 - e.w4)
			  - 3.2e-4*KMG.Math.dsin(3*e.l3 - 7*e.l4 + 2*e.p3 + 2*e.p4) + 3.0e-4*KMG.Math.dsin(4*(e.l3 - e.l4))
			  + 2.9e-4*KMG.Math.dsin(e.l3 + e.p4 - 2*LPEJ - 2*e.G) - 2.8e-4*KMG.Math.dsin(e.w3 + psi - 2*LPEJ - 2*e.G)
			  + 2.6e-4*KMG.Math.dsin(e.l3 - LPEJ - e.G) + 2.4e-4*KMG.Math.dsin(e.l2 - 3*e.l3 + 2*e.l4)
			  + 2.1e-4*KMG.Math.dsin(2*(e.l3 - LPEJ - e.G)) - 2.1e-4*KMG.Math.dsin(e.l3 - e.p2)
			  + 1.7e-4*KMG.Math.dsin(e.l3 - e.p3);
		Σ1 = KMG.Math.clamp(Σ1, 360.0);
		//sie.Gma = dee.GToRad(sie.Gma);
		var L = e.l3 + Σ1;

		//Calculate periodic terms for the tane.Gent of the latitude
		var B = 3.2402e-3*KMG.Math.dsin(L - e.w3) - 1.6911e-3*KMG.Math.dsin(L - psi)
		  + 6.847e-4*KMG.Math.dsin(L - e.w4) - 2.797e-4*KMG.Math.dsin(L - e.w2)
		  + 3.21e-5*KMG.Math.dsin(L + psi - 2*LPEJ - 2*e.G) + 5.1e-6*KMG.Math.dsin(L - psi + e.G)
		  - 4.5e-6*KMG.Math.dsin(L - psi - e.G) - 4.5e-6*KMG.Math.dsin(L + psi - 2*LPEJ)
		  + 3.7e-6*KMG.Math.dsin(L + psi - 2*LPEJ - 3*e.G) + 3.0e-6*KMG.Math.dsin(2*e.l2 - 3*L + 4.03*Σ1 + e.w2)
		  - 2.1e-6*KMG.Math.dsin(2*e.l2 - 3*L + 4.03*Σ1 + e.w3);
		B = KMG.Math.datan(B);

		//Calculate the periodic terms for distance
		var R = -1.4388e-3*KMG.Math.dcos(e.l3 - e.p3) - 7.919e-4*KMG.Math.dcos(e.l3 - e.p4)
		  + 6.342e-4*KMG.Math.dcos(e.l2 - e.l3) - 1.761e-4*KMG.Math.dcos(2*(e.l3 - e.l4))
		  + 2.94e-5*KMG.Math.dcos(e.l3 - e.l4) - 1.56e-5*KMG.Math.dcos(3*(e.l3 - e.l4))
		  + 1.56e-5*KMG.Math.dcos(e.l1 - e.l3) - 1.53e-5*KMG.Math.dcos(e.l1 - e.l2)
		  + 7.0e-6*KMG.Math.dcos(2*e.l2 - 3*e.l3 + e.p3) - 5.1e-6*KMG.Math.dcos(e.l3 + e.p3 - 2*LPEJ - 2*e.G);
		R = 14.98832 * KMG.Jupiter.radius * (1 + R) / KMG.AU_TO_KM;
		
		var T = (jd - 2433282.423) / 36525.0;
		var P = 1.3966626*T + 3.088e-4*T*T;
		L += P;
		//L += dee.GToRad(P);

		//L += JupAscendingNode;
		  
		
		//console.info([L, B, R]);
		
		L = L * KMG.PI_BY_180;
		B = B * KMG.PI_BY_180;
		
		B -= Math.PI / 2;
        L += Math.PI;
                   
		var x = Math.cos(L) * Math.sin(B) * R;
		var y = Math.cos(B) * R;
		var z = -Math.sin(L) * Math.sin(B) * R;

		var position = new THREE.Vector3(x, y, z);
		position.l = L;
		position.b = B;
		position.r = R;
		return position;
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 3.5511810791,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomGanymedeOrbit.prototype = Object.create( KMG.Orbit.prototype );




KMG.CustomCallistoOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	}
	
	function degToRad(v) {
		return v * KMG.PI_BY_180;
	}
	
	function positionAtTime(jd) {
		//var t = (jd - 2451545) / 36525;
		var t = (jd - 2443000.5);// / 36525;
		var e = KMG.Jupiter.computeElements(t);
		
		var LPEJ = e.Π;
		var psi = e.Ψ;
		var phi = e.Φ;

		  
		  
		//Calculate periodic terms for lone.Gitude
		var Σ1 =
			0.84287*KMG.Math.dsin(e.l4 - e.p4)
			+ 0.03431*KMG.Math.dsin(e.p4 - e.p3)
			- 0.03305*KMG.Math.dsin(2*(psi - LPEJ))
			- 0.03211*KMG.Math.dsin(e.G)
			- 0.01862*KMG.Math.dsin(e.l4 - e.p3)
			+ 0.01186*KMG.Math.dsin(psi - e.w4)
			+ 6.23e-3*KMG.Math.dsin(e.l4 + e.p4 - 2*e.G - 2*LPEJ)
			+ 3.87e-3*KMG.Math.dsin(2*(e.l4 - e.p4))
			- 2.84e-3*KMG.Math.dsin(5*e.G_ - 2*e.G + 0.9115)
			- 2.34e-3*KMG.Math.dsin(2*(psi - e.p4))
			- 2.23e-3*KMG.Math.dsin(e.l3 - e.l4)
			- 2.08e-3*KMG.Math.dsin(e.l4 - LPEJ)
			+ 1.78e-3*KMG.Math.dsin(psi + e.w4 - 2*e.p4)
			+ 1.34e-3*KMG.Math.dsin(e.p4 - LPEJ)
			+ 1.25e-3*KMG.Math.dsin(2*(e.l4 - e.G - LPEJ))
			- 1.17e-3*KMG.Math.dsin(2*e.G)
			- 1.12e-3*KMG.Math.dsin(2*(e.l3 - e.l4))
			+ 1.07e-3*KMG.Math.dsin(3*e.l3 - 7*e.l4 + 4*e.p4)
			+ 1.02e-3*KMG.Math.dsin(e.l4 - e.G - LPEJ)
			+ 9.6e-4*KMG.Math.dsin(2*e.l4 - psi - e.w4)
			+ 8.7e-4*KMG.Math.dsin(2*(psi - e.w4))
			- 8.5e-4*KMG.Math.dsin(3*e.l3 - 7*e.l4 + e.p3 + 3*e.p4)
			+ 8.5e-4*KMG.Math.dsin(e.l3 - 2*e.l4 + e.p4)
			- 8.1e-4*KMG.Math.dsin(2*(e.l4 - psi))
			+ 7.1e-4*KMG.Math.dsin(e.l4 + e.p4 - 2*LPEJ - 3*e.G)
			+ 6.1e-4*KMG.Math.dsin(e.l1 - e.l4)
			- 5.6e-4*KMG.Math.dsin(psi - e.w3)
			- 5.4e-4*KMG.Math.dsin(e.l3 - 2*e.l4 + e.p3)
			+ 5.1e-4*KMG.Math.dsin(e.l2 - e.l4)
			+ 4.2e-4*KMG.Math.dsin(2*(psi - e.G - LPEJ))
			+ 3.9e-4*KMG.Math.dsin(2*(e.p4 - e.w4))
			+ 3.6e-4*KMG.Math.dsin(psi + LPEJ - e.p4 - e.w4)
			+ 3.5e-4*KMG.Math.dsin(2*e.G_ - e.G + 3.2877)
			- 3.5e-4*KMG.Math.dsin(e.l4 - e.p4 + 2*LPEJ - 2*psi)
			- 3.2e-4*KMG.Math.dsin(e.l4 + e.p4 - 2*LPEJ - e.G)
			+ 3.0e-4*KMG.Math.dsin(2*e.G_ - 2*e.G + 2.6032)
			+ 2.9e-4*KMG.Math.dsin(3*e.l3 - 7*e.l4 + 2*e.p3 + 2*e.p4)
			+ 2.8e-4*KMG.Math.dsin(e.l4 - e.p4 + 2*psi - 2*LPEJ)
			- 2.8e-4*KMG.Math.dsin(2*(e.l4 - e.w4))
			- 2.7e-4*KMG.Math.dsin(e.p3 - e.p4 + e.w3 - e.w4)
			- 2.6e-4*KMG.Math.dsin(5*e.G_ - 3*e.G + 3.2877)
			+ 2.5e-4*KMG.Math.dsin(e.w4 - e.w3)
			- 2.5e-4*KMG.Math.dsin(e.l2 - 3*e.l3 + 2*e.l4)
			- 2.3e-4*KMG.Math.dsin(3*(e.l3 - e.l4))
			+ 2.1e-4*KMG.Math.dsin(2*e.l4 - 2*LPEJ - 3*e.G)
			- 2.1e-4*KMG.Math.dsin(2*e.l3 - 3*e.l4 + e.p4)
			+ 1.9e-4*KMG.Math.dsin(e.l4 - e.p4 - e.G)
			- 1.9e-4*KMG.Math.dsin(2*e.l4 - e.p3 - e.p4)
			- 1.8e-4*KMG.Math.dsin(e.l4 - e.p4 + e.G)
			- 1.6e-4*KMG.Math.dsin(e.l4 + e.p3 - 2*LPEJ - 2*e.G);
		Σ1 = KMG.Math.clamp(Σ1, 360.0);
		//Σ1 = dee.GToRad(Σ1);
		var L = e.l4 + Σ1;

		//Calculate periodic terms for the tane.Gent of the latitude
		var B =
			- 7.6579e-3 * KMG.Math.dsin(L - psi)
			+ 4.4134e-3 * KMG.Math.dsin(L - e.w4)
			- 5.112e-4  * KMG.Math.dsin(L - e.w3)
			+ 7.73e-5   * KMG.Math.dsin(L + psi - 2*LPEJ - 2*e.G)
			+ 1.04e-5   * KMG.Math.dsin(L - psi + e.G)
			- 1.02e-5   * KMG.Math.dsin(L - psi - e.G)
			+ 8.8e-6    * KMG.Math.dsin(L + psi - 2*LPEJ - 3*e.G)
			- 3.8e-6    * KMG.Math.dsin(L + psi - 2*LPEJ - e.G);
		B = KMG.Math.datan(B);

		//Calculate the periodic terms for distance
		var R =
			- 7.3546e-3 * KMG.Math.dcos(e.l4 - e.p4)
			+ 1.621e-4  * KMG.Math.dcos(e.l4 - e.p3)
			+ 9.74e-5   * KMG.Math.dcos(e.l3 - e.l4)
			- 5.43e-5   * KMG.Math.dcos(e.l4 + e.p4 - 2*LPEJ - 2*e.G)
			- 2.71e-5   * KMG.Math.dcos(2*(e.l4 - e.p4))
			+ 1.82e-5   * KMG.Math.dcos(e.l4 - LPEJ)
			+ 1.77e-5   * KMG.Math.dcos(2*(e.l3 - e.l4))
			- 1.67e-5   * KMG.Math.dcos(2*e.l4 - psi - e.w4)
			+ 1.67e-5   * KMG.Math.dcos(psi - e.w4)
			- 1.55e-5   * KMG.Math.dcos(2*(e.l4 - LPEJ - e.G))
			+ 1.42e-5   * KMG.Math.dcos(2*(e.l4 - psi))
			+ 1.05e-5   * KMG.Math.dcos(e.l1 - e.l4)
			+ 9.2e-6    * KMG.Math.dcos(e.l2 - e.l4)
			- 8.9e-6    * KMG.Math.dcos(e.l4 - LPEJ -e.G)
			- 6.2e-6    * KMG.Math.dcos(e.l4 + e.p4 - 2*LPEJ - 3*e.G)
			+ 4.8e-6    * KMG.Math.dcos(2*(e.l4 - e.w4));

		R = 26.36273 * KMG.Jupiter.radius * (1 + R) / KMG.AU_TO_KM;
		var T = (jd - 2433282.423) / 36525.0;
		var P = 1.3966626*T + 3.088e-4*T*T;
		L += P;
		//L += degToRad(P);

		//L += JupAscendingNode;
		  
		
		//console.info([L, B, R]);
		
		L = L * KMG.PI_BY_180;
		B = B * KMG.PI_BY_180;
		
		B -= Math.PI / 2;
        L += Math.PI;
                   
		var x = Math.cos(L) * Math.sin(B) * R;
		var y = Math.cos(B) * R;
		var z = -Math.sin(L) * Math.sin(B) * R;

		var position = new THREE.Vector3(x, y, z);
		position.l = L;
		position.b = B;
		position.r = R;
		return position;
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 3.5511810791,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomCallistoOrbit.prototype = Object.create( KMG.Orbit.prototype );





KMG.Saturn = {};
KMG.Saturn.radius = 60330.0;
KMG.Saturn.ascendingNode = 168.8112;
KMG.Saturn.tilt = 28.0817;

KMG.Saturn.computeSaturnElements = function(t) {
	
	
};


KMG.CustomMimasOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomMimasOrbit.prototype = Object.create( KMG.Orbit.prototype );
	
	

KMG.CustomEnceladusOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomEnceladusOrbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.CustomTethysOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomTethysOrbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.CustomDioneOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomDioneOrbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.CustomRheaOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomRheaOrbit.prototype = Object.create( KMG.Orbit.prototype );


KMG.CustomTitanOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 15.94544758,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomTitanOrbit.prototype = Object.create( KMG.Orbit.prototype );
	
	
	
	
	

KMG.CustomHyperionOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomHyperionOrbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.CustomIapetusOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomIapetusOrbit.prototype = Object.create( KMG.Orbit.prototype );



KMG.CustomPhoebeOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomPhoebeOrbit.prototype = Object.create( KMG.Orbit.prototype );







// ftp://ftp.imcce.fr/pub/ephem/planets/pluto95/pluto.doc
KMG.CustomPlutoOrbit = function() {
	KMG.Orbit.call( this );
	
	function distanceAtTime(jd) {
		var pos = positionAtTime(jd);
		return pos.r;
	};
	
	
	function positionAtTime(jd) {
		
	};
	
	return {
		positionAtTime : positionAtTime,
		distanceAtTime : distanceAtTime,
		period : 0,
		epoch : KMG.Util.julianNow()
	};
	
	
};
KMG.CustomPlutoOrbit.prototype = Object.create( KMG.Orbit.prototype );





KMG.CustomOrbits = {};

KMG.CustomOrbits.earth = function() {
	return new KMG.VSOP87Orbit(KMG.VSOPSeries.earth_L, KMG.VSOPSeries.earth_B, KMG.VSOPSeries.earth_R, 365.25);
};

KMG.CustomOrbits.mars = function() {
	return new KMG.VSOP87Orbit(KMG.VSOPSeries.mars_L, KMG.VSOPSeries.mars_B, KMG.VSOPSeries.mars_R, 689.998725);
};


KMG.CustomOrbitProxy = function(orbit)
{
	KMG.Orbit.call( this );
	
	return {
		positionAtTime : orbit.positionAtTime,
		distanceAtTime : orbit.distanceAtTime,
		period : orbit.period,
		epoch : orbit.epoch
	};
	
};
KMG.CustomOrbitProxy.prototype = Object.create( KMG.Orbit.prototype );

KMG.CustomEarthOrbit = function()
{
	return KMG.CustomOrbitProxy.call( this, KMG.CustomOrbits.earth() );
};

KMG.CustomMarsOrbit = function()
{
	return KMG.CustomOrbitProxy.call( this, KMG.CustomOrbits.mars() );
};

/* File: IAURotation.js */

KMG.IAU_SECULAR_TERM_VALID_CENTURIES = 50.0;
KMG.P03LP_VALID_CENTURIES = 5000.0;

KMG.IAURotation = function() {

	var scope = this;
	
	function clampCenturies(t) {
        if (t < -KMG.IAU_SECULAR_TERM_VALID_CENTURIES)
            t = -KMG.IAU_SECULAR_TERM_VALID_CENTURIES;
        else if (t > KMG.IAU_SECULAR_TERM_VALID_CENTURIES)
            t = KMG.IAU_SECULAR_TERM_VALID_CENTURIES;
		return t;
    };
	
	this.julianCentury = function(jd) {
		return (jd - 2451545.0) / 36525.0;
	}
	
	// T = Julian Centuries of 36525 days from epoch
	// d = Julian Days from epoch
	this.calculateOrientation = function(jd) {
		var t = this.julianCentury(jd);
		t = clampCenturies(t);
		jd = jd - 2451545.0;
		
		var result = this.__calculateOrientation(jd, t);
		var ra = result.ra;
		var dec = result.dec;
		
		var node = ra + 90.0;
        var inclination = 90.0 - dec;

		return {
			ra : ra,
			dec : dec,
			node : node,
			inclination : inclination
		};
	
	};
	
	// T = Julian Centuries of 36525 days from epoch
	// d = Julian Days from epoch
	this.computeSiderealRotation = function(jd) {
		var t = this.julianCentury(jd);
		jd = jd - 2451545.0;

		return {
			meridian : this.__computeSiderealRotation(jd, t).meridian
		};
		
	};
	
	this.computeRotationalQuaternion = function(jd, skipMeridian) {
	
		var orientation = this.calculateOrientation(jd);
		var meridian = KMG.Math.clamp(this.computeSiderealRotation(jd).meridian, 360) + 90;
		var nodeAxis = new THREE.Vector3( 1, 0, 0 );
		nodeAxis.rotateY((-orientation.node + 90) * KMG.PI_BY_180);
		
		var inclinationQ = new THREE.Quaternion();
		inclinationQ.setFromAxisAngle( nodeAxis, -orientation.inclination * KMG.PI_BY_180 );
	
		var noMeridian = inclinationQ.clone();
		if (!skipMeridian) {
			var meridianQ = new THREE.Quaternion();
			meridianQ.setFromAxisAngle(new THREE.Vector3( 0, 1, 0 ), meridian * KMG.PI_BY_180);
			inclinationQ.multiply(meridianQ);
		}
		
		/*
		var nodeAxis = new THREE.Vector3( 1, 0, 0 );
		nodeAxis.rotateY((orientation.node) * KMG.PI_BY_180);
		var satelliteQ = new THREE.Quaternion();
		satelliteQ.setFromAxisAngle( nodeAxis, -orientation.inclination * KMG.PI_BY_180 );	
		var meridianQ = new THREE.Quaternion();
		meridianQ.setFromAxisAngle(new THREE.Vector3( 0, 1, 0 ), -meridian * KMG.PI_BY_180);
		satelliteQ.multiply(meridianQ);
		*/
		
		inclinationQ.meridian = meridian - 90;
		inclinationQ.ra = orientation.ra;// * KMG._180_BY_PI;
		inclinationQ.dec = orientation.dec;
		inclinationQ.inclination = orientation.inclination;
		inclinationQ.node = orientation.node;
		inclinationQ.noMeridian = noMeridian;
		//inclinationQ.satelliteQ = satelliteQ;
		return inclinationQ;
	};
	
};





/* File: BaseObject.js */
KMG.BaseObject = function ( ) {
	
	THREE.Object3D.call( this );
	var scope = this;
	
	this.setVisibility = function(visible) {
		this.traverse(function(obj) {
			obj.visible = visible;
		});
	};
	
	this.setShadowInteraction = function(enable) {
		this.traverse(function(obj) {
			obj.castShadow = enable;
			obj.receiveShadow = enable;
		});
	};
	
};
KMG.BaseObject.prototype = Object.create( THREE.Object3D.prototype );





/* File: SurfaceObject.js */

KMG.DefaultSurfaceObjectConfig = {
	texture : KMG.textures[1].name,
	surfaceDetail : 0.0,
	elevationScale : 0,
	shininess : 0,
	diffuseIntensity : 170,
	specularIntensity : 4,
	ambientIntensity : 0,
	emissiveIntensity : 180,
	flattening : 0.0033528,
	axialTilt : 0.0,
	surfaceColorMode : "Normal",
	surfaceHue : 0.5,
	surfaceSaturation : 0.0,
	surfaceLightness : 0.75,
	surfaceWrapRGB : 0.031,
	surfaceRotation : 0.0,
	scaleSurface : 1.0
};

KMG.SurfaceObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultSurfaceObjectConfig);
	this.context = context;
	var scope = this;
	
	//var geometry = new THREE.SphereGeometry( scope.config.radius, 256, 256 );
	var geometry = new THREE.IcosahedronGeometry( scope.config.radius, 6 );
	geometry.computeTangents();
	
	var ambient = 0x000000, diffuse = 0xFFFFFF, specular = 0x040404, shininess = 15;

	var shader = KMG.ExtendedNormalMapShader;
	var uniforms = THREE.UniformsUtils.clone( shader.uniforms );
	
	var parameters = { 
		fragmentShader: shader.fragmentShader
		, vertexShader: shader.vertexShader
		, uniforms: uniforms
		, lights: true
		, fog : true
		, shading : THREE.SmoothShading
		, alphaTest : 0.2
		//, wireframe : true
	};
	var material = new THREE.ShaderMaterial( parameters );
	material.wrapAround = true;

	var mesh = new THREE.Mesh( geometry, material );
	mesh.position = new THREE.Vector3( 0, 0, 0 );
	this.position = new THREE.Vector3( 0, 0, 0 ); 
	
	this.props = {
		mesh : mesh,
		material : material,
		geometry : geometry
	};
	this.rotation.y = 180 * (Math.PI / 180);
	
	this.add(mesh);
	
	var lastCapture = (new Date()).getTime();
	
	
	this.update = function()
	{
		var texDefinition = KMG.TextureMap.getTextureDefinitionByName(this.config.texture);
	
		if (!this.context.configChanged) {
			return;
		}
	
		var tDiffuse = (texDefinition.texture) ? KMG.TextureMap.loadTexture(texDefinition.texture) : null;
		var tNormal = (texDefinition.normalMap) ? KMG.TextureMap.loadTexture(texDefinition.normalMap) : null;
		var tSpecular = (texDefinition.specularMap) ? KMG.TextureMap.loadTexture(texDefinition.specularMap) : null;
		var tDisplacement = (texDefinition.bumpMap) ? KMG.TextureMap.loadTexture(texDefinition.bumpMap) : null;
		
		var diffuseIntensity = KMG.Util.intensityToWhiteColor(this.config.diffuseIntensity);
		var specularIntensity = KMG.Util.intensityToWhiteColor(this.config.specularIntensity);
		var ambientIntensity = KMG.Util.intensityToWhiteColor(this.config.ambientIntensity);
		

		var diffuse = KMG.Util.intensityToWhiteColor(this.config.diffuseIntensity);
		var specular = KMG.Util.intensityToWhiteColor(this.config.specularIntensity);
		var ambient = KMG.Util.intensityToWhiteColor(this.config.ambientIntensity);
		var emissive = KMG.Util.intensityToWhiteColor(this.config.emissiveIntensity);

		uniforms[ "enableDisplacement" ].value = true;
		
		
		uniforms[ "enableSpecular" ].value = (tSpecular != null);
		
		uniforms[ "tDiffuse" ].value = tDiffuse;
		uniforms[ "tNormal" ].value = tNormal;
		uniforms[ "tSpecular" ].value = tSpecular;
		uniforms["tDisplacement"].value = tDisplacement;
		
		if (tNormal) {
			uniforms[ "uNormalScale" ].value.set( this.config.surfaceDetail, this.config.surfaceDetail );
		} else {
			uniforms[ "uNormalScale" ].value.set( 0, 0 );
		}
		uniforms["uDisplacementScale"].value = this.config.elevationScale;
	
		
		var hslColor = new THREE.Color(0xFFFFFF);
		hslColor.setHSL(this.config.surfaceHue, this.config.surfaceSaturation, this.config.surfaceLightness);
		
		var hslEmissiveColor = new THREE.Color(0xFFFFFF);
		hslEmissiveColor.setHSL(this.config.surfaceHue, this.config.surfaceSaturation, this.config.surfaceLightness);
		hslEmissiveColor.multiplyScalar( this.config.emissiveIntensity / 255.0 );
		
		
		uniforms["usingDirectionalLighting"].value = (this.config.lightingType === "Directional");
		
		uniforms[ "uDiffuseColor" ].value = hslColor;
		uniforms[ "uSpecularColor" ].value = specular ;
		uniforms[ "uAmbientColor" ].value = ambient ;
		uniforms[ "uEmissiveColor" ].value = hslEmissiveColor;
		uniforms[ "uShininess" ].value = this.config.shininess;
		
		uniforms[ "wrapRGB" ].value = new THREE.Vector3(this.config.surfaceWrapRGB, this.config.surfaceWrapRGB, this.config.surfaceWrapRGB); 
		
		uniforms["enableFog"].value = this.config.displayAtmosphere;
		
		mesh.castShadow = scope.config.shadows;
		mesh.receiveShadow = scope.config.shadows;
		
		mesh.rotation.set(0, this.config.surfaceRotation*(Math.PI/180.0), 0.0);
		
		// Doesn't really work this way, but visually, it's close enough
		this.scale.set(this.config.scaleSurface, this.config.scaleSurface - this.config.flattening, this.config.scaleSurface);
		//this.scale.y = scaleSurface - this.config.flattening;
		
		//
		this.rotation.z = -this.config.axialTilt * (Math.PI/180);
	};

};
KMG.SurfaceObject.prototype = Object.create( KMG.BaseObject.prototype );


/* File: CatalogStarsObject.js */

KMG.StarUtil = {};
KMG.StarUtil.colorForSpectralClass = function(code) {
	if (code == "O")
		return 0x9db4ff;
	else if (code == "B")
		return 0xaabfff;
	else if (code == "A")
		return 0xcad8ff;
	else if (code == "F")
		return 0xfbf8ff;
	else if (code == "G")
		return 0xfff4e8;
	else if (code == "K")
		return 0xffddb4;
	else if (code == "M")
		return 0xffbd6f;
	else if (code == "L")
		return 0xf84235;
	else if (code == "T")
		return 0xba3059;
	else if (code == "Y")
		return 0x605170;
	else // Including D,C,W,S,N,P for now until I get some concrete colors on these
		return 0xFFFFFF;
}



KMG.StarParticlesShader = {

	uniforms: THREE.UniformsUtils.merge( [
		{
			"tParticle"	   : { type: "t", value: null },
			"vAlpha": { type: "f", value: 0.15 },
			"vSizeMultiplier": { type: "f", value:3.5 },
			"color" : { type: "v3", value: new THREE.Vector3( 1, 1, 1 ) }
		}]),

	vertexShader: [
		'attribute float alpha;',
		'varying float vAlpha;',
		'uniform float vSizeMultiplier;',
		'void main() {',
		'	vAlpha = alpha;',
		'	vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );',
		'	gl_PointSize = (vSizeMultiplier * alpha);',
		'	gl_Position = projectionMatrix * mvPosition;',
		'}'
	].join("\n"),

	fragmentShader: [
		"uniform sampler2D tParticle;",
		'uniform vec3 color;',
		'varying float vAlpha;',
		'varying float vSizeMultiplier;',
		'void main() {',
		'	gl_FragColor = vec4( color, vAlpha );',
		'	gl_FragColor = gl_FragColor * texture2D( tParticle, vec2( gl_PointCoord.x, 1.0 - gl_PointCoord.y ) );',
		'}'
	].join("\n")

};


KMG.DefaultSpectralTypeStarParticleSystemOptions = {
	texture : '/img/star_particle.png',
	radius : 90000.0,
	sizeMultiplier : 6.5
};

KMG.SpectralTypeStarParticleSystem = function(context, config, typeCode) {
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultSpectralTypeStarParticleSystemOptions);
	this.context = context;
	var scope = this;
	
	var particleTexture = KMG.TextureMap.loadTexture(config.texture);
	particleTexture.wrapS = particleTexture.wrapT = THREE.ClampToEdgeWrapping;
	particleTexture.format = THREE.RGBAFormat;
	particleTexture.needsUpdate = true;
	
	var geometry = new THREE.Geometry();

	var shader = KMG.StarParticlesShader;
	var uniforms = THREE.UniformsUtils.clone( shader.uniforms );

	uniforms[ "color" ].value = KMG.Util.rgbToArray(KMG.StarUtil.colorForSpectralClass(typeCode));
	uniforms[ "vSizeMultiplier" ].value = config.sizeMultiplier;
	uniforms[ "tParticle" ].value = particleTexture;
	
	var attributes = {
        alpha: { type: 'f', value: [] }
    };
	
	var material = new THREE.ShaderMaterial( {
        uniforms:       uniforms,
        attributes:     attributes,
        vertexShader: shader.vertexShader,
		fragmentShader: shader.fragmentShader,
        transparent:    true
    });
	

	this.addStar = function(vertex, visualMagnitude) {	
		
		
		var alpha = 1.0 - ((visualMagnitude + 1.46) / 9.42);
		alpha = (alpha * 0.85) + (0.15);
		
		attributes.alpha.value.push(alpha);
		
		
		vertex.magnitude = visualMagnitude;
		geometry.vertices.push( vertex );
	};
	
	this.build = function() {
		var particles = new THREE.ParticleSystem( geometry, material );
		scope.add(particles);
	
	};
	
	
	this.update = function()
	{
		if (!this.context.configChanged)
			return;
			
		uniforms[ "vSizeMultiplier" ].value = this.config.sizeMultiplier;
		
	};
};
KMG.SpectralTypeStarParticleSystem.prototype = Object.create( KMG.BaseObject.prototype );




KMG.DefaultCatalogStarsObjectOptions = {
	radius : 90000.0,
	namesVisible : false
};
KMG.DefaultCatalogStarsObjectOptions = KMG.Util.extend(KMG.DefaultCatalogStarsObjectOptions, KMG.DefaultSpectralTypeStarParticleSystemOptions);


KMG.CatalogStarsObject = function ( context, config, onLoaded ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultCatalogStarsObjectOptions);
	this.context = context;
	var scope = this;
	
	var particleTexture = KMG.TextureMap.loadTexture('/img/star_particle.png');
	particleTexture.wrapS = particleTexture.wrapT = THREE.ClampToEdgeWrapping;
	particleTexture.format = THREE.RGBAFormat;
	
	var names = [];
	
	var spectralTypes = {};
	this.add(spectralTypes["O"] = new KMG.SpectralTypeStarParticleSystem(context, config, "O"));
	this.add(spectralTypes["B"] = new KMG.SpectralTypeStarParticleSystem(context, config, "B"));
	this.add(spectralTypes["A"] = new KMG.SpectralTypeStarParticleSystem(context, config, "A"));
	this.add(spectralTypes["F"] = new KMG.SpectralTypeStarParticleSystem(context, config, "F"));
	this.add(spectralTypes["G"] = new KMG.SpectralTypeStarParticleSystem(context, config, "G"));
	this.add(spectralTypes["K"] = new KMG.SpectralTypeStarParticleSystem(context, config, "K"));
	this.add(spectralTypes["M"] = new KMG.SpectralTypeStarParticleSystem(context, config, "M"));
	this.add(spectralTypes["L"] = new KMG.SpectralTypeStarParticleSystem(context, config, "L"));
	this.add(spectralTypes["T"] = new KMG.SpectralTypeStarParticleSystem(context, config, "T"));
	this.add(spectralTypes["Y"] = new KMG.SpectralTypeStarParticleSystem(context, config, "Y"));
	
	
	this.add(spectralTypes["D"] = new KMG.SpectralTypeStarParticleSystem(context, config, "D"));
	this.add(spectralTypes["C"] = new KMG.SpectralTypeStarParticleSystem(context, config, "C"));
	this.add(spectralTypes["W"] = new KMG.SpectralTypeStarParticleSystem(context, config, "W"));
	this.add(spectralTypes["S"] = new KMG.SpectralTypeStarParticleSystem(context, config, "S"));
	this.add(spectralTypes["N"] = new KMG.SpectralTypeStarParticleSystem(context, config, "N"));
	this.add(spectralTypes["P"] = new KMG.SpectralTypeStarParticleSystem(context, config, "P"));

	function createSystemForSpectralClass(classCode) {
	
	
	
	};
	
	
	function createStarLabel(vertex, name) {
		
		var text = new KMG.BillBoardTextObject(context, KMG.Util.replaceWithGreekLettersAbbreviated(name), {font : "8px sans-serif", fillStyle : "rgba(200,200,200,0.95)"});
		text.position = vertex.clone();
		scope.add(text);
		names.push(text);
		
	};
	
	
	$.ajax({
		url: "/api/stars/list/",
		dataType: "json",
		error: function( jqXHR, textStatus, errorThrown ) {
			console.warn("Error: " + errorThrown);
		
		},
		success: function(data, textStatus, jqxhr) {
			
		}
	}).done(function(data) {

		
		var lbls = 0;
		for (var i = 0; i < data.length; i++) {
			
			var vMag = data[i].Vmag;
			var l = data[i].eclLon;
			var b = data[i].eclLat;
			var specClass = data[i].SpClass.toUpperCase();
			var name = data[i].name;
			
			if (spectralTypes[specClass]) {

				var vertex = KMG.Math.getPoint3D(l, b, config.radius);
				
				spectralTypes[specClass].addStar(vertex, vMag);
				
				if (name.length > 0 && vMag < 6.0) {
					
					name = name.replace(/^[0-9]+/, "");
					if (vMag > 3) {
						name = name.replace(/[ ][A-Za-z]{3}$/i, "");
					}
					
					if (name.length > 0) {
						createStarLabel(vertex, name);
						lbls++;
					}
				}
				
				
			} else {
				console.warn("No particle system for spectral type " + specClass);
			}
		}
		
		console.info("Added " + lbls + " star labels");
		for (var key in spectralTypes) {
			spectralTypes[key].build();
		}
		
		scope.setTextVisibility(config.namesVisible);
		
		if (onLoaded) {
			onLoaded(scope);
		}
		
	});
	
	this.setTextVisibility = function(visible) {
		
		for (var i = 0; i < names.length; i++) {
			names[i].setVisibility(visible);
		}
		
	};
	
	this.update = function()
	{
		if (!this.context.configChanged)
			return;
		
		this.setTextVisibility(config.namesVisible);
		
		for (var key in spectralTypes) {
			spectralTypes[key].update();
		}
	};
};
KMG.CatalogStarsObject.prototype = Object.create( KMG.BaseObject.prototype );

/* File: ConstellationLines.js */
KMG.DefaultConstellationLinesConfig = {
	
	color : 0x12110C,
	radius : 90000.0,
	opacity : 0.45,
	lineThickness : 1.0
	
};

KMG.ConstellationLines = function ( context, config, onLoaded ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultConstellationLinesConfig);
	this.context = context;
	var scope = this;

	
	function buildPoint(point) {
		
		var coords = KMG.Math.convertEquatorialToEcliptic(point[0] * 15, point[1]);
		var vertex = KMG.Math.getPoint3D(coords.l, coords.b, config.radius);
		
		return vertex;
	}
	
	function buildPath(path) {
		
		var geometry = new THREE.Geometry();
		
		for (var i = 0; i < path.length; i++) {
			var point = path[i];
			var vertex = buildPoint(point);
			
			geometry.vertices.push( vertex );
			
		}
		
		var material = new THREE.LineBasicMaterial( { opacity: config.opacity, transparent : true, fog : false, color : config.color, linewidth:config.lineThickness } );
		var line = new THREE.Line( geometry,  material);

		return line;
	}
	
	function buildConstellation(constellation) {
		
		var const3d = new THREE.Object3D();
		
		for (var i = 0; i < constellation.length; i++) {
			var path = constellation[i];
			var path3d = buildPath(path);
			const3d.add(path3d);
		}
		
		return const3d;
		
	}
	
	
	$.ajax({
		url: "/api/constellations/list/",
		dataType: "json",
		error: function( jqXHR, textStatus, errorThrown ) {
			console.warn("Error: " + errorThrown);
		
		},
		success: function(data, textStatus, jqxhr) {
			
		}
	}).done(function(data) {

		console.info("Adding " + data.length + " constellations");
		
		for (var i = 0; i < data.length; i++) {
			
			scope.add(buildConstellation(data[i]));

		}
		
		if (onLoaded) {
			onLoaded(scope);
		}
	
	});
	
	

	
	this.update = function() {
		
		
	};
};
KMG.ConstellationLines.prototype = Object.create( KMG.BaseObject.prototype );





/* File: BackgroundImageSphereObject.js */

KMG.BackgroundImageSphereObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	this.config = config;
	this.context = context;
	var scope = this;
	
	var geometry = new THREE.SphereGeometry( 10000000, 64, 32 );
	geometry.computeTangents();
	
	
	var texDefinition = KMG.TextureMap.getBackgroundDefinitionByName(scope.config.backgroundImage);
	var tDiffuse = (texDefinition.texture != null) ? KMG.TextureMap.loadTexture(texDefinition.texture) : null;

	var material = new THREE.MeshBasicMaterial({
								shading		: THREE.SmoothShading,
								map			: tDiffuse,
								fog			: false,
								side		: THREE.BackSide,
								depthWrite  : false
							});
	
	
	var mesh = new THREE.Mesh( geometry, material );
	
	mesh.position = new THREE.Vector3( 0, 0, 0 );

	// Create Stuff

	this.mesh = mesh;
	this.material = material;
	this.geometry = geometry;
	
	this.add(mesh);
	
	this.update = function()
	{
		if (!this.context.configChanged) 
			return;
		
		
		var texDefinition = KMG.TextureMap.getBackgroundDefinitionByName(scope.config.backgroundImage);
		var tDiffuse = (texDefinition.texture) ? KMG.TextureMap.loadTexture(texDefinition.texture) : null;

		this.material.map = tDiffuse;
		
	};
};
KMG.BackgroundImageSphereObject.prototype = Object.create( KMG.BaseObject.prototype );


/* File: BackgroundImageSpriteObject.js */


KMG.BackgroundImageSpriteObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	this.config = config;
	this.context = context;
	var scope = this;

	
	var texDefinition = KMG.TextureMap.getBackgroundDefinitionByName(scope.config.backgroundImage);
	var tDiffuse = (texDefinition.texture != null) ? KMG.TextureMap.loadTexture(texDefinition.texture) : null;

	var material = new THREE.SpriteMaterial({
								shading		: THREE.SmoothShading,
								map			: tDiffuse,
								fog			: false,
								useScreenCoordinates : true,
								alignment: THREE.SpriteAlignment.topLeft,
								blending: THREE.AdditiveBlending
							});
	
	var sprite = new THREE.Sprite( material );
	sprite.position.set( 0, 0, -1000 );
	sprite.scale.set( 1000, 1000, 1 );

	this.sprite = sprite;
	this.material = material;
	
	this.add(sprite);
	
	this.update = function()
	{

		if (this.config.backgroundImageFitType === "stretch") {
			this.sprite.scale.set(this.context.containerWidth, this.context.containerHeight, 1);
		}
		if (!this.context.configChanged) {
			return;
		}
		var texDefinition = KMG.TextureMap.getBackgroundDefinitionByName(scope.config.backgroundImage);
		
		
		var tDiffuse = null;
		if (texDefinition.name = texDefinition.texture){
			 KMG.TextureMap.loadTexture(texDefinition.texture);
		}
		
		this.material.map = tDiffuse;
	};
};
KMG.BackgroundImageSpriteObject.prototype = Object.create( KMG.BaseObject.prototype );

/* File: BackgroundObject.js */

KMG.DefaultBackgroundConfig = {
	backgroundType : 'stars',
	backgroundImage : 'Starfield',
	backgroundImageType : 'flat',
	backgroundImageFitType : 'stretch',
	starQuantity : 3.5, // 0 - 10
	noStars : false
};

KMG.BackgroundObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultBackgroundConfig);
	this.context = context;
	var scope = this;
	
	
	this.starsConfig = {
		alphaMultiplier : 2.0
	};
	// Create Stuff
	if (!config.noStars) {
		this.stars = new KMG.CatalogStarsObject(context, this.starsConfig);
		this.add(this.stars);
	}
	
	this.imageSphere = new KMG.BackgroundImageSphereObject(context, config);
	this.imageFlat = new KMG.BackgroundImageSpriteObject(context, config);
	
	
	this.add(this.imageSphere);
	this.add(this.imageFlat);
	
	
	this.setShadowInteraction(false);
	
	this.update = function()
	{
	
		this.starsConfig.sizeMultiplier = scope.config.starQuantity;
	
		//if (scope.config.backgroundType === "stars") {
		if (this.stars) {
			this.stars.update();
		}
		//}
		
		if (scope.config.backgroundType === "image" && scope.config.backgroundImageType === "sphere") {
			this.imageSphere.update();
		}
		
		if (scope.config.backgroundType === "image" && scope.config.backgroundImageType === "flat") {
			this.imageFlat.update();
		}
		
		if (!this.context.configChanged) 
			return;
		
		// Handled within StarGroupObject object
		//this.stars.traverse(function(obj) {
		//
		//});
		
		this.imageSphere.traverse(function(obj) {
			obj.visible = (scope.config.backgroundType === "image" && scope.config.backgroundImageType === "sphere");
		});
		
		this.imageFlat.traverse(function(obj) {
			obj.visible = (scope.config.backgroundType === "image" && scope.config.backgroundImageType === "flat");
		});
	};
};
KMG.BackgroundObject.prototype = Object.create( KMG.BaseObject.prototype );

/* File: LocalStarObject.js */

KMG.DefaultLocalStarConfig = {
	localStarDistance : 1.0,
	displayLocalStar : true,
	localStarTexture : KMG.starFlares[0].name,
	localStarColor : [ 255, 255, 255 ],
	starColorAffectsPlanetLighting : true
};

KMG.LocalStarObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultLocalStarConfig);
	this.context = context;
	var scope = this;
	
	
	var star = new THREE.Object3D();

	var material = new THREE.SpriteMaterial({fog : false
											, color : 0xFFFFFF
											, sizeAttenuation : false
											, transparent : true
											, blending : THREE.AdditiveBlending
											, useScreenCoordinates: false
											, depthWrite: false
											, depthTest: true
											});

	var sprite = new THREE.Sprite(material);

	star.add(sprite);
	
	this.add(star);
		
	function getPosition(type, direction) {
		var position = null;
		
		if (type === "Directional") {
			position = new THREE.Vector3(-5000.0, 0, 0);
			position.rotateY(direction*(Math.PI/180.0));
		} else {
			position = new THREE.Vector3(0, 0, 0);
		}

		return position;
	}
	
	this.update = function()
	{
		//if (!this.context.configChanged) 
		//	return;
		
		var texDefinition = KMG.TextureMap.getFlareDefinitionByName(scope.config.localStarTexture);
		var tDiffuse = (texDefinition.texture != null) ? KMG.TextureMap.loadTexture(texDefinition.texture) : null;
		tDiffuse.wrapS = tDiffuse.wrapT = THREE.ClampToEdgeWrapping;
		tDiffuse.format = THREE.RGBAFormat;
		tDiffuse.needsUpdate = true;

		var starColor = KMG.Util.arrayToColor(config.localStarColor);
		
		
		material.map = tDiffuse;
		material.color = starColor;
		
		sprite.position = getPosition(this.config.lightingType, this.config.sunlightDirection);

		sprite.scale.set( 100 * this.config.localStarDistance, 100 * this.config.localStarDistance, 1 );
		
		
		this.traverse(function(obj) {
			obj.visible = scope.config.displayLocalStar;
		});
		
		sprite.updateMatrix();
		star.updateMatrix();
	};
	this.update();
};
KMG.LocalStarObject.prototype = Object.create( KMG.BaseObject.prototype );



/* File: LensFlareObject.js */

KMG.DefaultLensFlareConfig = {
	lensFlareEnabled : false
};

/**
 * Shamelessly copied from http://threejs.org/examples/webgl_lensflares.html
 */
KMG.LensFlareObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultLensFlareConfig);
	this.context = context;
	var scope = this;

	
	// Create Stuff
	// See http://planetmaker.wthr.us/img/LICENSE.txt
	var textureFlare1 = KMG.TextureMap.loadTexture( "/img/lensflare1.png" );
	var textureFlare2 = KMG.TextureMap.loadTexture( "/img/lensflare2.png" );
	var textureFlare3 = KMG.TextureMap.loadTexture( "/img/lensflare3.png" );
	
	var lensFlare = new THREE.LensFlare( textureFlare1, 1000, 1.0, THREE.AdditiveBlending, 0xFFFFFF );
	
	lensFlare.add( textureFlare2, 512, 0.0, THREE.AdditiveBlending );
	lensFlare.add( textureFlare3, 60, 0.6, THREE.AdditiveBlending );
	lensFlare.add( textureFlare3, 70, 0.7, THREE.AdditiveBlending );
	lensFlare.add( textureFlare3, 120, 0.9, THREE.AdditiveBlending );
	lensFlare.add( textureFlare3, 70, 1.0, THREE.AdditiveBlending );

	this.lensFlare = lensFlare;
	
	this.add(lensFlare);
	
	this.update = function()
	{
		if (!this.context.configChanged) 
			return;
		
		
		this.lensFlare.traverse(function(obj) {
			obj.visible = scope.config.lensFlareEnabled;
		});
		
		if (!this.config.lensFlareEnabled) {
			this.remove(this.lensFlare);
			return;
		} else {
			this.add(lensFlare);
		}

		
		if (this.config.lightingType === "Directional") {
			this.lensFlare.position = new THREE.Vector3(-5000.0, 0, 0);
			this.lensFlare.position.rotateY(scope.config.sunlightDirection*(Math.PI/180.0));
		} else {
			this.lensFlare.position = new THREE.Vector3(0, 0, 0);
		}
		
		this.lensFlare.updateMatrix();

	};
};
KMG.LensFlareObject.prototype = Object.create( KMG.BaseObject.prototype );

/* File: CometObject.js */

KMG.DefaultCometObjectConfig = {
	lookingTowards : null
};

KMG.CometObject = function ( context, config, ephemeris, tickController, centerObject, lookingTowards ) {
	
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultCometObjectConfig);
	var scope = this;
	
	this.lookingTowards = (config.lookingTowards) ? config.lookingTowards : { position: new THREE.Vector3(0, 0, 0) };
	
	var cometTexture = KMG.TextureMap.loadTexture("/img/basic-comet-3-trans.png");
	cometTexture.wrapS = cometTexture.wrapT = THREE.ClampToEdgeWrapping;
	cometTexture.format = THREE.RGBAFormat;
	cometTexture.needsUpdate = true;
	
	function createFaceMesh(rotate) {

		var geometry = new THREE.PlaneGeometry(400, 100, 10, 10);
		var materialOptions = { color: 0xFFFFFF
							, ambient : 0xFFFFFF
							, emissive : 0xAAAAAA
							, shading : THREE.NoShading
							, map : cometTexture
							, transparent : true
							, side: THREE.DoubleSide
							, blending : THREE.AdditiveBlending
							, depthWrite: false
							, depthTest: false
							};
		var mesh = new THREE.Mesh(geometry, new THREE.MeshLambertMaterial( materialOptions ));
		
		mesh.rotation.set(rotate, KMG.RAD_90, 0, 'YXZ');
		mesh.position.z -= 160;

		return mesh;
	}

	this.add(createFaceMesh(0));
	this.add(createFaceMesh(-KMG.RAD_90));
	this.add(createFaceMesh(-KMG.RAD_45));
	this.add(createFaceMesh(KMG.RAD_45));

	this.update = function() {

		if (this.lookingTowards && this.lookingTowards.position) {
			var lookAt = this.lookingTowards.position.clone();
			this.lookAt( lookAt );
		}
		
	};
};
KMG.CometObject.prototype = Object.create( KMG.BaseObject.prototype );
	

/* File: BillboardTextObject.js */


KMG.DefaultBillBoardTextConfig = {
	textColor : 0xFFFFFF,
	fillStyle : "rgba(255,255,255,0.95)",
	font : "10px sans-serif"

};

KMG.BillBoardTextObject = function(context, text, config) {
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultBillBoardTextConfig);
	this.context = context;
	var _text = text;
	var scope = this;
	
	function createTextTexture(text, width, height) {
		var canvas1 = document.createElement('canvas');
		var context1 = canvas1.getContext('2d');
		canvas1.width = width;
		canvas1.height = height;
		context1.font = scope.config.font;
		context1.fillStyle = scope.config.fillStyle;

		context1.textAlign="center";				
		context1.fillText(text
							, canvas1.width / 2
							, canvas1.height / 2 + 20);	
		context1.fill();

		var texture1 = new THREE.Texture(canvas1) 
		texture1.needsUpdate = true;

		return texture1;
	}
	
	
	
	var geometry = new THREE.Geometry();

	var material = new THREE.ParticleBasicMaterial( { map: createTextTexture(_text, 100, 100)
													, color: this.config.textColor
													, size: 100
													, fog : false
													, sizeAttenuation : false
													, transparent : true
													, opacity : 1.0
													, blending : THREE.AdditiveBlending
													, depthWrite: false
													, depthTest: false
													} );
	var vertex = new THREE.Vector3(0, 0, 0);
	geometry.vertices.push( vertex );
	
	var particles = new THREE.ParticleSystem( geometry, material );
	this.add(particles);
	
	
	this.setText = function(text) {
		_text = text;
		material.map = createTextTexture(_text, 100, 100);
	};

	this.update = function()
	{
		if (!this.context.configChanged)
			return;
			

	};
};
KMG.BillBoardTextObject.prototype = Object.create( KMG.BaseObject.prototype );
	
	
/* File: TexturedSphereObject.js */

KMG.MaterialPhong = 1;
KMG.MaterialLambert = 2;

KMG.DefaultTexturedSphereOptions = {
	texture : "Earth - Blue Marble",
	scale : 1,
	radius : 200,
	flattening : 0,
	ambient : 0x888888,
	color : 0xDDDDDD,
	emissive : 0x000000,
	material : KMG.MaterialLambert,
	specular : 0x444444,
	shadows : true,
	slices : 32,
	shading : true,
	transparent : false

};

/** A simpler sphere for small planets or moons
 *
 */
KMG.TexturedSphereObject = function(context, config) {
	KMG.BaseObject.call( this );
	this.config = config = KMG.Util.extend(config, KMG.DefaultTexturedSphereOptions);
	this.context = context;
	var scope = this;
	
	var geometry = new THREE.EllipsoidGeometry( this.config.radius, this.config.flattening, this.config.slices, this.config.slices );
	
	var texDefinition = KMG.TextureMap.getTextureDefinitionByName(this.config.texture);
	
	var tDiffuse = (texDefinition.texture) ? KMG.TextureMap.loadTexture(texDefinition.texture) : null;
	if (config.transparent) {
		tDiffuse.format = THREE.RGBAFormat;
	}
	
	var material;
	
	var shading = (config.shading) ? THREE.SmoothShading : THREE.NoShading;
	
	if (this.config.material == KMG.MaterialLambert) {
		material = new THREE.MeshLambertMaterial({
									ambient		: new THREE.Color(this.config.ambient),
									color		: new THREE.Color(this.config.color),
									emissive	: new THREE.Color(this.config.emissive),
									shading		: shading,
									map			: tDiffuse,
									fog			: this.config.fog,
									transparent : this.config.transparent
								});
	} else if (this.config.material == KMG.MaterialPhong) {
		var tSpecular = (texDefinition.specularMap) ? KMG.TextureMap.loadTexture(texDefinition.specularMap) : null;
		material = new THREE.MeshPhongMaterial({
									ambient		: new THREE.Color(this.config.ambient),
									color		: new THREE.Color(this.config.color),
									emissive	: new THREE.Color(this.config.emissive),
									specular	: new THREE.Color(this.config.specular),
									shading		: shading,
									map			: tDiffuse,
									specularMap	: tSpecular,
									fog			: this.config.fog,
									transparent : this.config.transparent
								});
	} 
								
	var mesh = new THREE.Mesh( geometry, material );
	mesh.position = new THREE.Vector3( 0, 0, 0 );
	
	mesh.castShadow = config.shadows;
	mesh.receiveShadow = config.shadows;
	
	
	this.add(mesh);
		
	
	this.sphereMesh = mesh;
	
	this.update = function()
	{
		if (!this.context.configChanged)
			return;
		
		this.scale.set(this.config.scale, this.config.scale, this.config.scale);
		mesh.castShadow = scope.config.shadows;
		mesh.receiveShadow = scope.config.shadows;
	};
};
KMG.TexturedSphereObject.prototype = Object.create( KMG.BaseObject.prototype );
	

/* File: BasicAsteroidBeltObject.js */

KMG.DefaultBasicAsteroidBeltConfig = {

	innerRadius : 260.0,
	outterRadius : 400.0,
	quantity : 2000,
	sizeAttenuation : true,
	size : 2,
	color : 0xFFFFFF,
	hue : 0.5,
	saturation : 0.0,
	lightness : 0.75,
	targetObject : 0
};


KMG.BasicAsteroidBeltObject = function ( context, config ) {
	
	KMG.BaseObject.call( this );
	
	
	this.config = config = KMG.Util.extend(config, KMG.DefaultBasicAsteroidBeltConfig);
	this.context = context;
	var scope = this;
	
	var lastBeltObject = null;
	
	function createBelt3dObject()
	{
		var innerRadius = config.innerRadius;
		var outterRadius = config.outterRadius;
		var quantity = config.quantity;
		
		var belt = new THREE.Object3D();
		
		var particleTexture = KMG.TextureMap.loadTexture('img/star_particle.png');
		particleTexture.wrapS = particleTexture.wrapT = THREE.ClampToEdgeWrapping;
		particleTexture.format = THREE.RGBAFormat;
		
		var particles, geometry, material, parameters, i, h, color;
		
		var hslColor = new THREE.Color(config.color);
		hslColor.setHSL(config.hue, config.saturation, config.lightness);
		
		geometry = new THREE.Geometry();
		material = new THREE.ParticleBasicMaterial( { color: hslColor
														, size: config.size
														, map: particleTexture
														, fog : false
														, sizeAttenuation : config.sizeAttenuation
														, transparent : true
														, blending : THREE.AdditiveBlending
														} );

		for ( i = 0; i < quantity; i ++ ) {
			var u = Math.random() * 360.0;
			var v = Math.random() * 180.0;
			
			var radius = innerRadius + (outterRadius - innerRadius) * Math.random();
			
			var vertex = new THREE.Vector3(0, 0, radius);
			vertex.rotateY(Math.random() * 360.0 *(Math.PI/180.0));
			geometry.vertices.push( vertex );
		
		}
		
		particles = new THREE.ParticleSystem( geometry, material );
			
		belt.add( particles );
	
		return belt;
	}
	
	
	
	
	
	this.update = function()
	{
		if (!this.context.configChanged) 
			return;
		
		if (lastBeltObject) {
			this.remove(lastBeltObject);
		}
		lastBeltObject = createBelt3dObject();
		this.add(lastBeltObject);	
		
		if (scope.config.targetObject) {
			this.position = scope.config.targetObject.position.clone();
		}
	};
	this.update();
};
KMG.BasicAsteroidBeltObject.prototype = Object.create( KMG.BaseObject.prototype );

/* File: StandardComposer.js */

KMG.StraightThroughRenderPass = function(scene, camera) {

	this.enabled = true;
	this.clear = true;
	this.needsSwap = false;
	
	this.render = function( renderer, writeBuffer, readBuffer, delta ) {
		renderer.render( scene, camera );
	}
};

KMG.StandardComposer = function(context, scene, camera, renderer) {

	THREE.EffectComposer.call( this, renderer );
	
	var scope = this;
	
	var renderPass = new KMG.StraightThroughRenderPass( scene, camera );
	this.addPass( renderPass );

	
	this.update = function() {
	
	
	};


};
KMG.StandardComposer.prototype = Object.create( THREE.EffectComposer.prototype );

/* File: FilmPassComposer.js */

KMG.DefaultFilmPassOptions = {
	noiseIntensity : 0.35,
	scanlinesIntensity : 0.75,
	scanlinesCount : 2048,
};


KMG.FilmPassComposer = function(context, scene, camera, renderer, config) {

	THREE.EffectComposer.call( this, renderer );
	this.config = config = KMG.Util.extend((config) ? config : {}, KMG.DefaultFilmPassOptions);
	var scope = this;
	
	var renderPass = new THREE.RenderPass( scene, camera );
	var effectFilm = new THREE.FilmPass( config.noiseIntensity, config.scanlinesIntensity, config.scanlinesCount );
	effectFilm.renderToScreen = true;

	this.addPass( renderPass );
	this.addPass( effectFilm );
	
	
	this.update = function() {
	
	
	};
	

};

KMG.FilmPassComposer.prototype = Object.create( THREE.EffectComposer.prototype );
/* File: BlurPassComposer.js */

KMG.DefaultBlurPassOptions = {
	hBlurAmount : 0.5,
	vBlurAmount : 0.5
};

KMG.BlurPassComposer = function(context, scene, camera, renderer, config) {

	THREE.EffectComposer.call( this, renderer );
	this.config = config = KMG.Util.extend((config) ? config : {}, KMG.DefaultBlurPassOptions);
	var scope = this;
	
	var renderPass = new THREE.RenderPass( scene, camera );
	
	var width = window.innerWidth || 2;
	var height = window.innerHeight || 2;
	
	var effectHBlur = new THREE.ShaderPass( THREE.HorizontalBlurShader );
	var effectVBlur = new THREE.ShaderPass( THREE.VerticalBlurShader );
	effectHBlur.uniforms[ 'h' ].value = config.hBlurAmount / ( width / 2 );
	effectVBlur.uniforms[ 'v' ].value = config.vBlurAmount / ( height / 2 );
	effectVBlur.renderToScreen = true;
	
	this.addPass( renderPass );
	this.addPass( effectHBlur );
	this.addPass( effectVBlur );
	
	this.update = function() {
		
		effectHBlur.uniforms[ 'h' ].value = config.hBlurAmount / ( width / 2 );
		effectVBlur.uniforms[ 'v' ].value = config.vBlurAmount / ( height / 2 );
	
	};
	

};

KMG.BlurPassComposer.prototype = Object.create( THREE.EffectComposer.prototype );
/* File: DynamicEffectsComposer.js */

KMG.DefaultEffectsConfig = {

	enableGodRays : false,
	godRaysIntensity : 0.75,
	
	enableBlur : false,
	blurAmount : 0.5,
	
	enableBloom : false,
	bloomStrength : 0.5,
	
	enableBleach : false,
	bleachAmount : 0.95,
	
	enableSepia : false,
	sepiaAmount : 0.9,
	
	enableFilm : false,
	noiseIntensity : 0.35,
	scanlinesIntensity : 0.75,
	scanlinesCount : 2048,
};


KMG.DynamicEffectsComposer = function(context, scene, camera, secondaryCamera, renderer, config) {

	THREE.EffectComposer.call( this, renderer );
	this.context = context;
	this.config = config = KMG.Util.extend(config, KMG.DefaultEffectsConfig);
	var scope = this;
	
	//var fxaaShader = THREE.FXAAShader;
	//var effectFXAA = new THREE.ShaderPass( fxaaShader );
	
	var shaderBleach = THREE.BleachBypassShader;
	var effectBleach = new THREE.ShaderPass( shaderBleach );
	
	var shaderSepia = THREE.SepiaShader;
	var effectSepia = new THREE.ShaderPass( shaderSepia );
	effectSepia.uniforms[ "amount" ].value = 0.9;
	
	//var godRaysPass = new THREE.GodRaysPass(config, context.primaryScene, context.secondaryScene, camera);
	var renderBackgroundPass = new THREE.RenderPass( [context.secondaryScene], secondaryCamera );
	var renderPass = new THREE.RenderPass( context.primaryScene, camera);
	renderPass.clear = false;
	
	var effectHBlur = new THREE.ShaderPass( THREE.HorizontalBlurShader );
	var effectVBlur = new THREE.ShaderPass( THREE.VerticalBlurShader );
	//var effectDotScreen = new THREE.DotScreenPass( new THREE.Vector2( 0, 0 ), 0.5, 0.8 );
	var effectBloom = new THREE.BloomPass(0.5);
	var effectFilm = new THREE.FilmPass( config.noiseIntensity, config.scanlinesIntensity, config.scanlinesCount);
	effectFilm.renderToScreen = true;
	//effectFXAA.renderToScreen = true;
	
	this.addPass( renderBackgroundPass );
	this.addPass( renderPass );
	//this.addPass( effectFXAA );
	//this.addPass( godRaysPass );

	this.addPass( effectBleach );
	this.addPass( effectBloom );
	this.addPass( effectHBlur );
	this.addPass( effectVBlur );
	this.addPass( effectSepia );
	this.addPass( effectFilm );
	
	
	//this.addPass( effectDotScreen );
	
	this.update = function() {
	
		if (!this.context.configChanged) 
			return;
		
		var width = window.innerWidth || 2;
		var height = window.innerHeight || 2;
		
		//renderPass.enabled = !this.config.enableGodRays;
		//godRaysPass.enabled = this.config.enableGodRays;
		
		//effectFXAA.uniforms["resolution"].value.x = 1 / width;
		//effectFXAA.uniforms["resolution"].value.y = 1 / height;
		
		//godRaysPass.setIntensity(this.config.godRaysIntensity);
		
		effectHBlur.enabled = this.config.enableBlur;
		effectVBlur.enabled = this.config.enableBlur;
		var blurAmount = (this.config.enableBlur) ? this.config.blurAmount : 0;
		effectHBlur.uniforms[ 'h' ].value = blurAmount / ( width / 2 );
		effectVBlur.uniforms[ 'v' ].value = blurAmount / ( height / 2 );
		
		
		effectBloom.enabled = this.config.enableBloom;
		var bloomStrength = (this.config.enableBloom) ? this.config.bloomStrength : 0;
		effectBloom.copyUniforms["opacity"].value = bloomStrength;
		
		effectBleach.enabled = this.config.enableBleach;
		var bleachAmount = (this.config.enableBleach) ? this.config.bleachAmount : 0;
		effectBleach.uniforms[ "opacity" ].value = bleachAmount;
	
		
		effectSepia.enabled = this.config.enableSepia;
		effectSepia.uniforms[ "amount" ].value = (this.config.enableSepia) ? this.config.sepiaAmount : 0;
		
		//effectFilm.enabled = this.config.enableFilm;
		effectFilm.uniforms.nIntensity.value = (this.config.enableFilm) ? this.config.noiseIntensity : 0;
		effectFilm.uniforms.sIntensity.value = (this.config.enableFilm) ? this.config.scanlinesIntensity : 0;
		effectFilm.uniforms.sCount.value = (this.config.enableFilm) ? this.config.scanlinesCount : 0;
	};


};
KMG.DynamicEffectsComposer.prototype = Object.create( THREE.EffectComposer.prototype );

/* File: Engine.js */


KMG.SCENE = { NONE : -1, PRIMARY : 0, SECONDARY : 1 };

KMG.Engine = function ( domElement, config, sceneCallbacks, cameraConfig, view ) {

	this.config = config;
	this.sceneCallbacks = sceneCallbacks;
	this.cameraConfig = cameraConfig;
	this.domElement = ( domElement !== undefined ) ? domElement : document;
	
	
	
	this.context = {
		composer : null,
		container : null, 
		stats : null,
		camera : null, 
		secondaryCamera : null,
		primaryScene : null,
		secondaryScene : null,
		renderer : null,
		controls : null,
		windowHalfX : window.innerWidth / 2,
		windowHalfY : window.innerHeight / 2,
		containerWidth : 0,
		containerHeight : 0,
		objects : [],
		configChanged : true,
		animationID : 0,
		script : null,
		background : null,
		lights : {
			ambient : null,
			primaryDirectional : null,
			secondaryDirectional : null,
			primaryPoint : null,
			secondaryPoint : null
		}
	};
	
	this.animators = [];
	
	this.stopAnimation = false;
	
	// Internals
	var scope = this;
	
	this.reset = function()
	{
	
	};
	
	
	
	function onWindowResize() 
	{
		scope.context.windowHalfX = $("#container").width() / 2;
		scope.context.windowHalfY = $("#container").height() / 2;
		
		scope.context.containerWidth = $("#container").width();
		scope.context.containerHeight = $("#container").height();
		
		scope.context.camera.aspect = $("#container").width() / $("#container").height();
		scope.context.camera.updateProjectionMatrix();

		scope.context.renderer.setSize( $("#container").width(), $("#container").height() );
	
		//scope.context.composer.reset();
	}
	
	function onDocumentMouseMove( event ) 
	{
		scope.context.mouseX = ( event.clientX - scope.context.windowHalfX );
		scope.context.mouseY = ( event.clientY - scope.context.windowHalfY );
	}
	
	this.start = function()
	{
		animate();
	};
	
	this.stop = function()
	{
		this.stopAnimation = true;
	
	};
	
	function checkGlError()
	{
		var err = scope.context.renderer.getContext().getError();
		if (err !== 0) {
			onWebGlException("WebGL Error Code " + err);
		}
	}
	
	function onWebGlException(ex)
	{
		if (scope.sceneCallbacks.webGlErrorCallback) {
			scope.sceneCallbacks.webGlErrorCallback(ex);
		}
	}
	

	
	function animate() 
	{	

		if (scope.stopAnimation) {
			scope.stopAnimation = false;
			
			if (scope.sceneCallbacks.animationStoppedCallback) {
				scope.sceneCallbacks.animationStoppedCallback();
			}
			
			return;
		}
		

		scope.context.animationID = requestAnimationFrame( animate );
		
		if (scope.config.useScript && fireOnFrameHandler()) {
			scope.context.configChanged = true;
		}
		
		
		for (var i = 0; i < scope.animators.length; i++) {
			var animator = scope.animators[i];
			animator.next();
		}

		scope.context.controls.update();
		
		
		if (scope.config.enableFps && scope.context.stats) {
			scope.context.stats.update();
		}
		

		render();
	}
	
	
	function updateDirectionalLight(configChanged, light)
	{
		
		if (scope.config.realtimeSunlight) {
			var positioner = new KMG.SunlightPositioning();
			
			// I'm not sure this is 100% correct... Need more validation...
			var sunlightDateTime = new Date(scope.config.sunlightDate);
			
			var position = positioner.getSunPositionOnDate(sunlightDateTime);
			light.position = position;
		}
	
		if (configChanged) {
			
			if (scope.config.lightingType === "Point") {
				light.intensity = 0.0;
				return;
			} else {
				light.intensity = 2.0;
			}
			
			
		
			var localStarLightColor = null;
			
			if (config.starColorAffectsPlanetLighting) {
				localStarLightColor = KMG.Util.arrayToColor(config.localStarColor);
			} else {
				localStarLightColor = new THREE.Color(0xFFFFFF);
			}
			light.color = localStarLightColor;
			
			if (!scope.config.realtimeSunlight) {
				light.position = new THREE.Vector3(-10000.0, 0, 0);
				light.position.rotateY(scope.config.sunlightDirection*(Math.PI/180.0)).normalize();
			}
				
			light.castShadow = scope.config.shadows;
				
			light.shadowDarkness = scope.config.shadowDarkness;
			light.updateMatrix();
		}
	
		
	}
	
	function updatePointLight(configChanged, light)
	{
		if (configChanged) {
			
			if (scope.config.lightingType === "Directional") {
				light.intensity = 0.0;
				return;
			} else {
				light.intensity = 2.0;
			}
			
			
			
			var localStarLightColor = null;
			
			if (config.starColorAffectsPlanetLighting) {
				localStarLightColor = KMG.Util.arrayToColor(config.localStarColor);
			} else {
				localStarLightColor = new THREE.Color(0xFFFFFF);
			}
			light.color = localStarLightColor;
		}

	}
	
	function updateLights(configChanged)
	{
		updateDirectionalLight(configChanged, scope.context.lights.primaryDirectional);
		updateDirectionalLight(configChanged, scope.context.lights.secondaryDirectional);
		updatePointLight(configChanged, scope.context.lights.primaryPoint);
		updatePointLight(configChanged, scope.context.lights.secondaryPoint);
	}
	
	function updateShadows()
	{
		if (!scope.context.configChanged) {
			return;
		}
		
		scope.context.renderer.shadowMapEnabled = scope.config.shadows;
		scope.context.renderer.shadowMapAutoUpdate = scope.config.shadows;
		
		if (!scope.config.shadows) {
			for (var i = 0; i < scope.context.lights.length; i++) {
				scope.context.renderer.clearTarget( scope.context.lights[i].shadowMap );
			}
		}
	}
	
	function areEffectsEnabled() {
		return (scope.config.enableBlur
			|| scope.config.enableBloom
			|| scope.config.enableBleach
			|| scope.config.enableFilm
			|| scope.config.enableSepia);
	}
	
	function render() 
	{	
	
		if (scope.config.ringAngle + scope.config.axialTilt > 0.0) {
			renderer.shadowMapCullFace = THREE.CullFaceFront;
		} else {
			renderer.shadowMapCullFace = THREE.CullFaceBack;
		}
		
		updateShadows();
		
		scope.context.containerWidth = $("#container").width();
		scope.context.containerHeight = $("#container").height();
		KMG.TextureMap.textureResolution = scope.config.textureResolution;

		updateLights(scope.context.configChanged);
		
		for (var i = 0; i < scope.context.objects.length; i++) {
			scope.context.objects[i].update();
		}
		
		scope.context.configChanged = false;

		var time = Date.now() * 0.0004;
		
		if (scope.config.useScript) {
			fireOnRenderHandler();
		}
		
		if (areEffectsEnabled() && scope.config.postprocessingEnabled) {
			scope.context.composer.render( time );
		} else {
			scope.context.renderer.clear();
			scope.context.renderer.render(scope.context.secondaryScene, scope.context.secondaryCamera);
			scope.context.renderer.render(scope.context.primaryScene, scope.context.camera);
		}
		
		
		

	}
	
	function fireOnFrameHandler()
	{
		if (scope.config.useScript && scope.context.script && scope.context.script.onFrameHandler) {
			return scope.context.script.onFrameHandler(scope, scope.config, scope.context);
		} else {
			return false;
		}
	}
	
	function fireOnRenderHandler()
	{
		if (scope.config.useScript && scope.context.script && scope.context.script.onRenderHandler) {
			return scope.context.script.onRenderHandler(scope, scope.config, scope.context);
		} else {
			return false;
		}
	}
	
	this.applySceneScriptInstance = function(scriptInstance) {
		
		var sceneChanged = false;
		if (this.context.script && this.context.script.onScriptDestroy(this, this.config, this.context)) {
			sceneChanged = true;
		}
		
		this.context.script = scriptInstance;
		
		// Make sure we weren't passed a null script instance (which effectively 
		// disables the script interface)
		if (this.context.script && this.context.script.onScriptInitialize(this, this.config, this.context)) {
			sceneChanged = true;
		}
		
		if (sceneChanged) {
			this.context.configChanged = true;
		}
	};
	
	
	this.context.container = this.domElement;
	
	
	KMG.TextureMap.sceneReadyCallback = this.sceneCallbacks.sceneReadyCallback;
	KMG.TextureMap.resourceLoadingStart = this.sceneCallbacks.resourceLoadingStart;
	KMG.TextureMap.resourceLoadingFinish = this.sceneCallbacks.resourceLoadingFinish;
	KMG.TextureMap.renderCallback = render;

	this.context.camera = new THREE.PerspectiveCamera( this.config.camera.fieldOfView, $("#container").width() / $("#container").height(), this.config.camera.near, this.config.camera.far );
	this.context.camera.forceDefaultDistance = false;
	
	if (this.config.camera.useSecondaryParameters) {
		this.context.secondaryCamera = new THREE.PerspectiveCamera( this.config.camera.fieldOfViewSecondary, $("#container").width() / $("#container").height(), this.config.camera.nearSecondary, this.config.camera.farSecondary );
	} else {
		this.context.secondaryCamera = new THREE.PerspectiveCamera( this.config.camera.fieldOfView, $("#container").width() / $("#container").height(), this.config.camera.near, this.config.camera.far );
	}
	this.context.secondaryCamera.forceDefaultDistance = true;
	
	
	var renderer = new THREE.WebGLRenderer( { antialias: true, alpha: true, preserveDrawingBuffer : true} );
	renderer.setSize( window.innerWidth, window.innerHeight );
	
	if (this.config.initWithShadows) {
		renderer.shadowMapEnabled = true;
		renderer.shadowMapSoft = true;
		renderer.shadowMapType = THREE.PCFSoftShadowMap;
		renderer.shadowMapCullFace = THREE.CullFaceBack;
	}
	
	renderer.autoClear = false;
	renderer.setClearColor(0x000000);
	renderer.gammaInput = true;
	renderer.gammaOutput = true;
	renderer.physicallyBasedShading = true;
	
	this.context.renderer = renderer;
	renderer.context.canvas.addEventListener("webglcontextlost", function(event) {
		event.preventDefault();
		console.error("WebGL Context has been lost!");
		if (scope.context.animationID) {
			//cancelAnimationFrame(scope.context.animationID); 
		}
		
		if (sceneCallbacks.contextLostCallback) {
			sceneCallbacks.contextLostCallback(event);
		}
		
	}, false);

	renderer.context.canvas.addEventListener("webglcontextrestored", function(event) {
		animate();
	}, false);

	
	this.context.container.appendChild( renderer.domElement );
	
	if (this.config.enableFps) {
		var stats = new Stats();
		stats.domElement.style.position = 'absolute';
		stats.domElement.style.bottom = '0px';
		stats.domElement.style.right = '0px';
		this.context.container.appendChild( stats.domElement );
		this.context.stats = stats;
	}
	
	document.addEventListener( 'mousemove', onDocumentMouseMove, false );
	window.addEventListener( 'resize', onWindowResize, false );


	this.context.controls = new KMG.ExamineControls( renderer.getContext(), [this.context.camera, this.context.secondaryCamera], this.context.container );
	if (view) {
		this.context.controls.fromConfig(view);
	}
	this.context.controls.update(true);
	
	this.context.controls.addEventListener( 'change', render );
	if (this.cameraConfig != null && this.cameraConfig.controlCenter !== undefined) {
		this.context.controls.center.set(this.cameraConfig.controlCenter.x, this.cameraConfig.controlCenter.y, this.cameraConfig.controlCenter.z);
	}
	
	
	this.initializePostProcessing = function() {
		this.context.composer = new KMG.DynamicEffectsComposer(this.context, this.context.primaryScene, this.context.camera, this.context.secondaryCamera, this.context.renderer, this.config);
		this.context.objects.push(this.context.composer);
	};
	
	
	
	
	
	
};

/* File: Planet.js */

KMG.Planet = function ( domElement, config, sceneCallbacks, cameraConfig, view) {
	
	KMG.Engine.call( this, domElement, config, sceneCallbacks, cameraConfig, view );

	// Internals
	var scope = this;
	

	
	function addObjectToScene(object, scene)
	{
		scope.context.objects.push(object);
		
		if (!scene || scene === KMG.SCENE.PRIMARY) {
			scope.context.primaryScene.add( object );
		} else if (scene === KMG.SCENE.SECONDARY) {
			scope.context.secondaryScene.add( object );
		}
		
		return object;
	}
	
	function initShadows(light) {
		light.castShadow = true;
		light.shadowCameraVisible = false;
		light.shadowMapWidth = 2048;
		light.shadowMapHeight = 2048;
		light.shadowCameraNear = -4000;
		light.shadowCameraFar = 100000;
		light.shadowCameraFov = 45.0;
		
		light.shadowCameraRight     =  700;
		light.shadowCameraLeft     = -700;
		light.shadowCameraTop      =  700;
		light.shadowCameraBottom   = -700;
		//light0.shadowBias = 0.03001;
		//light0.shadowDarkness = 0.5;
	}
	
	function buildScene()
	{
		scope.context.primaryScene = new THREE.Scene();

		scope.context.lights.primaryDirectional = new THREE.DirectionalLight( 0xFFFFFF, 2.0, 100);
		scope.context.lights.primaryDirectional.position.set( -10000, 0, 0 ).normalize();
		scope.context.primaryScene.add(scope.context.lights.primaryDirectional);

			
		
		scope.context.lights.primaryPoint = new THREE.PointLight( 0xFFFFFF, 0.0);
		scope.context.lights.primaryPoint.position.set( 0, 0, 0 );
		scope.context.primaryScene.add(scope.context.lights.primaryPoint);

		scope.context.lights.ambient = new THREE.AmbientLight( 0x888888 );
		scope.context.primaryScene.add( scope.context.lights.ambient );
		
		scope.context.lights.ambient = new THREE.AmbientLight( 0x888888 )
		
		if (scope.config.initWithShadows) {
			initShadows(scope.context.lights.primaryDirectional);
			//initShadows(scope.context.lights.primaryPoint);
		}
		
		
		if (!scope.config.noPlanet) {
			scope.context.surfaceObject = addObjectToScene(new KMG.SurfaceObject(scope.context, scope.config));
		}
		
		scope.context.secondaryScene = new THREE.Scene();
		
		scope.context.lights.secondaryDirectional = scope.context.lights.primaryDirectional.clone();
		scope.context.lights.secondaryPoint = scope.context.lights.primaryPoint.clone();

		scope.context.secondaryScene.add(scope.context.lights.secondaryDirectional);
		scope.context.secondaryScene.add(scope.context.lights.secondaryPoint);
		
		
		
		addObjectToScene(new KMG.LocalStarObject(scope.context, scope.config), KMG.SCENE.SECONDARY);
		addObjectToScene(new KMG.LensFlareObject(scope.context, scope.config), KMG.SCENE.SECONDARY);
		
		if (!scope.config.noBackground) {
			scope.context.background = new KMG.BackgroundObject(scope.context, scope.config);
			addObjectToScene(scope.context.background, KMG.SCENE.SECONDARY);
		}
		

		
		scope.context.primaryScene.updateMatrix();
		scope.context.secondaryScene.updateMatrix();
	}
	
	
	
	buildScene();
	
	this.initializePostProcessing();
	
	
};
KMG.Planet.prototype = Object.create( KMG.Engine.prototype );

/* File: RingGeometry2.js */

/** A modification of the standard three.js RingGeometry class, but with changes to support 
 * Celestia-like ring textures.
 */
THREE.RingGeometry2 = function ( innerRadius, outerRadius, thetaSegments, phiSegments, thetaStart, thetaLength ) {

    THREE.Geometry.call( this );

    innerRadius = innerRadius || 0;
    outerRadius = outerRadius || 50;

    thetaStart = thetaStart !== undefined ? thetaStart : 0;
    thetaLength = thetaLength !== undefined ? thetaLength : Math.PI * 2;

    thetaSegments = thetaSegments !== undefined ? Math.max( 3, thetaSegments ) : 8;
    phiSegments = phiSegments !== undefined ? Math.max( 3, phiSegments ) : 8;
    
    var i, o, uvs = [], radius = innerRadius, radiusStep = ( ( outerRadius - innerRadius ) / phiSegments);
	

	
    for( i = 0; i <= phiSegments; i++) {//concentric circles inside ring

        for( o = 0; o <= thetaSegments; o++) {//number of segments per circle

            var vertex = new THREE.Vector3();
            
            vertex.x = radius * Math.cos( thetaStart + o / thetaSegments * thetaLength );
            vertex.y = radius * Math.sin( thetaStart + o / thetaSegments * thetaLength );
            
            this.vertices.push( vertex );
			uvs.push( new THREE.Vector2((i / phiSegments), ( vertex.y / radius + 1 ) / 2));
        }
        
        radius += radiusStep;

    }
	
	
    var n = new THREE.Vector3( 0, 0, 1 );
    
    for( i = 0; i < phiSegments; i++) {//concentric circles inside ring

        for( o = 0; o <= thetaSegments; o++) {//number of segments per circle
            
            var v1, v2, v3;

            v1 = o + (thetaSegments * i) + i;
            v2 = o + (thetaSegments * i) + thetaSegments + i;
            v3 = o + (thetaSegments * i) + thetaSegments + 1 + i;
            
            this.faces.push( new THREE.Face3( v1, v2, v3, [ n, n, n ] ) );
            this.faceVertexUvs[ 0 ].push( [ uvs[ v1 ], uvs[ v2 ], uvs[ v3 ] ]);
            
            v1 = o + (thetaSegments * i) + i;
            v2 = o + (thetaSegments * i) + thetaSegments + 1 + i;
            v3 = o + (thetaSegments * i) + 1 + i;
            
            this.faces.push( new THREE.Face3( v1, v2, v3, [ n, n, n ] ) );
            this.faceVertexUvs[ 0 ].push( [ uvs[ v1 ], uvs[ v2 ], uvs[ v3 ] ]);

        }
    }
    
    this.computeCentroids();
    this.computeFaceNormals();

    this.boundingSphere = new THREE.Sphere( new THREE.Vector3(), radius ); 

};

THREE.RingGeometry2.prototype = Object.create( THREE.Geometry.prototype );

/* File: EllipsoidGeometry.js */
/** Modification of SphereGeometry to implements an Ellipsoid Oblate Spheriod
 * @author mrdoob / http://mrdoob.com/
 */

THREE.EllipsoidGeometry = function ( radius, flattening, widthSegments, heightSegments, phiStart, phiLength, thetaStart, thetaLength ) {

	THREE.Geometry.call( this );

	this.radius = radius = radius || 50;

	this.widthSegments = widthSegments = Math.max( 3, Math.floor( widthSegments ) || 8 );
	this.heightSegments = heightSegments = Math.max( 2, Math.floor( heightSegments ) || 6 );

	this.phiStart = phiStart = phiStart !== undefined ? phiStart : 0;
	this.phiLength = phiLength = phiLength !== undefined ? phiLength : Math.PI * 2;

	this.thetaStart = thetaStart = thetaStart !== undefined ? thetaStart : 0;
	this.thetaLength = thetaLength = thetaLength !== undefined ? thetaLength : Math.PI;

	var x, y, vertices = [], uvs = [];

	for ( y = 0; y <= heightSegments; y ++ ) {

		var verticesRow = [];
		var uvsRow = [];

		for ( x = 0; x <= widthSegments; x ++ ) {

			var u = x / widthSegments;
			var v = y / heightSegments;
			

			var lat =  thetaStart + v * thetaLength;
			var r = KMG.Math.radiusAtGeocentricLatitude(radius, lat - KMG.RAD_90, flattening);
			
			var vertex = new THREE.Vector3();
			vertex.x = - r * Math.cos( phiStart + u * phiLength ) * Math.sin( thetaStart + v * thetaLength );
			vertex.y = r * Math.cos( thetaStart + v * thetaLength );
			vertex.z = r * Math.sin( phiStart + u * phiLength ) * Math.sin( thetaStart + v * thetaLength );

			this.vertices.push( vertex );

			verticesRow.push( this.vertices.length - 1 );
			uvsRow.push( new THREE.Vector2( u, 1 - v ) );

		}

		vertices.push( verticesRow );
		uvs.push( uvsRow );

	}

	for ( y = 0; y < this.heightSegments; y ++ ) {

		for ( x = 0; x < this.widthSegments; x ++ ) {

			var v1 = vertices[ y ][ x + 1 ];
			var v2 = vertices[ y ][ x ];
			var v3 = vertices[ y + 1 ][ x ];
			var v4 = vertices[ y + 1 ][ x + 1 ];

			var n1 = this.vertices[ v1 ].clone().normalize();
			var n2 = this.vertices[ v2 ].clone().normalize();
			var n3 = this.vertices[ v3 ].clone().normalize();
			var n4 = this.vertices[ v4 ].clone().normalize();

			var uv1 = uvs[ y ][ x + 1 ].clone();
			var uv2 = uvs[ y ][ x ].clone();
			var uv3 = uvs[ y + 1 ][ x ].clone();
			var uv4 = uvs[ y + 1 ][ x + 1 ].clone();

			if ( Math.abs( this.vertices[ v1 ].y ) === this.radius ) {

				uv1.x = ( uv1.x + uv2.x ) / 2;
				this.faces.push( new THREE.Face3( v1, v3, v4, [ n1, n3, n4 ] ) );
				this.faceVertexUvs[ 0 ].push( [ uv1, uv3, uv4 ] );

			} else if ( Math.abs( this.vertices[ v3 ].y ) === this.radius ) {

				uv3.x = ( uv3.x + uv4.x ) / 2;
				this.faces.push( new THREE.Face3( v1, v2, v3, [ n1, n2, n3 ] ) );
				this.faceVertexUvs[ 0 ].push( [ uv1, uv2, uv3 ] );

			} else {

				this.faces.push( new THREE.Face3( v1, v2, v4, [ n1, n2, n4 ] ) );
				this.faceVertexUvs[ 0 ].push( [ uv1, uv2, uv4 ] );

				this.faces.push( new THREE.Face3( v2, v3, v4, [ n2.clone(), n3, n4.clone() ] ) );
				this.faceVertexUvs[ 0 ].push( [ uv2.clone(), uv3, uv4.clone() ] );

			}

		}

	}

	this.computeCentroids();
	this.computeFaceNormals();

	this.boundingSphere = new THREE.Sphere( new THREE.Vector3(), radius );

};

THREE.EllipsoidGeometry.prototype = Object.create( THREE.Geometry.prototype );

/* File: Vector3.js */

THREE.Vector3.prototype.rotateX = function (angle) {

	var cosX = Math.cos(angle);
	var sinX = Math.sin(angle);
			
	var ry = cosX * this.y + -sinX * this.z;
	var rz = sinX * this.y + cosX * this.z;
			
	this.y = ry;
	this.z = rz;
	
	return this;
};

THREE.Vector3.prototype.rotateY = function (angle) {

	var cosY = Math.cos(angle);
	var sinY = Math.sin(angle);
	
	var rx = cosY * this.x + sinY * this.z;
	var rz = -sinY * this.x + cosY * this.z;
	
	this.x = rx;
	this.z = rz;
	
	return this;
};

THREE.Vector3.prototype.rotateZ = function (angle) {

	var cosZ = Math.cos(angle);
	var sinZ = Math.sin(angle);
	
	var rx = cosZ * this.x + -sinZ * this.y;
	var ry = sinZ * this.x + cosZ * this.y;

	this.x = rx;
	this.y = ry;
	
	return this;
};


THREE.Vector3.prototype.rotate = function (angle, axis) {
	if (axis === undefined || axis === 'X') {
		return this.rotateX(angle);
	} else if (axis === 'Y') {
		return this.rotateY(angle);
	} else if (axis === 'Z') {
		return this.rotateZ(angle);
	}
	return this;
};

/* File: ConfigWrapper.js */

/** A very simple structure for creating property name aliases.
 *
 */
KMG.ConfigWrapper = function ( config ) {
	
	var scope = this;
	var config = config;
	var nameToPropMap = {};
	var propToNameMap = {};
	this.changeListener = null;
	
	this.add = function ( name, prop ) {
		propToNameMap[prop] = name;
		nameToPropMap[name] = prop;
		
		Object.defineProperty(this, name, {
			get : function()  { return config[nameToPropMap[name]]; },
			set : function(v) { 
				config[nameToPropMap[name]] = v;
				if (scope.changeListener) {
					scope.changeListener(nameToPropMap[name], v);
				}
			}
		});
		
	};
};

/* File: AnimateController.js */

KMG.AnimateController = function ( config, object ) {
	
	this.config = config;
	this.object = object;
	var scope = this;
	
	this.current = config.integralStart;
	
	this.isActive = false;

	function nextDistance()
	{
		var rawDistance = scope.config.distanceIntegral(scope.current);
		var distance = ((scope.config.startDistance - scope.config.minDistance) * ((rawDistance - scope.config.integralMin) / (scope.config.integralMax - scope.config.integralMin))) + scope.config.minDistance;
		return distance;
	}
	
	function nextRotation()
	{
		var rawRotation = scope.config.rotateIntegral(scope.current);
		return rawRotation*(Math.PI/180.0);
	}
	
	function onNext(direction)
	{
		if (!scope.isActive) {
			return;
		}
		
		if (!direction) {
			direction = 1.0;
		}
		
		scope.current += (scope.config.speed * direction);
		var distance = nextDistance();
		var rotation = nextRotation() * direction;
		
		var position = scope.object.position;
		var offset = position.clone();
		
		var theta = Math.atan2( offset.x, offset.z );
		var phi = Math.atan2( Math.sqrt( offset.x * offset.x + offset.z * offset.z ), offset.y );

		theta += rotation;
		phi += 0.0;
		
		offset.x = distance * Math.sin( phi ) * Math.sin( theta );
		offset.y = distance * Math.cos( phi );
		offset.z = distance * Math.sin( phi ) * Math.cos( theta );
	
		scope.object.position = offset;

		if (scope.current >= scope.config.integralEnd || scope.current <= scope.config.integralStart) {
			scope.stop();
		}
	}
	
	this.next = function() {
		onNext();
	};
	
	this.start = function() {
		this.isActive = true;
	};
	
	this.stop = function() {
		this.isActive = false;
	};
	
	
	this.rewind = function() {
		this.isActive = true;
		this.current = this.config.integralEnd;
		for (var i = this.config.integralEnd; i >= this.config.integralStart; i-=this.config.speed) {
			onNext(-1);
		}
		this.current = this.config.integralStart;
		this.isActive = false;
	};
};

/* File: ExtendedNormalMapShader.js */

KMG.ExtendedNormalMapShader = {

	uniforms: THREE.UniformsUtils.merge( [

		THREE.UniformsLib[ "fog" ],
		THREE.UniformsLib[ "lights" ],
		THREE.UniformsLib[ "shadowmap" ],

		{

		"enableSpecular"  : { type: "i", value: 0 },
		"enableReflection": { type: "i", value: 0 },
		"enableDisplacement": { type: "i", value: 0 },
		"usingDirectionalLighting" : { type : "i", value: 1 },
		"enableFog"		  : { type: "i", value: 1 },
		
		"tDisplacement": { type: "t", value: null }, // must go first as this is vertex texture
		"tDiffuse"	   : { type: "t", value: null },
		"tCube"		   : { type: "t", value: null },
		"tNormal"	   : { type: "t", value: null },
		"tSpecular"	   : { type: "t", value: null },
		
		"uNormalScale": { type: "v2", value: new THREE.Vector2( 1, 1 ) },

		"uDisplacementBias": { type: "f", value: 0.0 },
		"uDisplacementScale": { type: "f", value: 1.0 },

		"uDiffuseColor": { type: "c", value: new THREE.Color( 0xffffff ) },
		"uSpecularColor": { type: "c", value: new THREE.Color( 0x111111 ) },
		"uAmbientColor": { type: "c", value: new THREE.Color( 0xffffff ) },
		"uEmissiveColor": { type: "c", value: new THREE.Color( 0x000000 ) },
		"uShininess": { type: "f", value: 30 },
		"uOpacity": { type: "f", value: 1 },
		
		"useRefract": { type: "i", value: 0 },
		"uRefractionRatio": { type: "f", value: 0.98 },
		"uReflectivity": { type: "f", value: 0.5 },
	
		"uAOLevel": { type: "f", value: 1.0 },
		"uAlphaMulitplier": { type: "f", value: 1.0 },
		
		"uOffset" : { type: "v2", value: new THREE.Vector2( 0.001, 0.001 ) },
		"uRepeat" : { type: "v2", value: new THREE.Vector2( 0.998, 0.998 ) },

		"wrapRGB"  : { type: "v3", value: new THREE.Vector3( 0, 0, 0 ) },

		}

	] ),

	fragmentShader: [

		"uniform vec3 uAmbientColor;",
		"uniform vec3 uDiffuseColor;",
		"uniform vec3 uSpecularColor;",
		"uniform vec3 uEmissiveColor;",
		"uniform float uShininess;",
		"uniform float uAOLevel;",
		"uniform float uAlphaMulitplier;",
		
		"uniform bool enableSpecular;",
		"uniform bool enableReflection;",
		"uniform bool usingDirectionalLighting;",
		"uniform bool enableFog;",
		
		"uniform sampler2D tDiffuse;",
		"uniform sampler2D tNormal;",
		"uniform sampler2D tSpecular;",
		
		"uniform samplerCube tCube;",

		"uniform vec2 uNormalScale;",

		"uniform bool useRefract;",
		"uniform float uRefractionRatio;",
		"uniform float uReflectivity;",

		"varying vec3 vTangent;",
		"varying vec3 vBinormal;",
		"varying vec3 vNormal;",
		"varying vec2 vUv;",

		"uniform vec3 ambientLightColor;",
		
		"#if MAX_DIR_LIGHTS > 0",

			"uniform vec3 directionalLightColor[ MAX_DIR_LIGHTS ];",
			"uniform vec3 directionalLightDirection[ MAX_DIR_LIGHTS ];",

		"#endif",

		"#if MAX_HEMI_LIGHTS > 0",

			"uniform vec3 hemisphereLightSkyColor[ MAX_HEMI_LIGHTS ];",
			"uniform vec3 hemisphereLightGroundColor[ MAX_HEMI_LIGHTS ];",
			"uniform vec3 hemisphereLightDirection[ MAX_HEMI_LIGHTS ];",

		"#endif",

		"#if MAX_POINT_LIGHTS > 0",

			"uniform vec3 pointLightColor[ MAX_POINT_LIGHTS ];",
			"uniform vec3 pointLightPosition[ MAX_POINT_LIGHTS ];",
			"uniform float pointLightDistance[ MAX_POINT_LIGHTS ];",

		"#endif",

		"#if MAX_SPOT_LIGHTS > 0",

			"uniform vec3 spotLightColor[ MAX_SPOT_LIGHTS ];",
			"uniform vec3 spotLightPosition[ MAX_SPOT_LIGHTS ];",
			"uniform vec3 spotLightDirection[ MAX_SPOT_LIGHTS ];",
			"uniform float spotLightAngleCos[ MAX_SPOT_LIGHTS ];",
			"uniform float spotLightExponent[ MAX_SPOT_LIGHTS ];",
			"uniform float spotLightDistance[ MAX_SPOT_LIGHTS ];",

		"#endif",

		"#ifdef WRAP_AROUND",

			"uniform vec3 wrapRGB;",

		"#endif",

		"varying vec3 vWorldPosition;",
		"varying vec3 vViewPosition;",

		THREE.ShaderChunk[ "shadowmap_pars_fragment" ],
		THREE.ShaderChunk[ "fog_pars_fragment" ],

		"void main() {",

			"gl_FragColor = vec4(1.0 );",
			
			"vec3 specularTex = vec3( 1.0 );",
			
			"vec3 normalTex = texture2D( tNormal, vUv ).xyz * 2.0 - 1.0;",
			"normalTex.xy *= uNormalScale;",
			"normalTex = normalize( normalTex );",
			
			
			"vec4 texelColor = texture2D( tDiffuse, vUv );",
			
			"texelColor = vec4(texelColor.xyz, texelColor.w * uAlphaMulitplier);",
			"gl_FragColor = gl_FragColor * (texelColor * texelColor);",
			
			"vec3 aoColor;",

			"if( enableSpecular) {",
				"specularTex = texture2D( tSpecular, vUv ).xyz;",
			"}",
			"mat3 tsb = mat3( normalize( vTangent ), normalize( vBinormal ), normalize( vNormal ) );",
			"vec3 finalNormal = tsb * normalTex;",

			"#ifdef FLIP_SIDED",

				"finalNormal = -finalNormal;",

			"#endif",

			"vec3 normal = normalize( finalNormal );",
			"vec3 viewPosition = normalize( vViewPosition );",

			// point lights

			"#if MAX_POINT_LIGHTS > 0",

				"vec3 pointDiffuse = vec3( 0.0 );",
				"vec3 pointSpecular = vec3( 0.0 );",

				"for ( int i = 0; i < MAX_POINT_LIGHTS; i ++ ) {",

					"vec4 lPosition = viewMatrix * vec4( pointLightPosition[ i ], 1.0 );",
					"vec3 pointVector = lPosition.xyz + vViewPosition.xyz;",

					"float pointDistance = 1.0;",
					"if ( pointLightDistance[ i ] > 0.0 )",
						"pointDistance = 1.0 - min( ( length( pointVector ) / pointLightDistance[ i ] ), 1.0 );",

					"pointVector = normalize( pointVector );",

					// diffuse

					"#ifdef WRAP_AROUND",

						"float pointDiffuseWeightFull = max( dot( normal, pointVector ), 0.0 );",
						"float pointDiffuseWeightHalf = max( 0.5 * dot( normal, pointVector ) + 0.5, 0.0 );",

						"vec3 pointDiffuseWeight = mix( vec3 ( pointDiffuseWeightFull ), vec3( pointDiffuseWeightHalf ), wrapRGB );",

					"#else",

						"float pointDiffuseWeight = max( dot( normal, pointVector ), 0.0 );",

					"#endif",

					"pointDiffuse += pointDistance * pointLightColor[ i ] * uDiffuseColor * pointDiffuseWeight;",

					// specular

					"vec3 pointHalfVector = normalize( pointVector + viewPosition );",
					"float pointDotNormalHalf = max( dot( normal, pointHalfVector ), 0.0 );",
					"float pointSpecularWeight = specularTex.r * max( pow( pointDotNormalHalf, uShininess ), 0.0 );",

					"#ifdef PHYSICALLY_BASED_SHADING",

						// 2.0 => 2.0001 is hack to work around ANGLE bug

						"float specularNormalization = ( uShininess + 2.0001 ) / 8.0;",

						"vec3 schlick = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( pointVector, pointHalfVector ), 5.0 );",
						"pointSpecular += schlick * pointLightColor[ i ] * pointSpecularWeight * pointDiffuseWeight * pointDistance * specularNormalization;",

					"#else",

						"pointSpecular += pointDistance * pointLightColor[ i ] * uSpecularColor * pointSpecularWeight * pointDiffuseWeight;",

					"#endif",

				"}",

			"#endif",

			// spot lights

			"#if MAX_SPOT_LIGHTS > 0",

				"vec3 spotDiffuse = vec3( 0.0 );",
				"vec3 spotSpecular = vec3( 0.0 );",

				"for ( int i = 0; i < MAX_SPOT_LIGHTS; i ++ ) {",

					"vec4 lPosition = viewMatrix * vec4( spotLightPosition[ i ], 1.0 );",
					"vec3 spotVector = lPosition.xyz + vViewPosition.xyz;",

					"float spotDistance = 1.0;",
					"if ( spotLightDistance[ i ] > 0.0 )",
						"spotDistance = 1.0 - min( ( length( spotVector ) / spotLightDistance[ i ] ), 1.0 );",

					"spotVector = normalize( spotVector );",

					"float spotEffect = dot( spotLightDirection[ i ], normalize( spotLightPosition[ i ] - vWorldPosition ) );",

					"if ( spotEffect > spotLightAngleCos[ i ] ) {",

						"spotEffect = max( pow( spotEffect, spotLightExponent[ i ] ), 0.0 );",

						// diffuse

						"#ifdef WRAP_AROUND",

							"float spotDiffuseWeightFull = max( dot( normal, spotVector ), 0.0 );",
							"float spotDiffuseWeightHalf = max( 0.5 * dot( normal, spotVector ) + 0.5, 0.0 );",

							"vec3 spotDiffuseWeight = mix( vec3 ( spotDiffuseWeightFull ), vec3( spotDiffuseWeightHalf ), wrapRGB );",

						"#else",

							"float spotDiffuseWeight = max( dot( normal, spotVector ), 0.0 );",

						"#endif",

						"spotDiffuse += spotDistance * spotLightColor[ i ] * uDiffuseColor * spotDiffuseWeight * spotEffect;",

						// specular

						"vec3 spotHalfVector = normalize( spotVector + viewPosition );",
						"float spotDotNormalHalf = max( dot( normal, spotHalfVector ), 0.0 );",
						"float spotSpecularWeight = specularTex.r * max( pow( spotDotNormalHalf, uShininess ), 0.0 );",

						"#ifdef PHYSICALLY_BASED_SHADING",

							// 2.0 => 2.0001 is hack to work around ANGLE bug

							"float specularNormalization = ( uShininess + 2.0001 ) / 8.0;",

							"vec3 schlick = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( spotVector, spotHalfVector ), 5.0 );",
							"spotSpecular += schlick * spotLightColor[ i ] * spotSpecularWeight * spotDiffuseWeight * spotDistance * specularNormalization * spotEffect;",

						"#else",

							"spotSpecular += spotDistance * spotLightColor[ i ] * uSpecularColor * spotSpecularWeight * spotDiffuseWeight * spotEffect;",

						"#endif",

					"}",

				"}",

			"#endif",

			// directional lights

			"#if MAX_DIR_LIGHTS > 0",

				"vec3 dirDiffuse = vec3( 0.0 );",
				"vec3 dirSpecular = vec3( 0.0 );",

				"for( int i = 0; i < MAX_DIR_LIGHTS; i++ ) {",

					"vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ i ], 0.0 );",
					"vec3 dirVector = normalize( lDirection.xyz );",

					// diffuse

					"#ifdef WRAP_AROUND",

						"float directionalLightWeightingFull = max( dot( normal, dirVector ), 0.0 );",
						"float directionalLightWeightingHalf = max( 0.5 * dot( normal, dirVector ) + 0.5, 0.0 );",

						"vec3 dirDiffuseWeight = mix( vec3( directionalLightWeightingFull ), vec3( directionalLightWeightingHalf ), wrapRGB );",

					"#else",

						"float dirDiffuseWeight = max( dot( normal, dirVector ), 0.0 );",

					"#endif",

					"dirDiffuse += directionalLightColor[ i ] * uDiffuseColor * dirDiffuseWeight;",

					// specular

					"vec3 dirHalfVector = normalize( dirVector + viewPosition );",
					"float dirDotNormalHalf = max( dot( normal, dirHalfVector ), 0.0 );",
					"float dirSpecularWeight = specularTex.r * max( pow( dirDotNormalHalf, uShininess ), 0.0 );",

					"#ifdef PHYSICALLY_BASED_SHADING",

						// 2.0 => 2.0001 is hack to work around ANGLE bug

						"float specularNormalization = ( uShininess + 2.0001 ) / 8.0;",

						"vec3 schlick = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( dirVector, dirHalfVector ), 5.0 );",
						"dirSpecular += schlick * directionalLightColor[ i ] * dirSpecularWeight * dirDiffuseWeight * specularNormalization;",

					"#else",

						"dirSpecular += directionalLightColor[ i ] * uSpecularColor * dirSpecularWeight * dirDiffuseWeight;",

					"#endif",

				"}",

			"#endif",

			// hemisphere lights

			"#if MAX_HEMI_LIGHTS > 0",

				"vec3 hemiDiffuse  = vec3( 0.0 );",
				"vec3 hemiSpecular = vec3( 0.0 );" ,

				"for( int i = 0; i < MAX_HEMI_LIGHTS; i ++ ) {",

					"vec4 lDirection = viewMatrix * vec4( hemisphereLightDirection[ i ], 0.0 );",
					"vec3 lVector = normalize( lDirection.xyz );",

					// diffuse

					"float dotProduct = dot( normal, lVector );",
					"float hemiDiffuseWeight = 0.5 * dotProduct + 0.5;",

					"vec3 hemiColor = mix( hemisphereLightGroundColor[ i ], hemisphereLightSkyColor[ i ], hemiDiffuseWeight );",

					"hemiDiffuse += uDiffuseColor * hemiColor;",

					// specular (sky light)


					"vec3 hemiHalfVectorSky = normalize( lVector + viewPosition );",
					"float hemiDotNormalHalfSky = 0.5 * dot( normal, hemiHalfVectorSky ) + 0.5;",
					"float hemiSpecularWeightSky = specularTex.r * max( pow( hemiDotNormalHalfSky, uShininess ), 0.0 );",

					// specular (ground light)

					"vec3 lVectorGround = -lVector;",

					"vec3 hemiHalfVectorGround = normalize( lVectorGround + viewPosition );",
					"float hemiDotNormalHalfGround = 0.5 * dot( normal, hemiHalfVectorGround ) + 0.5;",
					"float hemiSpecularWeightGround = specularTex.r * max( pow( hemiDotNormalHalfGround, uShininess ), 0.0 );",

					"#ifdef PHYSICALLY_BASED_SHADING",

						"float dotProductGround = dot( normal, lVectorGround );",

						// 2.0 => 2.0001 is hack to work around ANGLE bug

						"float specularNormalization = ( uShininess + 2.0001 ) / 8.0;",

						"vec3 schlickSky = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( lVector, hemiHalfVectorSky ), 5.0 );",
						"vec3 schlickGround = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( lVectorGround, hemiHalfVectorGround ), 5.0 );",
						"hemiSpecular += hemiColor * specularNormalization * ( schlickSky * hemiSpecularWeightSky * max( dotProduct, 0.0 ) + schlickGround * hemiSpecularWeightGround * max( dotProductGround, 0.0 ) );",

					"#else",

						"hemiSpecular += uSpecularColor * hemiColor * ( hemiSpecularWeightSky + hemiSpecularWeightGround ) * hemiDiffuseWeight;",

					"#endif",

				"}",

			"#endif",

			// all lights contribution summation

			"vec3 totalDiffuse = vec3( 0.0 );",
			"vec3 totalSpecular = vec3( 0.0 );",

			"#if MAX_DIR_LIGHTS > 0",

				"totalDiffuse += dirDiffuse;",
				"totalSpecular += dirSpecular;",

			"#endif",

			"#if MAX_HEMI_LIGHTS > 0",

				"totalDiffuse += hemiDiffuse;",
				"totalSpecular += hemiSpecular;",

			"#endif",

			"#if MAX_POINT_LIGHTS > 0",

				"totalDiffuse += pointDiffuse;",
				"totalSpecular += pointSpecular;",

			"#endif",

			"#if MAX_SPOT_LIGHTS > 0",

				"totalDiffuse += spotDiffuse;",
				"totalSpecular += spotSpecular;",

			"#endif",
			

			THREE.ShaderChunk[ "shadowmap_fragment" ],

			"vec3 emissiveColor = uEmissiveColor * (texelColor.xyz + texelColor.xyz * uEmissiveColor);",
			
			"#ifdef METAL",

				"gl_FragColor.xyz = gl_FragColor.xyz * ( totalDiffuse + ambientLightColor * uAmbientColor + emissiveColor + totalSpecular );",

			"#else",
				
				"gl_FragColor.xyz = gl_FragColor.xyz * ( totalDiffuse + ambientLightColor * uAmbientColor) + totalSpecular + emissiveColor;",
				
			"#endif",
			
			
			
			
			
			THREE.ShaderChunk[ "linear_to_gamma_fragment" ],
			
			"if (enableFog) { ",
				"#ifdef USE_FOG",
				"float depth = gl_FragCoord.z / gl_FragCoord.w;",
				
				"vec3 _fogColor = fogColor;",
				
				"if (usingDirectionalLighting) { ",
					"_fogColor = _fogColor * dirDiffuse;",
				"}",
				
				"#ifdef FOG_EXP2",

					"const float LOG2 = 1.442695;",
					"float fogFactor = exp2( - fogDensity * fogDensity * depth * depth * LOG2 );",
					"fogFactor = 1.0 - clamp( fogFactor, 0.0, 1.0 );",

				"#else",

					"float fogFactor = smoothstep( fogNear, fogFar, depth );",

				"#endif",
				"gl_FragColor = mix( gl_FragColor, vec4( _fogColor, gl_FragColor.w ), fogFactor );",
				
				"#endif",
			"}",
			
			
			//THREE.ShaderChunk[ "fog_fragment" ],


			

			
			
			
			
		"}"

	].join("\n"),

	vertexShader: [

		"attribute vec4 tangent;",

		"uniform vec2 uOffset;",
		"uniform vec2 uRepeat;",

		"uniform bool enableDisplacement;",

		"#ifdef VERTEX_TEXTURES",

			"uniform sampler2D tDisplacement;",
			"uniform float uDisplacementScale;",
			"uniform float uDisplacementBias;",

		"#endif",

		"varying vec3 vTangent;",
		"varying vec3 vBinormal;",
		"varying vec3 vNormal;",
		"varying vec2 vUv;",

		"varying vec3 vWorldPosition;",
		"varying vec3 vViewPosition;",

		THREE.ShaderChunk[ "skinning_pars_vertex" ],
		THREE.ShaderChunk[ "shadowmap_pars_vertex" ],

		"void main() {",

			THREE.ShaderChunk[ "skinbase_vertex" ],
			THREE.ShaderChunk[ "skinnormal_vertex" ],

			// normal, tangent and binormal vectors

			"#ifdef USE_SKINNING",

				"vNormal = normalize( normalMatrix * skinnedNormal.xyz );",

				"vec4 skinnedTangent = skinMatrix * vec4( tangent.xyz, 0.0 );",
				"vTangent = normalize( normalMatrix * skinnedTangent.xyz );",

			"#else",

				"vNormal = normalize( normalMatrix * normal );",
				"vTangent = normalize( normalMatrix * tangent.xyz );",

			"#endif",

			"vBinormal = normalize( cross( vNormal, vTangent ) * tangent.w );",

			"vUv = uv * uRepeat + uOffset;",

			// displacement mapping

			"vec3 displacedPosition;",

			"#ifdef VERTEX_TEXTURES",

				"if ( enableDisplacement ) {",

					"vec3 dv = texture2D( tDisplacement, uv ).xyz;",
					"float df = uDisplacementScale * dv.x + uDisplacementBias;",
					"displacedPosition = position + normalize( normal ) * df;",

				"} else {",

					"#ifdef USE_SKINNING",

						"vec4 skinVertex = vec4( position, 1.0 );",

						"vec4 skinned  = boneMatX * skinVertex * skinWeight.x;",
						"skinned 	  += boneMatY * skinVertex * skinWeight.y;",

						"displacedPosition  = skinned.xyz;",

					"#else",

						"displacedPosition = position;",

					"#endif",

				"}",

			"#else",

				"#ifdef USE_SKINNING",

					"vec4 skinVertex = vec4( position, 1.0 );",

					"vec4 skinned  = boneMatX * skinVertex * skinWeight.x;",
					"skinned 	  += boneMatY * skinVertex * skinWeight.y;",

					"displacedPosition  = skinned.xyz;",

				"#else",

					"displacedPosition = position;",

				"#endif",

			"#endif",

			//

			"vec4 mvPosition = modelViewMatrix * vec4( displacedPosition, 1.0 );",
			"vec4 worldPosition = modelMatrix * vec4( displacedPosition, 1.0 );",

			"gl_Position = projectionMatrix * mvPosition;",

			//

			"vWorldPosition = worldPosition.xyz;",
			"vViewPosition = -mvPosition.xyz;",

			// shadows

			"#ifdef USE_SHADOWMAP",

				"for( int i = 0; i < MAX_SHADOWS; i ++ ) {",

					"vShadowCoord[ i ] = shadowMatrix[ i ] * worldPosition;",

				"}",

			"#endif",

		"}"

	].join("\n")

}




/* File: ExamineControls.js */

KMG.CenterPivot = 1;
KMG.SurfacePivot = 2;

KMG.Distance = 1;
KMG.Scale = 2;
KMG.FoV = 3;

KMG.ExamineControls = function ( gl, object, domElement, onChange ) {

	this.gl = gl;
	this.object = object;
	this.domElement = ( domElement !== undefined ) ? domElement : document;
	
	this.orientation = new THREE.Quaternion();
	
	console.info("Initializing examine controls");
	
	// API
	this.enabled = true;
	
	this.keys = { LEFT: 37, UP: 38, RIGHT: 39, BOTTOM: 40 };
	
	this.center = new THREE.Vector3(0, 0, 0);
	this.lastPosition = new THREE.Vector3(0, 0, 0);
	
	this.translate = new THREE.Vector3(0, 0, 0);
	
	this.pitchType = KMG.SurfacePivot;
	
	this.zoomType = KMG.Distance;
	
	this.maxDistance = 10000;
	this.minDistance = 210;
	this.distance = 700;
	this.defaultDistance = 700;
		
	this.maxScale = 10000.0;
	this.minScale = 0.0001;
	this.scale = 1.0;
	this.defaultScale = 1.0;
	
	this.maxFov = 90;
	this.minFov = 0.0001;
	this.fov = 45;
	this.defaultFov = 45;
	
	
	this.maxPitch = 90.0 * (Math.PI / 180.0);
	this.minPitch = 0.0;
	this._pitch = 0.0;
	this._yaw = 0.0;
	this._roll = 0.0;
	
	this.panVertical = 0;
	this.panHorizontal = 0;
	
	this.radius = 200;
	
	this.distanceMoveSpeed = 0.2;
	this.zoomSpeed = 0.001;
	this.rotateSpeed = 0.5;
	
	this.modelView = new THREE.Matrix4();
	
	var matrixRoll = new THREE.Matrix4();
	var matrixPitch = new THREE.Matrix4();
	var matrixYaw = new THREE.Matrix4();
	

	var lastX = -1;
	var lastY = -1;
	var mouseDownX = -1;
	var mouseDownY = -1;
	
	var scope = this;
	var changeEvent = { type: 'change' };
	
	var lastRotateX = 0;
	var lastRotateY = 0;
	
	var STATE = { NONE: -1, ROTATE: 0, ZOOM: 1, PAN: 2, TOUCH_ROTATE : 3, TOUCH_ZOOM : 4, TOUCH_PAN : 5, PITCH_ROLL_YAW : 6, ZOOM_SMOOTH : 7 };
	var state = STATE.NONE;
	
	var DIRECTION = { VERTICAL : 0, HORIZONTAL : 1 };
	
	this.toConfig = function() {
		var config = {
			type : "ExamineControls",
			pitch : this._pitch,
			roll : this._roll,
			yaw : this._yaw,
			orientation : this.orientation.toArray(),
			distance : this.distance,
			scale : this.scale
		};
		
		return config;
	};
	
	function isValidConfigNumber(v) {
		return v !== undefined && v !== null && v !== NaN;	
	}
	
	this.fromConfig = function(view) {
		if (!view) {
			return;
		}
		
		this._pitch = (isValidConfigNumber(view.pitch)) ? view.pitch : this._pitch;
		this._roll = (isValidConfigNumber(view.roll)) ? view.roll : this._roll;
		this._yaw = (isValidConfigNumber(view.yaw)) ? view.yaw : this._yaw;
		this.distance = (isValidConfigNumber(view.distance)) ? view.distance : this.distance;
		this.scale = (isValidConfigNumber(view.scale)) ? view.scale : this.scale;
		if (view.orientation) {
			this.orientation.fromArray(view.orientation);
		}
		this._update();
	};
	
	this.reset = function()
	{
		this._pitch = 0;
		this._roll = 0;
		this._yaw = 0;
		this.orientation = new THREE.Quaternion();
		
		this._update();
	};
	
	
	this.rotate = function( rotateX, rotateY ) {
		lastRotateX = rotateX;
		lastRotateY = rotateY;
		
		rotateX *= this.rotateSpeed;
		rotateY *= this.rotateSpeed;
	
		var xAxis = new THREE.Vector3(1, 0, 0);
		var yAxis = new THREE.Vector3(0, 1, 0);
		
		xAxis.rotateX(this._pitch);
		xAxis.rotateZ(this._roll);
		yAxis.rotateX(this._pitch);
		yAxis.rotateZ(this._roll);
	
		var xRot = new THREE.Quaternion();
		xRot.setFromAxisAngle(xAxis, -rotateX);
		
		var yRot = new THREE.Quaternion();
		yRot.setFromAxisAngle(yAxis, -rotateY);
		
		var newRot = yRot.multiply(xRot);
		this.orientation = this.orientation.multiply(newRot);
		
		this._update();
	};
	
	this.update = function(force) {
		if (force) {
			this._update();
		} else {
			if (state == STATE.NONE && (lastRotateX || lastRotateY) ) {
			//	this.rotate(lastRotateX, lastRotateY);
			}
		}
	};
	
	
	this._update = function(skipEventDispatch) {
		
		if (this.object instanceof Array) {
			for (var i = 0; i < this.object.length; i++) {
				this._updateObject(this.object[i]);
			}
		} else {
			this._updateObject(this.object);
		}
		
		if (!skipEventDispatch) {
			this.dispatchEvent( changeEvent );
		}
	}
	
	this._updateObject = function(object) {
	
		var translateMatrix = new THREE.Matrix4();

		matrixPitch.identity().makeRotationX(this._pitch);
		matrixYaw.identity().makeRotationY(this._yaw);
		matrixRoll.identity().makeRotationZ(this._roll);
		matrixRoll.multiply(matrixPitch);
		this.modelView.identity();
		
		
		var m = new THREE.Matrix4();
		m.identity();
		m.makeRotationFromQuaternion(this.orientation);
		this.modelView.multiply(m);
		
		if (this.pitchType == KMG.SurfacePivot) {
			translateMatrix.makeTranslation(0, 0, this.radius);
			this.modelView.multiply( translateMatrix );
		}
		
		this.modelView.multiply( matrixYaw );
		this.modelView.multiply( matrixRoll );
		
		if (this.pitchType == KMG.SurfacePivot) {
			translateMatrix.makeTranslation(0, 0, -this.radius);
			this.modelView.multiply( translateMatrix );
		}
		
		if (!object.forceDefaultDistance) {
			translateMatrix.makeTranslation(0, 0, this.distance);
			this.modelView.multiply( translateMatrix );
		} else {
			translateMatrix.makeTranslation(0, 0, this.defaultDistance);
			this.modelView.multiply( translateMatrix );
		}

		translateMatrix.makeTranslation(0, this.panVertical, 0);
		this.modelView.multiply( translateMatrix );
		
		translateMatrix.makeTranslation(this.panHorizontal, 0, 0);
		this.modelView.multiply( translateMatrix );
		
		//
		if (this.translate) {
			translateMatrix.makeTranslation(this.translate.x, this.translate.y, this.translate.z);
			this.modelView.multiply( translateMatrix );
		}
		
		object.matrix.identity();
		object.applyMatrix(this.modelView);

		
	};
	
	this.eyeDistanceToCenter = function() {
		var position = new THREE.Vector3(0.0, 0.0, this.radius);
		position.rotate(this.pitch, "X");
		position.negate();
		
		var a = this.distance + position.z;
		
		var distanceToCenter = Math.sqrt((position.y * position.y) + (a * a));
		return distanceToCenter;
	};
	
	this.eyeDistanceToSurface = function() {
		var distanceToSurface = this.eyeDistanceToCenter();// - (this.radius * this.scale);
		return distanceToSurface;
	};
	
	this.setScale = function( scale ) {
		if (scale > this.maxScale)
			scale = this.maxScale;
		if (scale < this.minScale)
			scale = this.minScale;
		this.scale = scale;
	};
	
	this.setMinScale = function( minScale ) {
		this.minScale = minScale;
		if (this.scale < minScale)
			this.scale = minScale;
	};
	
	this.setMaxScale = function( maxScale ) {
		this.maxScale = maxScale;
		if (this.scale > maxScale)
			this.scale = maxScale;
	};
	
	this.setFov = function( fov ) {
		if (this.maxFov > fov)
			fov = this.maxFov;
		if (this.minFov < fov)
			fov = this.minFov;
		this.fov = fov;
	};
	
	this.pitch = function ( pitch ) {
		this.setPitch(this._pitch + (pitch * this.rotateSpeed));
	};
	
	this.setPitch = function( pitch ) {
		if (pitch > this.maxPitch) 
			pitch = this.maxPitch;
		if (pitch < this.minPitch)
			pitch = this.minPitch;
		this._pitch = pitch;
		this._update();
	};
	
	this.setMinPitch = function( minPitch ) {
		this.minPitch = minPitch;
		if (this._pitch < minPitch)
			this._pitch = minPitch;
	};
	
	this.setMaxPitch = function( maxPitch ) {
		this.maxPitch = maxPitch;
		if (this._pitch > maxPitch)
			this._pitch = maxPitch;
	};
	
	this.roll = function ( roll ) {
		this.setRoll(this._roll + (roll * this.rotateSpeed));
	};
	
	this.setRoll = function( roll ) {
		this._roll = roll;
		this._update();
	};
	
	
	this.setDistance = function( distance ) {
		if (distance > this.maxDistance)
			distance = this.maxDistance;
		if (distance < this.minDistance)
			distance = this.minDistance;
		this.distance = distance;
		this._update();
	};
	
	this.setMinDistance = function( minDistance ) {
		this.minDistance = minDistance;
		if (this.distance < minDistance)
			this.distance = minDistance;
	};
	
	this.setMaxDistance = function( maxDistance ) {
		this.maxDistance = maxDistance;
		if (this.distance > maxDistance)
			this.distance = maxDistance;
	};
	
	
	this.pan = function(amount, direction) {
		if (!direction || direction === DIRECTION.VERTICAL) {
			this.panVertical = this.panVertical + amount;
		} else if (direction === DIRECTION.HORIZONTAL) {
			this.panHorizontal = this.panHorizontal + amount;
		}
		this._update();
	};
	
	
	function _adjustDistanceByDelta(delta) {
		
		if (scope.zoomType == KMG.Distance) {
		
			// Adjust the distance by a proportional amount every time
			var ratio = scope.distanceMoveSpeed / scope.defaultDistance;
			var distanceMove = ratio * scope.distance;
			
			scope.setDistance(scope.distance + (-delta * distanceMove));
			
		} else if (scope.zoomType == KMG.Scale) {
			
			// Adjust the distance by a proportional amount every time
			var ratio = scope.zoomSpeed / scope.defaultScale;
			var scaleMove = ratio * scope.scale;
			
			scope.setScale(scope.scale + (delta * scaleMove));
		} else if (scope.zoomType == KMG.FoV) {
			
		}
	}
	
	// Events
	function onMouseDown( event ) {
		if (!scope.enabled) return;
		
		lastX = event.clientX;
		lastY = event.clientY;
		
		mouseDownX = event.clientX;
		mouseDownY = event.clientY;
		
		lastRotateX = 0;
		lastRotateY = 0;
		
		// Should be:
		// Left Button -> Rotate
		// Shift + Left Button -> Pitch/Roll
		// Ctrl + Left Button -> Pan
		// Middle Button -> Pitch/Roll
		// Right Button -> Zoom

		if ( event.button === 0 && !event.ctrlKey && !event.shiftKey) {
			state = STATE.ROTATE;
		} else if (event.button == 0 && event.ctrlKey) {
			state = STATE.PAN;
		} else if (event.button == 0 && event.shiftKey) {
			state = STATE.PITCH_ROLL_YAW;
		} else if ( event.button === 1 ) {
			state = STATE.PITCH_ROLL_YAW;
		} else if ( event.button === 2) {
			state = STATE.ZOOM_SMOOTH;
		} 

		document.addEventListener( 'mousemove', onMouseMove, false );
		document.addEventListener( 'mouseup', onMouseUp, false );
	}
	
	function onMouseMove( event ) {
		if (!scope.enabled) return;
	
		event.preventDefault();
		
		if (state === STATE.NONE) {
			return;
		}
		
		var xDelta = event.clientX - lastX;
		var yDelta = event.clientY - lastY;
		
		if ( state === STATE.ROTATE ) {
			scope.rotate( (yDelta * (Math.PI / 180)), (xDelta * (Math.PI / 180)) );
		} else if ( state === STATE.ZOOM ) {
			_adjustDistanceByDelta(-yDelta);
		} else if ( state === STATE.ZOOM_SMOOTH) {
			_adjustDistanceByDelta(event.clientY - mouseDownY);
		} else if ( state === STATE.PAN ) {
			scope.pan(yDelta, DIRECTION.VERTICAL);
			scope.pan(-xDelta, DIRECTION.HORIZONTAL);
		} else if ( state === STATE.PITCH_ROLL_YAW ) {
			scope.pitch(yDelta * (Math.PI / 180));
			scope.roll(xDelta * (Math.PI / 180));
		}
		
		lastX = event.clientX;
		lastY = event.clientY;
		
		if (onChange) {
			onChange(scope);
		}
	}
	
	function onMouseUp( event ) {
		if (!scope.enabled) return;
		
		document.removeEventListener( 'mousemove', onMouseMove, false );
		document.removeEventListener( 'mouseup', onMouseUp, false );
		
		
		mouseDownX = -1;
		mouseDownY = -1;
		lastX = -1;
		lastY = -1;
		state = STATE.NONE;
		
		if (onChange) {
			onChange(scope);
		}
	}
	
	function onMouseWheel( event ) {
		if ( scope.enabled === false ) return;
		
		var delta = 0;

		if ( event.wheelDelta ) { // WebKit / Opera / Explorer 9
			delta = event.wheelDelta;
		} else if ( event.detail ) { // Firefox
			delta = -event.detail * 20;
		}
		
		_adjustDistanceByDelta(delta);
		
		if (onChange) {
			onChange(scope);
		}
	}
	
	function onKeyDown( event ) {
		
	}
	
	
	
	
	this.domElement.addEventListener( 'contextmenu', function ( event ) { event.preventDefault(); }, false );
	this.domElement.addEventListener( 'mousedown', onMouseDown, false );
	this.domElement.addEventListener( 'mousewheel', onMouseWheel, false );
	this.domElement.addEventListener( 'DOMMouseScroll', onMouseWheel, false ); // firefox
	this.domElement.addEventListener( 'keydown', onKeyDown, false );

};

KMG.ExamineControls.prototype = Object.create( THREE.EventDispatcher.prototype );

/* File: GUI.js */
/**
 * PlanetMaker JavaScript Controller API
 * http://planetmaker.wthr.us
 * 
 * Copyright 2013 Kevin M. Gill
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */
 
KMG.TopLeft = "TL";
KMG.TopRight = "TR";
KMG.BottomLeft = "BL";
KMG.BottomRight = "BR";

KMG.Opened = 0;
KMG.Closed = 1;

//http://stackoverflow.com/questions/105034/how-to-create-a-guid-uuid-in-javascript
KMG.GUID = {

	guid : function() {
		return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
			var r = Math.random()*16|0, v = c == 'x' ? r : (r&0x3|0x8);
			return v.toString(16);
		});
	}

};


KMG.OptionController = function(property, title, setValuesInterface, updateInterface) {

	var changeListeners = [];
	
	this.addChangeListener = function(changeListener) {
		changeListeners.push(changeListener);
		return this;
	};
	
	this.setValues = function(values) {
		if (setValuesInterface) {
			setValuesInterface(values);
		}
		return this;
	};
	
	this.update = function() {
		if (updateInterface) {
			updateInterface();
		}
		return this;	
	};
	
	this.onChange = function(oldValue, newValue) {
		
		for (var i = 0; i < changeListeners.length; i++) {
			changeListeners[i](property, title, oldValue, newValue);
		}
	
	};
	
	
	

};

/**
 * @class Provides an individual block of controllers. This is the base component as the user
 * will see it (though not the base in terms of the DOM). Individual controls will
 * be added to instances of this object.
 *
 * @param {string} title Title of the block
 * @param {Object} config Property object backing the controls
 * @param {function} changeListener Callback function fired when a control's value is modified by the user
 *
 * @member KMG
 */
KMG.Block = function ( title, config, changeListener ) {
	
	this.config = config;
	this.changeListener = changeListener;
	var scope = this;
	
	var expandState = KMG.Opened;
	
	var container = $("<div/>").addClass("control-outter-container");
	
	var id = KMG.GUID.guid();
	var element = $("<div/>");
	element.attr("id", id);
	element.addClass("control-container");
	element.appendTo(container);
	
	var expandIcon = $("<a/>").addClass("expand-icon")
							.attr("href", "#")
							.text("-")
							.attr("title", "Click to expand or retract")
							.click(function() { 
								scope.setExpandedState(
									(expandState == KMG.Opened) ? KMG.Closed : KMG.Opened
								);
							}).appendTo(element);
	
	element.append($("<div class='header-title'>" + title + "</div>"));

	var controlContainer = $("<div/>").addClass("inner-container").appendTo(element);

	var controlList = $("<ul/>");
	controlList.appendTo(controlContainer);
	
	var scrollBar = $("<div/>").addClass("control-container-scrollbar-vertical").appendTo(controlContainer);
	$("<hr/>").appendTo(scrollBar);
	
	var lastY = -1;
	var mouseDown = false;
	scrollBar.mousedown(function(e) {
		mouseDown = true;
		lastY = e.clientY;
		scrollBar.addClass("unselectable");
		//console.debug('Mouse Down');
	});
	
	$(document).mouseup(function(e) {
		mouseDown = false;
		scrollBar.removeClass("unselectable");
		//console.debug('Mouse Up');
	});
	
	controlList.hasScrollBar = function() {
        return this.get(0).scrollHeight > this.height();
    };

	
	$(document).mousemove(function(e) {
		if (!mouseDown || (lastY <= e.clientY && !controlList.hasScrollBar())) {
			return;
		}
		e.preventDefault();
		var height = e.clientY - controlList.offset().top - 5;
		//console.debug("Setting height to " + height);
		controlList.css("height", height);
		
		lastY = e.clientY;
	});
	
	function fireChangeListener() {
		if (changeListener != null) {
			changeListener();
		}
	}
	
	this.getElement = function()
	{
		return container;
	}
	
	function isValueFloatingPoint(n) {
		return !isNaN(parseFloat(n)) && isFinite(n);
	}
	
	/**
     * 
     * @param property
	 * @param title
	 * @param min
	 * @param max
	 * @param step
     */
	this.addRange = function(property, title, min, max, step) {
		if (!min)
			min = 0;
		if (!max) 
			max = 100;
		if (!step) 
			step = 1;
	
		var slider, text;
		
		var controller = new KMG.OptionController(property,
													title, 
													null, // Set Values
													function() { // Update
			
			text.val(config[property]);
			slider.slider({value:config[property]});
		});
		
		var onSlide = function( event, ui ) {
			var oldValue = config[property];
			config[property] = ui.value;
			text.val(ui.value);
			controller.onChange(oldValue, config[property]);
			fireChangeListener();
		};
	
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		$("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);
		slider = $("<div/>").attr("id", id)
				.addClass("slider-control")
				.slider({
					value : config[property],
					min : min,
					max : max,
					step : step,
					range : "min",
					slide : onSlide,
					change : onSlide
				})
				.appendTo(li);
		text = $("<input/>").attr("type", "text")
					.addClass("value-text")
					.css("width", "40px")
					.css("margin-left", "10px")
					.on('input', function(e) {

						var value = $(this).val();
						
						// If the user blanked the input, let them enter something new before we start messing with it
						if (value === "") {
							return;
						}
						
						// Validate that the input value is a number
						if (!isValueFloatingPoint(value)) {
							$(this).val(config[property]);
							return;
						}

						// Validate that the input value falls within proper ranges
						if (value < min) {
							value = min;
							$(this).val(value);
						} else if (value > max) {
							value = max;
							$(this).val(value);
						}
						
						slider.slider("value", value);
						fireChangeListener();
					}).val(config[property])
					.appendTo(li);
		li.appendTo(controlList);
		
		return controller;
	};
	

	
	/**
     * 
     * @param property
	 * @param title
	 * @param options
     */
	this.addSelect = function(property, title, options) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		$("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);
		var select = $("<select id='" + id + "'></select>").appendTo(li);
		
		var controller = new KMG.OptionController(property, title, function(values) {
			$(select).empty();
			$.each(values, function(key, value) {
				var option = $("<option/>").attr("value", value).text(value);
				if (value === config[property]) {
					option.attr("selected", "true");
				}
				option.appendTo(select);
			});
		});
		
		select.change(function(e) {
			var oldValue = config[property];
			config[property] = $(select).val();
			controller.onChange(oldValue, config[property]);
			fireChangeListener();
		});
		
		if (options) {
			controller.setValues(options);
		}
		li.appendTo(controlList);
		
		return controller;
	};
	
	/**
     * 
     * @param property
	 * @param title
     */
	this.addToggle = function(property, title) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		var oldValue;
		
		$("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);
		var check = $("<input type='checkbox'/>").attr("id", id);
		if (config[property] === true) {
			check.attr("checked", "true");
		}

		var controller = new KMG.OptionController(property, title);
		
		check.click(function(e) {
			oldValue = config[property];
			config[property] = check.prop('checked');
			controller.onChange(oldValue, config[property]);
			//console.debug("Setting '" + property + "' to '" + config[property] + "'");
			fireChangeListener();
		});
		check.appendTo(li);
		
		li.appendTo(controlList);
		return controller;
	};
	
	/**
     * 
     * @param array
     */
	function arrayToColor(array)
	{
		var r = parseInt(array[0]);
		var g = parseInt(array[1]);
		var b = parseInt(array[2]);
		var rgb = "rgb("+r+","+g+","+b+")";
		return rgb;
	}
	
	// https://github.com/vanderlee/colorpicker
	/**
     * 
     * @param property
	 * @param title
     */
	this.addColor = function(property, title) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		$("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);
		$("<input/>").attr("id", id)
					.attr("type", "text")
					.val(arrayToColor(config[property]))
					.colorpicker({
						colorFormat : "RGB",
						color : arrayToColor(config[property]),
						ok: function(event, color) {
							config[property] = [color.rgb.r*255, color.rgb.g*255, color.rgb.b*255];
							//console.debug("Set '" + property + "' to '" + config[property] + "'");
							fireChangeListener();
						}
					}).appendTo(li);

		li.appendTo(controlList);
	};
	
	/**
     * 
     */
	function getLocalTimeZoneOffsetMillis()
	{
		var dt = new Date();
		return dt.getTimezoneOffset() * 60000;
	}
	
	// http://trentrichardson.com/examples/timepicker/
	/**
     * 
     * @param property
	 * @param title
     */
	this.addDateTime = function(property, title) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		
		$("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);
		
		var picker = $("<input/>").attr("type", "text")
								.attr("id", id)
								.datetimepicker({
									showButtonPanel: true,
									changeMonth: true,
									changeYear: true,
									yearRange : "c-25:c+25",
									onSelect : function(dateText, el) {
										config[property] = picker.datetimepicker('getDate').getTime();// - getLocalTimeZoneOffsetMillis();
										//console.debug("Setting '" + property + "' to '" + (new Date(config[property] - getLocalTimeZoneOffsetMillis())) + "'");
										fireChangeListener();
									}
								}).appendTo(li);
		if (config[property]) {
			picker.datetimepicker('setDate', (new Date(config[property])) );
		}
		
		li.appendTo(controlList);
	};
	
	
	
	/**
     * 
     * @param property
	 * @param title
     */
	this.addDate = function(property, title) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		
		$("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);
		
		var picker = $("<input/>").attr("type", "text")
								.attr("id", id)
								.datepicker({
									showButtonPanel: true,
									changeMonth: true,
									changeYear: true,
									yearRange : "c-25:c+25",
									onSelect : function(dateText, el) {
										config[property] = picker.datepicker('getDate').getTime() - getLocalTimeZoneOffsetMillis();
										//console.debug("Setting '" + property + "' to '" + (new Date(config[property] - getLocalTimeZoneOffsetMillis())) + "'");
										fireChangeListener();
									}
								}).appendTo(li);
		if (config[property]) {
			picker.datetimepicker('setDate', (new Date(config[property] + getLocalTimeZoneOffsetMillis())) );
		}
		
		li.appendTo(controlList);
	};
	
	/**
     * 
     * @param input
	 * @param label
	 * @param property
	 * @param validation
     */
	function onTextInput(input, label, property, validation) {
		var valid = true;
		if (validation && typeof validation === 'function' ) {
			valid = validation(input.val());
		} else if (validation && validation instanceof RegExp) {
			valid = validation.test(input.val());
		}
		if (valid) {
			input.removeClass('invalid-value-input');
			label.removeClass('invalid-value-label');
			config[property] = input.val();
		} else {
			input.addClass('invalid-value-input');
			label.addClass('invalid-value-label');
		}
	}
	
	/**
     * 
     * @param property
	 * @param title
	 * @param validation
     */
	this.addText = function(property, title, validation) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		
		var label = $("<label for='" + id + "'>"+title+"</label>").addClass("control-label").appendTo(li);

		var text = $("<input/>").attr("type", "text")
								.on('input', function(e) {
									onTextInput($(this), label, property, validation);
								}).val(config[property]).appendTo(li);
		
		li.appendTo(controlList);
		
		onTextInput(text, label, property, validation);
	};
	
	/**
     * 
     * @param title
	 * @param callback
     */
	this.addAction = function(title, callback) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");

		$("<button/>").text(title)
					.attr("id", id)
					.button()
					.click(function(e) {
						callback(e, $(this));
					}).appendTo(li);
		
		li.appendTo(controlList);
	};
	
	/**
     * 
     * @param href
	 * @param text
	 * @param title
     */
	this.addLink = function(href, text, title, target) {
		if (!title) {
			title = text;
		}
		if (!target) {
			target = "_self";
		}
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		$("<a/>").text(text)
				.attr("id", id)
				.attr("href", href)
				.attr("title", title)
				.attr("target", target)
				.appendTo(li);
		
		li.appendTo(controlList);
	};

	/**
     * 
     * @param el
     */
	this.addElement = function(el) {
		var id = KMG.GUID.guid();
		var li = $("<li/>");
		$(el).appendTo(li);
		li.appendTo(controlList);
	};
	
	/**
     * 
     * @param anchor
	 * @param x
	 * @param y
     */
	this.setPosition = function(anchor, x, y) {
		switch(anchor[0]) {
		case "T":
			element.css("top", y);
			break;
		case "B":
			element.css("bottom", y);
			break;
		};
		
		switch(anchor[1]) {
		case "L":
			element.css("left", x);
			break;
		case "R":
			element.css("right", x);
			break;
		};
	};
	
	/**
     * 
     * @param visible
     */
	this.setVisible = function(visible) {
		element.css("display", (visible ? "inline-block" : "none"));
	};
	
	/**
     * 
     * @param state
     */
	this.setExpandedState = function(state) {
		expandState = state;
		if (state == KMG.Opened) {
			controlContainer.css("display", "inline-block");
			element.removeClass().addClass("control-container");
			controlContainer.addClass("inner-container");
			expandIcon.text("-");
		} else {
			controlContainer.css("display", "none");
			expandIcon.text("+");
		}
		
	};
	
};



/**
 * @class 
 *
 * @param {Object} config 
 * @param {function} changeListener
 *
 * @member KMG
 */
KMG.SideBar = function(config, changeListener, addClasses) {
	this.defaultConfig = config;
	this.defaultChangeListener = changeListener;
	var scope = this;
	
	var id = KMG.GUID.guid();
	var element = $("<div/>");
	element.attr("id", id);
	element.addClass("control-sidebar");
	//element.css("display", "inline");
	//element.css("position", "absolute");
	element.appendTo('body');
	
	for (var i = 0; i < addClasses.length; i++) {
		element.addClass(addClasses[i]);
	}
	
	/**
     * 
     * @param block
     */
	this.addBlock = function(block) {
		element.append(block.getElement());
	};
	
	this.removeBlock = function(block) {
		block.getElement().remove();
	};
	
	/**
     * 
     * @param title
	 * @param config
	 * @param changeListener
     */
	this.createBlock = function( title, config, changeListener ) {
		if (!config) {
			config = this.defaultConfig;
		}
		
		if (!changeListener) {
			changeListener = this.defaultChangeListener;
		}
		
		var block = new KMG.Block(title, config, changeListener);
		this.addBlock(block);
		return block;
	};
	

	
	
	/**
     * 
     * @param anchor
	 * @param x
	 * @param y
     */
	this.setPosition = function(anchor, x, y) {
		
		
		switch(anchor[0]) {
		case "T":
			//element.css("top", y);
			break;
		case "B":
			//element.css("bottom", y);
			break;
		};
		
		switch(anchor[1]) {
		case "L":
			element.addClass("control-sidebar-left");
			break;
		case "R":
			element.addClass("control-sidebar-right");
			break;
		};
		
	};
	
	/**
     * 
     * @param opacity
     */
	this.setOpacity = function(opacity) {
		element.css("opacity", opacity);
	};
	
	/**
     * 
     * @param visible
     */
	this.setVisible = function(visible) {
		element.css("display", (visible ? "inline-block" : "none"));
	};
	
	this.setHeight = function(height) {
		element.css("height", height);
	};
};


/**
 * @class 
 *
 * @param {Object} config 
 * @param {function} changeListener
 *
 * @member KMG
 */
KMG.GUI = function(config, changeListener) {

	
	var scope = this;
	
	this.onVisibilityChanged = null;
	
	this.left = new KMG.SideBar(config, changeListener, ["control-sidebar-left"]);
	this.right = new KMG.SideBar(config, changeListener, ["control-sidebar-right"]);
	
	this.left.setPosition("TL", "10px", "100px");
	this.right.setPosition("TR", "10px", "100px");
	
	var showGui = $("<div/>").addClass('control-show-gui')
							.css('display', 'none')
							.appendTo('#container');
	$("<a/>").attr("href", "#")
			.text('Show Controls')
			.on('click', function(e) {
				scope.setVisible(true);
			}).appendTo(showGui);
	
	function updateHeight() {
		scope.left.setHeight(($(window).height() - 100) + "px");
		scope.right.setHeight(($(window).height() - 100) + "px");
	}
	
	$(window).resize(function() {
		updateHeight();
	});
	updateHeight();
	
	
	/**
     * 
     * @param opacity
     */
	this.setOpacity = function(opacity) {
		this.left.setOpacity(opacity);
		this.right.setOpacity(opacity);
	};
	
	/**
     * 
     * @param visible
     */
	this.setVisible = function(visible, suppressShowBlock) {
		this.left.setVisible(visible);
		this.right.setVisible(visible);
		
		if (!suppressShowBlock) {
			showGui.css('display', (visible ? 'none' : 'inline-block'));
		}
		
		if (this.onVisibilityChanged) {
			this.onVisibilityChanged(visible);
		}
	};

	
};

/* File: AppEnv.js */

var AppEnv = {

	config : {
		devMode : false,
		embedded : false,
		disablePlanet : false,
		noAnalytics : false
	},
	
	_getConfig : function(urlParam, defaultValue) {
		var param = AppEnv.getUrlVar(urlParam);
		if (param) {
			return true;
		} else {
			return defaultValue;
		}
	
	},
	
	isDevMode : function() {
		return AppEnv._getConfig("devMode", AppEnv.config.devMode);
	},
	
	isEmbedded : function() {
		return AppEnv._getConfig("embedded", AppEnv.config.embedded);
	},
	
	isPlanetDisabled : function() {
		return AppEnv._getConfig("disablePlanet", AppEnv.config.disablePlanet);
	},
	
	noAnalytics : function() {
		return AppEnv._getConfig("noga", AppEnv.config.noAnalytics);
	},
	

	getUrlVars: function() {
		var vars = [], hash;
		var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
		for(var i = 0; i < hashes.length; i++) {
			hash = hashes[i].split('=');
			vars.push(hash[0]);
			vars[hash[0]] = hash[1];
		}
		return vars;
	},
	getUrlVar: function(name){
		return AppEnv.getUrlVars()[name];
	}

};

	