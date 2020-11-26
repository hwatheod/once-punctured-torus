const tolerance = 0.0001;

class Complex {
  constructor(re1, im1=0) {
    this.re = re1;
    this.im = im1;
  }

  static polar(r, theta) {
    return new Complex(r*Math.cos(theta), r*Math.sin(theta));
  }

  cadd(z) {
    return new Complex(this.re + z.re, this.im + z.im);
  }

  csubt(z) {
    return new Complex(this.re - z.re, this.im - z.im);
  }

  cmult(z) {
    return new Complex(this.re * z.re - this.im * z.im, this.re * z.im + this.im * z.re);
  }

  cmult_scalar(x) {
    return new Complex(this.re * x, this.im * x);
  }

  cdiv(z) {
    // (a+bi)(c-di)
    const conjMult = this.cmult(z.cconj());
    const norm_squared = z.re * z.re + z.im * z.im;
    return new Complex(conjMult.re / norm_squared, conjMult.im / norm_squared);
  }

  cconj() {
    return new Complex(this.re, -this.im);
  }

  cnorm() {
    return Math.sqrt(this.re * this.re + this.im * this.im);
  }

  cnorm_squared() {
    return this.re * this.re + this.im * this.im;
  }

  cneg() {
    return new Complex(-this.re, -this.im);
  }

  crecip () {
    const norm_squared = this.re * this.re + this.im * this.im;
    return new Complex(this.re / norm_squared, -this.im / norm_squared);
  }

  cnegrecip () {
    const norm_squared = this.re * this.re + this.im * this.im;
    return new Complex(-this.re / norm_squared, this.im / norm_squared);
  }

  cunitInvert () {
    /* The inverse of z with respect to the unit circle, i.e. 1/z' */

    const norm_squared = this.re * this.re + this.im * this.im;
    return new Complex(this.re / norm_squared, this.im / norm_squared);
  }

  carg() {
    return Math.atan2(this.im, this.re);
  }

  csqrt() {
    const r = this.cnorm();
    const theta = this.carg();
    return Complex.polar(Math.sqrt(r), theta/2);
  }

  equals(z2) {
    return (this.re == z2.re && this.im == z2.im);
  }

  equalsReal(x) {
    return (this.re == x && this.im == 0);
  }

  toString() {
    return this.re + " + " + this.im + " i";
  }
}

class MobiusTransformation {
  constructor(_a = new Complex(1), _b = new Complex(0), _c = new Complex(0), _d = new Complex(1)) {
    this.a = _a;
    this.b = _b;
    this.c = _c;
    this.d = _d;
  }

  normalize() {
    const _a = this.a;
    const _b = this.b;
    const _c = this.c;
    const _d = this.d;
    const det = _a.cmult(_d).csubt(_b.cmult(_c));
    const sqrt_det = det.csqrt();
    this.a = _a.cdiv(sqrt_det);
    this.b = _b.cdiv(sqrt_det);
    this.c = _c.cdiv(sqrt_det);
    this.d = _d.cdiv(sqrt_det);
  }

  apply(z) {
    return this.a.cmult(z).cadd(this.b).cdiv(this.c.cmult(z).cadd(this.d));
  }

  invert_normalized() {
    return new MobiusTransformation(this.d, this.b.cneg(), this.c.cneg(), this.a);
  }

  rightMultiply(m) {
    return new MobiusTransformation(this.a.cmult(m.a).cadd(this.b.cmult(m.c)), this.a.cmult(m.b).cadd(this.b.cmult(m.d)), this.c.cmult(m.a).cadd(this.d.cmult(m.c)), this.c.cmult(m.b).cadd(this.d.cmult(m.d)));
  }

  parabolicFixpoint() {
    // we don't verify parabolicity here because in some cases, we need to specify nearly parabolic words as special words.
    // See Indra's Pearls, p. 265.

    return this.a.csubt(this.d).cdiv(this.c.cmult_scalar(2));  /* (a-d)/(2c), returns NaN for point at infinity */
  }

  toString() {
    return (this.a.toString() + "," + this.b.toString() + "," + this.c.toString() + "," + this.d.toString());
  }
}

class Plane {
  constructor(_xmin, _xmax, _ymin, _ymax) {
    this.xmin = _xmin;
    this.ymin = _ymin;
    this.xmax = _xmax;
    this.ymax = _ymax;
    this.planeXtoScreenXFactor = width / (this.xmax - this.xmin);
    this.planeYtoScreenYFactor = height / (this.ymax - this.ymin);
  }

  toPoint(z) {
    return [Math.round((z.re - this.xmin) * this.planeXtoScreenXFactor), Math.round(height - (z.im - this.ymin) * this.planeYtoScreenYFactor)];
    // y-coordinate sign is reversed because graphics coordinates decrease going up, but standard Cartesian coordinates increase going up
  }

  fromPoint(p) {
    return new Complex(this.xmin + (p[0] / this.planeXtoScreenXFactor), this.ymax - (p[1] / this.planeYtoScreenYFactor));
    // y-coordinate sign is reversed because graphics coordinates decrease going up, but standard Cartesian coordinates increase going up
  }

  plot(z) {
    const p = this.toPoint(z);
    point(p[0], p[1]);
  }

  plotLine(oldPoint, newPoint) {
     const oldScreenPoint = this.toPoint(oldPoint);
     const newScreenPoint = this.toPoint(newPoint);

     /* Do nothing if we're going from a point to itself, and if new point is adjacent to old point, just plot the new point
        rather than drawing a line. This increases performance by a lot, since the graphics rendering forms a
        significant part of the run time of the algorithm. */
     if (!(oldScreenPoint[0] == newScreenPoint[0] && oldScreenPoint[1] == newScreenPoint[1])) {
       if (Math.abs(oldScreenPoint[0] - newScreenPoint[0]) <= 1 && Math.abs(oldScreenPoint[1] - newScreenPoint[1]) <= 1) {
         point(newScreenPoint[0], newScreenPoint[1]);
       }
       else {
         line(oldScreenPoint[0], oldScreenPoint[1], newScreenPoint[0], newScreenPoint[1]);
       }
       return 1;
     }
     return 0;
  }
}

class QuasiFuchsianPlotRoutine {
  constructor (a, b, terminationThreshold, maxDepth, specialWords) {
    this.transformList = [a, b, a.invert_normalized(), b.invert_normalized()];
    this.terminationThreshold = terminationThreshold;
    this.maxDepth = maxDepth;
    this.specialFixedPointList = new Array(4);
    const ok = this.setSpecialWords(specialWords); // also sets this.oldPoint to the beginning value
    if (!ok) { // parabolic fixed point at infinity
      this.terminateRun = true;
    } else {
      this.count = 0;
      this.wordLengthReached = 0;
      this.terminateRun = false;
    }
  }

  isRunTerminated() {
    return this.terminateRun;
  }

  static addAllCyclicShifts(specialWordsSplit, word) {
    let curWord = word;
    for (let i=0; i<curWord.length; i++) {
      const lastLetter = curWord.charAt(curWord.length - 1);
      const index = "abAB".indexOf(lastLetter);
      if (specialWordsSplit[index].indexOf(curWord) == -1)
        specialWordsSplit[index].push(curWord);
      curWord = curWord.substring(1) + curWord.charAt(0);  // cyclic shift, abc -> bca
    }
  }

  static inverse(s) {
     /* returns the inverse word: letters reversed with upper/lower case flipped too. */
    const letterList = "abAB";
    let result = "";
    for (let i = s.length - 1; i>=0; i--) {
       result += letterList.charAt((letterList.indexOf(s.charAt(i)) + 2) % 4);
    }
    return result;
  }

  plotFunction(plane, i, depth, curTransform) {
    let terminateBranch = true;
    let first = true;
    let prevPoint = this.oldPoint;
    const pointList = [];
    for (const fixedPoint of this.specialFixedPointList[i]) {
      if (first) {
        first = false; // skip first element since it should be the same as the last point plotted (Indra's Pearls, p. 185)
        continue;
      }
      const newPoint = curTransform.apply(fixedPoint);

      if (newPoint.csubt(prevPoint).cnorm() > this.terminationThreshold) {
        terminateBranch = false;
        break;
      }
      pointList.push(newPoint);

      prevPoint = newPoint;
    }

    if (terminateBranch) {
      const endPoint = prevPoint; // last point found in previous loop
      prevPoint = this.oldPoint;
      const curThis = this; // can't use "this" inside the forEach closure directly
      pointList.forEach(function(newPoint) {
        curThis.count += plane.plotLine(prevPoint, newPoint);
        prevPoint = newPoint;
      });
      this.oldPoint = endPoint;

      if (depth + 1 > this.wordLengthReached)
        this.wordLengthReached = depth + 1;
    }

    if (depth == this.maxDepth - 1) {
      if (!terminateBranch) { /* depth too long, we are in a non-discrete group or stuck at a parabolic fixed point */
        console.log("Terminating at matrix: " + curTransform);
        this.terminateRun = true;
      }
    }

    return terminateBranch;
  }

  setSpecialWords(specialWords) {
    this.specialWordsSplit = new Array(4);
    for (let i=0; i<4; i++) {
      this.specialWordsSplit[i] = [];  /* ith element is the list of special words ending in the ith symbol */
      this.specialFixedPointList[i] = []; /* ith element is the list of fixed points for the special words in specialWordsSplit[i] */
    }

    for (const word of specialWords) {
      QuasiFuchsianPlotRoutine.addAllCyclicShifts(this.specialWordsSplit, word);
      QuasiFuchsianPlotRoutine.addAllCyclicShifts(this.specialWordsSplit, QuasiFuchsianPlotRoutine.inverse(word));
    }

    QuasiFuchsianPlotRoutine.addAllCyclicShifts(this.specialWordsSplit, "abAB");
    QuasiFuchsianPlotRoutine.addAllCyclicShifts(this.specialWordsSplit, "BAba");
       
    const letterList = "abAB";
    for (let i=0; i<4; i++) {
       let lastIndex = i;
       this.specialWordsSplit[i].sort(function(s1, s2) {
         let prevIndex = lastIndex;
         for (let j=0; j<Math.min(s1.length, s2.length); j++) {
           const c1 = s1.charAt(j);
           const c2 = s2.charAt(j);
           if (c1 == c2) {
             prevIndex = letterList.indexOf(c1);
             continue;
           }
           const index1 = letterList.indexOf(c1);
           const index2 = letterList.indexOf(c2);

           /* ordering is: start at (prevIndex + 1) % 4 and count *backwards*. */
           let first = true;
           for (let k = (prevIndex + 1) % 4; k != (prevIndex + 1) % 4 || first; k = (k == 0 ? 3 : k-1)) {
              first = false;
              if (k == index1)
                 return -1;
              if (k == index2)
                 return 1;
           }
           throw ("Not supposed to reach here in special words comparator " + s1 + " and " + s2);
         }
         return s1.length - s2.length;
       });

       for (const word of this.specialWordsSplit[i]) {
         let transformationForWord = new MobiusTransformation();
         for (let j=0; j < word.length; j++) {
            const c = word.charAt(j);
            const index = letterList.indexOf(c);
            transformationForWord = transformationForWord.rightMultiply(this.transformList[index]);
         }
         const fixedPoint = transformationForWord.parabolicFixpoint();
         if (isNaN(fixedPoint.re)) {
           alert("Fixed point infinity for word: " + word + ". Changing the sign may help.");
           return false;
         }
         this.specialFixedPointList[i].push(fixedPoint);
       }
    }

    this.oldPoint = this.specialFixedPointList[0][0]; // This is the starting point.
    console.log(this.oldPoint);
    return true;
  }

  postPlotFunction() {
    if (this.terminateRun) {
      console.log("Run terminated early because maximum depth " + this.maxDepth + " was reached at oldPoint = " + this.oldPoint);
    }
    console.log("Points plotted:" + this.count);
    console.log("Maximum word length: " + this.wordLengthReached);
  }
}

dfsPlot = function(plane, plotRoutine, colorList) {
  let depth = 0;
  let word = new Array(plotRoutine.maxDepth);
  let i = -1;
  const curTransformList = new Array(plotRoutine.maxDepth);
  const id = new MobiusTransformation();
  let curTransform = id;
  let terminateBranch;
  let prevValue, doneValue;

  /*
     depth-first search, without explicit recursive calls.
     Invariants at the beginning of each iteration of the loop:
      word[j] for 0 <= j < depth: the jth transformation of the current word
      curTransformationList[j] for 0 <= j < depth: product (in order) from word[0] through word[j], inclusive.
      curTransform = product (in order) from word[0] to word[depth - 1] (identity if depth == 0).
      i = last transformation (0,1,2,3) tried at current depth, -1 if none.

      We need the curTransformationList[] array instead of multiplying by the inverse to undo the previous transform, because
      the latter process introduces numerical errors that eventually accumulate and cause problems with the drawing.
   */
  while (depth >= 0 && !plotRoutine.isRunTerminated()) {
    if (depth == 0) {
      prevValue = 0;
      doneValue = 1;
    } else {
      prevValue = word[depth - 1];
      doneValue = (prevValue + 2) % 4;
    }

    if (i == -1)
      i = (prevValue + 1) % 4;
    else {
      i = i-1;
      if (i<0) i=3;
      if (i == doneValue) {
        if (depth > 0) {
          i = word[depth - 1];
          if (depth > 1)
            curTransform = curTransformList[depth - 2];
          else curTransform = id;
        }
        depth--;
        continue;
      }
    }

    word[depth] = i;
    if (depth == 0) {
      const newColor = color(colorList[i]);
      /* Set both stroke and fill to the same color so that the p5.js code doesn't switch back and forth
         between them when plotting points.  See the prototype.point function in the p5.js code:
            https://github.com/processing/p5.js/blob/master/src/core/p5.Renderer2D.js
       */
      stroke(newColor);
      fill(newColor);
    }

    curTransform = curTransform.rightMultiply(plotRoutine.transformList[i]);
    curTransformList[depth] = curTransform;

    terminateBranch = plotRoutine.plotFunction(plane, i, depth, curTransform);

    if (depth == plotRoutine.maxDepth - 1) {
      terminateBranch = true;
    }

    if (terminateBranch) {
      if (depth > 0) {
        curTransform = curTransformList[depth - 1];
      } else {
        curTransform = id;
      }
    } else { /* next level */
      depth++;
      i = -1;
    }
  }

  plotRoutine.postPlotFunction();
};

oncePuncturedTorusLimitSetPlotter = function (ta, tb, sign, terminationThreshold, maxDepth, specialWords) {
  /* Reference, Indra's Pearls, p. 229 (Box 21, "Grandma's special parabolic commutator groups") */

  const ta2 = ta.cmult(ta);
  const tb2 = tb.cmult(tb);
  const disc = ta2.cmult(tb2).csubt(ta2.cadd(tb2).cmult_scalar(4));
  const tab = ta.cmult(tb).cadd(disc.csqrt().cmult_scalar(sign)).cmult_scalar(0.5); // need the right choice of sign in the quadratic formula in order to get the same pictures as in Indra's Pearls
  console.log("t_a = " + ta + ", t_b = " + tb + ", t_ab = " + tab);
  if (tab.csubt(new Complex(2)).cnorm() < tolerance) {
    throw ("t_ab is very close to 2, algorithm fails. Try the opposite sign, or replace one of the traces with 2.");
  }

  const z0_num = (tab.csubt(new Complex(2))).cmult(tb);
  const z0_den = (tb.cmult(tab)).csubt(ta.cmult_scalar(2)).cadd(tab.cmult(new Complex(0,2)));
  const z0 = z0_num.cdiv(z0_den);

  const ma = ta.cmult_scalar(0.5);
  const mb_num = ta.cmult(tab).csubt(tb.cmult_scalar(2)).cadd(new Complex(0,4));
  const mb_den = tab.cmult_scalar(2).cadd(new Complex(4)).cmult(z0);
  const mb = mb_num.cdiv(mb_den);
  const mc_num = ta.cmult(tab).csubt(tb.cmult_scalar(2)).csubt(new Complex(0,4)).cmult(z0);
  const mc_den = tab.cmult_scalar(2).csubt(new Complex(4));
  const mc = mc_num.cdiv(mc_den);

  const m = new MobiusTransformation(ma, mb, mc, ma);

  const m2a = tb.csubt(new Complex(0,2)).cmult_scalar(0.5);
  const m2b = tb.cmult_scalar(0.5);
  const m2d = tb.cadd(new Complex(0,2)).cmult_scalar(0.5);

  const m2 = new MobiusTransformation(m2a, m2b, m2b, m2d);

  m.normalize();
  m2.normalize();

  console.log("a = " + m);
  console.log("b = " + m2);
  const a0 = ta.cdiv(tab.cmult(tb));
  const a1 = tab.cdiv(tb.cmult(ta));
  const a2 = tb.cdiv(ta.cmult(tab));
  console.log("Complex probabilities = " + a0 + ", " + a1 + ", " + a2);
  console.log("z1, z2 = " + a0 + ", " + (a0.cadd(a1)));
  console.log("abAB = " + m.rightMultiply(m2).rightMultiply(m.invert_normalized()).rightMultiply(m2.invert_normalized()));

  const qfpr = new QuasiFuchsianPlotRoutine(m, m2, terminationThreshold, maxDepth, specialWords);
  if (!qfpr.terminateRun) { // parabolic fixed point at infinity
    return qfpr;
  }
  return null;
};

let plane = null;
let plotter = null;
/*
TODO

allow saving parameter sets and naming them

 */
function setup() {
  createCanvas(800, 800);
  background(255);
  noLoop();
}

function draw() {
  if (plane != null && plotter != null) {
    clear();
    background(255);
    const t0 = performance.now();
    dfsPlot(plane, plotter, ['#0000ff', '#c000c0', '#ff0000', '#00c000']);
    const t1 = performance.now();
    console.log("Time used: " + (t1 - t0)/1000 + " seconds");
    console.log("======");
  }
}

function validateNumeric(elementId, defaultValue = NaN) {
  const elt = document.getElementById(elementId);
  let result;
  if (elt.value.trim().length == 0) {
    result = defaultValue;
  } else {
    if (/^[-+]?[0-9]*(\.[0-9]*)?([eE][+-]?[0-9]+)?$/.test(elt.value)){
      result = parseFloat(elt.value);
    } else {
      result = NaN;
    }
  }

  if (isNaN(result)) {
    elt.style.color = "red";
  } else {
    elt.style.color = "black";
  }
  return result;
}

function setPredefinedDrawing() {
  const predefinedDrawings = {
    "None": {
      "taRe": "", "taIm": "", "tbRe": "", "tbIm": "", "sign": "-", "specialWords": "",
      "xMin": "", "xMax": "", "yMin": "", "yMax": "",
      "terminationThreshold": "", "maxDepth": ""
    },
    "fig_7_1": {
      "taRe": 2, "taIm": 0, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": "a,b",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_1": {
      "taRe": 1.87, "taIm": 0.1, "tbRe": 1.87, "tbIm": -0.1, "sign": "-", "specialWords": "",
      "xMin": -2.5, "xMax": 2.5, "yMin": -2.5, "yMax": 2.5,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_2_5": {
      "taRe": 2.2, "taIm": 0, "tbRe": 2.2, "tbIm": 0, "sign": "-", "specialWords": "",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_5_3": {
      "taRe": 2, "taIm": 0.5, "tbRe": 2, "tbIm": 0.5, "sign": "+", "specialWords": "",
      "xMin": -1.6, "xMax": 1.6, "yMin": -1.6, "yMax": 1.6,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_11": {
      "taRe": 1.91, "taIm": 0.05, "tbRe": 3, "tbIm": 0, "sign": "-", "specialWords": "",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_14": {
      "taRe": 3, "taIm": 0, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": "b",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_15_1": {
      "taRe": 2.0, "taIm": 0.05, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": "a,b",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_8_15_4": {
      "taRe": 1.887, "taIm": 0.05, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": "b",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.0005, "maxDepth": 5000
    },
    "fig_8_23": {
      "taRe": 1.96556, "taIm": -0.99536, "tbRe": 3, "tbIm": 0, "sign": "+", "specialWords": "aaB",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.001, "maxDepth": 2500
    },
    "fig_9_1": {
      "taRe": 1.95859, "taIm": -0.01128, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": `b,${"a".repeat(15)}B`,
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.005, "maxDepth": 2500
    },
    "fig_9_3": {
      "taRe": 1.64214, "taIm": -0.76659, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": "b,aaaBaaB",
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.005, "maxDepth": 2500
    },
    "fig_9_15_4": {
      "taRe": 1.90378, "taIm": -0.03958, "tbRe": 2, "tbIm": 0, "sign": "-", "specialWords": `b,${"a".repeat(10)}B${"a".repeat(9)}B`,
      "xMin": -1.1, "xMax": 1.1, "yMin": -1.1, "yMax": 1.1,
      "terminationThreshold": 0.005, "maxDepth": 2500
    }
  }

  const selection = document.getElementById('drawing').value;
  const drawing = predefinedDrawings[selection];
  for (key in drawing) {
    if (drawing.hasOwnProperty(key) && key != "sign") {
      const elt = document.getElementById(key);
      elt.value = drawing[key];
      elt.style.color = "black"; // in case it was previously red due to invalid input
    }
  }
  if (drawing["sign"] == "-") {
    document.getElementById("minusSign").checked = true;
  } else if (drawing["sign"] == "+") {
    document.getElementById("plusSign").checked = true;
  }
}

function plotLimitSet() {
  const taRe = validateNumeric("taRe", 0);
  const taIm = validateNumeric("taIm", 0);
  const tbRe = validateNumeric("tbRe", 0);
  const tbIm = validateNumeric("tbIm", 0);
  const xMin = validateNumeric("xMin");
  const xMax = validateNumeric("xMax");
  const yMin = validateNumeric("yMin");
  const yMax = validateNumeric("yMax");

  let sign = 1;
  if (document.getElementById('minusSign').checked) {
    sign = -1;
  }
  const terminationThreshold = validateNumeric("terminationThreshold");
  const maxDepth = validateNumeric("maxDepth");
  const specialWords = document.getElementById('specialWords').value.split(",").map(x => x.trim());
    
  if (isNaN(taRe) ||
  isNaN(taIm) ||
  isNaN(tbRe) ||
  isNaN(tbIm) ||
  isNaN(xMin) ||
  isNaN(xMax) ||
  isNaN(yMin) ||
  isNaN(yMax) ||
  isNaN(terminationThreshold) ||
  isNaN(maxDepth)) return;

  plane = new Plane(xMin, xMax, yMin, yMax);
  plotter = oncePuncturedTorusLimitSetPlotter(
    new Complex(taRe, taIm),
    new Complex(tbRe, tbIm),
    sign,
    terminationThreshold,
    maxDepth,
    specialWords);

  redraw();
}