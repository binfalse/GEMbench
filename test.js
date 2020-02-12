
const sources = ["Cell Line","Patient Data"]
const types = ["Microarray","RNA Seq","MS Proteomics"]
const imethods = ["GIMME","FASTCORE","INIT","iMAT"];

Object.values = Object.values || function(o){return Object.keys(o).map(function(k){return o[k]})};

function boxColor (str) {
  if (str.includes ("GIMME")) {
    return "#377eb8"
  } else if (str.includes ("FASTCORE")) {
    return "#4daf4a"
  } else if (str.includes ("INIT")) {
    return "#984ea3"
  }
  return "#ff7f00"
}
function datacolor (str) {
  if (str.includes ("Microarray")) {
    return "#fef0d9"
  } else if (str.includes ("MS Proteomics")) {
    return "#fdcc8a"
  }
  return "#fc8d59"
}

function dataset2sourcetype (dataset) {
  switch (dataset) {
      case "EMTAB-37":
        return [sources[0], types[0]]
      case "HPA":
        return [sources[0], types[1]]
      case "ProteomeNCI60":
        return [sources[0], types[2]]
      case "GSE2109":
        return [sources[1], types[0]]
      case "TCGA":
        return [sources[1], types[1]]
      case "ProteomePatients":
        return [sources[1], types[2]]
    }
  return ["unknown", "unknown"]
}

function sortSelect(selDom) {
    var options = [];
    for (var i=0;i<selDom.options.length;i++)
        options.push (selDom.options[i]);
    
    options.sort(function (a, b) {return a.value < b.value ? -1 : 1});
    for (var i=0;i<options.length;i++) {
        selDom.appendChild(options[i]);
    }
    return;
}

class Sample {
	constructor (name, source, type) {
		this.name = name
		this.scores = {}
    this.source = source
    this.type = type
	}
	setScore (measure, method, value) {
		this.scores[measure + "_" + method] = value
	}
}
function SampleSorter (scoreId) {
 return (a, b) => {return a.scores[scoreId]-b.scores[scoreId]}	
}

function sum (x) {
  var s = 0
  for (i = 0; i < x.length; i++)
    s += x[i]
  return s
}

function pearson (samples, metric1, metric2) {
  x = []
  y = []
  
  for (var i = 0; i < samples.length; i++){
    x.push (samples[i].scores[metric1])
    y.push (samples[i].scores[metric2])
  }
  let sumX = 0,
    sumY = 0,
    sumXY = 0,
    sumX2 = 0,
    sumY2 = 0;
  const minLength = x.length = y.length = Math.min(x.length, y.length),
    reduce = (xi, idx) => {
      const yi = y[idx];
      sumX += xi;
      sumY += yi;
      sumXY += xi * yi;
      sumX2 += xi * xi;
      sumY2 += yi * yi;
    }
  x.forEach(reduce);
  return (minLength * sumXY - sumX * sumY) / Math.sqrt((minLength * sumX2 - sumX * sumX) * (minLength * sumY2 - sumY * sumY));
}


function getMeth () {
  marr = {};
  for (i = 0; i < imethods.length; i++)
    marr[imethods[i]] = {
      "samples": [],
      "q1": 0,
      "median": 0,
      "q3": 0,
      "interQuantileRange": 0,
      "min": 0,
      "max": 0,
    };
  return marr;
};
function getType () {
  t = {}
  for (var i = 0; i < types.length; i++)
    t[types[i]] = getMeth ()
  return t
  //return {
    //"RNA Seq": getMeth (),
    //"Microarray": getMeth (),
    //"MS Proteomics": getMeth ()
  //};
};

function getSourceData () {
  s = {}
  for (var i = 0; i < sources.length; i++)
    s[sources[i]] = getType ()
  return s
	//return {
    //"Patient Data": getType (),
    //"Cell Line": getType ()
  //}
};




var margin = {top: 50, right: 0, bottom: 60, left: 40},
    width = d3.select("#my_dataviz").node().clientWidth - margin.left - margin.right,
    height = 800 - margin.top - margin.bottom;


boxplots = {}

//var xdomain = new Set ();
drawn = false
samples = {}
for (var i = 0; i < sources.length; i++)
for (var j = 0; j < types.length; j++) {
  samples[sources[i] + "__" + types[j] + "__single"] = {}
  samples[sources[i] + "__" + types[j] + "__double"] = {}
}
//MEASURE = "AFR"
//MEASURE = "EOR"
//MEASURE = "Hallmark"

//var slider = document.getElementById("myRange");
//console.log (slider.value)

// Update the current slider value (each time you drag the slider handle)
/*slider.oninput = function() {
  console.log (this.value);
} */

function draw_boxplots2 () {
drawn = true
  //sqrtscale = document.getElementById('sqrt').checked
  infinites_p = false
  infinites_m = false
  minY = 10000000;
  maxY = 0;
  MEASURE = document.getElementById("metric")
  MEASURE = MEASURE.options[MEASURE.selectedIndex].value
  xdomain = new Set ();
      
      d3.select("#my_dataviz svg").remove();
  
  
  
  
var marginWhole = {top: 10, right: 10, bottom: 10, left: 10},
    sizeWhole = 640 - marginWhole.left - marginWhole.right

// Create the svg area
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", sizeWhole  + marginWhole.left + marginWhole.right)
    .attr("height", sizeWhole  + marginWhole.top + marginWhole.bottom)
  .append("g")
    .attr("transform", "translate(" + marginWhole.left + "," + marginWhole.top + ")");
    
    
    var allVar = imethods
  var numVar = allVar.length
    mar = 20
  size = sizeWhole / numVar
  
  var position = d3.scalePoint()
    .domain(allVar)
    .range([0, sizeWhole-size])

  // Color scale: give me a specie name, I return a color
  //var color = d3.scaleOrdinal()
    //.domain(["setosa", "versicolor", "virginica" ])
    //.range([ "#402D54", "#D18975", "#8FD175"])
    
  
  //var data = []
  //for (const [type1key, type1value] of Object.entries(boxplots[MEASURE])) {
    //for (const [type2key, type2value] of Object.entries(type1value)) {
      //d = {
          
        //}
      //for (const [methkey, methvalue] of Object.entries(type2value)) {
        
      //}
      //data.push (d)
        
        
         //xdomain.add (type1key+"-"+type2key+"-"+methkey);
         //scoreId = MEASURE + "_" + methkey;
         //vals = []
         //inf_p = []
         //inf_m = []
         //for (i = 0; i< methvalue["samples"].length; i++)
         //{
           //v= methvalue["samples"][i].scores[scoreId]
           //if (v == 1000) {
             //infinites_p = true
             //inf_p.push (methvalue["samples"][i])
           //} else if (v == -1000) {
             //infinites_m = true
             //inf_m.push (methvalue["samples"][i])
           //} else {
             //vals.push (v)
           //}
         //}
         //vals = vals.sort(d3.ascending)
         //q1 = d3.quantile(vals,.25);
         //median = d3.quantile(vals,.5);
         //q3 = d3.quantile(vals,.75);
         //interQuantileRange = q3 - q1;
         //min = vals[0]
         //max = vals[vals.length - 1]
        //whiskersMin = Math.max(min, q1 - interQuantileRange * 1.5);
        //whiskersMax = Math.min(max, q3 + interQuantileRange * 1.5);
        //outliers = methvalue["samples"].filter (x => (x.scores[scoreId] < whiskersMin || x.scores[scoreId] > whiskersMax) && vals.includes (x.scores[scoreId]));
        //sumstat.push ({
          //"key": type1key+"-"+type2key+"-"+methkey,
          //"value": {q1: q1, median: median, q3: q3, interQuantileRange: interQuantileRange, min: min, max: max, whiskersMin: whiskersMin, whiskersMax: whiskersMax, outliers: outliers, scoreId: scoreId, inf_p: inf_p, inf_m: inf_m}});
        //if (minY > min && min != -1000 && min != 1000)
          //minY = min;
        //if (maxY < max && max != 1000 && max != -1000)
          //maxY = max;
      //}
    //}
  //}
  
  source = sources[0]
  type = types[0]
  cur_samples = undefined
  if (MEASURE == "BlandAltman" || MEASURE == "Jaccard")
    cur_samples = samples[source + "__" + type + "__double"]
  else
    cur_samples = samples[source + "__" + type + "__single"]
//console.log (samples[source + "__" + type])
  
  var color = d3.scaleLinear()
    .domain([-1, 0, 1])
    .range(["#fc8d59", "#ffffbf", "#91bfdb"]);
  
  console.log ("n samples: ", Object.values (cur_samples).length)
  
  // ------------------------------- //
  // Add charts
  // ------------------------------- //
  
  for (var_i in allVar){
    for (var_j in allVar){

      // Get current variable name
      var var1 = allVar[var_i]
      var var2 = allVar[var_j]

      scoreId_x = MEASURE + "_" + var1
      scoreId_y = MEASURE + "_" + var2

      // If var1 == var2 i'm on the diagonal, I skip that
      if (var1 === var2) { continue; }

console.log (source,type,var_i,var_j,var1,var2,scoreId_x,scoreId_y)

    if (var_i < var_j) {
      
      // Add X Scale of each graph
      xextent = d3.extent(Object.values (cur_samples), function(d) { return +d.scores[scoreId_x] })
console.log (xextent)
      var x = d3.scaleLinear()
        .domain(xextent).nice()
        .range([ 0, size-2*mar ]);

      // Add Y Scale of each graph
      yextent = d3.extent(Object.values (cur_samples), function(d) { return +d.scores[scoreId_y] })
      var y = d3.scaleLinear()
        .domain(yextent).nice()
        .range([ size-2*mar, 0 ]);

      // Add a 'g' at the right position
      var tmp = svg
        .append('g')
        .attr("transform", "translate(" + (position(var1)+mar) + "," + (position(var2)+mar) + ")");

      // Add X and Y axis in tmp
      tmp.append("g")
        .attr("transform", "translate(" + 0 + "," + (size-mar*2) + ")")
        .call(d3.axisBottom(x).ticks(3));
      tmp.append("g")
        .call(d3.axisLeft(y).ticks(3));
        
        
      if (false && Object.values (cur_samples).length < 1500) {

      // Add circle
      tmp
        .selectAll("myCircles")
        .data(Object.values (cur_samples))
        .enter()
        .append("circle")
          .attr("cx", function(d){ return x(+d.scores[scoreId_x]) })
          .attr("cy", function(d){ return y(+d.scores[scoreId_y]) })
          .attr("r", 3)
          .attr("fill", "#000")
        }
        else {
          
          var color = d3.scaleLinear()
      .domain([0, 1]) // Points per square pixel.
      .range(["white", "#69b3a2"])
      
      var densityData = d3.contourDensity()
    .x(function(d) { return x(+d.scores[scoreId_x]); })
    .y(function(d) { return y(+d.scores[scoreId_y]); })
    .size([size, size])
    .bandwidth(10)
    (Object.values (cur_samples))

  // show the shape!
  tmp.insert("g", "g")
        //.attr("transform", "translate(" + 0 + "," + (size-mar*2) + ")")
    .selectAll("path")
    .data(densityData)
    .enter().append("path")
      .attr("d", d3.geoPath())
      //.attr("fill", function(d) { return color(d.value); })
            .attr("fill", "none")
      .attr("stroke", "#69b3a2")
      .attr("stroke-linejoin", "round")
      //.attr("fill", "#69b3a2")
      
        }
    }
    else {
      
      corr = pearson (Object.values (cur_samples), scoreId_x, scoreId_y)
console.log (source,type,var_i,var_j,var1,var2,scoreId_x,scoreId_y,corr)
      
      var tmp = svg
        .append('g')
        .attr("transform", "translate(" + (position(var1)+mar) + "," + (position(var2)) + ")");
        
          tmp.append("rect")
             .attr("x", 0)
              .attr("y", mar)
             .attr("width", size-mar*2)
             .attr("height", size-mar*2)
             .style("fill", color(corr))
    tmp.append("text")
      .attr("y", size/2+5)
      .attr("x", (size-mar*2)/2)
      .text(d3.format(".2f")(corr))
      .style("font-size", 11)
      .style("text-align", "center")
    .style("text-anchor","middle")
      //.style("fill", color(corr));

    }
  }
  }
  
  
  
  
  // ------------------------------- //
  // Add histograms = diagonal
  // ------------------------------- //
  for (i in allVar){
    for (j in allVar){

      // variable names
      var var1 = allVar[i]
      var var2 = allVar[j]

      scoreId_x = MEASURE + "_" + var1
      scoreId_y = MEASURE + "_" + var2

      // If var1 == var2 i'm on the diagonal, otherwisee I skip
      if (i != j) { continue; }

      // create X Scale
      //xextent = d3.extent(data, function(d) { return +d[var1] })
      xextent = d3.extent(Object.values (cur_samples), function(d) { return +d.scores[scoreId_x] })
      var x = d3.scaleLinear()
        .domain(xextent).nice()
        .range([ 0, size-2*mar ]);

      // Add a 'g' at the right position
      var tmp = svg
        .append('g')
        .attr("transform", "translate(" + (position(var1)+mar) + "," + (position(var2)+mar) + ")");

      // Add x axis
      tmp.append("g")
        .attr("transform", "translate(" + 0 + "," + (size-mar*2) + ")")
        .call(d3.axisBottom(x).ticks(3));

      // set the parameters for the histogram
       var histogram = d3.histogram()
           .value(function(d) { return +d.scores[scoreId_x] })   // I need to give the vector of value
           .domain(x.domain())  // then the domain of the graphic
           .thresholds(x.ticks(15)); // then the numbers of bins

       // And apply this function to data to get the bins
       var bins = histogram(Object.values (cur_samples));

       // Y axis: scale and draw:
       var y = d3.scaleLinear()
            .range([ size-2*mar, 0 ])
            .domain([0, d3.max(bins, function(d) { return d.length; })]);   // d3.hist has to be called before the Y axis obviously

       // append the bar rectangles to the svg element
       tmp.append('g')
          .selectAll("rect")
          .data(bins)
          .enter()
          .append("rect")
             .attr("x", 1)
             .attr("transform", function(d) { return "translate(" + x(d.x0) + "," + y(d.length) + ")"; })
             .attr("width", function(d) { return x(d.x1) - x(d.x0)  ; })
             .attr("height", function(d) { return (size-2*mar) - y(d.length); })
             .style("fill", "#b8b8b8")
             .attr("stroke", "white")
    tmp.append("text")
      .attr("y", 0)
      .attr("x", (size-2*mar)/2)
      .text(var1)
      .style("font-size", 11)
      .style("text-align", "center")
    .style("text-anchor","middle")
    }
  }
  
  
    
    
}






















function draw_boxplots () {
drawn = true
  sqrtscale = document.getElementById('sqrt').checked
  infinites_p = false
  infinites_m = false
  minY = 10000000;
  maxY = 0;
  MEASURE = document.getElementById("metric")
  MEASURE = MEASURE.options[MEASURE.selectedIndex].value
  xdomain = new Set ();
      
      d3.select("#my_dataviz svg").remove();
      
// append sthe svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

  var sumstat = []
  for (const [type1key, type1value] of Object.entries(boxplots[MEASURE])) {
    for (const [type2key, type2value] of Object.entries(type1value)) {
      for (const [methkey, methvalue] of Object.entries(type2value)) {
         xdomain.add (type1key+"-"+type2key+"-"+methkey);
         scoreId = MEASURE + "_" + methkey;
         vals = []
         inf_p = []
         inf_m = []
         for (i = 0; i< methvalue["samples"].length; i++)
         {
           v= methvalue["samples"][i].scores[scoreId]
           if (v == 1000) {
             infinites_p = true
             inf_p.push (methvalue["samples"][i])
           } else if (v == -1000) {
             infinites_m = true
             inf_m.push (methvalue["samples"][i])
           } else {
             vals.push (v)
           }
         }
         vals = vals.sort(d3.ascending)
         q1 = d3.quantile(vals,.25);
         median = d3.quantile(vals,.5);
         q3 = d3.quantile(vals,.75);
         interQuantileRange = q3 - q1;
         min = vals[0]
         max = vals[vals.length - 1]
        whiskersMin = Math.max(min, q1 - interQuantileRange * 1.5);
        whiskersMax = Math.min(max, q3 + interQuantileRange * 1.5);
        outliers = methvalue["samples"].filter (x => (x.scores[scoreId] < whiskersMin || x.scores[scoreId] > whiskersMax) && vals.includes (x.scores[scoreId]));
        sumstat.push ({
          "key": type1key+"-"+type2key+"-"+methkey,
          "value": {q1: q1, median: median, q3: q3, interQuantileRange: interQuantileRange, min: min, max: max, whiskersMin: whiskersMin, whiskersMax: whiskersMax, outliers: outliers, scoreId: scoreId, inf_p: inf_p, inf_m: inf_m}});
        if (minY > min && min != -1000 && min != 1000)
          minY = min;
        if (maxY < max && max != 1000 && max != -1000)
          maxY = max;
      }
    }
  }
  console.log (minY, maxY);
  
  
  //var boxtip = d3.tip().attr('class', 'd3-tip').direction('e').offset([0,5])
  var boxtip = d3.tip().attr('class', 'd3-tip').direction(function(d) {if (x(d.key) < width/2) return 'e'; return 'w'}).offset(function(d) {if (x(d.key) < width/2) return [0,5]; return [0,-5]})
            .html(function(d) {
                var content = "<span style='margin-left: 2.5px;'><b>" + d.key + "</b></span><br>";
                content +=`
                    <table style="margin-top: 2.5px;">
                            <tr><td>Max: </td><td style="text-align: right">` + d3.format(".2f")(d.value.max) + `</td></tr>
                            <tr><td>Q3: </td><td style="text-align: right">` + d3.format(".2f")(d.value.q3) + `</td></tr>
                            <tr><td>Median: </td><td style="text-align: right">` + d3.format(".2f")(d.value.median) + `</td></tr>
                            <tr><td>Q1: </td><td style="text-align: right">` + d3.format(".2f")(d.value.q1) + `</td></tr>
                            <tr><td>Min: </td><td style="text-align: right">` + d3.format(".2f")(d.value.min) + `</td></tr>
                    </table>
                    `;
                return content;
            });
        svg.call(boxtip);
  var outliertip = d3.tip().attr('class', 'd3-tip').direction(function(d) {if (x(d.n) < width/2) return 'e'; return 'w'}).offset(function(d) {if (x(d.n) < width/2) return [0,5]; return [0,-5]})
            .html(function(d) {
                var content = "<span style='margin-left: 2.5px;'><b>" + d.n + "</b></span><br>";
                content += "<span style='margin-left: 2.5px;'>Sample: " + d.s + "</span><br>";
                content += "<span style='margin-left: 2.5px;'>Value: " + d3.format(".2f")(d.p) + "</span><br>";
                return content;
            });
        svg.call(outliertip);
  
  
  //var btm = undefined
  
  xdomain = Array.from(xdomain).sort(function(a, b) { return d3.ascending(a, b); })
  var x = d3.scaleBand()
    .range([ 0, width ])
    .domain(xdomain)
    .paddingInner(1)
    .paddingOuter(.5)
        
  chart_top = 0
  chart_bottom = height
  
  r_u = chart_top + 20
  r_b = chart_bottom
  
  if (infinites_p || infinites_m)
    r_u = r_u +30
  if (infinites_p && infinites_m)
    r_u = r_u + 30
  
  sqrtscale = true
  if (sqrtscale) {
    extra = 0
    //extra = Math.pow (.01*(maxY-minY), document.getElementById("myRange").value)
    //console.log (extra, .01*(maxY-min), document.getElementById("myRange").value)
  //var y = d3.scaleSqrt()
  var y = d3.scalePow().exponent(document.getElementById("myRange").value)
    .domain([minY - extra,maxY + extra])
    .range([r_b, r_u])
  } else {
  var y = d3.scaleLinear()
    .domain([minY - .1*(maxY-minY),maxY + .1*(maxY-minY)])
    .range([r_b, r_u])
    
  }
  
    
  
  var boxWidth = width / sumstat.length
  svg
    .selectAll("databackground")
    .data(sumstat)
    .enter()
    .append("rect")
        .attr("x", function(d){return(x(d.key)-boxWidth/2)})
        .attr("y", chart_top)
        .attr("height", chart_bottom + 5)
        .attr("width", boxWidth )
        .style("opacity", ".2")
        .style("fill", function(d){return datacolor(d.key)})
        
        
        
  console.log (sumstat);
  
  for (i = 0; i < sumstat.length; i++) {
    const arr = []
    for (j = 0; j < sumstat[i].value.outliers.length; j++)
      arr.push ({n: sumstat[i].key, s:sumstat[i].value.outliers[j].name, p: sumstat[i].value.outliers[j].scores[sumstat[i].value.scoreId]})
    console.log (arr)
    svg
      .selectAll("outliers")
      .data(arr)
      .enter()
      .append("circle")
        .attr("r", 1)
        .attr("cx", function(d){return(x(d.n) + 2 * (Math.random () - .5))})
        .attr("cy", function(d){return(y(d.p))})
        .on('mouseover', outliertip.show)
        .on('mouseout', outliertip.hide);
        //.attr("stroke", "black")
        //.style("width", 40)
  }
  
  inf_y = r_u
  if (infinites_p) {
    svg.append("text")
    .attr("x", 0)
    .attr("y", inf_y - 10)
    .attr("text-anchor","end")
    .attr("font-size","10px")
    .text("Infinity");
    svg
      .append("line")
        .attr("x1", 0)
        .attr("x2", width)
        .attr("y1", inf_y - 5)
        .attr("y2", inf_y - 5)
        .attr("stroke", "black")
        .style("width", 10)
          .style("opacity", ".2")
    
    for (i = 0; i < sumstat.length; i++) {
    const arr = []
    for (j = 0; j < sumstat[i].value.inf_p.length; j++)
      arr.push ({n: sumstat[i].key, s:sumstat[i].value.inf_p[j].name, p: Infinity})
    console.log (arr)
    svg
      .selectAll("inf_p")
      .data(arr)
      .enter()
      .append("circle")
        .attr("r", 1)
        .attr("cx", function(d){return(x(d.n) + 2 * (Math.random () - .5))})
        .attr("cy", inf_y - 13)
        .on('mouseover', outliertip.show)
        .on('mouseout', outliertip.hide);
        //.attr("stroke", "black")
        //.style("width", 40)
  }
    inf_y = inf_y - 30
  }
  
  if (infinites_m) {
    svg.append("text")
    .attr("x", 0)
    .attr("y", inf_y - 10)
    .attr("text-anchor","end")
    .attr("font-size","10px")
    .text("-Infinity");
    svg
      .append("line")
        .attr("x1", 0)
        .attr("x2", width)
        .attr("y1", inf_y - 5)
        .attr("y2", inf_y - 5)
        .attr("stroke", "black")
        .style("width", 10)
          .style("opacity", ".2")
    
    for (i = 0; i < sumstat.length; i++) {
    const arr = []
    for (j = 0; j < sumstat[i].value.inf_m.length; j++)
      arr.push ({n: sumstat[i].key, s:sumstat[i].value.inf_m[j].name, p: -Infinity})
    console.log (arr)
    svg
      .selectAll("inf_m")
      .data(arr)
      .enter()
      .append("circle")
        .attr("r", 1)
        .attr("cx", function(d){return(x(d.n) + 2 * (Math.random () - .5))})
        .attr("cy", inf_y - 13)
        .on('mouseover', outliertip.show)
        .on('mouseout', outliertip.hide);
        //.attr("stroke", "black")
        //.style("width", 40)
  }
  }
  
  
  svg
    .selectAll("vertLines")
    .data(sumstat)
    .enter()
    .append("line")
      .attr("x1", function(d){return(x(d.key))})
      .attr("x2", function(d){return(x(d.key))})
      .attr("y1", function(d){return(y(d.value.whiskersMin))})
      .attr("y2", function(d){return(y(d.value.whiskersMax))})
      .attr("stroke", "black")
      .style("width", 40)
  
  
  

  
  
  
  
  
  // rectangle for the main box
  var boxWidth = width / sumstat.length - 20
  svg
    .selectAll("boxes")
    .data(sumstat)
    .enter()
    .append("rect")
        .attr("x", function(d){return(x(d.key)-boxWidth/2)})
        .attr("y", function(d){return(y(d.value.q3))})
        .attr("height", function(d){return(y(d.value.q1)-y(d.value.q3))})
        .attr("width", boxWidth )
        //.attr("stroke", function(d){return boxColor(d.key)})
      .attr("stroke", "black")
        .style("fill", function(d){return boxColor(d.key)})
        .on('mouseover', boxtip.show)
        .on('mouseout', boxtip.hide);
  
  svg
    .selectAll("medianLines")
    .data(sumstat)
    .enter()
    .append("line")
      .attr("x1", function(d){return(x(d.key)-boxWidth/2) })
      .attr("x2", function(d){return(x(d.key)+boxWidth/2) })
      .attr("y1", function(d){return(y(d.value.median))})
      .attr("y2", function(d){return(y(d.value.median))})
      .attr("stroke", "black")
      .style("width", 80)
        .on('mouseover', boxtip.show)
        .on('mouseout', boxtip.hide);
  
  
  types1 = []
  types2 = []
  ctypes1 = {"name": undefined}
  ctypes2 = {"name": undefined}
  for (i =0; i < xdomain.length; i++) {
    t1 = xdomain[i].replace (/([a-zA-Z ]+)-.*-.*/g, '$1')
    t2 = xdomain[i].replace (/.*-([a-zA-Z ]+)-.*/g, '$1')
    
    if (t1 == ctypes1["name"])
      ctypes1["end"] = xdomain[i]
    else {
      if (ctypes1["name"])
        types1.push (ctypes1)
      ctypes1 = {
        "start": xdomain[i],
        "end": xdomain[i],
        "name": t1
      }
    }
    
    if (t2 == ctypes2["name"])
      ctypes2["end"] = xdomain[i]
    else {
      if (ctypes2["name"])
        types2.push (ctypes2)
      ctypes2 = {
        "start": xdomain[i],
        "end": xdomain[i],
        "name": t2
      }
    }
    
  }
      if (ctypes2["name"])
        types2.push (ctypes2)
      if (ctypes1["name"])
        types1.push (ctypes1)
  
  
  for (i =0; i < types1.length; i++) {
    console.log (types1[i]);
    txt = svg.append("text")
    .attr("x", (x(types1[i]["start"]) + x(types1[i]["end"]))/2)
    .attr("text-anchor","middle")
    .attr("y", -margin.top + 30)
    .text(types1[i]["name"]);
  }
  for (i =0; i < types2.length; i++) {
    console.log (types2[i]);
    txt = svg.append("text")
    .attr("x", (x(types2[i]["start"]) + x(types2[i]["end"]))/2)
    .attr("text-anchor","middle")
    .attr("y", 20)
    .text(types2[i]["name"]);
  }
  
  
  
  
  
  
  
  
  
  svg.append("g")
    .attr("transform", "translate(0," + (height + 5) + ")")
    .call(d3.axisBottom(x))
      .selectAll("text")
        .attr("transform", "translate(-12,10)rotate(-90)")
        .style("text-anchor", "end")
        .style("font-size", 28)
    .each(function(d, i){
      d3.select(this).text(d.replace (/.*-/g,''));
    })
  svg.append("g").call(d3.axisLeft(y))
  
  
  
  
  
  
  
};


d3.select("#download").on("click", function(){
  d3.select(this)
    .attr("href", 'data:application/octet-stream;base64,' + btoa(d3.select("#my_dataviz").html()))
    .attr("download", "plot.svg") 
})


d3.csv("data/combined-afr-eor-hallmark.csv").then (function(data) {
  metrics = new Set ()
  for (row=0; row < data.length; row++){
    metrics.add (data[row]["Score"])
    var [source, type] = dataset2sourcetype (data[row]["Dataset"])
    
		if (!samples[source + "__" + type + "__single"][data[row]["Sample"]]) {
			samples[source + "__" + type + "__single"][data[row]["Sample"]] = new Sample (data[row]["Sample"], source, type)
		}
		
    const sample = samples[source + "__" + type + "__single"][data[row]["Sample"]]
		if (!boxplots[data[row]["Score"]]) {
			boxplots[data[row]["Score"]] = getSourceData ();
    }
    
    for (i = 0; i < imethods.length; i++) {
      boxplots[data[row]["Score"]][source][type][imethods[i]]["samples"].push (sample);
      sample.setScore (data[row]["Score"], imethods[i], parseFloat (data[row][imethods[i]]))
    }
  }
  
  m_select = document.getElementById("metric")
  for (let metric of metrics) {
    var opt = document.createElement('option');
    opt.value = metric;
    opt.innerHTML = metric;
    m_select.appendChild(opt);
  }
  
  sortSelect (m_select)
  
  if (!drawn)
    draw_boxplots ();
});


d3.csv("data/combined-jaccard-ba.csv").then (function(data) {
  metrics = new Set ()
  for (row=0; row < data.length; row++){
    metrics.add (data[row]["Score"])
    var [source, type] = dataset2sourcetype (data[row]["Dataset"])
		
    sid = data[row]["Dx1"] + " -vs- " + data[row]["Dx2"]
    
		if (!samples[source + "__" + type + "__double"][sid]) {
			samples[source + "__" + type + "__double"][sid] = new Sample (sid, source, type)
		}
    
    const sample = samples[source + "__" + type + "__double"][sid]
		if (!boxplots[data[row]["Score"]]) {
			boxplots[data[row]["Score"]] = getSourceData ();
    }
    
    for (i = 0; i < imethods.length; i++) {
      boxplots[data[row]["Score"]][source][type][imethods[i]]["samples"].push (sample);
      sample.setScore (data[row]["Score"], imethods[i], parseFloat (data[row][imethods[i]]))
    }
  }
  
  m_select = document.getElementById("metric")
  for (let metric of metrics) {
    var opt = document.createElement('option');
    opt.value = metric;
    opt.innerHTML = metric;
    m_select.appendChild(opt);
  }
  
  sortSelect (m_select)
  
  if (!drawn)
    draw_boxplots ();
});

