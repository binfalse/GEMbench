
const sources = ["Cell Line","Patient Data"];
const types = ["Microarray","RNA Seq","MS Proteomics"];
const imethods = ["FASTCORE","GIMME","INIT","iMAT"];
const metrics = ["AFR","EOR","Hallmark","BlandAltman","Jaccard"]
const data_sources = ["EMTAB-37","GSE2109","HPA","ProteomeNCI60","ProteomePatients","TCGA"]
metric_dict = {
  "AFR": 0,
  "EOR": 1,
  "Hallmark": 2,
  "BlandAltman": 3,
  "Jaccard": 4
}
data_dict = {
  "EMTAB-37": 0,
  "GSE2109": 1,
  "HPA": 2,
  "ProteomeNCI60": 3,
  "ProteomePatients": 4,
  "TCGA": 5
}
imeth_dict = {
  "FASTCORE": 0,
  "GIMME": 1,
  "INIT": 2,
  "iMAT": 3,
}


Object.values = Object.values || function(o){return Object.keys(o).map(function(k){return o[k]})};

function boxColor (imeth) {
  //if (str.includes ("GIMME")) {
  if (imeth == 1) {
    return "#377eb8"
  //} else if (str.includes ("FASTCORE")) {
  } else if (imeth == 0) {
    return "#4daf4a"
  //} else if (str.includes ("INIT")) {
  } else if (imeth == 2) {
    return "#984ea3"
  }
  return "#ff7f00"
}
function datacolor (data) {
  if (data == data_dict["EMTAB-37"] || data == data_dict["GSE2109"]) {
  //if (str.includes ("Microarray")) {
    return "#fef0d9"
  } else if (data == data_dict["ProteomeNCI60"] || data == data_dict["ProteomePatients"]) {
  //} else if (str.includes ("MS Proteomics")) {
    return "#fdddb2"
  }
  return "#f3e793"
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

function getDataFromSource (source, type) {
  if (source == sources[0]) {
    if (type == types[0])
      return "EMTAB-37"
    else if (type == types[1])
      return "HPA"
    else if (type == types[2])
      return "ProteomeNCI60"
  }
  else if (source == sources[1]) {
    if (type == types[0])
      return "GSE2109"
    else if (type == types[1])
      return "TCGA"
    else if (type == types[2])
      return "ProteomePatients"
  }
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
function addSelect (metric, metricid) {
  
  var opt = document.createElement('option');
  opt.value = metricid;
  opt.innerHTML = metric;
  
  var m_select = document.getElementById("metric");
  m_select.appendChild(opt);
  
  sortSelect (m_select)
}
function get_x_for_boxplot (boxplot) {
    var tmp = ""
    switch (boxplot[11]) {
        case data_dict["EMTAB-37"]:
          tmp = sources[0] + "__" + types[0]
    break;
        case data_dict["HPA"]:
          tmp = sources[0] + "__" + types[1]
    break;
        case data_dict["ProteomeNCI60"]:
          tmp = sources[0] + "__" + types[2]
    break;
        case data_dict["GSE2109"]:
          tmp = sources[1] + "__" + types[0]
    break;
        case data_dict["TCGA"]:
          tmp = sources[1] + "__" + types[1]
    break;
        case data_dict["ProteomePatients"]:
          tmp = sources[1] + "__" + types[2]
    break;
    }
    return metrics[boxplot[10]] + "_" + tmp + "_" + imethods[boxplot[12]]
  }
function get_sample_name (item) {
  if (item.length == 3)
    return samples[item[1]] + " -vs- " + samples[item[2]]
  return samples[item[1]]
}
//class Sample {
	//constructor (name, source, type) {
		//this.name = name
		//this.scores = {}
    //this.source = source
    //this.type = type
	//}
	//setScore (measure, method, value) {
		//this.scores[measure + "_" + method] = value
	//}
//}
//function SampleSorter (scoreId) {
 //return (a, b) => {return a.scores[scoreId]-b.scores[scoreId]}	
//}

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
    x.push (samples[i].values[metric1])
    y.push (samples[i].values[metric2])
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


function update_slider_value () {
  slider = document.getElementById("slider")
  //console.log (slider)
  document.getElementById("slider_value").innerHTML = slider.value
}

function get_outlier_table_id (metric, data_id) {
  return "outlier__" + metric + "__" + data_id.replace(/[\W_]+/g,"_")
}

var margin = {top: 50, right: 0, bottom: 70, left: 40},
    width = d3.select("#my_dataviz").node().clientWidth - margin.left - margin.right,
    height = 800 - margin.top - margin.bottom;


//boxplots = {}

//var xdomain = new Set ();
drawn = false
samples = []
//samples = {}
//for (var i = 0; i < sources.length; i++)
//for (var j = 0; j < types.length; j++) {
  //samples[sources[i] + "__" + types[j] + "__single"] = {}
  //samples[sources[i] + "__" + types[j] + "__double"] = {}
//}
//MEASURE = "AFR"
//MEASURE = "EOR"
//MEASURE = "Hallmark"

//var slider = document.getElementById("myRange");
//console.log (slider.value)

// Update the current slider value (each time you drag the slider handle)
/*slider.oninput = function() {
  console.log (this.value);
} */

$(".metric_expl").hide ()
















function draw_correlation (measure, sumstat, source, type, domnode) {
  console.log (sumstat)
  console.log (measure, source, type)
  console.log (getDataFromSource (source, type))
  var data_id = data_dict[getDataFromSource (source, type)]
  var measure_id = metric_dict[measure]
  
  cor_samples = {}
  
  //if (measure == "BlandAltman" || measure == "Jaccard")
    //for (var i = 0; i < samples.length; i++)
    //for (var j = 0; j < samples.length; j++)
      //cor_samples[samples[i]] = {"id": i, "values": Array(imethods.length).fill(null)}
  //else
    //for (var i = 0; i < samples.length; i++)
      //cor_samples[samples[i]] = {"id": i, "values": Array(imethods.length).fill(null)}
  //console.log (cor_samples)
  
  for (var i = 0; i < imethods.length; i++) {
    var boxplot = sumstat[measure_id + "_" + data_id + "_" + i]
    
    var values = boxplot[15] // values in whisker range
        .concat (boxplot[8]) // outliers min
        .concat (boxplot[9]) // outliers max
    
    // values in whisker range
    for (var j = 0; j < values.length; j++) {
      const sid = get_sample_name (values[j])
      if (!(sid in cor_samples))
        cor_samples[sid] = {"id": sid, "values": Array(imethods.length).fill(null)}
      cor_samples[sid]["values"][i] = values[j][0]
    }
  }
  //console.log (cor_samples)
  
  
  //for (const [sumstat_key, boxplot] of Object.entries(sumstat)) {
    //console.log (sumstat_key,"---",measure_id,data_id)
  //} 
  
  $("#" + domnode).empty ();
  $("#" + domnode).append ("<p><a title='Correlation of "+measure+" scores for "+ source + " &mdash; " + type + "' id='"+domnode+"_a' href='#"+domnode+"_svg'>" + source + "<br>" + type + "</a></p>");
  
var marginWhole = {top: 10, right: 10, bottom: 10, left: 10},
    sizeWhole = 600 - marginWhole.left - marginWhole.right

// Create the svg area
var svg = d3.select("#" + domnode + " p a")
  .append("svg")
  //.attr ("preserveAspectRatio", "xMinYMin meet")
    .attr("width", 150)
    .attr("height", 150)
    .attr("id", domnode+"_svg")
    //.attr("width", sizeWhole  + marginWhole.left + marginWhole.right)
    //.attr("height", sizeWhole  + marginWhole.top + marginWhole.bottom)
  .attr ("viewBox", "0 0 " + (sizeWhole  + marginWhole.left + marginWhole.right) + " " + (sizeWhole  + marginWhole.left + marginWhole.right))
  //.classed("svg-content", true)
  .append("g")
    .attr("transform", "translate(" + marginWhole.left + "," + marginWhole.top + ")");
    
  $("#" + domnode+"_a").fancybox ({
    beforeLoad: function ()  {
      var e = $("#" + domnode+"_svg")
      e.attr("height", 600);
      e.attr("width", 600);
    },
    afterClose: function ()  {
      var e = $("#" + domnode+"_svg")
      e.attr("height", 150);
      e.attr("width", 150);
      e.show ()
      $("#" + domnode + " p a").append (e)
      console.log ("test", domnode)
    },
    title: "<h5>Correlation of "+measure+" scores for "+ source + " &mdash; " + type + "</h5>",
		helpers		: {
			title	: { type : 'inside',
                position : 'top' },
    }
	});
    
    var allVar = imethods
  var numVar = allVar.length
    var mar = 20
  var size = sizeWhole / numVar
  
  var position = d3.scalePoint()
    .domain(allVar)
    .range([0, sizeWhole-size])
  
  //var cur_samples = {}
  
  //if (MEASURE == "BlandAltman" || MEASURE == "Jaccard")
    //cur_samples = samples[source + "__" + type + "__double"]
  //else
    //cur_samples = samples[source + "__" + type + "__single"]
//console.log (samples[source + "__" + type])
  
  var color = d3.scaleLinear()
    .domain([-1, 0, 1])
    .range(["#2b9d67", "#fffbf0", "#2b9d67"])
    .unknown ("#fff");
  
  console.log ("n samples: ", Object.values (cor_samples).length)
  
  // ------------------------------- //
  // Add charts
  // ------------------------------- //
  
  for (var_i in allVar){
    for (var_j in allVar){

      // Get current variable name
      var var1 = allVar[var_i]
      var var2 = allVar[var_j]

      scoreId_x = measure + "_" + var1
      scoreId_y = measure + "_" + var2

      // If var1 == var2 i'm on the diagonal, I skip that
      if (var1 === var2) { continue; }

//console.log (source,type,var_i,var_j,var1,var2,scoreId_x,scoreId_y)

    if (var_i < var_j) {
      
      // Add X Scale of each graph
      xextent = d3.extent(Object.values (cor_samples), function(d) { return +d.values[var_i] })
//console.log (xextent)
      var x = d3.scalePow().exponent(document.getElementById("slider").value)
      //scaleLinear()
        .domain(xextent).nice()
        .range([ 0, size-2*mar ]);

      // Add Y Scale of each graph
      yextent = d3.extent(Object.values (cor_samples), function(d) { return +d.values[var_j] })
      var y = d3.scalePow().exponent(document.getElementById("slider").value)
      //scaleLinear()
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
        
        
      
          
          //var color = d3.scaleLinear()
      //.domain([0, 1]) // Points per square pixel.
      //.range(["white", "#69b3a2"])
      
      var densityData = d3.contourDensity()
    .x(function(d) { return x(+d.values[var_i]); })
    .y(function(d) { return y(+d.values[var_j]); })
    .size([size, size])
    .bandwidth(10)
    (Object.values (cor_samples))

  // show the shape!
  tmp.insert("g", "g")
        //.attr("transform", "translate(" + 0 + "," + (size-mar*2) + ")")
    .selectAll("path")
    .data(densityData)
    .enter().append("path")
      .attr("d", d3.geoPath())
      //.attr("fill", function(d) { return color(d.value); })
            .attr("fill", "none")
      .attr("stroke", "#2b9d67")
      .attr("stroke-linejoin", "round")
      //.attr("fill", "#69b3a2")
      
        
    }
    else {
      
      corr = pearson (Object.values (cor_samples), var_i, var_j)
//console.log (source,type,var_i,var_j,var1,var2,scoreId_x,scoreId_y,corr)
      
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
  for (var i in allVar){
    for (var j in allVar){

      // variable names
      var var1 = allVar[i]
      var var2 = allVar[j]

      scoreId_x = measure + "_" + var1
      scoreId_y = measure + "_" + var2

      // If var1 == var2 i'm on the diagonal, otherwisee I skip
      if (i != j) { continue; }

      // create X Scale
      //xextent = d3.extent(data, function(d) { return +d[var1] })
      xextent = d3.extent(Object.values (cor_samples), function(d) { return +d.values[var_i] })
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
           .value(function(d) { return +d.values[var_i] })   // I need to give the vector of value
           .domain(x.domain())  // then the domain of the graphic
           .thresholds(x.ticks(15)); // then the numbers of bins

       // And apply this function to data to get the bins
       var bins = histogram(Object.values (cor_samples));

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
};














global_sumstat = {}






function draw_boxplots () {
  update_slider_value ()
  drawn = true
  //sqrtscale = document.getElementById('sqrt').checked
  var MEASURE = document.getElementById("metric")
  MEASURE = metrics[MEASURE.options[MEASURE.selectedIndex].value]
$(".metric_expl").hide ()
  $("#" + MEASURE + "_expl").show ();
      
      d3.select("#my_dataviz svg").remove();
      
// append sthe svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

  var sumstat = global_sumstat[MEASURE]["boxplots"];
  var n_boxplots = Object.keys(sumstat).length;
  var minY  = global_sumstat[MEASURE]["minY"];
  var maxY  = global_sumstat[MEASURE]["maxY"];
  var infinites_p = global_sumstat[MEASURE]["infinites_p"];
  var infinites_m = global_sumstat[MEASURE]["infinites_m"];
  var xdomain = global_sumstat[MEASURE]["xdomain"];
  //for (var i = 0; i < xdomain.length; i++) {
    //var [t1, t2] = dataset2sourcetype (xdomain[i].replace (/.*_([a-zA-Z0-9 -]+)_.*/g, '$1'))
    //xdomain[i] = t1 + "__" + t2 + "__" + xdomain[i]
  //}
  console.log (minY, maxY);
  var columnWidth = width / n_boxplots;
  // rectangle for the main box
  var boxWidth = width / n_boxplots - 20;
  console.log (xdomain)
  console.log (sumstat)
  
  //var boxtip = d3.tip().attr('class', 'd3-tip').direction('e').offset([0,5])
  var boxtip = d3.tip().attr('class', 'd3-tip').direction(function(d) {if (x(get_x_for_boxplot(d)) < width/2) return 'e'; return 'w'}).offset(function(d) {if (x(get_x_for_boxplot(d)) < width/2) return [0,5]; return [0,-5]})
            .html(function(d) {
                var content = "<span style='margin-left: 2.5px;'><b>" + get_x_for_boxplot(d) + "</b></span><br>";
                content +=`
                    <table style="margin-top: 2.5px;">
                            <tr><td>Max: </td><td style="text-align: right">` + d3.format(".2f")(d[5]) + `</td></tr>
                            <tr><td>Q3: </td><td style="text-align: right">` + d3.format(".2f")(d[2]) + `</td></tr>
                            <tr><td>Median: </td><td style="text-align: right">` + d3.format(".2f")(d[1]) + `</td></tr>
                            <tr><td>Q1: </td><td style="text-align: right">` + d3.format(".2f")(d[0]) + `</td></tr>
                            <tr><td>Min: </td><td style="text-align: right">` + d3.format(".2f")(d[4]) + `</td></tr>
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
  
  xdomain = xdomain.sort(function(a, b) { return d3.ascending(a, b); })
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
  
  var y = d3.scalePow().exponent(document.getElementById("slider").value)
    .domain([minY,maxY])
    .range([r_b, r_u])
  
  
  
  
  svg
    .selectAll("databackground")
    .data(Object.values (sumstat))
    .enter()
    .append("rect")
        .attr("x", function(d){console.log (d); return(x(get_x_for_boxplot(d))-columnWidth/2)})
        .attr("y", chart_top)
        .attr("height", chart_bottom + 5)
        .attr("width", columnWidth )
        .style("opacity", ".3")
        .style("fill", function(d){return datacolor(d[11])})
  svg
    .append("rect")
        .attr("x", 0)
        .attr("y", chart_top - 50)
        .attr("height", 50)
        .attr("width", width )
        //.style("opacity", ".1")
        .style("fill", "#fafafa")
  svg
    .append("line")
        .attr("x1", 0)
        .attr("x2", width)
        .attr("y1", chart_top - 50)
        .attr("y2", chart_top - 50)
        .attr("stroke", "black")
        
        
        
  //console.log (sumstat);
  
    //for (j = 0; j < boxplot[13].length; j++)
      //arr.push ({n: get_x_for_boxplot(boxplot), s:boxplot[13][j], p: Infinity})
  //for (i = 0; i < sumstat.length; i++) {
  for (const [sumstat_key, boxplot] of Object.entries(sumstat)) {
    const arr_min = []
    const arr_max = []
    for (j = 0; j < boxplot[8].length; j++)
      arr_min.push ({n: get_x_for_boxplot(boxplot), s:get_sample_name (boxplot[8][j]), p: boxplot[8][j][0]})
    for (j = 0; j < boxplot[9].length; j++)
      arr_max.push ({n: get_x_for_boxplot(boxplot), s:get_sample_name (boxplot[9][j]), p: boxplot[9][j][0]})
    
    arr_min.sort(function (a, b) {return a.p < b.p ? -1 : 1});
    arr_max.sort(function (a, b) {return a.p < b.p ? -1 : 1});
    
    const outlier_table = boxplot[16]
    
    function draw_outliers (outliers) {
      
      if (outliers.length < 70) {
        svg
          .selectAll("outliers")
          .data(outliers)
          .enter()
          .append("circle")
            .attr("r", 1)
            .attr("cx", function(d){return(x(d.n) + 3 * (Math.random () - .5))})
            .attr("cy", function(d){return(y(d.p))})
            .on('mouseover', outliertip.show)
            .on('mouseout', outliertip.hide)
            .on("click", function(){
				console.log ("not too many outliers to draw")
				//console.log (outlier_table)
              $.fancybox( $(outlier_table) );
            });
      } else {
        console.log ("too many outliers to draw", outliers.length)
        
        xval = x(outliers[0].n)
        min = y(outliers[outliers.length - 1].p)
        max = y(outliers[0].p)
        
        points = [[xval,min]]
        
        left = true
        //cur = min
        
        //console.log (min, min + 10, max)
        for (var cur = min + 10; cur < max; cur += Math.floor(3+7*Math.random())) {
          //console.log (cur)
          points.push ([left ? xval - boxWidth/10 : xval + boxWidth/10, cur])
          left = !left;
        }
        points.push ([xval,max])
        
        //console.log (points)
        
        var lineGenerator = d3.line()
                              .curve(d3.curveCardinal);
        var pathData = lineGenerator(points);

        svg.append('path')
          .attr('d', pathData)
                .attr("stroke", "#666")
                //.style("opacity", ".2")
                .style("fill", "none")    
          .on("click", function(){
				console.log ("too many outliers to draw")
				//console.log (outlier_table)
              $.fancybox( $(outlier_table) );
            //~ $.fancybox( "#" + outlier_table );
          });
  
      }
    }
    draw_outliers (arr_min)
    draw_outliers (arr_max)
  }
  
  inf_y = r_u
  if (infinites_p) {
    svg.append("text")
    .attr("x", -5)
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
    
    for (const [sumstat_key, boxplot] of Object.entries(sumstat)) {
    const arr = []
    for (j = 0; j < boxplot[13].length; j++)
      arr.push ({n: get_x_for_boxplot(boxplot), s:get_sample_name (boxplot[13][j]), p: Infinity})
    //console.log (arr)
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
        ////.attr("stroke", "black")
        ////.style("width", 40)
  }
    inf_y = inf_y - 30
  }
  
  if (infinites_m) {
    svg.append("text")
    .attr("x", -5)
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
    
    //for (i = 0; i < sumstat.length; i++) {
    for (const [sumstat_key, boxplot] of Object.entries(sumstat)) {
    const arr = []
    for (j = 0; j < boxplot[14].length; j++)
      //arr.push ({n: sumstat[i].key, s:sumstat[i].value.inf_m[j].name, p: -Infinity})
      arr.push ({n: get_x_for_boxplot(boxplot), s:get_sample_name (boxplot[14][j]), p: -Infinity})
    //console.log (arr)
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
        ////.attr("stroke", "black")
        ////.style("width", 40)
  }
  }
  
  
  svg
    .selectAll("vertLines")
    .data(Object.values (sumstat))
    .enter()
    .append("line")
      .attr("x1", function(d){return(x(get_x_for_boxplot(d)))})
      .attr("x2", function(d){return(x(get_x_for_boxplot(d)))})
      .attr("y1", function(d){return(y(d[6]))})
      .attr("y2", function(d){return(y(d[7]))})
      .attr("stroke", "black")
      .style("width", 40)
  
  
  

  
  
  
  
  
  svg
    .selectAll("boxes")
    .data(Object.values (sumstat))
    .enter()
    .append("rect")
        .attr("x", function(d){return(x(get_x_for_boxplot(d))-boxWidth/2)})
        .attr("y", function(d){return(y(d[2]))})
        .attr("height", function(d){return(y(d[0])-y(d[2]))})
        .attr("width", boxWidth )
        //.attr("stroke", function(d){return boxColor(d.key)})
      .attr("stroke", "black")
        .style("fill", function(d){return boxColor(d[12])})
        .on('mouseover', boxtip.show)
        .on('mouseout', boxtip.hide)
        // TODO:
            .on("click", function(d){
              $.fancybox( $(d[16]));
            });
  
  svg
    .selectAll("medianLines")
    .data(Object.values (sumstat))
    .enter()
    .append("line")
      .attr("x1", function(d){return(x(get_x_for_boxplot(d))-boxWidth/2) })
      .attr("x2", function(d){return(x(get_x_for_boxplot(d))+boxWidth/2) })
      .attr("y1", function(d){return(y(d[1]))})
      .attr("y2", function(d){return(y(d[1]))})
      .attr("stroke", "black")
      .style("width", 80)
        // TODO:
        .on('mouseover', boxtip.show)
        .on('mouseout', boxtip.hide)
            .on("click", function(d){
              $.fancybox( $(d[16]));
            });
  
  
  types1 = []
  types2 = []
  ctypes1 = {"name": undefined}
  ctypes2 = {"name": undefined}
  corId = 1
  for (var x_i =0; x_i < xdomain.length; x_i++) {
    const t1 = xdomain[x_i].replace (/.*_([a-zA-Z -]+)__.*/g, '$1')
    const t2 = xdomain[x_i].replace (/.*__([a-zA-Z -]+)_.*/g, '$1')
    
    //var [t2, t1] = dataset2sourcetype (xdomain[x_i].replace (/.*_([a-zA-Z0-9 -]+)_.*/g, '$1'))
    //console.log (xdomain[x_i], t1, t2)
    
    if (t1 == ctypes1["name"])
      ctypes1["end"] = xdomain[x_i]
    else {
      if (ctypes1["name"])
        types1.push (ctypes1)
      ctypes1 = {
        "start": xdomain[x_i],
        "end": xdomain[x_i],
        "name": t1
      }
    }
    
    if (t2 == ctypes2["name"])
      ctypes2["end"] = xdomain[x_i]
    else {
      if (ctypes2["name"])
        types2.push (ctypes2)
      ctypes2 = {
        "start": xdomain[x_i],
        "end": xdomain[x_i],
        "name": t2
      }
      
      //console.log ("drawing:", t1, t2, "cor" + corId)
      
      
      const tmpCorId = corId
      //if (x_i == 0)
      setTimeout (function () {
        draw_correlation (MEASURE, sumstat, t1, t2, "cor" + tmpCorId)
        }, 1000* Math.random ());
      corId = corId + 1
      //console.log ("drawn")
    }
    
  }
      if (ctypes2["name"])
        types2.push (ctypes2)
      if (ctypes1["name"])
        types1.push (ctypes1)
  
  
  for (i =0; i < types1.length; i++) {
    //console.log (types1[i]);
    //svg.append("rect")
        //.attr("x", x(types1[i]["start"]))
        //.attr("y", -margin.top + 30)
        //.attr("height", 30)
        //.attr("width",  (x(types1[i]["start"]) + x(types1[i]["end"])))
        ////.attr("stroke", function(d){return boxColor(d.key)})
      //.attr("stroke", "black")
        ////.style("fill", function(d){return boxColor(d.key)})
        //.style("fill", "red")
        
    txt = svg.append("text")
    .attr("x", (x(types1[i]["start"]) + x(types1[i]["end"]))/2)
    .attr("text-anchor","middle")
    .attr("y", -margin.top + 30)
    .text(types1[i]["name"])
  svg
    .append("line")
      .attr("x1", x(types1[i]["end"])+columnWidth/2)
      .attr("x2", x(types1[i]["end"])+columnWidth/2)
      .attr("y1", chart_top - 50)
      .attr("y2", chart_bottom + 5)
      .attr("stroke", "black")
    svg.append("line")
      .attr("x1", x(types1[i]["start"])-columnWidth/2)
      .attr("x2", x(types1[i]["start"])-columnWidth/2)
      .attr("y1", chart_top - 50)
      .attr("y2", chart_bottom + 5)
      .attr("stroke", "black")
      .style("width", 10)
  }
  for (i =0; i < types2.length; i++) {
    //console.log (types2[i]);
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
        //.style("font-size", 28)
    .each(function(d, i){
      d3.select(this).text(d.replace (/.*_/g,''));
    })
  svg.append("g").call(d3.axisLeft(y))
  
  
  
  
  
  //for (var i = 0; i < 6; i++) {
    
    //var svg = d3.select("#cor" + (i+1) + " svg").node ()//svg.node ()//document.getElementById(domnode + "svg");
    //console.log (svg)
    //svg.setAttribute("viewBox", "0 0 600 600")
    //svg.setAttribute("width",  "150")
    //svg.setAttribute("height", "150")
  //}
  
};


d3.select("#download").on("click", function(){
  d3.select(this)
    .attr("href", 'data:application/octet-stream;base64,' + btoa(d3.select("#my_dataviz").html()))
    .attr("download", "plot.svg") 
})




$.getJSON( "data/data.json", function( data ) {  
  console.log( "success" );
  console.log( data );
  samples = data["samples"]
  global_sumstat = data["sumstat"]
  
  for (const [metric_key, metric_value] of Object.entries(metric_dict)) {
  for (const [sumstat_key, sumstat] of Object.entries(global_sumstat[metric_key])) {
    for (const [key, boxplot] of Object.entries(sumstat)) {
      //console.log (boxplot)
      const outliers_min = boxplot[8];
      const outliers_max = boxplot[9];
      const inf_p = boxplot[13];
      const inf_m = boxplot[14];
      
      //outliers_min.sort (function (a, b) {return a[0] < b[0] ? -1 : 1});
      //outliers_max.sort (function (a, b) {return a[0] < b[0] ? -1 : 1});
      
      var [source, type] = dataset2sourcetype (data_sources[boxplot[11]]);
      var imeth = imethods[boxplot[12]]
      
			var outlier_table = "<div><h3>"+(outliers_min.length+outliers_max.length)+" outliers for "+ metric_key + " of " + source + " -- " + type + " using " + imeth +"</h3><table class='outliers table'><thead><tr><th>Sample</th><th>Value</th></tr></thead><tbody>";

			for (var o = 0; o < inf_m.length; o++) {
			  outlier_table += "<tr><td>"+get_sample_name (inf_m[o])+"</td><td>-Infinity</td></tr>";
			}
      
			for (var o = 0; o < outliers_min.length; o++) {
			  outlier_table += "<tr><td>"+get_sample_name (outliers_min[o])+"</td><td>"+outliers_min[o][0]+"</td></tr>";
			}
			outlier_table += "<tr><th>--- MEDIAN ---</td><th>"+boxplot[1]+"</th></tr>";
			for (var o = 0; o < outliers_max.length; o++) {
			  outlier_table += "<tr><td>"+get_sample_name (outliers_max[o])+"</td><td>"+outliers_max[o][0]+"</td></tr>";
			}
			for (var o = 0; o < inf_p.length; o++) {
			  outlier_table += "<tr><td>"+get_sample_name (inf_p[o])+"</td><td>Infinity</td></tr>";
			}
			outlier_table += "</tbody></table></div>"
      boxplot[16] =  (outlier_table)
    }
  }
    addSelect (metric_key, metric_value)
	}
  
  if (!drawn)
    draw_boxplots ();
})
  .done(function(data) {
  })
  .fail(function(d, textStatus, error) {
    // TODO !!!!
    console.log( "error" );
    console.log( d );
        console.error("getJSON failed, status: " + textStatus + ", error: "+error)
  });





//function do_sumstat (metric) {
	
    
    
    
    
    
  //console.log ("starting sumstat")
  //var sumstat = []
  //var minY = 10000000;
  //var maxY = 0;
  //var infinites_p = false
  //var infinites_m = false
  //var xdomain = new Set ();
  
  
  
  //if (typeof(Storage) !== "undefined") {
	  
	  //cache = localStorage.getItem("sumstat_" + metric)
	  ////console.log ("found in cache:")
	  ////console.log (global_sumstat[metric])
	  ///*if (cached !== null) {
		  //sumstat = cached["sumstat"]
		  //minY = cached["minY"]
		  //maxY = cached["maxY"]
		  //infinites_p = cached["infinites_p"]
		  //infinites_m = cached["infinites_m"]
		  //xdomain = cached["xdomain"]
	  //}*/
  //}
  //if (cache === null) {
	  //for (const [type1key, type1value] of Object.entries(boxplots[metric])) {
		//for (const [type2key, type2value] of Object.entries(type1value)) {
		  //for (const [methkey, methvalue] of Object.entries(type2value)) {
			//var key = type1key+"-"+type2key+"-"+methkey
			 //xdomain.add (key);
			 //scoreId = metric + "_" + methkey;
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
			//outliers_min = methvalue["samples"].filter (x => (x.scores[scoreId] < whiskersMin) && vals.includes (x.scores[scoreId]));
			//outliers_max = methvalue["samples"].filter (x => (x.scores[scoreId] > whiskersMax) && vals.includes (x.scores[scoreId]));
			//outliers_min.sort (function (a, b) {return a.scores[scoreId] < b.scores[scoreId] ? -1 : 1});
			//outliers_max.sort (function (a, b) {return a.scores[scoreId] < b.scores[scoreId] ? -1 : 1});
			
			//var outlier_table = "<div><h3>"+(outliers_min.length+outliers_max.length)+" outliers for "+ metric + " of " + key +"</h3><table class='outliers table'><thead><tr><th>Sample</th><th>Value</th></tr></thead><tbody>";

			//for (var o = 0; o < outliers_min.length; o++) {
			  //outlier_table += "<tr><td>"+outliers_min[o].name+"</td><td>"+outliers_min[o].scores[scoreId]+"</td></tr>";
			//}
			//outlier_table += "<tr><th>--- MEDIAN ---</td><th>"+median+"</th></tr>";
			//for (var o = 0; o < outliers_max.length; o++) {
			  //outlier_table += "<tr><td>"+outliers_max[o].name+"</td><td>"+outliers_max[o].scores[scoreId]+"</td></tr>";
			//}
			//outlier_table += "</tbody></table></div>"
			
			//sumstat.push ({
			  //"key": key,
			  //"value": {q1: q1, median: median, q3: q3, interQuantileRange: interQuantileRange, min: min, max: max, whiskersMin: whiskersMin, whiskersMax: whiskersMax, outliers_min: outliers_min, outliers_max: outliers_max, scoreId: scoreId, inf_p: inf_p, inf_m: inf_m, outlier_table: outlier_table}});
			//if (minY > min && min != -1000 && min != 1000)
			  //minY = min;
			//if (maxY < max && max != 1000 && max != -1000)
			  //maxY = max;
		  //}
		//}
	  //}
	  //global_sumstat[metric] = {
		//"sumstat": sumstat,
		//"minY": minY,
		//"maxY": maxY,
		//"infinites_p": infinites_p,
		//"infinites_m": infinites_m,
		//"xdomain": xdomain
	  //}
	  //console.log (metric)
	  //console.log (global_sumstat[metric])
	  ////if (typeof(Storage) !== "undefined") {
		  ////localStorage.setItem("sumstat_" + metric, JSON.stringify(global_sumstat[metric]))
	  ////}
  //} else {
	  //global_sumstat[metric] = JSON.parse(cache)
  //}
  //console.log ("done sumstat")
  
  
  //addSelect (metric)
  
  //if (!drawn)
    //draw_boxplots ();
  
//}







//d3.csv("data/combined-afr-eor-hallmark.csv").then (function(data) {
  //var metrics = new Set ()
  //for (row=0; row < data.length; row++){
    //metrics.add (data[row]["Score"])
    //var [source, type] = dataset2sourcetype (data[row]["Dataset"])
    
		//if (!samples[source + "__" + type + "__single"][data[row]["Sample"]]) {
			//samples[source + "__" + type + "__single"][data[row]["Sample"]] = new Sample (data[row]["Sample"], source, type)
		//}
		
    //const sample = samples[source + "__" + type + "__single"][data[row]["Sample"]]
		//if (!boxplots[data[row]["Score"]]) {
			//boxplots[data[row]["Score"]] = getSourceData ();
    //}
    
    //for (i = 0; i < imethods.length; i++) {
      //boxplots[data[row]["Score"]][source][type][imethods[i]]["samples"].push (sample);
      //sample.setScore (data[row]["Score"], imethods[i], parseFloat (data[row][imethods[i]]))
    //}
  //}
  //for (let metric of metrics)
    //setTimeout (function () {do_sumstat (metric)}, 10);
//});


//d3.csv("data/combined-jaccard-ba.csv").then (function(data) {
  //var metrics = new Set ()
  //for (row=0; row < data.length; row++){
    //metrics.add (data[row]["Score"])
    //var [source, type] = dataset2sourcetype (data[row]["Dataset"])
		
    //sid = data[row]["Dx1"] + " -vs- " + data[row]["Dx2"]
    
		//if (!samples[source + "__" + type + "__double"][sid]) {
			//samples[source + "__" + type + "__double"][sid] = new Sample (sid, source, type)
		//}
    
    //const sample = samples[source + "__" + type + "__double"][sid]
		//if (!boxplots[data[row]["Score"]]) {
			//boxplots[data[row]["Score"]] = getSourceData ();
    //}
    
    //for (i = 0; i < imethods.length; i++) {
      //boxplots[data[row]["Score"]][source][type][imethods[i]]["samples"].push (sample);
      //sample.setScore (data[row]["Score"], imethods[i], parseFloat (data[row][imethods[i]]))
    //}
  //}
  
  //for (let metric of metrics)
    //setTimeout (function () {do_sumstat (metric)}, 1000*Math.random ());
  
  ///*m_select = document.getElementById("metric")
  //for (let metric of metrics) {
    //var opt = document.createElement('option');
    //opt.value = metric;
    //opt.innerHTML = metric;
    //m_select.appendChild(opt);
  //}
  
  //sortSelect (m_select)
  
  //if (!drawn)
    //draw_boxplots ();*/
//});

