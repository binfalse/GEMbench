
const sources = ["Cell Line","Patient Data"];
const types = ["Microarray","RNA Seq","MS Proteomics"];
const imethods = ["FASTCORE","GIMME","INIT","iMAT"];
const metrics = ["AFR","EOR","Hallmark","BlandAltman","Jaccard","Clusterability"]
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
};

function getSourceData () {
  s = {}
  for (var i = 0; i < sources.length; i++)
    s[sources[i]] = getType ()
  return s
};


function update_slider_value () {
  slider = document.getElementById("slider")
  document.getElementById("slider_value").innerHTML = slider.value
}

function get_outlier_table_id (metric, data_id) {
  return "outlier__" + metric + "__" + data_id.replace(/[\W_]+/g,"_")
}

var margin = {top: 50, right: 0, bottom: 70, left: 40},
    width = d3.select("#my_dataviz").node().clientWidth - margin.left - margin.right,
    height = 900 - margin.top - margin.bottom;












var drawn = false
var samples = []
const pca = {}

$(".metric_expl").hide ()

const subsystems = []
const affected_subsystems = {}









function draw_pca () {
  //console.log (pca);
        
        
      d3.select("#my_dataviz svg").remove();
      
// append sthe svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    
  
  
var pattern = svg.append("defs")
	.append("pattern")
  .attr ("id", "fillpattern")
  .attr ("width", "8")
  .attr ("height", "8")
  .attr ("patternUnits", "userSpaceOnUse")
  .attr ("patternTransform", "rotate(60)")
		//.attr({ id:"fillpattern", width:"8", height:"8", patternUnits:"userSpaceOnUse", patternTransform:"rotate(60)"})
	.append("rect")
  .attr ("width", "4")
  .attr ("height", "8")
  .attr ("transform", "translate(0,0)")
  .attr ("fill", "#f0f0f0")
  
  svg = svg.append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");
  
  var size_y = height / imethods.length;
  var size_x = width / data_sources.length;
  
  var pcatip = d3.tip().attr('class', 'd3-tip').direction("e").offset([0,-5])
            .html(function(d) {
                var content = "<span style='margin-left: 2.5px;'><b>" + d[0] + "</b></span><br>";
                content += "<span style='margin-left: 2.5px;'>PC1: " + d3.format(".2f")(d[1]) + "</span><br>";
                content += "<span style='margin-left: 2.5px;'>PC2: " + d3.format(".2f")(d[2]) + "</span><br>";
                return content;
            });
        svg.call(pcatip);
  
  var position_x = d3.scalePoint()
    .domain(data_sources)
    .range([0, width-size_x]);
  var position_y = d3.scalePoint()
    .domain(imethods)
    .range([0, height-size_y]);
  
  var nBarcodeColors = 5
  
  var color = d3.scaleLinear()
    .domain([...Array(nBarcodeColors).keys()])
    .range(["#e41a1c","#377eb8","#4daf4a","#984ea3", "#ff7f00"])
  var colorBackground = d3.scaleLinear()
    .domain([...Array(nBarcodeColors).keys()])
    .range(["#aaa", "#eee", "#ccc", "#ddd", "#bbb"])
  
  var padding_x = 30
  var padding_y = 80
  
  for (var s = 0; s < data_sources.length; s++) {
    for (var i = 0; i < imethods.length; i++) {
      const source = data_sources[s];
      const imeth = imethods[i];
      const pca_samples = pca[source + "_" + imeth]
      
      //console.log (source, imeth, position_x(imeth), size_x, position_x(imeth) + size_x)
      
      var tmp = svg
        .append('g')
        .attr("transform", "translate(" + (position_x(source)) + "," + (position_y(imeth)) + ")");
        
        
          //tmp.append("rect")
             //.attr("x", 0)
              //.attr("y", 0)
             //.attr("width", size_x)
             //.attr("height", size_y)
             //.style("fill", function(d,i){return color(Math.random ())})
        
        
        //ax = Math.min (size_x, size_y) - 20;
        
        xextent = d3.extent(pca_samples, function(d) { return d[1] })
        //console.log (xextent)
      var x = d3.scalePow().exponent(document.getElementById("slider").value)
        .domain(xextent).nice()
        .range([ padding_x, size_x-10 ]);

      // Add Y Scale of each graph
      yextent = d3.extent(pca_samples, function(d) { return d[2] })
        //console.log (yextent)
      var y = d3.scalePow().exponent(document.getElementById("slider").value)
        .domain(yextent).nice()
        .range([ 0, size_y-padding_y ]);

      // Add a 'g' at the right position
      //var tmp = svg
        //.append('g')
        //.attr("transform", "translate(" + (position(var1)+mar) + "," + (position(var2)+mar) + ")");

      // Add X and Y axis in tmp
      tmp.append("g")
        .attr("transform", "translate(" + 0 + "," + (size_y-padding_y) + ")")
        .call(d3.axisBottom(x).ticks(3));
      tmp.append("g")
        .attr("transform", "translate(" + padding_x + "," + 0 + ")")
        .call(d3.axisLeft(y).ticks(3));
             
             tmp
             .selectAll("myCircles")
             .data(pca_samples)
             .enter()
             .append("circle")
             .attr("cx", function(d){ return x(d[1]) })
             .attr("cy", function(d){ return y(d[2]) })
             .attr("r", 2)
             .attr("fill", "#000")
            .on('mouseover', pcatip.show)
            .on('mouseout', pcatip.hide)
             
      //tmp.append("text")
      //.attr("y", 0)
      //.attr("x", 0)
      //.text(source+"_"+imeth)
      //.style("font-size", 11)
      //.style("text-align", "center")
    //.style("text-anchor","middle");
    
    
    // draw barcode
    
      console.log (source + "_" + imeth)
      
      if (affected_subsystems[source + "_" + imeth]) {
        console.log (affected_subsystems[source + "_" + imeth])
        const unitl = (size_x - 10) / affected_subsystems[source + "_" + imeth].reduce (function(acc, val) { return acc + val; }, subsystems.length)
        var cum = 10;
        for (var sub = 0; sub < subsystems.length; sub++) {
          var c = "#fff"
          if (affected_subsystems[source + "_" + imeth][sub] > 0)
            c = color (sub % nBarcodeColors)
          else
            c = colorBackground (sub % nBarcodeColors)
          tmp.append("rect")
             .attr("x", cum)
              .attr("y", size_y - padding_y + 25)
             .attr("width", (1 + affected_subsystems[source + "_" + imeth][sub]) * unitl)
             .attr("height", 20)
             //.style("fill", color(Math.random ()))
             .style("fill", c)
          cum += (1 + affected_subsystems[source + "_" + imeth][sub]) * unitl
        }
        
      } else {
        
          tmp.append("rect")
             .attr("x", 10)
              .attr("y", size_y - padding_y + 25)
             .attr("width", size_x - 10)
             .attr("height", 20)
             //.style("fill", color(Math.random ()))
             .style("fill", "url(#fillpattern)")
      }
    
    
    }
  }
  
  for (var s = 0; s < data_sources.length; s++) {
    svg.append("text")
      .attr("y", 0)
      .attr("x", position_x(data_sources[s]) + size_x/2 + padding_x)
      .text(data_sources[s])
      .style("font-size", 11)
      .style("text-align", "center")
    .style("text-anchor","middle");
  }
  for (var i = 0; i < imethods.length; i++) {
    svg.append("text")
      .attr("x", -(position_y(imethods[i]) + size_y/2 - padding_y/2))
      .attr("y", 0)
      .text(imethods[i])
      .style("font-size", 11)
      .style("text-align", "center")
    .style("text-anchor","middle")
            .attr("transform", function(d) {
                return "rotate(-90)" 
                });
  }
  
  
  
}





function draw_correlation (measure, sumstat, source, type, domnode) {
  //console.log (sumstat)
  //console.log (measure, source, type)
  //console.log (getDataFromSource (source, type))
  var data_id = data_dict[getDataFromSource (source, type)]
  var measure_id = metric_dict[measure]
  
  var cor_samples = {}
  
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
  
  $("#" + domnode).empty ();
  $("#" + domnode).append ("<p><a title='Correlation of "+measure+" scores for "+ source + " &mdash; " + type + "' id='"+domnode+"_a' href='#"+domnode+"_svg'>" + source + "<br>" + type + "</a></p>");
  
var marginWhole = {top: 10, right: 10, bottom: 10, left: 10},
    sizeWhole = 600 - marginWhole.left - marginWhole.right

// Create the svg area
var svg = d3.select("#" + domnode + " p a")
  .append("svg")
    .attr("width", 150)
    .attr("height", 150)
    .attr("id", domnode+"_svg")
  .attr ("viewBox", "0 0 " + (sizeWhole  + marginWhole.left + marginWhole.right) + " " + (sizeWhole  + marginWhole.left + marginWhole.right))
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
      //console.log ("test", domnode)
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
  
  var color = d3.scaleLinear()
    .domain([-1, 0, 1])
    .range(["#2b9d67", "#fffbf0", "#2b9d67"])
    .unknown ("#fff");
  
  //console.log ("n samples: ", Object.values (cor_samples).length)
  
  
  
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





function draw_boxplots (measure) {
  update_slider_value ()
  drawn = true
  //sqrtscale = document.getElementById('sqrt').checked
      
      d3.select("#my_dataviz svg").remove();
      
// append sthe svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

  var sumstat = global_sumstat[measure]["boxplots"];
  var n_boxplots = Object.keys(sumstat).length;
  var minY  = global_sumstat[measure]["minY"];
  var maxY  = global_sumstat[measure]["maxY"];
  var infinites_p = global_sumstat[measure]["infinites_p"];
  var infinites_m = global_sumstat[measure]["infinites_m"];
  var xdomain = global_sumstat[measure]["xdomain"];
  //console.log (minY, maxY);
  var columnWidth = width / n_boxplots;
  // rectangle for the main box
  var boxWidth = width / n_boxplots - 20;
  //console.log (xdomain)
  //console.log (sumstat)
  
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
        .attr("x", function(d){return(x(get_x_for_boxplot(d))-columnWidth/2)})
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
				//console.log ("not too many outliers to draw")
				//console.log (outlier_table)
              $.fancybox( $(outlier_table) );
            });
      } else {
        //console.log ("too many outliers to draw", outliers.length)
        
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
				//console.log ("too many outliers to draw")
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
        draw_correlation (measure, sumstat, t1, t2, "cor" + tmpCorId)
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
  
  
  
};





















function draw_analysis () {
  var measure = document.getElementById("metric")
  //console.log (measure)
  measure = metrics[measure.options[measure.selectedIndex].value]
$(".metric_expl").hide ()
  $("#" + measure + "_expl").show ();
  //console.log (measure)
  if (measure == "Clusterability") {
    draw_pca ()
  } else {
    draw_boxplots (measure)
  }
}











d3.select("#download").on("click", function(){
  d3.select(this)
    .attr("href", 'data:application/octet-stream;base64,' + btoa(d3.select("#my_dataviz").html()))
    .attr("download", "plot.svg") 
})















$.getJSON("data/data.json", function( data ) {
  //console.log( "success" );
  //console.log( data );
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
      
			var outlier_table = "<div><h3>"+(outliers_min.length+outliers_max.length)+" outliers for "+ metric_key + " of " + source + ": " + type + " using " + imeth +"</h3><table class='outliers table'><thead><tr><th>Sample</th><th>Value</th></tr></thead><tbody>";

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
    draw_analysis ();
})
  .fail(function(d, textStatus, error) {
    // TODO !!!!
    console.log( "error" );
    console.log( d );
        console.error("getJSON failed, status: " + textStatus + ", error: "+error)
  });





d3.csv("data/affects.csv").then (function( data ) {
  for (row=0; row < data.length; row++){
    if (!subsystems.includes (data[row]["subsystem"]))
      subsystems.push (data[row]["subsystem"]);
  }
  subsystems.sort();
  
  for (row=0; row < data.length; row++){
    const id = data[row]["Data"] + "_" + data[row]["Method"];
    if (!(id in affected_subsystems))
      affected_subsystems[id] = Array(subsystems.length).fill(0)
    affected_subsystems[id][subsystems.indexOf(data[row]["subsystem"])]++
  }
  
  //console.log (subsystems);
  //console.log (affected_subsystems);
});



var pcas_expected = 0;
for (var s = 0; s < data_sources.length; s++) {
  for (var i = 0; i < imethods.length; i++) {
    pcas_expected++;
    const source = data_sources[s]
    const imeth = imethods[i]
    //console.log ("data/pca/" + source + "_" + imeth + ".csv")
    pca[source + "_" + imeth] = undefined
    
    d3.csv("data/pca/" + source + "_" + imeth + ".csv").then (function( data ) {
      //console.log (data)
      pca[source + "_" + imeth] = []
      for (row=0; row < data.length; row++){
        pca[source + "_" + imeth].push ([
          data[row][""],
          parseFloat(data[row]["PC1"]),
          parseFloat(data[row]["PC2"])
        ])
      }
      if (pcas_expected-- == 1) {
        addSelect ("Clusterability", 5)
      }
      //console.log (pcas_expected)
    });
  }
}

