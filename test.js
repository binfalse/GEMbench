

//run ()



  
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

imethods = ["GIMME","FASTCORE","INIT","iMAT"];
function getMeth () {
  marr = {};
  for (i = 0; i < imethods.length; i++)
    marr[imethods[i]] = {
      "points": [],
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
  return {
  "RNA Seq": getMeth (),
  "Microarray": getMeth (),
  "MS Proteomics": getMeth ()
  };
};

boxplots = {
"Patient Data": getType (),
"Cell Line": getType ()
};
minY = 10000000;
maxY = 0;

var xdomain = [];

d3.csv("data/combined.csv").then (function(data) {
  // console.log(data)
  for (row=0; row < data.length; row++){
    t=undefined;
    b=undefined;
    switch (data[row]["DS"]) {
      case "EMTAB37":
        t = "Microarray";
        b = "Cell Line";
        break;
      case "HPA":
        t = "RNA Seq";
        b = "Cell Line";
        break;
      case "ProteomeNCI60":
        t = "MS Proteomics";
        b = "Cell Line";
        break;
      case "GSE2019":
        t = "Microarray";
        b = "Patient Data";
        break;
      case "TCGA":
        t = "RNA Seq";
        b = "Patient Data";
        break;
      case "ProteomePatients":
        t = "MS Proteomics";
        b = "Patient Data";
        break;
    }
    for (i = 0; i < imethods.length; i++) {
      boxplots[b][t][imethods[i]]["points"].push (parseInt (data[row][imethods[i]]));
      xdomain.push (b+"-"+t+"-"+imethods[i]);
    }
  }
  
  console.log (minY,maxY)
  
 // set the dimensions and margins of the graph
var margin = {top: 50, right: 30, bottom: 60, left: 40},
    width = 1000 - margin.left - margin.right,
    height = 800 - margin.top - margin.bottom;

// append sthe svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

// Read the data and compute summary statistics for each specie
// d3.csv("iris.csv").then (function(data) {

    //width = 160/4,
    //height = 400;
  var sumstat = []
  for (const [type1key, type1value] of Object.entries(boxplots)) {
    for (const [type2key, type2value] of Object.entries(type1value)) {
      for (const [methkey, methvalue] of Object.entries(type2value)) {
        q1 = d3.quantile(methvalue["points"].sort(d3.ascending),.25);
        median = d3.quantile(methvalue["points"].sort(d3.ascending),.5);
        q3 = d3.quantile(methvalue["points"].sort(d3.ascending),.75);
        interQuantileRange = q3 - q1;
        min = Math.min(...methvalue["points"])
        max = Math.max(...methvalue["points"])
        //min = q1 - 1.5 * interQuantileRange;
        //max = q3 + 1.5 * interQuantileRange;
        whiskersMin = Math.max(min, q1 - interQuantileRange * 1.5);
        whiskersMax = Math.min(max, q3 + interQuantileRange * 1.5);
        outliers = methvalue["points"].filter (x => x < whiskersMin || x > whiskersMax);
        sumstat.push ({
          "key": type1key+"-"+type2key+"-"+methkey,
          "value": {q1: q1, median: median, q3: q3, interQuantileRange: interQuantileRange, min: min, max: max, whiskersMin: whiskersMin, whiskersMax: whiskersMax, outliers: outliers}});
        if (minY > min)
          minY = min;
        if (maxY < max)
          maxY = max;
      }
    }
  }
  console.log (minY, maxY);
  
  
  var boxtip = d3.tip().attr('class', 'd3-tip').direction('e').offset([0,5])
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
  var outliertip = d3.tip().attr('class', 'd3-tip').direction('e').offset([0,5])
            .html(function(d) {
                var content = "<span style='margin-left: 2.5px;'><b>" + d.n + "</b></span><br>";
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
        
  
  var y = d3.scaleLinear()
    .domain([.9*minY,maxY*1.1])
    .range([height, 0])
  
  
  
  var boxWidth = width / sumstat.length
  svg
    .selectAll("databackground")
    .data(sumstat)
    .enter()
    .append("rect")
        .attr("x", function(d){return(x(d.key)-boxWidth/2)})
        .attr("y", y(maxY*1.1))
        .attr("height", y(.9*minY))
        .attr("width", boxWidth )
        .style("opacity", ".2")
        .style("fill", function(d){return datacolor(d.key)})
        
        
        
  console.log (sumstat);
  
  for (i = 0; i < sumstat.length; i++) {
    const arr = []
    for (j = 0; j < sumstat[i].value.outliers.length; j++)
      arr.push ({n: sumstat[i].key, p: sumstat[i].value.outliers[j]})
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
  
  //var boxWidth = width / sumstat.length - 10
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
  //svg
    //.selectAll("medianCircles")
    //.data(sumstat)
    //.enter()
    //.append("circle")
      //.attr("cx", function(d){return(x(d.key)) })
      //.attr("cy", function(d){return(y(d.value.median))})
      //.attr("r", boxWidth/2)
      //.attr("stroke", "black")
      //.style("fill", function(d){return boxColor(d.key)})
  
  
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
  
  
  //console.log (xg)
  //console.log (xg.node ().getBBox ())
  for (i =0; i < types1.length; i++) {
    console.log (types1[i]);
    txt = svg.append("text")
    .attr("x", (x(types1[i]["start"]) + x(types1[i]["end"]))/2)
    .attr("text-anchor","middle")
    .attr("y", -margin.top + 30)
    //.attr("y", height + xg.node ().getBBox ().width + 80)
    //.attr("stroke","black")
    .text(types1[i]["name"]);
    //bb = txt.node ().getBBox ()
    //oben = bb.y - 10
    //unten = 0//bb.y + bb.height + 20
    //left = x (types1[i]["start"]) - width / (2*sumstat.length)
    //right = x (types1[i]["end"]) + width / (2*sumstat.length)
    //svg.append("path")
    //.attr("d", d3.line()([[left, unten], [left, oben], [right,oben], [right,unten]]))
    //.attr("stroke","black")
    //.style("fill", "none")
    //.style("stroke-opacity", ".3")
    //.attr("stroke-width", ".3")
  }
  for (i =0; i < types2.length; i++) {
    console.log (types2[i]);
    txt = svg.append("text")
    .attr("x", (x(types2[i]["start"]) + x(types2[i]["end"]))/2)
    .attr("text-anchor","middle")
    //.attr("y", y (minY*.9)-10)
    .attr("y", 20)
    //.attr("y", height + xg.node ().getBBox ().width + 40)
    //.attr("stroke","black")
    .text(types2[i]["name"]);
    //bb = txt.node ().getBBox ()
    //oben = bb.y - 10
    //unten = bb.y + bb.height + 5
    //left = x (types2[i]["start"]) - width / (2*sumstat.length)
    //right = x (types2[i]["end"]) + width / (2*sumstat.length)
    //svg.append("path")
    //.attr("d", d3.line()([[left, oben], [left, unten], [right,unten], [right,oben]]))
    //.attr("stroke","black")
    //.style("fill", "none")
    //.style("stroke-opacity", "1")
  }
  
  
  
  
  
  
  
  
  
  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(x))
      .selectAll("text")
        .attr("transform", "translate(-10,10)rotate(-90)")
        .style("text-anchor", "end")
        .style("font-size", 28)
    .each(function(d, i){
      d3.select(this).text(d.replace (/.*-/g,''));
    })
  svg.append("g").call(d3.axisLeft(y))
  
  
  
  
  
  
  
  
console.log(xdomain)
});


    d3.select("#download").on("click", function(){
      d3.select(this)
        .attr("href", 'data:application/octet-stream;base64,' + btoa(d3.select("#my_dataviz").html()))
        .attr("download", "plot.svg") 
    })
//var type1 = d3.select("#my_dataviz").append ("div").attr("class", "row");
//var type2 = d3.select("#my_dataviz").append ("div").attr("class", "row");
//var method = d3.select("#my_dataviz").append ("div").attr("class", "row");
//var bp = d3.select("#my_dataviz").append ("div").attr("class", "row");
////for (b = 0; b < boxplots.length; b++) {
//for (const [type1key, type1value] of Object.entries(boxplots)) {
  //type1.append ("div").attr("class", "text-center col-sm-" + (12/Object.entries(boxplots).length)).text (type1key);
  //for (const [type2key, type2value] of Object.entries(type1value)) {
    //type2.append ("div").attr("class", "text-center col-sm-" + (6/Object.entries(type1value).length)).text (type2key);
    //meth = method.append ("div").attr("class", "text-center col-sm-" + (6/Object.entries(type1value).length)).append ("div").attr("class", "row");
    //bp = method.append ("div").attr("class", "text-center col-sm-" + (6/Object.entries(type1value).length)).append ("div").attr("class", "row");
    //for (const [methkey, methvalue] of Object.entries(type2value)) {
      //meth.append ("div").attr("class", "method col-sm-" + (12/Object.entries(type2value).length)).text (methkey);
      //var svg = bp.append("svg")
                  //.attr("width", width + margin.left + margin.right)
                  //.attr("height", height + margin.top + margin.bottom)
                //.append("g")
                  //.attr("transform",
                        //"translate(" + margin.left + "," + margin.top + ")");
      
      
      
      
      
    //}
  //}
//}























