// based on https://bl.ocks.org/kerryrodden/766f8f6d31f645c39f488a0befa1e3c8


function hash (node) {
  var s = node.data.name;
  if (node.parent)
    s += node.parent.data.name;
    
    var h = 0, i, chr;
    for (i = 0; i < s.length; i++) {
      chr   = s.charCodeAt(i);
      h  = ((h << 5) - h) + chr;
      h |= 0; // Convert to 32bit integer
    }
    return h;
}

function shortenLabel (label, arcsize) {
  if (label == "root")
    return "";
  if (label.includes ("("))
    label = label.replace (/ \(.*/, '')
  if (label.length < arcsize*60)
    return label;
  return label.substring (0, Math.min (label.length, arcsize*65)) + "..."
  //return arcsize + label;
}

// Dimensions of sunburst.
var width = 750;
var height = 600;
var radius = Math.min(width, height) / 2;

// Breadcrumb dimensions: width, height, spacing, width of tip/tail.
var b = {
  w: 10, h: 30, s: 3, t: 10
};

// Mapping of step names to colors.
var colors = {
  "Patient": "rgb(22,79,134)",
  
  "Cell-line": "rgb(142,0,0)",
  
  "Microarray (EMTAB-37)": "rgb(218,65,65)",
  "RNA-Seq (HPA)": "rgb(244,91,91)",
  "Mass Spec (ProteomeNCI60)": "rgb(255,161,161)",
  
  "Microarray (GSE2109)": "rgb(73,130,185)",
  "RNA-Seq (TCGA)": "rgb(99,156,211)",
  "Mass Spec (ProteomePatients)": "rgb(166,223,255)"
};
var fgColors = {
  "Patient": "#fff",
  
  "Cell-line": "#fff",
  
  "Microarray (EMTAB-37)": "#fff",
  "RNA-Seq (HPA)": "#fff",
  "Mass Spec (ProteomeNCI60)": "#000",
  
  "Microarray (GSE2109)": "#fff",
  "RNA-Seq (TCGA)": "#000",
  "Mass Spec (ProteomePatients)": "#000"
};



function colorMap (node) {
  if (node.data.name == "root")
    return "#000";
  if (node.data.name == "Proteome") {
    if (node.parent.data.name == "Patient")
      return "rgb(124,181,236)"
    else
      return "rgb(255,119,119)"
  }
  if (node.data.name == "Transcriptome") {
    if (node.parent.data.name == "Patient")
      return "rgb(48,105,160)"
    else
      return "rgb(168,15,15)"
  }
  if (colors[node.data.name])
    return colors[node.data.name];
  return colors[node.parent.data.name];
}
function fgColorMap (node) {
  if (node.data.name == "root")
    return "#000";
  if (node.data.name == "Proteome") {
    if (node.parent.data.name == "Patient")
      return "#000"
    else
      return "#000"
  }
  if (node.data.name == "Transcriptome") {
    if (node.parent.data.name == "Patient")
      return "#fff"
    else
      return "#fff"
  }
  if (fgColors[node.data.name])
    return fgColors[node.data.name];
  return fgColors[node.parent.data.name];
}


// Total size of all segments; we set this later, after loading the data.
var totalSize = 0; 

var vis = d3.select("#sunburst").append("svg:svg")
  .attr ("viewBox", "0 0 " + width + " " + height)
    .append("svg:g")
    .attr("id", "container")
    .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

var partition = d3.partition()
    .size([2 * Math.PI, radius * radius]);

var arc = d3.arc()
    .startAngle(function(d) { return d.x0; })
    .endAngle(function(d) { return d.x1 - 0.002; })
    .innerRadius(function(d) { return Math.sqrt(d.y0) + 1; })
    .outerRadius(function(d) { return Math.sqrt(d.y1) - 1; });

  console.log ("test");
d3.dsv(";","/data/data_viz_sunburst.csv", function(d) {
        return [
          d.Source,
          d.Location,
          d.Type,
          d.Cancer,
          d.Samples
        ];
   }).then (function(text) {
  console.log ("test2");
  console.log (text);
  var json = buildHierarchy(text);
  console.log (json);
  createVisualization(json);
})
  .catch(error => console.log(error));;

// Main function to draw and set up the visualization, once we have the data.
function createVisualization(json) {

  vis.append("svg:circle")
      .attr("r", radius)
      .style("opacity", 0);

  // Turn the data into a d3 hierarchy and calculate the sums.
  var root = d3.hierarchy(json)
      .sum(function(d) { return d.size; })
      .sort(function(a, b) { return b.value - a.value; });
  
  // For efficiency, filter nodes to keep only those large enough to see.
  var nodes = partition(root).descendants();
  var labelled_nodes = partition(root).descendants()
      .filter(function(d) {
          return (d.x1 - d.x0 > 0.1); // 0.005 radians = 0.29 degrees
      });
      
      console.log (nodes);


  var path = vis.data([json]).selectAll("path")
      .data(nodes)
      .enter().append("svg:path")
      .attr("display", function(d) { return d.depth ? null : "none"; })
      .attr("d", arc)
      .style("fill", function(d) { return colorMap (d); })
      .style("opacity", 1)
      .on("mouseover", mouseover)
    .each(function(d,i) {
      // based on https://www.visualcinnamon.com/2015/09/placing-text-on-arcs.html
        var firstArcSection = /(^.+?)L/;
        var newArc = firstArcSection.exec( d3.select(this).attr("d") )[1];
        newArc = newArc.replace(/,/g , " ");

        //Create a new invisible arc that the text can flow along
        vis.append("path")
            .attr("class", "hiddenarcs")
            .attr("id", "arc"+hash (d))
            .attr("d", newArc)
            .style("fill", "none");
    });
      
  var path2 = vis.data([json]).selectAll("names")
      .data(labelled_nodes)
      .enter().append("text")
      .attr("dy", 20)
      .append("textPath")
      .attr("xlink:href",function(d) { return "#arc" + hash (d)})
    .style("text-anchor","middle") 
    .attr("startOffset", "50%")
    .attr("font-size","8px")
      .attr("fill",function(d) { return fgColorMap (d); })
      .text(function(d) { return shortenLabel (d.data.name, d.x1 - d.x0)});

  // Add the mouseleave handler to the bounding circle.
  //d3.select("#container").on("mouseleave", mouseleave);

  // Get total size of the tree = value of root node from partition.
  totalSize = path.datum().value;
 };


// Fade all but the current sequence, and show it in the breadcrumb trail.
function mouseover(d) {

  var percentage = (100 * d.value / totalSize).toPrecision(3);
  var percentageString = percentage + "%";
  if (percentage < 0.1) {
    percentageString = "< 0.1%";
  }

  d3.select("#percentage")
      .text(percentageString);

  d3.select("#sunburst-explanation")
      .style("visibility", "");

  var sequenceArray = d.ancestors().reverse();
  sequenceArray.shift(); // remove root node from the array
  updateBreadcrumbs(sequenceArray, d.value);

  // Fade all the segments.
  d3.selectAll("path")
      .style("opacity", 0.3);

  // Then highlight only those that are an ancestor of the current segment.
  vis.selectAll("path")
      .filter(function(node) {
                return (sequenceArray.indexOf(node) >= 0);
              })
      .style("opacity", 1);
}


// Update the breadcrumb trail to show the current sequence and percentage.
function updateBreadcrumbs(nodeArray, nSamples) {
console.log (nodeArray);

$("#sunburst-bc").empty ();
for (var i = 0; i < nodeArray.length; i++) {
  var arrow = "<span class='sbbcnoarrow' ></span>";
  var ml = "";
  var numbers = "";
  
  if (i < nodeArray.length - 1)
    arrow = "<span class='sbbcarrow' style='border-left-color:"+colorMap(nodeArray[i])+";' ></span>";
  if (i > 0)
    ml = ";margin-left:-1em;"
  if (i == nodeArray.length - 1)
    numbers = "<div class='sbbcpc'>" + nSamples + " samples</div>";
  
  $("#sunburst-bc").append ("<div class='my-2 sbbc' style='background-color:"+colorMap(nodeArray[i])+";color:"+fgColorMap(nodeArray[i])+ml+"'>" + nodeArray[i].data.name + arrow + "</div>" + numbers);
}



}


function buildHierarchy(csv) {
  var root = {"name": "root", "children": []};
  for (var i = 0; i < csv.length; i++) {
    var l = csv[i].length;
    var parts = csv[i].slice (0, l - 1);
    
    var size = +csv[i][l - 1];
    if (isNaN(size)) { // e.g. if this is a header row
      continue;
    }
    //var parts = sequence.split("-");
    var currentNode = root;
    for (var j = 0; j < parts.length; j++) {
      var children = currentNode["children"];
      var nodeName = parts[j];
      var childNode;
      if (j + 1 < parts.length) {
   // Not yet at the end of the sequence; move down the tree.
 	var foundChild = false;
 	for (var k = 0; k < children.length; k++) {
 	  if (children[k]["name"] == nodeName) {
 	    childNode = children[k];
 	    foundChild = true;
 	    break;
 	  }
 	}
  // If we don't already have a child node for this branch, create it.
 	if (!foundChild) {
 	  childNode = {"name": nodeName, "children": []};
 	  children.push(childNode);
 	}
 	currentNode = childNode;
      } else {
 	// Reached the end of the sequence; create a leaf node.
 	childNode = {"name": nodeName, "size": size};
 	children.push(childNode);
      }
    }
  }
  return root;
};
