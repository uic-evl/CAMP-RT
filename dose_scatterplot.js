//function to move svg elements  to the front
//from https://github.com/wbkd/d3-extended
d3.selection.prototype.moveToFront = function() {  
  return this.each(function(){
	this.parentNode.appendChild(this);
  });
};

function DoseScatterPlot(data){
	this.data = data;
	this.selectedColor = '#e41a1c';
}

DoseScatterPlot.prototype.draw = function(target, selectedPerson = null){
	div = document.getElementById(target);
	this.width = div.clientWidth;
	this.height = .5*this.width;
	this.xMargin = .04*this.width;
	this.yMargin = .03*this.width;
	d3.select("#"+target).selectAll('.scatterSvg').remove();
	this.svg = d3.select("#"+target).insert('svg',':first-child')
					.attr('class', 'scatterSvg')
					.attr('width', this.width)
					.attr('height', this.height)
					.append('g');
	this.tooltip = d3.select("div.tooltip")
			.attr('class','tooltip')
			.style('visibility','hidden');
	this.getColor = this.getColorMap();
	this.drawCircles();
	if (selectedPatient != null){
		this.highlightSelectedPatients(selectedPerson);
	}
	this.setupTooltip();
	this.setupSwitchButtons();
}

DoseScatterPlot.prototype.setupSwitchButtons = function(){
	var setSelected = function(element){
		element.style.opacity = 1;
		element.style.background = 'hsl(240, 90%, 50%)';
	}
	var setUnselected = function(element){
		element.style.opacity = .4;
		element.style.background = 'hsl(240, 50%, 60%)';
	}
	var doseButton = document.getElementById('doseScatterButton');
	var distanceButton = document.getElementById('distanceScatterButton');
	setSelected(doseButton);
	setUnselected(distanceButton);
	var self = this;
	doseButton.addEventListener('click', function(){
		self.switchAxisVariable('dose');
		setSelected(doseButton);
		setUnselected(distanceButton);
	});
	distanceButton.addEventListener('click', function(){
		self.switchAxisVariable('distance');
		setSelected(distanceButton);
		setUnselected(doseButton);
	});
}

DoseScatterPlot.prototype.getAxisScales = function(getAxis){
	var xDomain = d3.extent(this.data, function(d){return getAxis(d)[0];});
	var yDomain = d3.extent(this.data, function(d){return getAxis(d)[1];})
	
	var xScale = d3.scaleLinear()
		.domain(xDomain)
		.range([this.xMargin, this.width - this.xMargin]);
	var yScale = d3.scaleLinear()
		.domain(yDomain)
		.range([this.height - this.yMargin, this.yMargin]);
	return([xScale, yScale])
}

DoseScatterPlot.prototype.getFeatureScales = function(){
	var sizeScale = d3.scalePow().exponent(2)
		.domain( d3.extent(this.data, function(d){ return d.mean_error; }) )
		.range([3, 8]);
	return sizeScale;
}

DoseScatterPlot.prototype.getColorMap = function(){
	//get a cluster value based color map.  assumes clusters index 1-n_clusters
	//fir
	var classColors = ['#ffff33', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628', '#f781bf', '#999999', '#eeeeee'];
	var errorPrinted = false;
	//scale luminosity based on best score?
	var luminosityScale = d3.scaleLinear()
		.domain(d3.extent(this.data, function(d){ return d.scores_ssim[1]; }))
		.range([.3, .8]);
	var colorMap = function(dataPoint){
		var cluster = dataPoint.cluster;
		//print error if more classes than colors, defaults excess to first color
		if( cluster > classColors.length && !errorPrinted){
			console.log('number of classes out of range we have colors for');
			errorPrinted = true;
			cluster = 1
		}
		var color = d3.hsl(classColors[cluster - 1]);
		color.l = luminosityScale(dataPoint.scores_ssim[1]);
		return color;
	}
	return colorMap;
}

DoseScatterPlot.prototype.drawCircles = function(){
	var getAxis = function(x){ return x.dose_pca; };
	var [xScale, yScale] = this.getAxisScales(getAxis);
	var sizeScale = this.getFeatureScales();
	var getColor = this.getColor;
	this.circles = this.svg.selectAll('.dot')
		.data(this.data).enter()
		.append('circle')
		.attr('class', 'dot')
		.attr('id', function(d){ return 'scatterDot' + d.ID;})
		.attr('r', function(d){
			return sizeScale(d.mean_error)})
		.attr('cx', function(d) {
			return xScale(getAxis(d)[0]);})
		.attr('cy', function(d) {
			return yScale(getAxis(d)[1]);});
	this.circles
		.attr('fill', function(d){ 
			return getColor(d);})
		.attr('stroke', 'black')
		.attr('opacity', .5)
		.attr('stroke-width', 1)
		.on('click', function(d){
			switchPatient(d.ID_internal);//from camprt.js
		});
}

DoseScatterPlot.prototype.switchAxisVariable = function(type){
	var getAxis;
	if(type == 'distance'){
		var getAxis = function(d){ return d.distance_pca; };
	} else{
		var getAxis = function(d){ return d.dose_pca; };
	}
	var [xScale, yScale] = this.getAxisScales(getAxis);
	console.log(xScale);
	this.circles.transition().duration(800)
		.attr('cx', function(d) {
			return xScale(getAxis(d)[0]);})
		.attr('cy', function(d) {
			return yScale(getAxis(d)[1]);}); 
}

DoseScatterPlot.prototype.highlightSelectedPatients = function(selectedPerson){
	//recolor everything first
	var getColor = this.getColor;
	this.circles
		.attr('fill', function(d){ 
			return getColor(d);})
		.attr('stroke', 'black')
		.attr('stroke-width', 1)
		.attr('opacity', .5);
	//get entries linked to the person
	var selectedData = this.data[selectedPerson - 1];
	this.subsetData = new Array();
	selectedData.similarity_ssim.forEach(function(x){
		this.subsetData.push(this.data[x - 1]);
	}, this);
	//recolor people matched with selected patient, make them darker and colorfuler
	this.subsetData.forEach(function(x){
		var dataPoint = d3.select('#scatterDot' + x.ID)
			.attr('opacity', 1)
			.attr('stroke-width', 2)
			.moveToFront();
		console.log(dataPoint);
	});
	//make main patient red and stuff
	d3.select('#scatterDot' + selectedData.ID)
		.attr('opacity', 1)
		.attr('stroke-width', 2)
		.attr('fill', 'FireBrick')
		.moveToFront();
	console.log(this.subsetData);
}

DoseScatterPlot.prototype.setupTooltip = function(d){
	var tooltip = this.tooltip;
	this.circles.on('mouseover', function(d){
		tooltip.html(d.name + '</br>' 
			+ 'Error: ' + Math.round(d.mean_error*100)/100 + '</br>')
			.style('left', d3.event.pageX + 10 + 'px')
			.style('top', d3.event.pageY - 30 + 'px');
		tooltip.transition().duration(50).style('visibility','visible');
	}).on('mouseout', function(d){
		tooltip.transition().duration(50).style('visibility', 'hidden');
	});
}