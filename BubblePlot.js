var OrganBubblePlot = (function(){
	"use strict";
	var instance;
	var num_organs;
	var self = this;
	function Graph(){
		if (instance){
			return instance;
		}
		instance = this
	}
	Graph.getInstance = function(){
		return instance || new Graph();
	}
	Graph.init = function(target, patientInternalId, patientData){
		console.log('drawing');
		div = document.getElementById(target);
		this.data = patientData;
		this.setPatient(patientInternalId, patientData);
		
		this.width = div.clientWidth;
		this.height = div.clientHeight
		this.xMargin = .025*this.width;
		this.yMargin = .1*this.height;
		this.binWidth = (this.width - 2*this.xMargin)/num_organs;
		this.xAxisSize = 80;
		
		d3.select("#"+target).selectAll('.bubbleSvg').remove();
		this.svg = d3.select("#"+target).append('svg')
						.attr('class', 'bubbleSvg')
						.attr('width', this.width)
						.attr('height', this.height)
						.append('g');
		this.tooltip = d3.select("div.bubbletooltip")
				.attr('class','tooltip')
				.style('visibility','hidden');
		console.log(this.organList);
		this.drawAxes();
		this.drawOrgans();
	}
	
	Graph.switchPatient = function(patientInternalId){
		this.setPatient(patientInternalId);
		this.drawOrgans();
	}
	
	Graph.setPatient = function(patientInternalId){
		this.patient = this.data[patientInternalId - 1];
		
		var organList = Object.keys(this.patient.organData)
		num_organs = organList.length;
		if('GTVp' in this.patient.organData){ num_organs = num_organs - 1; }
		if('GTVn' in this.patient.organData){ num_organs = num_organs - 1; }
		this.organList = organList.slice(0, num_organs);
		
		this.patientMatches = [];
		for(var i = 1; i < this.patient.similarity_ssim.length; i++){
			let match = this.patient.similarity_ssim[i];
			this.patientMatches.push( this.data[match - 1] );
		}
	}
	
	Graph.drawOrgans = function(){
		var [xAxis, yAxis] = this.getScales();
		var rAxis = this.getRScale(this.patient, this.binWidth);
		this.organList.forEach(function(d){
			this.drawOrganBubbles(d, xAxis, yAxis, rAxis);
		},this);
		this.drawDoseRects(xAxis, yAxis);	
	}
	
	Graph.drawAxes = function(){
		var[xScale, yScale] = this.getScales();
		var organList = this.organList;
		var self = this;
		var xAxis = d3.axisBottom(xScale)
			.ticks(this.organList.length)
			.tickFormat(function(d){
				return organList[d];
			});
		var yAxis = d3.axisLeft(yScale).ticks(5);
		var translate = "translate( 0," + (this.height - .5*this.yMargin- this.xAxisSize).toString() + ")";
		this.svg.append('g')
			.attr('class', 'axis')
			.attr('transform', translate)
			.call(xAxis)
			.selectAll('text')
			.attr('transform', 'rotate(-45)')
			.style('text-anchor', 'end')
		this.svg.append('g')
			.attr('class', 'axis')
			.attr('transform', 'translate(' + .8*this.xMargin + ', ' + .5*this.yMargin + ')' )
			.call(yAxis);
			
		var tickSpaces = [];
		for(var i = 0; i < organList.length; i++){
			tickSpaces.push(xScale(i));
		}
		this.svg.selectAll('.axisLine')
			.data(tickSpaces).enter()
			.append('line')
			.attr('class', 'axisLine')
			.attr('x1', function(d){return d;})
			.attr('x2', function(d){return d;})
			.attr('y2', 1.5*self.yMargin)
			.attr('y1', self.height - .5*self.yMargin- self.xAxisSize)
			.attr('stroke', 'silver')
			.attr('stroke-width', .05*self.binWidth);
	}
	
	Graph.drawDoseRects = function(xScale, yScale){
		var patient = this.patient;
		var width = .8*this.binWidth;
		var height = .15*this.binWidth;
		var self = this;
		var drawDose = function(variable, color){
			d3.selectAll('.' + variable + 'Rect').remove();
			self.svg.selectAll('.'+ variable + 'Rect')
				.data(self.organList).enter()
				.append('rect')
				.attr('class', 'meanDoseRect')
				.attr('x', function(o,index){ return xScale(index) - width/2;})
				.attr('y', function(o) { return yScale(patient.organData[o][variable]) - height/2; })
				.attr('fill', color)
				.attr('width', width)
				.attr('height', height)
				.attr('opacity', .8)
				.attr('stroke-width', .5)
				.attr('stroke', 'black');
		}
		drawDose('meanDose', '#8b0000');
		drawDose('estimatedDose', '#3d32ff');
	}
	
	Graph.drawOrganBubbles = function(organName, xScale, yScale, rScale){
		var position = this.organList.indexOf(organName);
		var x = xScale(position);
		var patient = this.patient;
		var tooltip = this.tooltip;
		var getRadius = function(x){
			var index = patient.similarity_ssim.indexOf(+x.ID_internal);
			if(index == 1){ console.log("bad stuff");}
			if(index != -1){
				return rScale(patient.scores_ssim[index]);
			} else{
				return 0;
				console.log('error in getRadius?');
			}
		}
		this.svg.selectAll('.organDot').filter('.' + organName).remove();
		var bubbles = this.svg.selectAll('.organDot').filter('.' + organName)
			.data(this.patientMatches).enter()
			.append('circle')
			.attr('class', 'organDot ' + organName)
			.attr('cx', x)
			.attr('cy', function(d){ return yScale(d.organData[organName].meanDose);})
			.attr('r', getRadius)
			.attr('opacity', .5)
			.attr('fill', '#acae56')
			.attr('stroke', 'black')
			.attr('stroke-width', .1)
			.on('mouseover', function(d,i){
				tooltip.html(d.name + '</br>'
				+ 'dose: ' + d.organData[organName].meanDose + ' Gy')
				.style('left', d3.event.pageX + 10 + 'px')
				.style('top', d3.event.pageY - 30 + 'px');
				tooltip.transition().duration(50).style('visibility','visible');
			}).on('mouseout', function(d){
				tooltip.transition().duration(50).style('visibility', 'hidden');
			});
	}
	
	Graph.getRScale = function(patient, width){
		var ratio = .4;
		var scoreExtent = d3.extent( patient.scores_ssim.slice(1), function(d) {return d;});
		console.log([ratio*width, width*(1-ratio)]);
		var rScale = d3.scalePow().exponent(2)
			.domain(scoreExtent)
			.range([ratio*width/2, width*(1-ratio)/2]);
		return(rScale)
	}
	
	Graph.getScales = function(){
		var xScale = d3.scaleLinear()
			.domain([0, this.organList.length])
			.range([this.binWidth/2 + this.xMargin, this.width - this.xMargin - this.binWidth/2]);
		//convoluted way of getting ths max and min dose for everything
		var allPatients = [...this.patientMatches];
		allPatients.push(this.patient);
		let maxDose = 0;
		let minDose = 0;
		allPatients.forEach(function(d){
			this.organList.forEach(function(x){
				var newDose = d.organData[x].meanDose;
				maxDose = (maxDose > newDose)? maxDose: newDose;
				minDose = (minDose < newDose)? minDose: newDose;
			});
		}, this);
		var yScale = d3.scaleLinear()
			.domain( [minDose, maxDose])
			.range([this.height - this.yMargin - this.xAxisSize, this.yMargin]);
		return [xScale, yScale];
	}
	return Graph;
}());