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
		this.xMargin = .04*this.width;
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
			tickSpaces.push(i);
		}
		this.svg.selectAll('.axisLine')
			.data(tickSpaces).enter()
			.append('line')
			.attr('class', 'axisLine')
			.attr('id', function(d) {return organList[d] + 'axisLine';})
			.attr('x1', function(d){return xScale(d);})
			.attr('x2', function(d){return xScale(d);})
			.attr('y2', self.yMargin)
			.attr('y1', self.height - .5*self.yMargin- self.xAxisSize)
			.attr('stroke', 'silver')
			.attr('stroke-width', .05*self.binWidth)
		this.svg.selectAll('.bins')
			.data(tickSpaces).enter()
			.append('rect')
			.attr('class', '.bins')
			.attr('x', function(d){ return xScale(d) - self.binWidth/2;})
			.attr('y', self.yMargin)
			.attr('width', self.binWidth)
			.attr('height', self.height - self.yMargin - self.xAxisSize)
			.attr('opacity', 0)
			.on('mouseover', function(d){
				var axisLine = d3.select( "#" + organList[d] + 'axisLine' )
				axisLine.attr('stroke', 'white')
					.attr('stroke-width', .08*self.binWidth);
				var organ = self.patient.organData[self.organList[d]]
				self.tooltip.html(self.organList[d] + '</br>'
				+ 'Predicted: ' + organ.estimatedDose + ' Gy </br>'
				+ 'Actual: ' + organ.meanDose + ' Gy')
				.style('left', d3.event.pageX + .5*self.binWidth + 'px')
				.style('top', d3.event.pageY - 50 + 'px');
				self.tooltip.transition().duration(50).style('visibility','visible');
				d3.select(this).on('mousemove', function(d){
					self.tooltip.style('left', d3.event.pageX + .5*self.binWidth + 'px')
						.style('top', d3.event.pageY - 50 + 'px');
				});
			}).on('mouseout', function(d){
				var axisLine = d3.select( "#" + organList[d] + 'axisLine' )
				axisLine.attr('stroke', 'silver')
					.attr('stroke-width', .05*self.binWidth);
				self.tooltip.transition().duration(50).style('visibility', 'hidden');
				d3.select(this).on('mousemove', null);
			});;
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
			.attr('fill', 'hsl(60, 60%, 50%)')
			.attr('stroke', 'black')
			.attr('stroke-width', .1)
			.on('mouseover', function(d,i){
				tooltip.html(d.name + '</br>'
				+ 'dose: ' + d.organData[organName].meanDose + ' Gy')
				.style('left', d3.event.pageX + 10 + 'px')
				.style('top', d3.event.pageY - 30 + 'px');
				tooltip.transition().duration(50).style('visibility','visible');
				d3.select(this)
					.moveToFront()
					.transition().duration(50)
					.attr('opacity', 1)
					.attr('stroke-width', .3)
					.attr('stroke', 'white');
			}).on('mouseout', function(d){
				tooltip.transition().duration(50).style('visibility', 'hidden');
				d3.select(this).transition().duration(50)
					.attr('opacity', .5)
					.attr('stroke', 'black')
					.attr('stroke-width', .1);
				d3.selectAll('.meanDoseRect').moveToFront();
				d3.selectAll('.estimatedDoseRect').moveToFront();
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