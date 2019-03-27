var Controller = (function(){
	var scatterData = {};
	var bubbleData = {};
	var enableBrush = false;
	var doseColor = d3.scaleLinear()
		.domain([0,70])
		.range(['#ffffe0','#8b0000']);
	var doseErrorColor = d3.scaleLinear()
		.domain([0, 20])
		.range(['#999999','#3d32ff']);

	return {
		
		brushOrgan: function(organ){
			try{
				var axisLine = d3.select( "#" + organ + 'axisLine' );
				axisLine.attr('stroke', 'white')
					.attr('stroke-width', .1*OrganBubblePlot.binWidth);
				scenes.forEach(function(scene){
					var node = scene.getObjectByName(organ);
					var model = scene.getObjectByName(organ + '_model');
					//model.material.color.set('white');
					model.material.opacity = 1;
				});
			} catch{}
		},
		
		unbrushOrgan: function(organ){
			try{
				var axisLine = d3.select( "#" + organ + 'axisLine' )
				var currentWidth = axisLine.attr('stroke-width');
				axisLine.attr('stroke', 'silver')
					.attr('stroke-width', .05*OrganBubblePlot.binWidth);
				scenes.forEach(function(scene){
					var node = scene.getObjectByName(organ);
					var model = scene.getObjectByName(organ + '_model');
					var organColor = doseColor(node.userData.meanDose);
					//model.material.color.set(organColor);
					var currentOpacity = document.getElementById("opacSlider").value/100.0;
					model.material.opacity = currentOpacity;
				});
			} catch{}
		},
		
		switchScene: function(scene, type, data){
			var id = scene.userData.element.parentElement.id;
			scene.children.forEach(function(d){
				if(d.userData.type == 'node' && d.name != 'GTVp' && d.name != 'GTVn'){
					var organName = d.name;
					if(type.toLowerCase() == 'predict'){
						var dose = data.getEstimatedDose(id, organName);
						var color = doseColor(dose);
					} else if(type.toLowerCase() == 'error'){
						var dose = data.getEstimationError(id, organName);
						var color = doseErrorColor(dose);
					} else{
						var dose = data.getMeanDose(id, organName);
						var color = doseColor(dose);
					}
					d.material.color.set(color);

					d.userData.meanDose = dose.toFixed(3);
					d.userData.dosePerVolume = (dose/data.getOrganVolume(id, organName)).toFixed(3);
					if(type.toLowerCase() == 'actual'){
						d.userData.maxDose = data.getMaxDose(id, organName).toFixed(3);
						d.userData.minDose = data.getMinDose(id, organName).toFixed(3);
					} else{
						d.userData.maxDose = '';
						d.userData.minDose = '';
					}
					var model = scene.getObjectByName(organName + '_model');
					if(model != undefined){
						model.material.color.set(color);
					}
				}
			});
		},
		
		getDoseColor: function(d){ return doseColor(d); },
		
		getDoseErrorColor: function(d){ return doseErrorColor(d); },
		
		toggleBrush: function(isEnabled){
			enableBrush = isEnabled;
		},
		
		brush: function(id){
			if(enableBrush){
				//highlight patient matches in other views on mousover
				var patientView = document.getElementById(id);
				if(patientView != null){
					var description = patientView.querySelector('.description');
					description.style.color = 'white';;
				}
				
				var bubble = d3.selectAll('#organBubble'+id);
				if(!bubble.empty()){
					bubbleData.fill = bubble.attr('fill');
					bubbleData.stroke = bubble.attr('stroke');
					bubbleData.opacity = bubble.attr('opacity');
					bubble.attr('fill', 'white')
						.attr('stroke', 'white')
						.attr('opacity', 1)
						.moveToFront();
				}
				var scatterSelection =  d3.selectAll('path').filter('#scatterDot' + id);
				if(!scatterSelection.empty()){
					scatterData.stroke = scatterSelection.attr('stroke');
					scatterSelection.attr('stroke', 'white')
						.moveToFront();
				}
			}
		},
		
		unbrush: function(id){
			if(enableBrush){
				var patientView = document.getElementById(id);
				if(patientView != null){
					var description = patientView.querySelector('.description');
					description.style.color = 'black';
				}
				if(!d3.selectAll('#organBubble'+id).empty()){
					d3.selectAll('#organBubble'+id)
						.attr('fill', bubbleData.fill)
						.attr('stroke', bubbleData.stroke)
						.attr('opacity', bubbleData.opacity);
				}
				if(!d3.selectAll('path').filter('#scatterDot' + id).empty()){
					d3.selectAll('path').filter('#scatterDot' + id)
						.attr('stroke', scatterData.stroke);
				}
				if( !d3.selectAll('.doseRect').empty()){
					d3.selectAll('.doseRect').moveToFront();
				}
			}
		},
		 
		setup: function(){
			var self = this;
			d3.selectAll('.description').on('mouseover', function(d){
				var id = this.parentNode.id;
				self.brush(id);
			}).on('mouseout', function(d){
				var id = this.parentNode.id;
				self.unbrush(id);
			}).on('click', function(d){
				var id = this.parentNode.id;
				switchPatient(id);
			});
		}
	}
})();