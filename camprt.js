'use strict';

if (!Detector.webgl) {
    Detector.addGetWebGLMessage();
}
//storing everything in like 100 global variables.  The best programming practice.
var parent = document.getElementById("content");
var nodeDetails = document.getElementById("details");

var detailsOffsetX = 15;
var detailsOffsetY = 15;

var organName = document.getElementById("details_organName"),
    dosePerVolume = document.getElementById("details_dosePerVolume"),
    lineSeparator = document.getElementById("details_line"),
    volumeVal = document.getElementById("details_volume_val"),
    meanDoseVal = document.getElementById("details_meanDose_val"),
    minDoseVal = document.getElementById("details_minDose_val"),
    maxDoseVal = document.getElementById("details_maxDose_val"),
    pDisplayed = document.getElementById("pDisplayed");


var scenes = [],
    renderer;

	
var selectedPatient = 1;
//patients shown on load screen?
var patientsToShow = 11;

var totalModelCount;

var syncCameras = true,
    syncCamerasInterval,
    detailsOnRotate = true;

var raycaster;

var mouse = new THREE.Vector2(-500, -500);

var mouseNorm = new THREE.Vector2(-500, -500),
    INTERSECTED = null,
    nodeHover = null;

var cameraDistZ = 500;
// data
var organs, oAtlas, links;
var organModels = new THREE.Group();

var partitions = ["Oral Cavity & Jaw", "Throat", "Salivary Glands", "Eyes", "Brainstem & Spinal Cord", "Other"];

var organRef = [];

var currScene;

var master = document.getElementById("masterList");

var materialArray;

var canvas = document.getElementById("c");
var template = document.getElementById("template").text;

var manager = new THREE.LoadingManager();

manager.onStart = function(url, itemsLoaded, itemsTotal){
	document.getElementById("loadScreen").style.display = "block";
}

manager.onLoad = function () {
	//this may break this because I moved it to the front
	initializeRiskPrediction(selectedPatient);
    document.getElementById("loadScreen").style.display = "none";
	Controller.toggleBrush(true);
};

manager.onProgress = function (url, itemsLoaded, itemsTotal) {
    document.getElementById("loadProgress").innerHTML = parseInt(itemsLoaded / itemsTotal * 100) + " %"
};

var scatter;
var bubbleChart;

var files = ["data/organAtlas.json", "PYTHON/data/patient_dataset_v23.json"];
var promises = [];
var data; 

files.forEach(function (url) {
    promises.push(d3.json(url));
});

Promise.all(promises).then(function (values) {
    start(values[0], values[1]);
});

function start(organAtlas, patientsData) {
    oAtlas = organAtlas[0];
	data = Data(patientsData, oAtlas);
    selectedPatient = populateDropDownMenu();

    init(); // initialize

    populateOrganMasterList();

    currScene = scenes[0];

    document.addEventListener("mousedown", onMouseDown, false);
    document.addEventListener("mouseup", onMouseUp, false);

    document.addEventListener("touchstart", onTouchStart, false);
    document.addEventListener("touchend", onTouchEnd, false);

    document.addEventListener("mousemove", onDocumentMouseMove, false);

    animate(); // render
	scatter = new DoseScatterPlot(data); //ok If I don't make this a global I have to like, rewrite half this code
	scatter.draw('organErrorViz', selectedPatient);
	OrganBubblePlot.init('bubbleChart', selectedPatient, data);
	Controller.setup();
	ColorScale.draw();
	window.addEventListener('resize', function(d){
		scatter.draw('organErrorViz', selectedPatient);
		OrganBubblePlot.init('bubbleChart', selectedPatient, data);
		Controller.setup();
	});
}

// ----------------------------------------------------------------

function populateDropDownMenu() {
	//holds an array of patient internal ids
    var menu = document.getElementById("patientMenu");
    // copy of patients sorted
    var patients_sorted = data.getSortedPatients();
    patients_sorted.forEach(function (patient, index) {
        var tempOption = document.createElement("option");
        tempOption.value = patient.ID_internal;
        tempOption.innerHTML = patient.name;
        menu.appendChild(tempOption);
    });

    // first patient 
    var firstPatient = patients_sorted[0].ID_internal;
    //THis appears to look at the url to see if a patient is selected there
	//if so, sets "first patient" to this guy? otherwise, uses the lowest id (above)
    var patientURL = getQueryVariable("id");

    if (patientURL != false) {
        var convertedPatient = getInternalID(patientURL);
        if (convertedPatient != false) {
            firstPatient = convertedPatient;
            menu.value = convertedPatient;
        }
    }
    patients_sorted.length = 0;
    return firstPatient;
}

function getInternalID(searchString) {
	var patient;
    for (var i = 1; i <= data.PatientCount; i++) {
		patient = data.getPatient(i);
        if (patient.ID == searchString)
            return patient.ID_internal;

    }

    return false;
}

function getQueryVariable(variable) {
    var query = window.location.search.substring(1);
    var vars = query.split("&");
    for (var i = 0; i < vars.length; i++) {
        var pair = vars[i].split("=");
        if (pair[0] == variable) {
            return pair[1];
        }
    }
    return (false);
}


function handleCheckBoxSingle(event) {
    if (event.checked) {

        scenes.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);
            var model = scene.getObjectByName(String(event.value) + "_model");
            if (node && model) {
                node.visible = true;
                model.visible = true;
            }
        });

        //d3.select("#line_" + event.value).style("opacity", 1.0);
        d3.select("#line_" + event.value).attr("display", null);

    } else {

        scenes.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);
            var model = scene.getObjectByName(String(event.value) + "_model");

            if (node && model) {
                node.visible = false;
                model.visible = false;
            }
        });

        d3.select("#line_" + event.value).attr("display", "none");
    }
}

function handleCheckBoxGroup(event) {
    var children = master.getElementsByClassName(event.id[0] + "_GroupChildren");

    if (event.checked) {

        for (var i = 0; i < children.length; i++) {

            children[i].checked = true;

            d3.select("#line_" + children[i].value).attr("display", null);
        }

        scenes.forEach(function (scene, index) {

            for (var i = 0; i < children.length; i++) {

                var node = scene.getObjectByName(children[i].value);
                var model = scene.getObjectByName(String(children[i].value) + "_model");

                if (node && model) {
                    node.visible = true;
                    model.visible = true;
                }

            }
        });

    } else {

        for (var i = 0; i < children.length; i++) {

            children[i].checked = false;
            
            d3.select("#line_" + children[i].value).attr("display", "none");
        }

        scenes.forEach(function (scene, index) {

            for (var i = 0; i < children.length; i++) {

                var node = scene.getObjectByName(children[i].value);
                var model = scene.getObjectByName(String(children[i].value) + "_model");
				
                if (node && model) {
                    node.visible = false;
                    model.visible = false;
                }

            }
        });
    }


}

function populateOrganMasterList() {
    // make partition input first
    partitions.forEach(function (group, i) {

        var tempDiv = document.createElement("div");

        tempDiv.setAttribute("class", "checkbox_group");
        tempDiv.setAttribute("id", String(i + 1) + "_group_container");

        var tempInput = document.createElement("INPUT");

        tempInput.setAttribute("type", "checkbox");
        tempInput.setAttribute("id", String(i + 1) + "_title");
        tempInput.setAttribute("value", group);
 
        tempInput.setAttribute("onchange", "handleCheckBoxGroup(this)");

        tempInput.setAttribute("checked", true);

        var tempLabel = document.createElement("label");
        tempLabel.setAttribute("for", String(i + 1) + "_title");

        tempLabel.style.fontSize = "16px";
        tempLabel.style.fontWeight = "bold";
        tempLabel.style.paddingLeft = "15px";
        tempLabel.innerHTML = group;
        // ----------
        tempDiv.appendChild(tempInput);
        tempDiv.appendChild(tempLabel);
        // ----------
        var tempDiv2 = document.createElement("div");

        tempDiv2.setAttribute("id", String(i + 1) + "_single_container");

        master.appendChild(tempDiv);
        master.appendChild(tempDiv2);

        var tempDiv3 = document.createElement("div");

        tempDiv3.setAttribute("class", "dummy");

        master.appendChild(tempDiv3);

    });

	////unhelpful documentation:
    // for loop bad, iterates in an unspecified order!!!!!!!!!!!!!!!!!!!!!!

    // individual organs
    for (var organ in oAtlas) {

        if (organ != "GTVn" && organ != "GTVp") {

            var parent = document.getElementById(String(oAtlas[organ].partition) + "_single_container");

            var tempDiv = document.createElement("div");
            tempDiv.style.marginRight = "20px";
            tempDiv.style.marginLeft = "50px";
            tempDiv.style.paddingLeft = "10px";
            tempDiv.style.paddingTop = "1px";
            tempDiv.style.paddingBottom = "1px";
            tempDiv.setAttribute("class", "GroupChildren");

            var tempInput = document.createElement("INPUT");

            tempInput.setAttribute("type", "checkbox");
            tempInput.setAttribute("id", organ + "_checkList");
            tempInput.setAttribute("class", String(oAtlas[organ].partition) + "_GroupChildren");
            tempInput.setAttribute("value", organ);
            tempInput.setAttribute("onchange", "handleCheckBoxSingle(this)");

            tempInput.setAttribute("checked", true);

            var tempLabel = document.createElement("label");
            tempLabel.setAttribute("for", organ + "_checkList");
            tempLabel.style.paddingLeft = "15px";

            tempLabel.innerHTML = organ;

            // ----------

            tempDiv.appendChild(tempInput);
            tempDiv.appendChild(tempLabel);

            parent.appendChild(tempDiv);

        }

    }
}

function checkOrganMasterList() {
    organs.forEach(function (organ, index) {

        var tempItem = document.getElementById(organ.name);

        if (scenes[selectedPatient - 1].getObjectByName(organ.name).userData.dosePerVolume != null) {

            if (tempItem.checked != true)
                tempItem.setAttribute("checked", true);

        } else {

            if (tempItem.checked != false)
                tempItem.setAttribute("checked", false);
        }
    });

}

function formatOrganMasterList() {
    var organList = master.children;

    for (var i = 0; i < organList.length; i++) {

        if (i % 2 == 0)
            organList[i].style["backgroundColor"] = "#3a3a3a";
        else
            organList[i].style["backgroundColor"] = "#444444";

    }
}

function init() {
	//renderer for main views?
	var getRenderer = function(canvas, isAlpha){
		var r = new THREE.WebGLRenderer({
			canvas: canvas,
			antialias: true,
			alpha: isAlpha
		});
		r.setClearColor(0x888888, 1);
		r.setPixelRatio(window.devicePixelRatio);
		r.sortObjects = true;
		return r
	}
    
	renderer = getRenderer(canvas, false);

    raycaster = new THREE.Raycaster();

    var maxAnisotropy = renderer.getMaxAnisotropy();

    var textureLoader = new THREE.TextureLoader();

    var texture0 = textureLoader.load('resources/anterior.png'), // xpos, Right
        texture1 = textureLoader.load('resources/posterior.png'), // xneg, Left
        texture2 = textureLoader.load('resources/superior.png'), // ypos, Top
        texture3 = textureLoader.load('resources/inferior.png'), // yneg, Bottom
        texture4 = textureLoader.load('resources/right.png'), // zpos, Back
        texture5 = textureLoader.load('resources/left.png'); // zneg, Front

    texture0.anisotropy = maxAnisotropy;
    texture1.anisotropy = maxAnisotropy;
    texture2.anisotropy = maxAnisotropy;
    texture3.anisotropy = maxAnisotropy;
    texture4.anisotropy = maxAnisotropy;
    texture5.anisotropy = maxAnisotropy;
    //getcontext2d. draw image
    materialArray = [
            new THREE.MeshBasicMaterial({
            map: texture0
        }),
            new THREE.MeshBasicMaterial({
            map: texture1
        }),
            new THREE.MeshBasicMaterial({
            map: texture2
        }),
            new THREE.MeshBasicMaterial({
            map: texture3
        }),
            new THREE.MeshBasicMaterial({
            map: texture4
        }),
            new THREE.MeshBasicMaterial({
            map: texture5
        })
    ];
    
	var patientObject = data.getPatient(selectedPatient);
	scenes = updateScenes(selectedPatient, materialArray);//populates required views	
	updateOrder(selectedPatient);
}

function updateScenes(selectedPatient, material){
	var scenes = []; //scenes is a wonderful global for now
	var matches = data.getPatientMatches(selectedPatient);
	for (var i = 0; i < patientsToShow && i < matches.length; i++) {
		var id = matches[i];
		var target = (i == 0)? "leftContent" : "content";
		var newScene = showPatient(material, id, target);
		scenes.push(newScene);
	}
	return scenes
}

function placeOrganModels(pOrgan, organProperties, scene, nodeColor) {
    let loader = new THREE.VTKLoader(manager);

    if (!(pOrgan == "GTVn" || pOrgan == "GTVp")) {

        loader.load('resources/models/' + pOrgan + '.vtk', function (geometry) {

            geometry.computeVertexNormals();
            geometry.center();

            let material = new THREE.MeshBasicMaterial({
                color: nodeColor,
                opacity: 0.2,
                transparent: true,
                depthTest: true,
                depthWrite: true,
                depthFunc: THREE.LessEqualDepth
            });

            let mesh = new THREE.Mesh(geometry, material);
            mesh.name = (String(pOrgan) + "_model");
            mesh.userData.type = "node_model";

            mesh.position.x = organProperties.x;
            mesh.position.y = organProperties.y;
            mesh.position.z = organProperties.z;

            mesh.rotation.x = -Math.PI / 2.0;
            mesh.rotation.z = -Math.PI / 2;


            // oral cavity
            if (pOrgan == "Tongue")
                mesh.renderOrder = -10;
            else if (pOrgan == "Genioglossus_M")
                mesh.renderOrder = -10;
            else if (pOrgan == "Lt_Ant_Digastric_M")
                mesh.renderOrder = -10;
            else if (pOrgan == "Mylogeniohyoid_M")
                mesh.renderOrder = -10;
            else if (pOrgan == "Rt_Ant_Digastric_M")
                mesh.renderOrder = -10;

            else if (pOrgan == "Extended_Oral_Cavity")
                mesh.renderOrder = -9;

            // throat
            else if (pOrgan == "Larynx")
                mesh.renderOrder = -10;
            else if (pOrgan == "Supraglottic_Larynx")
                mesh.renderOrder = -9;


            scene.add(mesh);
        });
    }
}

function showPatient(materialArray, id, parentDivId){
	var scene = new THREE.Scene();
	var patient = data.getPatient(id);
	var patientOrganList = patient.organData;

	// make a list item
	var element = document.createElement("div");
	element.className = "list-item";
	element.id = id;
	element.innerHTML = template.replace('$', patient.name);
	
	var totDoseElement = element.querySelector(".totDose");
	totDoseElement.innerHTML = "Total Dose: " + "<b>" + patient.total_Dose + "</b>" + " GY";

	var tVolumeElement = element.querySelector(".tVolume");
	tVolumeElement.innerHTML = "GTV: " + "<b>" + patient.tumorVolume + "</b>" + " cc";

	var lateralityElement = element.querySelector(".laterality");
	lateralityElement.innerHTML = "<b>(" + patient.laterality + ")</b> " + " " + patient.tumorSubsite;

	// Look up the element that represents the area
	// we want to render the scene
	scene.userData.element = element.querySelector(".scene");
	
	if(!document.getElementById( element.id )){
		document.getElementById(parentDivId).appendChild(element);
	}

	var scalarVal = 2.4; //4.1

	var camera = new THREE.OrthographicCamera(scene.userData.element.offsetWidth / -scalarVal, 
		scene.userData.element.offsetWidth / scalarVal, 
		scene.userData.element.offsetHeight / scalarVal, 
		scene.userData.element.offsetHeight / -scalarVal, 
		1, 100000);
		
	camera.position.z = cameraDistZ;

	camera.updateProjectionMatrix();
	scene.userData.camera = camera;

	// orientation marker, patient coordinate system
	var MovingCubeMat = new THREE.MultiMaterial(materialArray);
	var MovingCubeGeom = new THREE.CubeGeometry(25, 25, 25, 1, 1, 1, materialArray);
	var MovingCube = new THREE.Mesh(MovingCubeGeom, MovingCubeMat);

	camera.add(MovingCube);
	MovingCube.position.set(121, -121, -250);

	//
	var controls = new THREE.OrbitControls(scene.userData.camera, scene.userData.element);
	controls.minDistance = 2;
	controls.maxDistance = 5000;
	controls.enablePan = false;
	controls.enableZoom = false;

	scene.userData.controls = controls;

	var geometry = new THREE.SphereGeometry(4, 16, 16);

	var material = new THREE.MeshStandardMaterial({
		color: new THREE.Color().setHex(0xa0a0a0),
		roughness: 0.5,
		metalness: 0,
		shading: THREE.FlatShading
	});

	var outlineMaterial = new THREE.MeshBasicMaterial({
		color: 0x3d3d3d,
		side: THREE.BackSide
	});

	var linkMaterial = new THREE.LineBasicMaterial({
		color: 0x3d3d3d,
		opacity: 1,
		linewidth: 3
	});
	
	for (var pOrgan in patientOrganList) {
		//this looks like it draws the organs in each patient?
		
		var organSphere = new THREE.Mesh(geometry, material.clone());

		organSphere.position.x = (patientOrganList[pOrgan].x);
		organSphere.position.y = (patientOrganList[pOrgan].y);
		organSphere.position.z = (patientOrganList[pOrgan].z);

		organSphere.name = pOrgan;
		organSphere.userData.type = "node";

		// outline
		var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

		outlineMesh.name = pOrgan + "_outline";

		if (organSphere.name == "GTVp")
			outlineMesh.scale.multiplyScalar(1.6);
		else if (organSphere.name == "GTVn")
			outlineMesh.scale.multiplyScalar(1.5);
		else
			outlineMesh.scale.multiplyScalar(1.3);

		// color
		var nodeColor;

		organSphere.userData.volume = data.getOrganVolume(id, pOrgan);
		organSphere.userData.minDose = data.getMinDose(id, pOrgan);
		organSphere.userData.meanDose = data.getMeanDose(id, pOrgan);
		organSphere.userData.maxDose = data.getMaxDose(id, pOrgan);
		
		organSphere.userData.estimatedDose = data.getEstimatedDose(id, pOrgan);

		// do this in python script maybe
		//grays are already in joules per kilogram?!?!? I might want to delete this because it's misleading to users
		organSphere.userData.dosePerVolume = (data.getMeanDose(id, pOrgan) / data.getOrganVolume(id, pOrgan)).toFixed(3);

		if (organSphere.userData.meanDose >= 0.0) //null == -1 in json, pearson problems
			nodeColor = Controller.getDoseColor(organSphere.userData.meanDose);
		else {
			nodeColor = '#a0a0a0'; 
			organSphere.userData.meanDose = undefined;
			organSphere.userData.dosePerVolume = undefined;
		}

		organSphere.material.color.setStyle(nodeColor);

		scene.add(organSphere);
		organSphere.add(outlineMesh);

		placeOrganModels(pOrgan, patientOrganList[pOrgan], scene, nodeColor);
	}

	var tmp_geo = new THREE.Geometry();

	var source = scene.getObjectByName("GTVp");
	var target = scene.getObjectByName("GTVn");

	if (source != null && target != null) {
		//draws a line between the gtvp and gtvn is they are there
		tmp_geo.vertices.push(source.position);
		tmp_geo.vertices.push(target.position);

		var line = new THREE.LineSegments(tmp_geo, linkMaterial);
		line.scale.x = line.scale.y = line.scale.z = 1;
		line.originalScale = 1;

		scene.add(line);
	}

	// check for missing data
	for (var organ in oAtlas) {

		if (!patientOrganList.hasOwnProperty(organ)) {

			// node
			var organSphere = new THREE.Mesh(geometry, material.clone());

			organSphere.position.x = (oAtlas[organ].x);
			organSphere.position.y = (oAtlas[organ].y);
			organSphere.position.z = (oAtlas[organ].z);

			organSphere.name = organ;
			organSphere.userData.type = "node";

			// outline
			var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

			outlineMesh.name = organ + "_outline";

			if (organSphere.name == "GTVp")
				outlineMesh.scale.multiplyScalar(1.6);
			else if (organSphere.name == "GTVn")
				outlineMesh.scale.multiplyScalar(1.5);
			else
				outlineMesh.scale.multiplyScalar(1.3);

			// color
			var nodeColor = '#a0a0a0';

			organSphere.userData.volume = undefined;
			organSphere.userData.minDose = undefined;
			organSphere.userData.meanDose = undefined;
			organSphere.userData.maxDose = undefined;
			
			organSphere.userData.estimatedDose = undefined;
			
			organSphere.userData.dosePerVolume = undefined;

			organSphere.material.color.setStyle(nodeColor);

			scene.add(organSphere);
			organSphere.add(outlineMesh);

			placeOrganModels(organ, oAtlas[organ], scene, nodeColor);
		}
	}

	scene.add(camera);

	// light
	var light = new THREE.AmbientLight(0xffffff, 1.0); // white light   

	scene.add(light);
	return scene;
}

function removeOldViews(selectedPatientObject){
	//remove list-items not matched to the patient
	var matches = selectedPatientObject.similarity_ssim;
	var patientViews = document.getElementsByClassName('list-item');
	var element;
	for(var i = patientViews.length - 1; i >= 0; i--){
		element = patientViews[i];
		element.parentElement.removeChild(element);
	}
}

function switchPatient(updatedPatient){
	if(updatedPatient == selectedPatient){ 
		return;
	}
	selectedPatient = updatedPatient;
	document.getElementById("patientMenu").value = selectedPatient
	var patientObject = data.getPatient(updatedPatient);
	removeOldViews(patientObject); //removes old views
	scenes = updateScenes(selectedPatient, materialArray);//populates required views
	updateOrder(updatedPatient);
	scatter.highlightSelectedPatients(updatedPatient); 
	OrganBubblePlot.switchPatient(updatedPatient);
	Controller.toggleBrush(false);
	Controller.setup();
}

function formatFirstPatient(updatedPatient){
	var firstPatient = document.getElementById(updatedPatient);
	firstPatient.style.display = "none";
	firstPatient.parentElement.insertBefore(firstPatient, firstPatient.parentElement.childNodes[2] );
	firstPatient.style.zIndex = 1;
	var description = firstPatient.querySelector('.description');
	description.innerHTML += ' &#10010'
		.fontcolor(data.getClusterColor(updatedPatient));//add a colored cross by the selected patients name
	description.style.width = '131px';//320 is whole width
    // first patient always has score of 1, clear it
    firstPatient.querySelector(".pScore").remove();
	var buttonNames = ['Error', 'Predict', 'Actual'];
	buttonNames.forEach(function(name){
		var sceneElement = firstPatient.querySelector('.scene');
		var differenceButton = document.createElement('div')
		differenceButton.className = 'sceneToggleButton';
		differenceButton.innerHTML = name;
		firstPatient.insertBefore(differenceButton, firstPatient.children[1]);
		if(name == 'Actual'){
			differenceButton.style.opacity = 1;
		}
	});
	d3.selectAll('.sceneToggleButton').on('click', function(){
		var buttons = document.getElementsByClassName('sceneToggleButton');
		Array.prototype.forEach.call(buttons, function(e){
			e.style.opacity = '';
		});
		this.style.opacity = 1;
		Controller.switchScene(scenes[0], this.innerHTML, data);
	});
}

function updateOrder(updatedPatient) {
	//sorts the divs of list-items for the patients based on similarity score
    var lastPatient = document.getElementById(data.getPatientMatches(selectedPatient)[scenes.length - 1]);
	formatFirstPatient(updatedPatient);
    //insert last element from patientMatches in last place (before null)
    parent.insertBefore(lastPatient, null);
	
	var first;
	var second;
	var patientMatches = data.getPatientMatches(selectedPatient);
    for (var i = (scenes.length - 2); i > 0; i--) {

        first = document.getElementById(patientMatches[i]);
        second = document.getElementById(patientMatches[i + 1]);
        // order div elements
        parent.insertBefore(first, second);
		//updates the similarity score for the second patient
        pScoreElement = second.querySelector(".pScore");
        // update patient score
        pScoreElement.innerHTML = data.getPatientSimilarityScores(selectedPatient)[i+1].toFixed(5);
        // hide patients
        second.style.display = "none";
    }
	//update similarity for the first non-self match
	if(scenes.length > 1){
		var pScoreElement = first.querySelector(".pScore");
		pScoreElement.innerHTML = data.getPatientSimilarityScores(selectedPatient)[i+1].toFixed(5);
		first.style.display = "none";
	}

}

function initializeRiskPrediction(rank) {
	
    var simScores = data.getPatientSimilarityScores(rank);
    for (var j = 0; j < patientsToShow && j < simScores.length; j++) {
        var p = document.getElementById(data.getPatientMatches(selectedPatient)[j]);
        p.style.display = "inline-block";
    }

}

function animate() {
    render();
    requestAnimationFrame(animate);
}

function render() {

    updateSize();

    renderer.setClearColor(0xbbbbbb);//will be background color
    renderer.setScissorTest(false);
    renderer.clear();

    renderer.setClearColor(0x888888);//will be color in viewport?
    renderer.setScissorTest(true);

    updateMainView();
	
}

function updateMainView(rotMatrix) {
	var raycaster = new THREE.Raycaster();
    var rotMatrix = new THREE.Matrix4();
	//scenes = updateScenes(selectedPatient, materialArray);

	for (var index = 0; index < scenes.length; index++) {
		var scene = scenes[index];
		var controls = scene.userData.controls;
		var camera = scene.userData.camera;

		var orientMarkerCube = camera.children[0];

		// get the element that is a place holder for where we want to
		// draw the scene
		var element = scene.userData.element;

		// get its position relative to the page's viewport
		var rect = element.getBoundingClientRect();

		// check if it's offscreen. If so skip it
		if (rect.bottom < 0 || rect.top > renderer.domElement.clientHeight ||
			rect.right < 0 || rect.left > renderer.domElement.clientWidth) {
			continue; // it's off screen
		}

		// update orientation marker
		rotMatrix.extractRotation(controls.object.matrix);
		orientMarkerCube.rotation.setFromRotationMatrix(rotMatrix.transpose());

		// set the viewport
		var width = rect.right - rect.left;
		var height = rect.bottom - rect.top;
		var left = rect.left;
		var bottom = renderer.domElement.clientHeight - rect.bottom;

		renderer.setViewport(left, bottom, width, height);
		renderer.setScissor(left, bottom, width, height);

		// raycaster
		raycaster.setFromCamera(mouseNorm, currScene.userData.camera);

		var intersects = raycaster.intersectObjects(currScene.children);

		if (intersects.length >= 1 && detailsOnRotate) {

			for (var i = intersects.length - 1; i >= 0; i--) {

				if (intersects[i].object.userData.type == "node") {

					nodeHover = intersects[i].object;
					var tempObject = scene.getObjectByName(nodeHover.name + "_outline");

					if (INTERSECTED != tempObject) {

						if (INTERSECTED) {
							INTERSECTED.material.color.setHex(INTERSECTED.currentHex);
						}

						INTERSECTED = tempObject;

						if (INTERSECTED) {
							INTERSECTED.currentHex = INTERSECTED.material.color.getHex();
							INTERSECTED.material.color.setHex(0x00e4ff);
						}

						// details

						populateAndPlaceDetails("SHOW");

					}

					break;

				} else {
					populateAndPlaceDetails("HIDE");
				}


			}
		} else {

			if (INTERSECTED) {
				INTERSECTED.material.color.setHex(INTERSECTED.currentHex);
				// details
				populateAndPlaceDetails("HIDE");

			}

			INTERSECTED = null;
		}
		renderer.render(scene, camera);
	}
}

function updateSize() {

    var width = canvas.clientWidth;
    var height = canvas.clientHeight;

    if (canvas.width !== width || canvas.height != height)
        renderer.setSize(width, height, false);
}

function populateAndPlaceDetails(state) {

    if (state == "SHOW") {

        nodeDetails.style.display = "block";
        // PLACEMENT
        // check if details are offscreen, then shift appropriately
        // X, add 10 pixels for buffer, since width is dynamic
        if (mouse.x + detailsOffsetX + nodeDetails.offsetWidth + 5 >= canvas.clientWidth) {
            nodeDetails.style.left = (mouse.x - detailsOffsetX - nodeDetails.offsetWidth) + "px";
        } else {
            nodeDetails.style.left = (mouse.x + detailsOffsetX) + "px";
        }

        // Y
        if (mouse.y + detailsOffsetY + nodeDetails.offsetHeight + 5 >= canvas.clientHeight) {
            nodeDetails.style.top = (mouse.y - detailsOffsetY - nodeDetails.offsetHeight) + "px";
        } else {
            nodeDetails.style.top = (mouse.y + detailsOffsetY) + "px";
        }
       
        // POPULATE

        // Organ name
        organName.innerHTML = nodeHover.name;

        // Dose Per Volume
        dosePerVolume.innerHTML = nodeHover.userData.dosePerVolume;

        // line separator
        lineSeparator.style["borderColor"] = "#" + nodeHover.material.color.getHexString();

        // Volume
        volumeVal.innerHTML = nodeHover.userData.volume + "";

        // Mean Dose
        meanDoseVal.innerHTML = nodeHover.userData.meanDose + "  GY";

        // Min Dose
        minDoseVal.innerHTML = nodeHover.userData.minDose + "";

        // Max Dose
        maxDoseVal.innerHTML = nodeHover.userData.maxDose + "";

    } else if (state == "HIDE") {

        nodeDetails.style.display = "none";
        nodeDetails.style.top = -500 + "px";
        nodeDetails.style.left = -500 + "px";
    }
}

function onMouseDown(event) {

    if (event.target) {

        detailsOnRotate = false;
        handleInputRotate(event);
    }
}

function onTouchStart(event) {

    if (event.target) {

        detailsOnRotate = false;
        handleInputRotate(event);
    }
}

function getSceneIndex(internalId){
	//gets the index in the scene list from a given internal id
	var index = data.getPatientMatches(selectedPatient).indexOf( +internalId )
	return(index)
}

function handleInputRotate(event) {

    var targ = event.target,
        cameraToCopy;

    if (targ.className == "scene" && syncCameras == true) {
		var index;
        if (targ.parentNode.hasAttribute("id")) {
			index = getSceneIndex( +targ.parentNode.id );
            cameraToCopy = scenes[index].userData.camera;
        } 
		else{
			cameraToCopy = scenes[scenes.length - 1].userData.camera;
		}

        // 20 milliseconds interval => 50 FPS
        syncCamerasInterval = setInterval(syncAllCameras, 20, cameraToCopy);
    }
}

function syncAllCameras(cameraToCopy) {

    for( var i = 0; i < scenes.length; i++) {

        var scene = scenes[i];
        var camera = scene.userData.camera;
        var controls = scene.userData.controls;

        camera.position.subVectors(cameraToCopy.position, controls.target);
        camera.position.setLength(cameraDistZ);
        camera.lookAt(scene.position);
    };
}

function onMouseUp(event) {
    detailsOnRotate = true;
    clearInterval(syncCamerasInterval);
}

function onTouchEnd(event) {
    detailsOnRotate = true;
    clearInterval(syncCamerasInterval);
}

function onDocumentMouseMove(event) {
    if (event.target) {

        var targ = event.target;

        if (targ.className == "scene") {

			let index = getSceneIndex(+targ.parentNode.id);
			currScene = (index > -1)? scenes[index]: scenes[scenes.length-1];

            mouse.x = event.clientX;
            mouse.y = event.clientY;

            mouseNorm.x = (event.offsetX / targ.offsetWidth) * 2 - 1;
            mouseNorm.y = -(event.offsetY / targ.offsetHeight) * 2 + 1;
        }
    }
}

document.getElementById("opacSlider").oninput = function () {

    var opac = (this.value / 100.0);
	ColorScale.setOpacity(opac);
	console
    scenes.forEach(function (scene, index) {

        for (var pOrgan in oAtlas) {

            var tempObject = scene.getObjectByName(pOrgan + "_model");

            if (tempObject)
                tempObject.material.opacity = opac;
        }
    });
}
