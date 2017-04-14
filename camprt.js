'use strict';

// SOME QUICK, TEMPORARY NOTES:
//      update color scale/range
//      dynamic dropdown
//      top navbar slide down to show extended view of selected patient

if (!Detector.webgl) {
    Detector.addGetWebGLMessage();
}

var canvas;

var parent = document.getElementById("content");
var nodeDetails = document.getElementById("details");

var detailsOffsetX = 10;
var detailsOffsetY = 10;

var organName = document.getElementById("details_organName");
var dosePerVolume = document.getElementById("details_dosePerVolume");
var lineSeparator = document.getElementById("details_line");
var volumeVal = document.getElementById("details_volume_val");
var meanDoseVal = document.getElementById("details_meanDose_val");
var minDoseVal = document.getElementById("details_minDose_val");
var maxDoseVal = document.getElementById("details_maxDose_val");

var scenes = [],
    renderer;

var sceneCenter = [0, 0, 0];

var selectedPatient = 1;

var syncCameras = true,
    syncCamerasInterval;

var raycaster;

var mouse = new THREE.Vector2();

var mouseNorm = new THREE.Vector2(),
    INTERSECTED, nodeHover;

var width, height;

// 36 steps

var domainColorScale = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0];
var rangeColorScale = ['#ffffe0', '#fff8d2', '#fff0c4', '#ffe9b8', '#ffe2ae', '#ffdaa3', '#ffd39a', '#ffcb91', '#ffc389', '#ffbb82', '#ffb27c', '#ffab77', '#ffa272', '#ff986e', '#fe906a', '#fb8768', '#f98065', '#f67762', '#f26f60', '#ee675d', '#eb5f5b', '#e75758', '#e25055', '#dd4852', '#d8404e', '#d3394a', '#cc3146', '#c62a41', '#c0223b', '#b91c35', '#b3152f', '#ab0e28', '#a40820', '#9b0317', '#93010e', '#8b0000'];
var color = d3.scaleLinear()
    .domain(domainColorScale)
    .range(rangeColorScale);

// data
var organs = [
    {
        name: "Brainstem",
        x: 480,
        y: 210,
        z: 0
    },
    {
        name: "Cricopharyngeal_Muscle",
        x: 410,
        y: 700,
        z: 0
    },
    {
        name: "Esophagus",
        x: 450,
        y: 810,
        z: 0
    },
    {
        name: "Extended_Oral_Cavity",
        x: 225,
        y: 400,
        z: 0
    },
    {
        name: "Glottic_Area",
        x: 340,
        y: 715,
        z: 0
    },
    {
        name: "IPC",
        x: 405,
        y: 625,
        z: 0
    },
    {
        name: "Lower_Lip",
        x: 65,
        y: 445,
        z: 0
    },
    {
        name: "Lt_Anterior_Seg_Eyeball",
        x: 140,
        y: 75,
        z: 65
    },
    {
        name: "Lt_Parotid_Gland",
        x: 400,
        y: 340,
        z: 135
    },
    {
        name: "Lt_Posterior_Seg_Eyeball",
        x: 190,
        y: 75,
        z: 65
    },
    {
        name: "Lt_Submandibular_Gland",
        x: 325,
        y: 500,
        z: 90
    },
    {
        name: "Lt_thyroid_lobe",
        x: 385,
        y: 760,
        z: 45
    },
    {
        name: "MPC",
        x: 400,
        y: 560,
        z: 0
    },
    {
        name: "Pituitary_Gland",
        x: 395,
        y: 135,
        z: 0
    },
    {
        name: "Rt_Anterior_Seg_Eyeball",
        x: 140,
        y: 75,
        z: -65
    },
    {
        name: "Rt_Parotid_Gland",
        x: 400,
        y: 340,
        z: -135
    },
    {
        name: "Rt_Posterior_Seg_Eyeball",
        x: 190,
        y: 75,
        z: -65
    },
    {
        name: "Rt_Submandibular_Gland",
        x: 325,
        y: 500,
        z: -90
    },
    {
        name: "Rt_thyroid_lobe",
        x: 385,
        y: 760,
        z: -45
    },
    {
        name: "SPC",
        x: 375,
        y: 430,
        z: 0
    },
    {
        name: "Spinal_Cord",
        x: 505,
        y: 610,
        z: 0
    },
    {
        name: "Supraglottic_Larynx",
        x: 310,
        y: 590,
        z: 0
    },
    {
        name: "Thyroid_cartilage",
        x: 345,
        y: 655,
        z: 0
    },
    {
        name: "Upper_Lip",
        x: 70,
        y: 270,
        z: 0
    }
];

var links = [
    {
        source: 2,
        target: 1
    },
    {
        source: 4,
        target: 1
    },
    {
        source: 5,
        target: 1
    },
    {
        source: 5,
        target: 4
    },
    {
        source: 6,
        target: 3
    },
    {
        source: 9,
        target: 7
    },
    {
        source: 10,
        target: 3
    },
    {
        source: 10,
        target: 8
    },
    {
        source: 11,
        target: 1
    },
    {
        source: 11,
        target: 2
    },
    {
        source: 11,
        target: 4
    },
    {
        source: 12,
        target: 5
    },
    {
        source: 13,
        target: 0
    },
    {
        source: 13,
        target: 9
    },
    {
        source: 16,
        target: 9
    },
    {
        source: 16,
        target: 13
    },
    {
        source: 16,
        target: 14
    },
    {
        source: 17,
        target: 3
    },
    {
        source: 17,
        target: 15
    },
    {
        source: 18,
        target: 1
    },
    {
        source: 18,
        target: 2
    },
    {
        source: 18,
        target: 4
    },
    {
        source: 19,
        target: 3
    },
    {
        source: 19,
        target: 8
    },
    {
        source: 19,
        target: 10
    },
    {
        source: 19,
        target: 12
    },
    {
        source: 19,
        target: 15
    },
    {
        source: 19,
        target: 17
    },
    {
        source: 20,
        target: 0
    },
    {
        source: 20,
        target: 1
    },
    {
        source: 20,
        target: 2
    },
    {
        source: 20,
        target: 5
    },
    {
        source: 20,
        target: 12
    },
    {
        source: 20,
        target: 19
    },
    {
        source: 21,
        target: 3
    },
    {
        source: 21,
        target: 5
    },
    {
        source: 21,
        target: 10
    },
    {
        source: 21,
        target: 12
    },
    {
        source: 21,
        target: 17
    },
    {
        source: 22,
        target: 1
    },
    {
        source: 22,
        target: 4
    },
    {
        source: 22,
        target: 5
    },
    {
        source: 22,
        target: 11
    },
    {
        source: 22,
        target: 18
    },
    {
        source: 22,
        target: 21
    },
    {
        source: 23,
        target: 3
    },
    {
        source: 23,
        target: 6
    }
];

var patients = [
    {
        id: 1,
        organData: [0.558468835, 5.71615, 4.788434315, 0.43009521, 103.2081633, 23.0260989, 7.947596154, 4.276848214, 1.394489302, 0.159385569, 10.42983139, 8.29548, 54.33392857, 13.88121212, 4.539666667, 1.488854721, 0.20479439, 7.163132022, 9.0047346, 5.690875561, 1.156974432, 4.944950922, 6.040249516, 5.65825188],
        similarity: [1, 6, 9, 10, 15, 11, 13, 5, 7, 8, 3, 14, 12, 2, 4],
        scores: [1, 0.98852, 0.97781, 0.95906, 0.85954, 0.82487, 0.81777, 0.79799, 0.72638, 0.72515, 0.67101, 0.66705, 0.65221, 0.64747, 0.63087]
    },
    {
        id: 2,
        organData: [0.280421735, 8.626207165, 5.094619387, 0.5438902, 9.425348558, 8.438579137, 7.649184261, 2.663206081, 0.666099979, 0.107645727, 8.296629213, 7.869070513, 48.65405629, 4.625408654, 2.460351351, 0.783078038, 0.107504118, 6.705035634, 8.126458991, 5.329027952, 1.546326185, 3.755602981, 2.707105179, 1.416886062],
        similarity: [2, 12, 14, 3, 4, 7, 8, 5, 11, 13, 10, 1, 6, 9, 15],
        scores: [1, 0.99271, 0.9846, 0.98368, 0.981, 0.97851, 0.97527, 0.96433, 0.93783, 0.92883, 0.8219, 0.64747, 0.61675, 0.54321, 0.53921]
    },
    {
        id: 3,
        organData: [0.737508165, 4.18150134, 3.768583333, 0.440693347, 7.815465686, 7.676552463, 3.658615819, 2.524771739, 0.728994956, 0.128690249, 5.299900398, 7.332373978, 32.87736345, 7.134477273, 2.424077778, 1.09028581, 0.116227011, 7.373717159, 7.279570601, 4.96604118, 1.384190107, 3.406903907, 2.629715447, 1.646442203],
        similarity: [3, 12, 8, 2, 14, 4, 7, 5, 13, 11, 10, 1, 6, 15, 9],
        scores: [1, 0.99154, 0.9884, 0.98368, 0.98304, 0.98197, 0.96642, 0.96392, 0.95597, 0.9329, 0.83861, 0.67101, 0.62817, 0.56578, 0.56431]
    },
    {
        id: 4,
        organData: [0.620242425, 10.0390431, 4.605362289, 0.317255936, 7.979099431, 6.989797735, 5.127522732, 5.0, 1.157517642, 0.008829494, 8.038825449, 10.69167602, 44.39099993, 11.69622222, 3.5, 0.410665789, 0.017084337, 5.499135206, 7.779198376, 5.441099164, 2.248953763, 2.421176186, 2.421139198, 2.791346835],
        similarity: [4, 3, 2, 12, 8, 14, 7, 5, 13, 11, 10, 1, 6, 15, 9],
        scores: [1, 0.98197, 0.981, 0.98075, 0.97407, 0.96526, 0.95264, 0.95057, 0.92991, 0.90908, 0.80366, 0.63087, 0.59485, 0.55144, 0.53177]
    },
    {
        id: 5,
        organData: [0.655318663, 6.462222222, 3.877905973, 0.589977007, 30.5016282, 15.61439032, 8.356816072, 2.454563094, 0.741313635, 0.145199527, 6.438250475, 10.18152721, 69.46666667, 13.15310789, 2.384924343, 0.729239751, 0.148193886, 8.493031953, 8.820359111, 7.504534283, 1.2838198, 3.769131093, 3.820019524, 4.636100213],
        similarity: [5, 11, 7, 8, 14, 12, 13, 2, 3, 4, 10, 1, 6, 9, 15],
        scores: [1, 0.98795, 0.98563, 0.97639, 0.97226, 0.97034, 0.96867, 0.96433, 0.96392, 0.95057, 0.91762, 0.79799, 0.7638, 0.70176, 0.64824]
    },
    {
        id: 6,
        organData: [0.343556805, 20.3412069, 3.341508264, 0.277945811, 114.721558, 23.06959248, 4.171597813, 3.473519737, 0.679473426, 0.14261899, 5.852044189, 8.227883489, 55.84084967, 6.084488636, 3.372006579, 0.716827483, 0.141357072, 6.166703997, 6.096322899, 4.103134498, 1.919095722, 5.105682484, 5.900119485, 1.724078624],
        similarity: [6, 1, 9, 10, 15, 11, 13, 5, 7, 8, 3, 14, 12, 2, 4],
        scores: [1, 0.98852, 0.9877, 0.94409, 0.86195, 0.79693, 0.77674, 0.7638, 0.69279, 0.68934, 0.62817, 0.62793, 0.61798, 0.61675, 0.59485]
    },
    {
        id: 7,
        organData: [0.393492914, 4.189319765, 3.382826814, 0.283459525, 22.01029504, 15.71520422, 4.945201797, 3.262997821, 0.835493174, 0.111638493, 11.22516959, 8.712209537, 71.81075635, 6.705153996, 2.356240928, 0.758031311, 0.127919841, 9.618584075, 5.854392232, 5.289597016, 1.302177426, 4.256181153, 4.009335079, 0.949267228],
        similarity: [7, 14, 5, 11, 12, 2, 8, 3, 4, 13, 10, 1, 6, 9, 15],
        scores: [1, 0.98704, 0.98563, 0.98273, 0.98172, 0.97851, 0.97288, 0.96642, 0.95264, 0.94067, 0.87093, 0.72638, 0.69279, 0.61334, 0.58484]
    },
    {
        id: 8,
        organData: [0.365695988, 6.532828283, 3.811310358, 0.270691115, 12.12019036, 12.64056437, 3.585918909, 3.879217645, 0.474309432, 0.177343265, 5.371368389, 6.435659686, 38.15218864, 8.874735707, 5.004254511, 0.926379401, 0.201799338, 6.530889074, 6.647172765, 4.664850434, 1.355461327, 3.236156089, 2.918632967, 1.423110339],
        similarity: [8, 3, 14, 12, 5, 2, 4, 7, 13, 11, 10, 1, 6, 15, 9],
        scores: [1, 0.9884, 0.98479, 0.98459, 0.97639, 0.97527, 0.97407, 0.97288, 0.9652, 0.95848, 0.87998, 0.72515, 0.68934, 0.63871, 0.6315]
    },
    {
        id: 9,
        organData: [0.593020745, 29.99923106, 3.406970392, 0.64181476, 145.6687243, 34.84779372, 14.22666132, 3.220452311, 1.472626122, 0.197087583, 7.595372379, 10.76234157, 50.39661376, 15.76883867, 6.138791423, 1.93430534, 0.272766766, 7.527368547, 9.800206892, 5.85563747, 2.184065764, 4.870993562, 6.713010938, 9.77048433],
        similarity: [9, 6, 1, 10, 15, 13, 11, 5, 8, 7, 3, 14, 12, 2, 4],
        scores: [1, 0.9877, 0.97781, 0.9149, 0.86564, 0.7409, 0.72949, 0.70176, 0.6315, 0.61334, 0.56431, 0.55797, 0.54463, 0.54321, 0.53177]
    },
    {
        id: 10,
        organData: [0.502817027, 8.653924707, 3.088327985, 0.365751339, 45.56495726, 18.12276369, 4.479639861, 2.25950275, 1.284455705, 0.094597395, 8.445689412, 6.815160273, 41.49472606, 7.227896576, 2.065948525, 1.633247, 0.096893273, 5.988334002, 8.51438815, 4.936995698, 1.860516297, 4.002530134, 4.642835979, 2.765339622],
        similarity: [10, 1, 6, 11, 13, 5, 9, 8, 7, 3, 14, 12, 2, 15, 4],
        scores: [1, 0.95906, 0.94409, 0.9331, 0.93257, 0.91762, 0.9149, 0.87998, 0.87093, 0.83861, 0.83795, 0.82705, 0.8219, 0.81726, 0.80366]
    },
    {
        id: 11,
        organData: [0.276963108, 3.678883521, 1.413787674, 0.397958187, 37.8617284, 23.343195, 4.722454186, 3.134282756, 0.732981938, 0.130035589, 7.648657137, 7.692927973, 75.85784695, 7.601089825, 2.487824799, 0.64939162, 0.133592375, 8.008232949, 5.637612104, 5.144510236, 1.056667194, 4.479145227, 5.459451355, 2.242058538],
        similarity: [11, 5, 7, 14, 8, 13, 12, 2, 10, 3, 4, 1, 6, 9, 15],
        scores: [1, 0.98795, 0.98273, 0.96009, 0.95848, 0.94964, 0.94628, 0.93783, 0.9331, 0.9329, 0.90908, 0.82487, 0.79693, 0.72949, 0.66962]
    },
    {
        id: 12,
        organData: [0.428760143, 6.921007902, 4.181903115, 0.473941929, 9.32425281, 10.69307049, 3.749659464, 2.265828427, 0.604519304, 0.156541212, 5.739402209, 8.320531846, 48.20406213, 6.824231984, 2.227761675, 1.109763757, 0.200330527, 9.134059298, 8.518690966, 5.522966531, 1.506542402, 2.988626328, 2.745435244, 2.36697597],
        similarity: [12, 2, 14, 3, 8, 7, 4, 5, 11, 13, 10, 1, 6, 9, 15],
        scores: [1, 0.99271, 0.99234, 0.99154, 0.98459, 0.98172, 0.98075, 0.97034, 0.94628, 0.93967, 0.82705, 0.65221, 0.61798, 0.54463, 0.5356]
    },
    {
        id: 13,
        organData: [0.977529748, 6.153779807, 4.231566822, 0.560924555, 21.31210686, 14.50122483, 6.561851852, 0, 1.991365267, 0, 8.830564398, 9.762879074, 37.86170407, 12.63507985, 0, 1.10968098, 0, 9.690783948, 5.970156376, 6.813810501, 1.195058281, 4.36867367, 3.887660923, 4.939640837],
        similarity: [13, 5, 8, 3, 11, 14, 7, 12, 10, 4, 2, 1, 6, 9, 15],
        scores: [1, 0.96867, 0.9652, 0.95597, 0.94964, 0.9462, 0.94067, 0.93967, 0.93257, 0.92991, 0.92883, 0.81777, 0.77674, 0.7409, 0.65807]
    },
    {
        id: 14,
        organData: [0.325656505, 4.525185185, 2.982656908, 0.375619949, 10.63184271, 16.35188345, 5.960115856, 1.98013468, 0.66581052, 0.122201296, 7.240405214, 8.229796685, 51.29126222, 6.750709877, 2.541971014, 1.028787226, 0.114437862, 9.080851911, 6.121684057, 5.126579596, 1.600694251, 3.983301788, 2.857290525, 2.20016835],
        similarity: [14, 12, 7, 8, 2, 3, 5, 4, 11, 13, 10, 1, 6, 9, 15],
        scores: [1, 0.99234, 0.98704, 0.98479, 0.9846, 0.98304, 0.97226, 0.96526, 0.96009, 0.9462, 0.83795, 0.66705, 0.62793, 0.55797, 0.53908]
    },
    {
        id: 15,
        organData: [0.336028642, 22.17806026, 3.0915369, 0.475937091, 93.2597161, 18.53507778, 5.963149923, 41.95178451, 1.031331953, 0.276064505, 9.750421919, 10.30118219, 38.63746985, 7.622099694, 51.35733333, 1.220450746, 0.294246771, 9.190413943, 9.049147614, 4.665860977, 1.740893437, 3.157518832, 5.912535815, 5.269529543],
        similarity: [15, 9, 6, 1, 10, 11, 13, 5, 8, 7, 3, 4, 2, 14, 12],
        scores: [1, 0.86564, 0.86195, 0.85954, 0.81726, 0.66962, 0.65807, 0.64824, 0.63871, 0.58484, 0.56578, 0.55144, 0.53921, 0.53908, 0.5356]
    }
];


var pRankingOrder = patients[selectedPatient - 1].similarity;
var pScores = patients[selectedPatient - 1].scores;

populateColorScale();
flipGraph(); // fixes orientation of organs
computeCenterOfGraph(); // compute center of graph
shiftGraphToOrigin(); // center graph to origin
init(); // initialize

var listItems = document.getElementsByClassName("list-item"),
    arrayOfDivs = [];

for (var i = 0, ref = arrayOfDivs.length = listItems.length; i < ref; i++) {
    arrayOfDivs[i] = listItems[i];
}

var currScene = scenes[0];

document.addEventListener("mousedown", onMouseDown, false);
document.addEventListener("mouseup", onMouseUp, false);

document.addEventListener("touchstart", onTouchStart, false);
document.addEventListener("touchend", onTouchEnd, false);

document.addEventListener("mousemove", onDocumentMouseMove, false)

updateOrder(selectedPatient); // update order in GUI
animate(); // render

// ----------------------------------------------------------------

function populateColorScale() {

    var parentDiv = document.getElementById("colorScale");

    rangeColorScale.forEach(function (color, index) {

        var tempDiv = document.createElement("div");

        tempDiv.style.height = "100%";
        tempDiv.style.width = "9px";
        tempDiv.style["backgroundColor"] = color;
        tempDiv.style.display = "inline-block";
        tempDiv.style.transition = "0.3s";

        parentDiv.appendChild(tempDiv);
    });


}

function flipGraph() {

    organs.forEach(function (organ, index) {

        organ.x = (organ.x * -1);
        organ.y = (organ.y * -1);
        organ.z = (organ.z * -1);
    });
}

function computeCenterOfGraph() {

    var xyzMin = new Array(3);
    var xyzMax = new Array(3);

    xyzMin = getMin();
    xyzMax = getMax();

    sceneCenter = [
        Math.round((xyzMin[0] + xyzMax[0]) / 2),
        Math.round((xyzMin[1] + xyzMax[1]) / 2),
        Math.round((xyzMin[2] + xyzMax[2]) / 2)
    ];
}

function getMin() {

    var x = organs.reduce(function (min, obj) {
        return obj.x < min ? obj.x : min;
    }, Infinity);

    var y = organs.reduce(function (min, obj) {
        return obj.y < min ? obj.y : min;
    }, Infinity);

    var z = organs.reduce(function (min, obj) {
        return obj.z < min ? obj.z : min;
    }, Infinity);


    return [x, y, z];
}

function getMax() {

    var x = organs.reduce(function (max, obj) {
        return obj.x > max ? obj.x : max;
    }, -Infinity);

    var y = organs.reduce(function (max, obj) {
        return obj.y > max ? obj.y : max;
    }, -Infinity);

    var z = organs.reduce(function (max, obj) {
        return obj.z > max ? obj.z : max;
    }, -Infinity);


    return [x, y, z];
}

function shiftGraphToOrigin() {

    organs.forEach(function (organ, index) {

        organ.x = (organ.x - sceneCenter[0]);
        organ.y = (organ.y - sceneCenter[1]);
        organ.z = (organ.z - sceneCenter[2]);
    });
}

function init() {

    canvas = document.getElementById("c");

    var template = document.getElementById("template").text;

    renderer = new THREE.WebGLRenderer({
        canvas: canvas,
        antialias: true
    });

    renderer.setClearColor(0xffffff, 1);
    renderer.setPixelRatio(window.devicePixelRatio);
    //renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.sortObjects = false;

    raycaster = new THREE.Raycaster();

    //
    var maxAnisotropy = renderer.getMaxAnisotropy();

    var textureLoader = new THREE.TextureLoader();

    var texture0 = textureLoader.load('data/anterior.png'); // xpos, Right
    var texture1 = textureLoader.load('data/posterior.png'); // xneg, Left
    var texture2 = textureLoader.load('data/superior.png'); // ypos, Top
    var texture3 = textureLoader.load('data/inferior.png'); // yneg, Bottom
    var texture4 = textureLoader.load('data/right.png'); // zpos, Back
    var texture5 = textureLoader.load('data/left.png'); // zneg, Front

    texture0.anisotropy = maxAnisotropy;
    texture1.anisotropy = maxAnisotropy;
    texture2.anisotropy = maxAnisotropy;
    texture3.anisotropy = maxAnisotropy;
    texture4.anisotropy = maxAnisotropy;
    texture5.anisotropy = maxAnisotropy;

    var materialArray = [
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
    //

    for (var i = 0; i < patients.length; i++) {

        var scene = new THREE.Scene();

        // make a list item
        var element = document.createElement("div");
        element.className = "list-item";
        element.id = i + 1;
        element.innerHTML = template.replace('$', i + 1);

        // Look up the element that represents the area
        // we want to render the scene
        scene.userData.element = element.querySelector(".scene");

        parent.appendChild(element);

        var camera = new THREE.OrthographicCamera(scene.userData.element.offsetWidth / -.8, scene.userData.element.offsetWidth / .8, scene.userData.element.offsetHeight / .8, scene.userData.element.offsetHeight / -.8, 1, 100000);
        camera.position.z = 2000;

        camera.updateProjectionMatrix();
        scene.userData.camera = camera;

        var controls = new THREE.OrbitControls(scene.userData.camera, scene.userData.element);
        controls.minDistance = 2;
        controls.maxDistance = 1000;
        controls.enablePan = false;
        controls.enableZoom = false;

        scene.userData.controls = controls;

        var geometry = new THREE.SphereGeometry(20, 16, 16);

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
            linewidth: 1
        });

        // nodes (scene.children from index 0 to 23)
        organs.forEach(function (organ, index) {

            scene.add(new THREE.Mesh(geometry, material.clone()));

            scene.children[index].position.x = (organ.x);
            scene.children[index].position.y = (organ.y);
            scene.children[index].position.z = (organ.z);

            scene.children[index].name = organ.name;
            scene.children[index].userData.type = "node";
        });

        // node outline
        organs.forEach(function (organ, index) {

            var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

            outlineMesh.position.x = (organ.x);
            outlineMesh.position.y = (organ.y);
            outlineMesh.position.z = (organ.z);

            outlineMesh.name = organ.name + "_outline";

            outlineMesh.scale.multiplyScalar(1.1);
            scene.add(outlineMesh);
        });

        // color the nodes
        var pOrganData = patients[i].organData;

        organs.forEach(function (organ, index) {

            if (scene.children[index].name != organ.name)
                console.log("Something isn't right");

            scene.children[index].userData.dosePerVolume = pOrganData[index];

            var nodeColor = color(pOrganData[index]);

            scene.children[index].material.color.setStyle(nodeColor);
        });

        // links
        links.forEach(function (link, index) {

            var tmp_geo = new THREE.Geometry();

            tmp_geo.vertices.push(scene.children[link.source].position);
            tmp_geo.vertices.push(scene.children[link.target].position);

            var line = new THREE.LineSegments(tmp_geo, linkMaterial);
            line.scale.x = line.scale.y = line.scale.z = 1;
            line.originalScale = 1;

            //line.frustumCulled = false;

            scene.add(line);
        });


        scene.add(camera);

        // orientation marker, patient coordinate system
        var MovingCubeMat = new THREE.MultiMaterial(materialArray);
        var MovingCubeGeom = new THREE.CubeGeometry(65, 65, 65, 1, 1, 1, materialArray);
        var MovingCube = new THREE.Mesh(MovingCubeGeom, MovingCubeMat);

        camera.add(MovingCube);
        MovingCube.position.set(365, -365, -1000);

        var light = new THREE.AmbientLight(0xffffff, 1.0); // white light
        scene.add(light);


        scenes.push(scene);
    }

    //renderer.autoClear = false;
}

function updateOrder(updatedPatient) {

    selectedPatient = updatedPatient;
    pRankingOrder = patients[selectedPatient - 1].similarity;
    pScores = patients[selectedPatient - 1].scores;

    //insert last element from pRankingOrder in last place (before null)
    parent.insertBefore(arrayOfDivs[(pRankingOrder[pRankingOrder.length - 1] - 1)], null);

    var pScoreElement = arrayOfDivs[(pRankingOrder[0] - 1)].querySelector(".pScore");

    // first patient always has score of 1, clear it
    pScoreElement.innerHTML = "";

    for (var i = (pRankingOrder.length - 2); i >= 0; i--) {

        // order div elements
        parent.insertBefore(arrayOfDivs[(pRankingOrder[i] - 1)], arrayOfDivs[(pRankingOrder[i + 1] - 1)]);

        pScoreElement = arrayOfDivs[(pRankingOrder[i + 1] - 1)].querySelector(".pScore");

        // update patient score
        pScoreElement.innerHTML = pScores[i + 1];
    }
}

function animate() {

    render();
    requestAnimationFrame(animate);
}

function render() {

    updateSize();

    renderer.setClearColor(0xffffff);
    renderer.setScissorTest(false);
    renderer.clear();

    renderer.setClearColor(0xa5a5a5);
    renderer.setScissorTest(true);

    pRankingOrder = patients[selectedPatient - 1].similarity;

    var rotMatrix = new THREE.Matrix4();

    pRankingOrder.forEach(function (rank, index) {

        var scene = scenes[rank - 1];
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

            return; // it's off screen
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

        //controls.update();

        // raycaster
        //raycaster.setFromCamera(mouseNorm, camera);
        raycaster.setFromCamera(mouseNorm, currScene.userData.camera);

        //var intersects = raycaster.intersectObjects(scene.children);
        var intersects = raycaster.intersectObjects(currScene.children);

        if (intersects.length > 1 && intersects[0].object.userData.type == "node") {

            nodeHover = intersects[0].object;
            var tempObject = scene.getObjectByName(nodeHover.name + "_outline");

            if (INTERSECTED != tempObject) {

                if (INTERSECTED)
                    INTERSECTED.material.color.setHex(INTERSECTED.currentHex);

                INTERSECTED = tempObject;

                INTERSECTED.currentHex = INTERSECTED.material.color.getHex();
                INTERSECTED.material.color.setHex(0x00e4ff);

                // details

                // mouse
                // check if details are offscreen, then shift appropriately
                // X, add 10 pixels for buffer, since width is dynamic
                if (mouse.x + detailsOffsetX + nodeDetails.offsetWidth + 10 > canvas.clientWidth) {
                    nodeDetails.style.left = (mouse.x - detailsOffsetX - nodeDetails.offsetWidth) + "px";
                } else {
                    nodeDetails.style.left = (mouse.x + detailsOffsetX) + "px";
                }

                // Y
                if (mouse.y + detailsOffsetY + nodeDetails.offsetHeight > canvas.clientHeight) {
                    nodeDetails.style.top = (mouse.y - detailsOffsetY - nodeDetails.offsetHeight) + "px";
                } else {
                    nodeDetails.style.top = (mouse.y + detailsOffsetY) + "px";
                }

                nodeDetails.style.display = "block";
                //nodeDetails.style.opacity = .95;

                //nodeDetails.style["borderColor"] = "#" + nodeHover.material.color.getHexString();

                // Organ name
                organName.innerHTML = nodeHover.name;

                // Dose Per Volume
                dosePerVolume.innerHTML = nodeHover.userData.dosePerVolume + " (GY/cc)";

                // line separator
                lineSeparator.style["borderColor"] = "#" + nodeHover.material.color.getHexString();


                // Volume
                //volumeVal.innerHTML =

                // Mean Dose
                //meanDoseVal.innerHTML =

                // Min Dose
                //minDoseVal.innerHTML =

                // Max Dose
                //maxDoseVal.innerHTML =

                //console.log(mouse.x);
                //console.log(mouse.y);

            }
        } else {

            if (INTERSECTED) {
                INTERSECTED.material.color.setHex(INTERSECTED.currentHex);

                // details

                // mouse

                nodeDetails.style.display = "none";
                //nodeDetails.style.opacity = .3;

                nodeDetails.style.top = -500 + "px";
                nodeDetails.style.left = -500 + "px";

            }

            INTERSECTED = null;
        }
        //

        renderer.render(scene, camera);
    });

    //renderer.render( currScene, currScene.userData.camera );
}

function updateSize() {

    width = canvas.clientWidth;
    height = canvas.clientHeight;

    if (canvas.width !== width || canvas.height != height)
        renderer.setSize(width, height, false);
}

function onMouseDown(event) {

    if (event.target) {

        handleInputRotate(event);
    }
}

function onTouchStart(event) {

    if (event.target) {

        handleInputRotate(event);
    }
}

function handleInputRotate(event) {

    var targ = event.target,
        cameraToCopy;

    if (targ.className == "scene" && syncCameras == true) {

        cameraToCopy = scenes[targ.parentNode.id - 1].userData.camera;

        // 20 milliseconds interval => 50 FPS
        syncCamerasInterval = setInterval(syncAllCameras, 20, cameraToCopy);
    }
}

function syncAllCameras(cameraToCopy) {

    pRankingOrder.forEach(function (rank, index) {

        var scene = scenes[rank - 1];
        var camera = scene.userData.camera;
        var controls = scene.userData.controls;

        camera.position.subVectors(cameraToCopy.position, controls.target);
        camera.position.setLength(2000);
        camera.lookAt(scene.position);
    });
}

function onMouseUp(event) {

    clearInterval(syncCamerasInterval);
}

function onTouchEnd(event) {

    clearInterval(syncCamerasInterval);
}

function onDocumentMouseMove(event) {

    //event.preventDefault();

    if (event.target) {

        var targ = event.target;

        if (targ.className == "scene") {

            currScene = scenes[targ.parentNode.id - 1];

            mouse.x = event.clientX;
            mouse.y = event.clientY;

            mouseNorm.x = (event.offsetX / targ.offsetWidth) * 2 - 1;
            mouseNorm.y = -(event.offsetY / targ.offsetHeight) * 2 + 1;
        }
    }
}
