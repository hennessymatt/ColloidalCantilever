SetFactory("OpenCASCADE");

// Substrate dimensions (non-dimensionalised)
H = 10;
W = 8/22;
T = 75e-6 / 22e-3;

// Characteristic element lengths
lcBulk = H/20;
lcRect = W / 10000;

// -----------------------------------------------------------------------------
// Geometry
// -----------------------------------------------------------------------------
Point(101) = {0, 0, 0, lcRect};


// substrate
Box(1) = {0, 0, -H/2, H/2, H/2, H};

// box
Box(100) = {0, 0, -T, 0.5, W/2, T};

// Split the box so the rectangle becomes a separate surface
BooleanDifference(200) = { Volume{1}; Delete; }{ Volume{100}; Delete; };

// -----------------------------------------------------------------------------
// Mesh refinement around the contact line
// -----------------------------------------------------------------------------

// Distance from the contact line
Field[1] = Distance;
Field[1].CurvesList = {18, 19};   // line IDs of the contact line
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcRect;      // fine mesh
Field[2].SizeMax = lcBulk;      // coarse mesh
Field[2].DistMin = 0;
Field[2].DistMax = 3;

Background Field = 2;

// Disable automatic sizing
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

// create IDs for use in FEniCS
Physical Surface(999) = {5};
Physical Volume(1000) = {200};

