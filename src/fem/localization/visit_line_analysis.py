name = 'line_vector_1000000.vtu'

OpenDatabase("localhost:/mn/sarpanitu/eksterne-u1/mirok/Documents/Programming/FractalFlow/src/fem/localization/line_curve000000.vtu", 0)
AddPlot("Pseudocolor", "f", 1, 1)
SetActivePlots(0)
SetActivePlots(0)
AddOperator("Threshold", 0)
ThresholdAtts = ThresholdAttributes()
ThresholdAtts.outputMeshType = 0
ThresholdAtts.boundsInputType = 0
ThresholdAtts.listedVarNames = ("default")
ThresholdAtts.zonePortions = (1)
ThresholdAtts.lowerBounds = (1)
ThresholdAtts.upperBounds = (1)
ThresholdAtts.defaultVarName = "f"
ThresholdAtts.defaultVarIsScalar = 1
ThresholdAtts.boundsRange = ("1:1")
SetOperatorOptions(ThresholdAtts, 0)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "xray"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 2
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 0
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)

OpenDatabase("localhost:/mn/sarpanitu/eksterne-u1/mirok/Documents/Programming/FractalFlow/src/fem/localization/%s" % name, 0)
AddPlot("Pseudocolor", "f", 1, 0)
SetActivePlots(1)
SetActivePlots(1)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "RdYlBu"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.Ramp  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 0.74902
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)
DrawPlots()
AddPlot("Pseudocolor", "f", 1, 0)
SetActivePlots(2)
AddOperator("Transform", 0)
TransformAtts = TransformAttributes()
TransformAtts.doRotate = 0
TransformAtts.rotateOrigin = (0, 0, 0)
TransformAtts.rotateAxis = (0, 0, 1)
TransformAtts.rotateAmount = 0
TransformAtts.rotateType = TransformAtts.Deg  # Deg, Rad
TransformAtts.doScale = 0
TransformAtts.scaleOrigin = (0, 0, 0)
TransformAtts.scaleX = 1
TransformAtts.scaleY = 1
TransformAtts.scaleZ = 1
TransformAtts.doTranslate = 1
TransformAtts.translateX = 1.1
TransformAtts.translateY = 0
TransformAtts.translateZ = 0
TransformAtts.transformType = TransformAtts.Similarity  # Similarity, Coordinate, Linear
TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
TransformAtts.outputCoordSys = TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
TransformAtts.continuousPhi = 0
TransformAtts.m00 = 1
TransformAtts.m01 = 0
TransformAtts.m02 = 0
TransformAtts.m03 = 0
TransformAtts.m10 = 0
TransformAtts.m11 = 1
TransformAtts.m12 = 0
TransformAtts.m13 = 0
TransformAtts.m20 = 0
TransformAtts.m21 = 0
TransformAtts.m22 = 1
TransformAtts.m23 = 0
TransformAtts.m30 = 0
TransformAtts.m31 = 0
TransformAtts.m32 = 0
TransformAtts.m33 = 1
TransformAtts.invertLinearTransform = 0
TransformAtts.vectorTransformMethod = TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
TransformAtts.transformVectors = 1
SetOperatorOptions(TransformAtts, 0)
AddOperator("Isosurface", 0)
DrawPlots()
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.27193, -0.434311, 0.858736)
View3DAtts.focus = (1.05, 0.5, 0.5)
View3DAtts.viewUp = (0.204958, 0.84575, 0.492646)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1.2659
View3DAtts.nearPlane = -2.5318
View3DAtts.farPlane = 2.5318
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (1.05, 0.5, 0.5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "RdYlBu"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.142627, -0.67995, 0.719254)
View3DAtts.focus = (1.05, 0.5, 0.5)
View3DAtts.viewUp = (-0.522676, 0.565356, 0.638108)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1.2659
View3DAtts.nearPlane = -2.5318
View3DAtts.farPlane = 2.5318
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (1.05, 0.5, 0.5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.20997, -0.6438, 0.735822)
View3DAtts.focus = (1.05, 0.5, 0.5)
View3DAtts.viewUp = (-0.520894, 0.563233, 0.641434)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1.2659
View3DAtts.nearPlane = -2.5318
View3DAtts.farPlane = 2.5318
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (1.05, 0.5, 0.5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.247432, -0.691605, 0.678572)
View3DAtts.focus = (1.05, 0.5, 0.5)
View3DAtts.viewUp = (-0.453384, 0.536298, 0.711918)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1.2659
View3DAtts.nearPlane = -2.5318
View3DAtts.farPlane = 2.5318
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (1.05, 0.5, 0.5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 1
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "RdYlBu"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 0
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.247432, -0.691605, 0.678572)
View3DAtts.focus = (1.05, 0.5, 0.5)
View3DAtts.viewUp = (-0.453384, 0.536298, 0.711918)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1.2659
View3DAtts.nearPlane = -2.5318
View3DAtts.farPlane = 2.5318
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (1.05, 0.5, 0.5)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = name + 'png'
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 2048
SaveWindowAtts.height = 2048
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SaveWindowAtts.subWindowAtts.win1.position = (0, 0)
SaveWindowAtts.subWindowAtts.win1.size = (128, 128)
SaveWindowAtts.subWindowAtts.win1.layer = 0
SaveWindowAtts.subWindowAtts.win1.transparency = 0
SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
SaveWindowAtts.subWindowAtts.win2.position = (0, 0)
SaveWindowAtts.subWindowAtts.win2.size = (128, 128)
SaveWindowAtts.subWindowAtts.win2.layer = 0
SaveWindowAtts.subWindowAtts.win2.transparency = 0
SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
SaveWindowAtts.subWindowAtts.win3.position = (0, 0)
SaveWindowAtts.subWindowAtts.win3.size = (128, 128)
SaveWindowAtts.subWindowAtts.win3.layer = 0
SaveWindowAtts.subWindowAtts.win3.transparency = 0
SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
SaveWindowAtts.subWindowAtts.win4.position = (0, 0)
SaveWindowAtts.subWindowAtts.win4.size = (128, 128)
SaveWindowAtts.subWindowAtts.win4.layer = 0
SaveWindowAtts.subWindowAtts.win4.transparency = 0
SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
SaveWindowAtts.subWindowAtts.win5.position = (0, 0)
SaveWindowAtts.subWindowAtts.win5.size = (128, 128)
SaveWindowAtts.subWindowAtts.win5.layer = 0
SaveWindowAtts.subWindowAtts.win5.transparency = 0
SaveWindowAtts.subWindowAtts.win5.omitWindow = 0
SaveWindowAtts.subWindowAtts.win6.position = (0, 0)
SaveWindowAtts.subWindowAtts.win6.size = (128, 128)
SaveWindowAtts.subWindowAtts.win6.layer = 0
SaveWindowAtts.subWindowAtts.win6.transparency = 0
SaveWindowAtts.subWindowAtts.win6.omitWindow = 0
SaveWindowAtts.subWindowAtts.win7.position = (0, 0)
SaveWindowAtts.subWindowAtts.win7.size = (128, 128)
SaveWindowAtts.subWindowAtts.win7.layer = 0
SaveWindowAtts.subWindowAtts.win7.transparency = 0
SaveWindowAtts.subWindowAtts.win7.omitWindow = 0
SaveWindowAtts.subWindowAtts.win8.position = (0, 0)
SaveWindowAtts.subWindowAtts.win8.size = (128, 128)
SaveWindowAtts.subWindowAtts.win8.layer = 0
SaveWindowAtts.subWindowAtts.win8.transparency = 0
SaveWindowAtts.subWindowAtts.win8.omitWindow = 0
SaveWindowAtts.subWindowAtts.win9.position = (0, 0)
SaveWindowAtts.subWindowAtts.win9.size = (128, 128)
SaveWindowAtts.subWindowAtts.win9.layer = 0
SaveWindowAtts.subWindowAtts.win9.transparency = 0
SaveWindowAtts.subWindowAtts.win9.omitWindow = 0
SaveWindowAtts.subWindowAtts.win10.position = (0, 0)
SaveWindowAtts.subWindowAtts.win10.size = (128, 128)
SaveWindowAtts.subWindowAtts.win10.layer = 0
SaveWindowAtts.subWindowAtts.win10.transparency = 0
SaveWindowAtts.subWindowAtts.win10.omitWindow = 0
SaveWindowAtts.subWindowAtts.win11.position = (0, 0)
SaveWindowAtts.subWindowAtts.win11.size = (128, 128)
SaveWindowAtts.subWindowAtts.win11.layer = 0
SaveWindowAtts.subWindowAtts.win11.transparency = 0
SaveWindowAtts.subWindowAtts.win11.omitWindow = 0
SaveWindowAtts.subWindowAtts.win12.position = (0, 0)
SaveWindowAtts.subWindowAtts.win12.size = (128, 128)
SaveWindowAtts.subWindowAtts.win12.layer = 0
SaveWindowAtts.subWindowAtts.win12.transparency = 0
SaveWindowAtts.subWindowAtts.win12.omitWindow = 0
SaveWindowAtts.subWindowAtts.win13.position = (0, 0)
SaveWindowAtts.subWindowAtts.win13.size = (128, 128)
SaveWindowAtts.subWindowAtts.win13.layer = 0
SaveWindowAtts.subWindowAtts.win13.transparency = 0
SaveWindowAtts.subWindowAtts.win13.omitWindow = 0
SaveWindowAtts.subWindowAtts.win14.position = (0, 0)
SaveWindowAtts.subWindowAtts.win14.size = (128, 128)
SaveWindowAtts.subWindowAtts.win14.layer = 0
SaveWindowAtts.subWindowAtts.win14.transparency = 0
SaveWindowAtts.subWindowAtts.win14.omitWindow = 0
SaveWindowAtts.subWindowAtts.win15.position = (0, 0)
SaveWindowAtts.subWindowAtts.win15.size = (128, 128)
SaveWindowAtts.subWindowAtts.win15.layer = 0
SaveWindowAtts.subWindowAtts.win15.transparency = 0
SaveWindowAtts.subWindowAtts.win15.omitWindow = 0
SaveWindowAtts.subWindowAtts.win16.position = (0, 0)
SaveWindowAtts.subWindowAtts.win16.size = (128, 128)
SaveWindowAtts.subWindowAtts.win16.layer = 0
SaveWindowAtts.subWindowAtts.win16.transparency = 0
SaveWindowAtts.subWindowAtts.win16.omitWindow = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()