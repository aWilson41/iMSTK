/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details. 
*/

%module(directors="1") Utils
#pragma SWIG nowarn=302,314,317,401,476,501,503,505,516,844,
%{
/* Common */
#include "imstkMacros.h"
#include "imstkMath.h"
#include "imstkAbstractDataArray.h"
#include "imstkDataArray.h"
#include "imstkVecDataArray.h"
#include "imstkLogger.h"
#ifdef iMSTK_SYNCHRONOUS_LOGGING
#include "imstkLoggerSynchronous.h"
#else
#include "imstkLoggerG3.h"
#endif
#include "imstkModule.h"
#include "imstkModuleDriver.h"
#include "imstkColor.h"
#include "imstkEventObject.h"
#include "imstkTypes.h"
#include "imstkFactory.h"

/*
 * DataStructures
 */
#include "imstkNeighborSearch.h"

/* 
 * Geometry 
 */
#include "imstkGeometry.h"
#include "imstkGeometryUtilities.h"
#include "imstkPointSet.h"
#include "imstkAbstractCellMesh.h"
#include "imstkCellMesh.h"
#include "imstkSurfaceMesh.h"
#include "imstkLineMesh.h"
#include "imstkImageData.h"
#include "imstkVolumetricMesh.h"
#include "imstkTetrahedralMesh.h"
#include "imstkHexahedralMesh.h"
#include "imstkImplicitGeometry.h"
#include "imstkAnalyticalGeometry.h"
#include "imstkCompositeImplicitGeometry.h"
#include "imstkPlane.h"
#include "imstkSphere.h"
#include "imstkOrientedBox.h"
#include "imstkCapsule.h"
#include "imstkCylinder.h"
#include "imstkGeometryUtilities.h"
#include "imstkSignedDistanceField.h"
#include "imstkImplicitFunctionFiniteDifferenceFunctor.h"

/*
 * GeometryMappers
 */
#include "imstkGeometryMap.h"
#include "imstkPointwiseMap.h"
#include "imstkPointToTetMap.h"

/*
 * Filter
 */
#include "imstkGeometryAlgorithm.h"
#include "imstkImplicitGeometryToImageData.h"
#include "imstkQuadricDecimate.h"
#include "imstkSelectEnclosedPoints.h"
#include "imstkSurfaceMeshFlyingEdges.h"
#include "imstkSurfaceMeshSmoothen.h"
#include "imstkSurfaceMeshSubdivide.h"
#include "imstkSurfaceMeshTextureProject.h"

/* 
 * MeshIO 
 */
#include "imstkMeshIO.h"

/* 
 * DynamicalModels
 */
#include "imstkAbstractDynamicalModel.h"
#include "imstkDynamicalModel.h"
#include "imstkPbdModelConfig.h"
#include "imstkPbdModel.h"
#include "imstkPbdFemConstraint.h"
#include "imstkPbdCollisionConstraint.h"
#include "imstkSphBoundaryConditions.h"
#include "imstkRigidBodyState2.h"
#include "imstkRigidBodyModel2.h"
#include "imstkSphState.h"
#include "imstkSphModel.h"

/*
 * DynamicalModelsVegaFEM
 */
#ifdef iMSTK_USE_VegaFEM
#include "imstkBackwardEuler.h"
#include "imstkFemDeformableBodyModel.h"
#include "imstkInternalForceModelTypes.h"
#include "imstkTimeIntegrator.h"
#include "imstkVectorizedState.h"
#endif

/* 
 * Rendering
 */
#include "imstkRenderMaterial.h"
#include "imstkTexture.h"

/*
 * Constraints
 */
#include "imstkPbdBody.h"
#include "imstkPbdConstraint.h"
#include "imstkRbdConstraint.h"

/*
 * ComponentModel
 */
#include "imstkEntity.h"
#include "imstkComponent.h"

/* 
 * SceneEntities
 */
#include "imstkSceneObject.h"
#include "imstkCollidingObject.h"
#include "imstkDynamicObject.h"
#include "imstkPbdConnectiveTissueConstraintGenerator.h"
#include "imstkPbdObject.h"
#include "imstkVisualModel.h"
#include "imstkCamera.h"
#include "imstkLight.h"
#include "imstkDirectionalLight.h"
#include "imstkRigidObject2.h"
#include "imstkSphObject.h"
#ifdef iMSTK_USE_VegaFEM
#include "imstkFeDeformableObject.h"
#endif

/*
 * CollisionDetection
 */
#include "imstkClosedSurfaceMeshToMeshCD.h"
#include "imstkCollisionData.h"
#include "imstkCollisionDetectionAlgorithm.h"
#include "imstkBidirectionalPlaneToSphereCD.h"
#include "imstkCollisionDetectionAlgorithm.h"
#include "imstkCollisionUtils.h"
#include "imstkImplicitGeometryToPointSetCCD.h"
#include "imstkImplicitGeometryToPointSetCD.h"
#include "imstkPointSetToCapsuleCD.h"
#include "imstkPointSetToOrientedBoxCD.h"
#include "imstkPointSetToPlaneCD.h"
#include "imstkPointSetToSphereCD.h"
#include "imstkSphereToCylinderCD.h"
#include "imstkSphereToSphereCD.h"
#include "imstkSurfaceMeshToCapsuleCD.h"
#include "imstkSurfaceMeshToSphereCD.h"
#include "imstkSurfaceMeshToSurfaceMeshCD.h"
#include "imstkTetraToLineMeshCD.h"
#include "imstkTetraToPointSetCD.h"
#include "imstkUnidirectionalPlaneToSphereCD.h"

/*
 * CollisionHandling
 */
#include "imstkCollisionHandling.h"
#include "imstkRigidBodyCH.h"

/*
 * Controller
 */
#include "imstkDeviceControl.h"
#include "imstkMouseControl.h"
#include "imstkKeyboardControl.h"
#include "imstkTrackingDeviceControl.h"
#include "imstkSceneObjectController.h"
#include "imstkRigidObjectController.h"
#include "imstkPbdObjectController.h"

/*
 * Needle
 */
#include "imstkPuncture.h"
#include "imstkNeedle.h"
#include "imstkStraightNeedle.h"
#include "imstkArcNeedle.h"
#include "imstkPuncturable.h"

/*
 * Scene
 */
#include "imstkScene.h"
#include "imstkCollisionInteraction.h"
#include "imstkRigidObjectCollision.h"
#include "imstkPbdObjectCutting.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdRigidObjectCollision.h"
#include "imstkPbdRigidObjectGrasping.h"
#include "imstkSphObjectCollision.h"

/*
 * SimulationManager
 */
#include "imstkModule.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkMouseSceneControl.h"
#include "imstkKeyboardSceneControl.h"

/*
 * ViewerCore
 */
#include "imstkViewer.h"

#ifdef iMSTK_USE_RENDERING_VTK
/*
 * ViewerVTK
 */
#include "imstkAbstractVTKViewer.h"
#include "imstkVTKViewer.h"
#endif

/*
 * Devices
 */
#include "imstkDeviceClient.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkMouseDeviceClient.h"

#ifdef iMSTK_USE_HAPLY
#include "imstkHaplyDeviceManager.h"
#include "imstkHaplyDeviceClient.h"
#endif

#ifdef iMSTK_USE_OpenHaptics
#include "imstkOpenHapticDeviceManager.h"
#include "imstkOpenHapticDeviceClient.h"
#endif

#ifdef iMSTK_USE_VRPN
#include "imstkVRPNDeviceManager.h"
#include "imstkVRPNDeviceClient.h"
#endif

#include "imstkDeviceManager.h"
#include "imstkDeviceManagerFactory.h"

%} /* end of module */

/*
 * stl
 */
%include <stdint.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_pair.i>
namespace std
{
  %template(VectorInt) vector<int>;
  %template(VectorSizet) vector<std::size_t>;
  %template(VectorDouble) vector<double>;
  %template(VectorCollisionElement) vector<imstk::CollisionElement>;
  %template(VectorPbdBody) vector<imstk::PbdBody>;
}

%include <std_except.i>
%include <exception.i>


%include "shared_ptr_instantiation.i"
%include "weak_ptr.i"
%include "ignored.i"
%include "modifiers.i"
%include "type_cast.i"
%include "std_function.i"
%include "callback.i"
%include "except.i"

/* rename these operators to "compute" due to lack of operator overloading */
%rename(compute) imstk::ImplicitFunctionGradient::operator();
%rename(compute) imstk::ImplicitFunctionCentralGradient::operator();

/*
 * * * * * * * * * * * * * * * *
 * list C/C++ declarations
 * * * * * * * * * * * * * * * * *
 */
/*
 * Common
 */
%include "common.i"

/*
 * DataStructures
 */
%include "../../Source/DataStructures/imstkNeighborSearch.h"

/*
 * Geometry
 */
%include "../../Source/Geometry/imstkGeometry.h";
%include "../../Source/Geometry/Mesh/imstkPointSet.h"
%include "../../Source/Geometry/Mesh/imstkImageData.h"
%include "../../Source/Geometry/Mesh/imstkAbstractCellMesh.h"
%include "../../Source/Geometry/Mesh/imstkCellMesh.h"
%template(CellMesh2) imstk::CellMesh<2>;
%template(CellMesh3) imstk::CellMesh<3>;
%template(CellMesh4) imstk::CellMesh<4>;
%template(CellMesh8) imstk::CellMesh<8>;
%include "../../Source/Geometry/Mesh/imstkLineMesh.h"
%include "../../Source/Geometry/Mesh/imstkSurfaceMesh.h"
%include "../../Source/Geometry/Mesh/imstkVolumetricMesh.h"
%template(VolumetricMesh4) imstk::VolumetricMesh<4>;
%template(VolumetricMesh8) imstk::VolumetricMesh<8>;
%include "../../Source/Geometry/Mesh/imstkTetrahedralMesh.h"
%include "../../Source/Geometry/Mesh/imstkHexahedralMesh.h"
%include "../../Source/Geometry/Implicit/imstkImplicitGeometry.h"
%include "../../Source/Geometry/Implicit/imstkCompositeImplicitGeometry.h"
%include "../../Source/Geometry/Analytic/imstkAnalyticalGeometry.h"
%include "../../Source/Geometry/Analytic/imstkPlane.h"
%include "../../Source/Geometry/Analytic/imstkSphere.h"
%include "../../Source/Geometry/Analytic/imstkOrientedBox.h"
%include "../../Source/Geometry/Analytic/imstkCapsule.h"
%include "../../Source/Geometry/Analytic/imstkCylinder.h"
%include "../../Source/Geometry/imstkGeometryUtilities.h"
%include "../../Source/Geometry/Implicit/imstkSignedDistanceField.h"
%include "../../Source/Geometry/Implicit/imstkImplicitFunctionFiniteDifferenceFunctor.h"

/*
 * GeometryMap
 */
%include "../../Source/GeometryMappers/imstkGeometryMap.h"
%include "../../Source/GeometryMappers/imstkPointwiseMap.h"
%include "../../Source/GeometryMappers/imstkPointToTetMap.h"

/*
 * FilteringCore
 */
%include "../../Source/FilteringCore/imstkGeometryAlgorithm.h"

/*
 * Filtering
 */
%include "../../Source/Filtering/imstkImplicitGeometryToImageData.h"
%include "../../Source/Filtering/imstkQuadricDecimate.h"
%include "../../Source/Filtering/imstkSelectEnclosedPoints.h"
%include "../../Source/Filtering/imstkSurfaceMeshFlyingEdges.h"
%include "../../Source/Filtering/imstkSurfaceMeshSmoothen.h"
%include "../../Source/Filtering/imstkSurfaceMeshSubdivide.h"
%include "../../Source/Filtering/imstkSurfaceMeshTextureProject.h"

/*
 * MeshIO
 */
%include "../../Source/MeshIO/imstkMeshIO.h";
%template(readImageData) imstk::MeshIO::read<imstk::ImageData>;
%template(readPointSet) imstk::MeshIO::read<imstk::PointSet>;
%template(readSurfaceMesh) imstk::MeshIO::read<imstk::SurfaceMesh>;
%template(readTetrahedralMesh) imstk::MeshIO::read<imstk::TetrahedralMesh>;

/*
 * Constraint
 */
%include "../../Source/Constraint/PbdConstraints/imstkPbdBody.h"
%include "../../Source/Constraint/PbdConstraints/imstkPbdConstraint.h"
%include "../../Source/Constraint/PbdConstraints/imstkPbdCollisionConstraint.h"
%include "../../Source/Constraint/PbdConstraints/imstkPbdFemConstraint.h"
%include "../../Source/Constraint/RigidBodyConstraints/imstkRbdConstraint.h"

/*
 * DynamicalModels
 */
%include "../../Source/DynamicalModels/ObjectModels/imstkAbstractDynamicalModel.h"
%include "../../Source/DynamicalModels/ObjectModels/imstkDynamicalModel.h"
%include "../../Source/DynamicalModels/ObjectModels/imstkPbdModelConfig.h"
%include "../../Source/DynamicalModels/ObjectModels/imstkPbdModel.h"
%include "../../Source/DynamicalModels/ObjectModels/imstkSphBoundaryConditions.h"
%include "../../Source/DynamicalModels/ObjectStates/imstkRigidBodyState2.h"
%template(DynamicalModelRigidBodyState2) imstk::DynamicalModel<imstk::RigidBodyState2>;
%include "../../Source/DynamicalModels/ObjectModels/imstkRigidBodyModel2.h"
%include "../../Source/DynamicalModels/ObjectStates/imstkSphState.h"
%template(DynamicalModelSphState) imstk::DynamicalModel<imstk::SphState>;
%include "../../Source/DynamicalModels/ObjectModels/imstkSphModel.h"

/*
 * DynamicalModelsVegaFEM
 */
#ifdef iMSTK_USE_VegaFEM
%include "../../Source/DynamicalModelsVegaFEM/imstkVectorizedState.h"
%template(DynamicalModelFeDeformBodyState) imstk::DynamicalModel<imstk::FeDeformBodyState>;
%include "../../Source/DynamicalModelsVegaFEM/InternalForceModel/imstkInternalForceModelTypes.h"
%include "../../Source/DynamicalModelsVegaFEM/imstkFemDeformableBodyModel.h"
%include "../../Source/DynamicalModelsVegaFEM/TimeIntegrators/imstkTimeIntegrator.h"
%include "../../Source/DynamicalModelsVegaFEM/TimeIntegrators/imstkBackwardEuler.h"
#endif

/* 
 * Rendering 
 */
%include "../../Source/Materials/imstkRenderMaterial.h";
%include "../../Source/Materials/imstkTexture.h";

/*
 * ComponentModel
 */
%include "../../Source/ComponentModel/imstkEntity.h"
%include "../../Source/ComponentModel/imstkComponent.h"

/*
 * SceneEntities
 */
%include "../../Source/SceneEntities/Components/imstkVisualModel.h";
%include "../../Source/SceneEntities/Objects/imstkSceneObject.h";
%include "../../Source/SceneEntities/Objects/imstkCollidingObject.h";
%include "../../Source/SceneEntities/Objects/imstkDynamicObject.h";
%include "../../SceneEntities/Objects/imstkPbdConnectiveTissueConstraintGenerator.h";
%include "../../Source/SceneEntities/Objects/imstkPbdObject.h";
%include "../../Source/SceneEntities/Objects/imstkRigidObject2.h";
%include "../../Source/SceneEntities/Objects/imstkSphObject.h";
%include "../../Source/SceneEntities/Camera/imstkCamera.h";
%include "../../Source/SceneEntities/Lights/imstkLight.h";
%include "../../Source/SceneEntities/Lights/imstkDirectionalLight.h";
#ifdef iMSTK_USE_VegaFEM
%include "../../Source/SceneEntities/Objects/imstkFeDeformableObject.h";
#endif

/*
 * CollisionDetection
 */
%include "../../Source/CollisionDetection/imstkCollisionData.h"
%include "../../Source/CollisionDetection/imstkCollisionDetectionAlgorithm.h"
%include "../../Source/CollisionDetection/imstkCollisionUtils.h"

%include "../../Source/CollisionDetection/CollisionDetection/imstkBidirectionalPlaneToSphereCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkClosedSurfaceMeshToMeshCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkImplicitGeometryToPointSetCCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkImplicitGeometryToPointSetCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkPointSetToCapsuleCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkPointSetToOrientedBoxCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkPointSetToPlaneCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkPointSetToSphereCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkSphereToCylinderCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkSphereToSphereCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkSurfaceMeshToCapsuleCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkSurfaceMeshToSphereCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkSurfaceMeshToSurfaceMeshCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkTetraToLineMeshCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkTetraToPointSetCD.h"
%include "../../Source/CollisionDetection/CollisionDetection/imstkUnidirectionalPlaneToSphereCD.h"

/*
 * CollisionHandling
 */ 
%include "../../Source/CollisionHandling/imstkCollisionHandling.h";
%include "../../Source/CollisionHandling/imstkRigidBodyCH.h";

/* 
 * Controllers
 */
%include "../../Controllers/imstkDeviceControl.h"
%include "../../Controllers/imstkMouseControl.h"
%include "../../Controllers/imstkKeyboardControl.h"
%include "../../Controllers/imstkTrackingDeviceControl.h"
%include "../../Controllers/imstkSceneObjectController.h"
%include "../../Controllers/imstkRigidObjectController.h"
%include "../../Controllers/imstkPbdObjectController.h"

/*
 * Needle
 */
%include "../../Source/Needle/imstkPuncture.h"
%include "../../Source/Needle/imstkNeedle.h"
%include "../../Source/Needle/imstkStraightNeedle.h"
%include "../../Source/Needle/imstkArcNeedle.h"
%include "../../Source/Needle/imstkPuncturable.h"

/* 
 * Scene
 */
%include "../../Source/Scene/imstkScene.h";
%include "../../Source/Scene/imstkCollisionInteraction.h"
%include "../../Source/Scene/imstkRigidObjectCollision.h"
%include "../../Source/Scene/imstkPbdObjectCutting.h"
%include "../../Source/Scene/imstkPbdObjectGrasping.h"
%include "../../Source/Scene/imstkPbdObjectCollision.h"
%include "../../Source/Scene/imstkPbdRigidObjectCollision.h"
%include "../../Source/Scene/imstkPbdRigidObjectGrasping.h"
%include "../../Source/Scene/imstkSphObjectCollision.h"

/*
 * ViewerCore
 */
%include "../../Source/ViewerCore/imstkViewer.h";

#ifdef iMSTK_USE_RENDERING_VTK
/*
 * ViewerVTK
 */
%include "../../Source/ViewerVTK/imstkAbstractVTKViewer.h";
%include "../../Source/ViewerVTK/imstkVTKViewer.h";
#endif

/*
 * SimulationManager
 */
%include "../../Source/SimulationManager/imstkSceneManager.h"
%include "../../Source/SimulationManager/imstkSimulationManager.h"
%include "../../Source/SimulationManager/imstkMouseSceneControl.h"
%include "../../Source/SimulationManager/imstkKeyboardSceneControl.h"

/*
 * Devices
 */
%include "../../Source/Devices/imstkDeviceClient.h"
%include "../../Source/Devices/imstkKeyboardDeviceClient.h"
%include "../../Source/Devices/imstkMouseDeviceClient.h"
%include "../../Source/Devices/imstkDeviceManager.h"
%include "../../Source/Devices/imstkDeviceManagerFactory.h"

/*
 * The Superclass static functions don't seem to get exposed, 
 * this adds a "local" static function that just invokes the builtin
 * contains function
 */
%extend imstk::DeviceManagerFactory {
	static bool contains(const std::string& val) {
		return imstk::DeviceManagerFactory::contains(val);
	}
}

#ifdef iMSTK_USE_HAPLY
	%include "../../Source/Devices/imstkHaplyDeviceManager.h"
	%include "../../Source/Devices/imstkHaplyDeviceClient.h"
#endif

#ifdef iMSTK_USE_OpenHaptics
	#define HDCALLBACK
	%include "../../Source/Devices/imstkOpenHapticDeviceManager.h"
	%include "../../Source/Devices/imstkOpenHapticDeviceClient.h"
#endif

#ifdef iMSTK_USE_VRPN
	// The static calls in DeviceClient are getting ignored anyway define these
	// Rather than dealing with the correct includes for VRPN
	#define VRPN_CALLBACK
	#define _vrpn_TRACKERCB void*
	#define _vrpn_TRACKERVELCB void*
	#define _vrpn_ANALOGCB void*
	#define _vrpn_BUTTONCB void*

	// Swig things that VRPNDeviceManager is abstract, that will
	// mean no constructors are created
	%feature("notabstract") imstk::VRPNDeviceManager;

	%include "../../Source/Devices/imstkVRPNDeviceManager.h"
	%include "../../Source/Devices/imstkVRPNDeviceClient.h"
#endif
