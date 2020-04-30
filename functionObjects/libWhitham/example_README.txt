This is a POD function object directory.
The function object works with OpenFOAM-v1812.

Current version of POD is able to read only the T field on the given surface every <timeIndexPOD> iterations
and create a covariance matrix, acoeffs, timelist, basis and map of coordinates of this data after <nSnapshots> interpolations.

To let the function object work User should:
- Run "wmake" in current directory.
- include libs ("libsampledEnsightSurfaceMesh.so") to controlDict file.
- include functionPOD_Global instructions in the end of controlDict file.
- include functionPOD_* instructions int the end of controlDict file for each scalar field like:
	functions
	{
	        #include "functionUcmpts"
		#include "functionPOD_T"
		#include "functionPOD_p"
		#include "functionPOD_Ux"
		#include "functionPOD_Uy"
		#include "functionPOD_Uz"
	}
- Add interpolation surface (*.stl) in directory <casefoler>/constant/triSurface/ .
- Check settings in controlDict file.
- Double-check settings in conrolDict file.

PS:
In current case after "wmake" User copies the text of "control.txt" in the end of controlDict file,
then copies "outlet.stl" into /triSurface folder and that's all.

WARNING! Also read comments in control.txt to evade mistakes.

Good luck!
