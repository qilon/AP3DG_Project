#include "PolyMesh.h"
//=============================================================================
int loadMesh(PolyMesh& mesh, const char* filename, bool read_vertex_colors)
{
	OpenMesh::IO::Options ropt, wopt;

	if (read_vertex_colors)
	{
		mesh.request_vertex_colors();
		ropt += OpenMesh::IO::Options::VertexColor;
	}

	std::cout << "Reading file... " << std::endl;
	if (!OpenMesh::IO::read_mesh(mesh, filename))
	{
		std::cout << "Could not read file: " << filename << std::endl << std::endl;

		return -1;
	}

	return 0;
}
//=============================================================================