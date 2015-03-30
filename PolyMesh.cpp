#include "PolyMesh.h"
//=============================================================================
void setColors(PolyMesh* mesh, Eigen::ArrayXf values, float min_value,
	float max_value)
{
	std::cout << "Coloring mesh..." << std::endl;

	float mean = 0.f;
	float std_coeff_left = 1.f;
	float std_coeff_right = 1.f;
	float range;

	std::vector<float> l_values(values.size());
	for (int i = 0; i < values.size(); i++)
	{
		l_values[i] = values(i);
	}
	sort(l_values.begin(), l_values.end());
	float value1_10 = l_values[ceil(values.size() / 10)];
	float value9_10 = l_values[9 * ceil(values.size() / 10)];
	range = value9_10 - value1_10;

	if (min_value > 0.f)
	{
		min_value = -min_value;
	}

	range = std::min(max_value, std::max(min_value, range));

	float new_min_value = -range;
	float new_max_value = range;

	mesh->request_vertex_colors();
	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	int i_v = 0;
	for (v_it = mesh->vertices_begin(); v_it != v_end; v_it++)
	{
		unsigned char red = 0;
		unsigned char blue = 0;
		unsigned char green = 0;

		if (values(i_v) >= new_min_value && new_max_value >= values(i_v))
		{
			float ratio = 2 * (values(i_v) - new_min_value) / (new_max_value - new_min_value);
			red = unsigned char(std::max(0.0f, 255 * (ratio - 1)));
			green = unsigned char(std::max(0.0f, 255 * (1 - ratio)));
			blue = 255 - green - red;
		}
		else
		{
			if (values(i_v) < new_min_value)
			{
				red = 0;
				green = 255;
				blue = 0;
			}
			else //new_max_value < values(i_v)
			{
				red = 255;
				green = 0;
				blue = 0;
			}
		}

		PolyMesh::Color color(red, green, blue);

		mesh->set_color(*v_it, color);

		i_v++;
	}
}
//=============================================================================
Eigen::ArrayXf meanCurvatureUL(PolyMesh* mesh, PolyMesh* colored_mesh)
{
	std::cout << "Calculating mean curvature..." << std::endl;

	int numVertices = mesh->n_vertices();
	Eigen::ArrayXf H(numVertices);
	bool color_mesh = colored_mesh != nullptr;

	float minH = std::numeric_limits<float>::max();
	float maxH = std::numeric_limits<float>::min();
	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	int i_v = 0;

	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		Eigen::Vector3f e_p;
		e_p << p[0], p[1], p[2];
		Eigen::Vector3f e_p_laplacian = Eigen::Vector3f::Zero();
		int p_valence = 0;


		PolyMesh::VertexVertexIter vv_it;
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			PolyMesh::Point p_neighbour = mesh->point(*vv_it);
			Eigen::Vector3f e_p_neighbour;
			e_p_neighbour << p_neighbour[0], p_neighbour[1], p_neighbour[2];
			e_p_laplacian += (e_p_neighbour - e_p);
			p_valence++;
		}

		if (p_valence > 0)
		{
			e_p_laplacian = e_p_laplacian / float(p_valence); 
			H(i_v) = 0.5f * ( e_p_laplacian.norm() );
			minH = std::min(minH, H(i_v));
			maxH = std::max(maxH, H(i_v));
		}
		else
		{
			H(i_v) = 0.0f;
		}

		//std::cout << i_v << ": " << H[i_v] << std::endl;

		i_v++;
	}

	if (color_mesh)
	{
		setColors(colored_mesh, H, minH, maxH);
	}

	return H;
}
//=============================================================================
Eigen::ArrayXf meanCurvatureLB(PolyMesh* mesh, PolyMesh* colored_mesh)
{
	std::cout << "Calculating mean curvature..." << std::endl;

	int numVertices = mesh->n_vertices();
	Eigen::ArrayXf H(numVertices);
	bool color_mesh = colored_mesh != nullptr;

	float minH = std::numeric_limits<float>::max();
	float maxH = std::numeric_limits<float>::min();
	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	int i_v = 0;
	for (v_it = mesh->vertices_sbegin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		Eigen::Vector3f e_p;
		e_p << p[0], p[1], p[2];
		Eigen::Vector3f e_p_laplacian = Eigen::Vector3f::Zero();
		int p_valence = 0.0f;
		float area = 0.0f;

		std::vector<Eigen::Vector3f> l_e_p_neighbours;
		PolyMesh::VertexVertexIter vv_it;
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			PolyMesh::Point p_neighbour = mesh->point(*vv_it);
			Eigen::Vector3f e_p_neighbour;
			e_p_neighbour << p_neighbour[0], p_neighbour[1], p_neighbour[2];
			l_e_p_neighbours.push_back(e_p_neighbour); 
		}
		p_valence = l_e_p_neighbours.size();

		for (int i_neigh = 0; i_neigh < l_e_p_neighbours.size(); i_neigh++)
		{
			int i_prev = (i_neigh - 1)<0 ? l_e_p_neighbours.size() - 1 : i_neigh - 1;
			int i_next = (i_neigh + 1) % l_e_p_neighbours.size();

			Eigen::Vector3f e_prev_p = l_e_p_neighbours.at(i_prev);
			Eigen::Vector3f e_curr_p = l_e_p_neighbours.at(i_neigh);
			Eigen::Vector3f e_next_p = l_e_p_neighbours.at(i_next);

			Eigen::Vector3f e_prev_v_1 = (e_p - e_prev_p).normalized();
			Eigen::Vector3f e_prev_v_2 = (e_curr_p - e_prev_p).normalized();

			Eigen::Vector3f e_next_v_1 = (e_p - e_next_p).normalized();
			Eigen::Vector3f e_next_v_2 = (e_curr_p - e_next_p).normalized();

			float dot_prev = e_prev_v_1.dot(e_prev_v_2);
			float cot_a = dot_prev / (e_prev_v_1 - dot_prev*e_prev_v_2).norm();

			float dot_next = e_next_v_1.dot(e_next_v_2);
			float cot_b = dot_next / (e_next_v_1 - dot_next*e_next_v_2).norm();

			e_p_laplacian += (cot_a + cot_b) * (e_curr_p - e_p);

			area += (e_next_p - e_p).cross(e_curr_p - e_p).norm() / 6;
		}

		if (p_valence > 0 && area > 0.0f)
		{
			e_p_laplacian = e_p_laplacian / (2 * area);
			H(i_v) = 0.5f * ( e_p_laplacian.norm() );
			minH = std::min(minH, H(i_v));
			maxH = std::max(maxH, H(i_v));
		}
		else
		{
			H(i_v) = 0.0f;
		}

		i_v++;
	}

	if (color_mesh)
	{
		setColors(colored_mesh, H, minH, maxH);
	}

	return H;
}
//=============================================================================
Eigen::ArrayXf meanCurvature(PolyMesh* mesh, PolyMesh* colored_mesh,
	int laplacian_type)
{
	if (laplacian_type == UL)
	{
		return meanCurvatureUL(mesh, colored_mesh);
	}
	else //LB
	{
		return meanCurvatureLB(mesh, colored_mesh);
	}
}
//=============================================================================
Eigen::ArrayXf gaussianCurvature(PolyMesh* mesh, PolyMesh* colored_mesh)
{
	std::cout << "Calculating gaussian curvature..." << std::endl;

	int numVertices = mesh->n_vertices();
	Eigen::ArrayXf K(numVertices);
	bool color_mesh = colored_mesh != nullptr;

	float minK = std::numeric_limits<float>::max();
	float maxK = std::numeric_limits<float>::min();
	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	int i_v = 0;
	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		Eigen::Vector3f e_p;
		e_p << p[0], p[1], p[2];
		float angle = 2.0f * M_PI;
		float area = 0.0f;

		std::vector<Eigen::Vector3f> l_e_p_neighbours;
		PolyMesh::VertexVertexIter vv_it;
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			PolyMesh::Point p_neighbour = mesh->point(*vv_it);
			Eigen::Vector3f e_p_neighbour;
			e_p_neighbour << p_neighbour[0], p_neighbour[1], p_neighbour[2];
			l_e_p_neighbours.push_back(e_p_neighbour);
		}

		for (int i_neigh = 0; i_neigh < l_e_p_neighbours.size(); i_neigh++)
		{
			Eigen::Vector3f e_curr_p = l_e_p_neighbours.at(i_neigh);

			int i_next = (i_neigh + 1) % l_e_p_neighbours.size();
			Eigen::Vector3f e_next_p = l_e_p_neighbours.at(i_next);

			Eigen::Vector3f e_v_curr = e_curr_p - e_p;
			Eigen::Vector3f e_v_next = e_next_p - e_p;
			angle -= acos(e_v_curr.normalized().dot(e_v_next.normalized()));

			area += e_v_next.cross(e_v_curr).norm() / 6;
		}

		if (area > 0.0f)
		{
			K(i_v) = angle / area;
			minK = std::min(minK, K(i_v));
			maxK = std::max(maxK, K(i_v));
		}
		else
		{
			K(i_v) = 0.0f;
		}

		//std::cout << i_v << ": " << K(i_v) << std::endl;

		i_v++;
	}

	if (color_mesh)
	{
		setColors(colored_mesh, K, minK, maxK);
	}


	return K;
}
//=============================================================================
void principalK1K2(PolyMesh* mesh, Eigen::ArrayXf H,
	Eigen::ArrayXf K, Eigen::ArrayXf *k1, Eigen::ArrayXf *k2,
	PolyMesh* k1_mesh, PolyMesh* k2_mesh)
{
	std::cout << "Calculating principal curvatures..." << std::endl;

	Eigen::ArrayXf zeros = Eigen::ArrayXf::Zero(H.size());

	Eigen::ArrayXf diffH2K = (H.array().pow(2) - K).max(zeros);
	Eigen::ArrayXf srqtH2K = (diffH2K).pow(0.5);
	*k1 = H + srqtH2K;
	*k2 = H - srqtH2K;

	if (k1_mesh != nullptr)
	{
		setColors(k1_mesh, *k1, k1->minCoeff(), k1->maxCoeff());
	}

	if (k2_mesh != nullptr)
	{
		setColors(k2_mesh, *k2, k2->minCoeff(), k2->maxCoeff());
	}
}
//=============================================================================
PolyMesh explicitSmoothingUL(PolyMesh* mesh, float diffusion_const, 
	PolyMesh::Color color)
{
	std::cout << "Explicit smoothing... " << std::endl;
	
	PolyMesh smoothed_mesh = *mesh;
	smoothed_mesh.request_vertex_colors();

	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	PolyMesh::VertexIter v_it_smooth = smoothed_mesh.vertices_begin();
	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		Eigen::Vector3f e_p;
		e_p << p[0], p[1], p[2];
		Eigen::Vector3f e_p_laplacian = Eigen::Vector3f::Zero();
		int p_valence = 0;

		PolyMesh::VertexVertexIter vv_it;
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			PolyMesh::Point p_neighbour = mesh->point(*vv_it);
			Eigen::Vector3f e_p_neighbour;
			e_p_neighbour << p_neighbour[0], p_neighbour[1], p_neighbour[2];
			e_p_laplacian += (e_p_neighbour - e_p);
			p_valence++;
		}

		if (p_valence > 0) {
			e_p_laplacian = e_p_laplacian / float(p_valence);

			PolyMesh::Point smoothed_point = smoothed_mesh.point(*v_it_smooth);

			Eigen::Vector3f e_smoothed_point;
			e_smoothed_point << smoothed_point[0], smoothed_point[1], smoothed_point[2];

			e_smoothed_point += diffusion_const * e_p_laplacian;

			smoothed_mesh.set_point(*v_it_smooth, 
				PolyMesh::Point(e_smoothed_point(0), e_smoothed_point(1), 
				e_smoothed_point(2)));
		}

		smoothed_mesh.set_color(*v_it_smooth, color);

		++v_it_smooth;
	}

	return smoothed_mesh;
}
//=============================================================================
PolyMesh explicitSmoothingLB(PolyMesh* mesh, float diffusion_const,
	PolyMesh::Color color)
{
	std::cout << "Explicit smoothing... " << std::endl;

	PolyMesh smoothed_mesh = *mesh;
	smoothed_mesh.request_vertex_colors();

	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	PolyMesh::VertexIter v_it_smooth = smoothed_mesh.vertices_begin();
	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		Eigen::Vector3f e_p;
		e_p << p[0], p[1], p[2];
		Eigen::Vector3f e_p_laplacian = Eigen::Vector3f::Zero();
		int p_valence = 0.0f;
		float area = 0.0f;

		std::vector<Eigen::Vector3f> l_e_p_neighbours;
		PolyMesh::VertexVertexIter vv_it;
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			PolyMesh::Point p_neighbour = mesh->point(*vv_it);
			Eigen::Vector3f e_p_neighbour;
			e_p_neighbour << p_neighbour[0], p_neighbour[1], p_neighbour[2];
			l_e_p_neighbours.push_back(e_p_neighbour);
		}
		p_valence = l_e_p_neighbours.size();

		for (int i_neigh = 0; i_neigh < l_e_p_neighbours.size(); i_neigh++)
		{
			int i_prev = (i_neigh - 1)<0 ? l_e_p_neighbours.size() - 1 : i_neigh - 1;
			int i_next = (i_neigh + 1) % l_e_p_neighbours.size();

			Eigen::Vector3f e_prev_p = l_e_p_neighbours.at(i_prev);
			Eigen::Vector3f e_curr_p = l_e_p_neighbours.at(i_neigh);
			Eigen::Vector3f e_next_p = l_e_p_neighbours.at(i_next);

			Eigen::Vector3f e_prev_v_1 = (e_p - e_prev_p).normalized();
			Eigen::Vector3f e_prev_v_2 = (e_curr_p - e_prev_p).normalized();

			Eigen::Vector3f e_next_v_1 = (e_p - e_next_p).normalized();
			Eigen::Vector3f e_next_v_2 = (e_curr_p - e_next_p).normalized();

			float dot_prev = e_prev_v_1.dot(e_prev_v_2);
			float cot_a = dot_prev / (e_prev_v_1 - dot_prev*e_prev_v_2).norm();

			float dot_next = e_next_v_1.dot(e_next_v_2);
			float cot_b = dot_next / (e_next_v_1 - dot_next*e_next_v_2).norm();

			e_p_laplacian += (cot_a + cot_b) * (e_curr_p - e_p);

			area += (e_next_p - e_p).cross(e_curr_p - e_p).norm() / 6;
		}

		if (p_valence > 0 && area > 0.0f)
		{
			e_p_laplacian = e_p_laplacian / (2 * area);

			PolyMesh::Point smoothed_point = smoothed_mesh.point(*v_it_smooth);

			Eigen::Vector3f e_smoothed_point;
			e_smoothed_point << smoothed_point[0], smoothed_point[1], smoothed_point[2];

			e_smoothed_point += diffusion_const * e_p_laplacian;

			smoothed_mesh.set_point(*v_it_smooth,
				PolyMesh::Point(e_smoothed_point(0), e_smoothed_point(1),
				e_smoothed_point(2)));
		}

		smoothed_mesh.set_color(*v_it_smooth, color);

		++v_it_smooth;
	}

	return smoothed_mesh;
}
//=============================================================================
PolyMesh explicitSmoothing(PolyMesh* mesh, float diffusion_const,
	int num_iter, int laplacian_type, PolyMesh::Color color)
{
	PolyMesh smoothed_mesh = *mesh;
	if (laplacian_type == UL)
	{
		for (int i = 0; i < num_iter; i++)
		{
			smoothed_mesh = 
				explicitSmoothingUL(&smoothed_mesh, diffusion_const, color);
		}
	}
	else //LB
	{
		for (int i = 0; i < num_iter; i++)
		{
			smoothed_mesh = 
				explicitSmoothingLB(&smoothed_mesh, diffusion_const, color);
		}
	}

	return smoothed_mesh;
}
//=============================================================================
PolyMesh implicitSmoothing(PolyMesh* mesh, float diffusion_const,
	int num_iter, PolyMesh::Color color)
{
	std::cout << "Implicit smoothing... " << std::endl;

	PolyMesh smoothed_mesh = *mesh;
	smoothed_mesh.request_vertex_colors();

	int num_vertices = mesh->n_vertices();

	std::vector<Eigen::Triplet<float>> coefficients;
	Eigen::SparseMatrix<float> A(num_vertices, num_vertices);
	Eigen::VectorXf x(num_vertices), y(num_vertices), z(num_vertices);

	float p_coeff = 1.f + diffusion_const;

	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		int p_idx = v_it.handle().idx();

		coefficients.push_back(
			Eigen::Triplet<float>(p_idx, p_idx, p_coeff));

		int p_valence = 0;
		PolyMesh::VertexVertexIter vv_it;
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			p_valence++;
		}

		float neigh_coeff = -diffusion_const / float(p_valence);
		for (vv_it = mesh->vv_iter(*v_it); vv_it; ++vv_it)
		{
			int p_neigh_idx = vv_it.handle().idx();
			coefficients.push_back(
				Eigen::Triplet<float>(p_idx, p_neigh_idx, neigh_coeff));
		}
	}

	A.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;
	solver.compute(A);

	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		int p_idx = v_it.handle().idx();
		x(p_idx) = p[0];
		y(p_idx) = p[1];
		z(p_idx) = p[2];
	}

	for (int i_iter = 0; i_iter < num_iter; i_iter++)
	{
		x = solver.solve(x);
		y = solver.solve(y);
		z = solver.solve(z);
	}

	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		int p_idx = v_it.handle().idx();

		smoothed_mesh.set_point(*v_it, 
			PolyMesh::Point(x(p_idx), y(p_idx), z(p_idx)) );

		smoothed_mesh.set_color(*v_it, color);
	}

	return smoothed_mesh;
}
//=============================================================================
void addNoise(PolyMesh* mesh, float mean, float sigma, PolyMesh::Color color)
{
	mesh->request_vertex_colors();

	std::default_random_engine generator;
	std::normal_distribution<float> distribution(mean, sigma);

	PolyMesh::VertexIter v_it, v_end(mesh->vertices_end());
	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		PolyMesh::Point p = mesh->point(*v_it);
		int p_idx = v_it.handle().idx();

		float noise;
		noise = distribution(generator);
		p[0] += noise;
		noise = distribution(generator);
		p[1] += noise;
		noise = distribution(generator);
		p[2] += noise;

		mesh->set_point(*v_it, p);

		mesh->set_color(*v_it, color);
	}
}
//=============================================================================