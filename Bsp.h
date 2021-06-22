
// *****************************************************************
// BSP
// *****************************************************************
namespace BSP
{
	const real c_precision = 1e-5;
	const real c_scale = 10;

	struct triangle_t
	{
		vec3 p[3];
		bool valid = true;

		triangle_t* parent = 0;
		triangle_t* children[3] = {0};

		triangle_t(crvec p1, crvec p2, crvec p3)
		{
			p[0] = p1;
			p[1] = p2;
			p[2] = p3;
		}
		triangle_t() {}
	};
	struct tri_t
	{
		int ind;
		real d;
		vec3 p;
		vec3 n;
	};
	using convex_t = std::vector<tri_t>;
	std::vector<convex_t> convexlist;
	struct bsp_t
	{
		std::vector<triangle_t> triangles;
		tri_t cuttri;
		bsp_t* parent = 0;
		bsp_t* left = 0, * right = 0;
	};
	real check_front(crvec o, crvec n, crvec p)
	{
		if ((p - o).len() <= c_precision)
			return 0;

		real d = (p - o).dot(n);
		return d;
	}
	vector3 triangside_blend(vector3 v1, vector3 v2, real alpha)
	{
		if (alpha < c_precision || alpha > 1.0f - c_precision)
			PRINT("triangside_blend alpha = " << alpha);

		return v1 * (1 - alpha) + v2 * alpha;
	}
	int cuttriangle(const triangle_t& tri, crvec o, crvec n, triangle_t front[2], triangle_t back[2])
	{
		crvec p1 = tri.p[0];
		crvec p2 = tri.p[1];
		crvec p3 = tri.p[2];

		real a = check_front(o, n, p1);
		real b = check_front(o, n, p2);
		real c = check_front(o, n, p3);

		if (fabs(a) < c_precision && fabs(b) < c_precision && fabs(c) < c_precision)
		{
			//PRINT("!cuttriangle a=" << a << " b=" << b << " c=" << c);
			return 0;
		}

		if ((a >= 0 || fabs(a) < c_precision) && (b >= 0 || fabs(b) < c_precision) && (c >= 0 || fabs(c) < c_precision))
			return -1;

		if ((a <= 0 || fabs(a) < c_precision) && (b <= 0 || fabs(b) < c_precision) && (c <= 0 || fabs(c) < c_precision))
			return -2;
		{
			if (fabs(a) < c_precision && b > 0 && c < 0)
			{
				vec3 p23 = triangside_blend(p2, p3, b / (b - c));
				front[0] = triangle_t(p1, p2, p23);
				back[0] = triangle_t(p1, p23, p3);

				return 3;
			}
			if (fabs(a) < c_precision && b < 0 && c > 0)
			{
				vec3 p23 = triangside_blend(p2, p3, b / (b - c));
				back[0] = triangle_t(p1, p2, p23);
				front[0] = triangle_t(p1, p23, p3);
				
				return 3;
			}
		}
		{// bug fixed
			if (fabs(b) < c_precision && a > 0 && c < 0)
			{
				vec3 p13 = triangside_blend(p1, p3, a / (a - c));
				front[0] = triangle_t(p2, p13, p1);
				back[0] = triangle_t(p2, p3, p13); 
				
				return 3;
			}
			if (fabs(b) < c_precision && a < 0 && c > 0)
			{
				vec3 p13 = triangside_blend(p1, p3, a / (a - c));
				back[0] = triangle_t(p2, p13, p1);
				front[0] = triangle_t(p2, p3, p13); 
				
				return 3;
			}
		}
		{
			if (fabs(c) < c_precision && a > 0 && b < 0)
			{
				vec3 p12 = triangside_blend(p1, p2, a / (a - b));
				front[0] = triangle_t(p3, p1, p12);
				back[0] = triangle_t(p3, p12, p2);
				
				return 3;
			}
			if (fabs(c) < c_precision && a < 0 && b > 0)
			{
				vec3 p12 = triangside_blend(p1, p2, a / (a - b));
				back[0] = triangle_t(p3, p1, p12);
				front[0] = triangle_t(p3, p12, p2);

				return 3;
			}
		}
		{
			if (a > 0 && b < 0 && c < 0)
			{
				vec3 p12 = triangside_blend(p1, p2, a / (a - b));
				vec3 p13 = triangside_blend(p1, p3, a / (a - c));
				front[0] = triangle_t(p1, p12, p13);
				back[0] = triangle_t(p12, p2, p3);
				back[1] = triangle_t(p13, p12, p3);
				return 1;
			}
			if (a < 0 && b > 0 && c > 0)
			{
				vec3 p12 = triangside_blend(p1, p2, a / (a - b));
				vec3 p13 = triangside_blend(p1, p3, a / (a - c));
				back[0] = triangle_t(p1, p12, p13);
				front[0] = triangle_t(p12, p2, p3);
				front[1] = triangle_t(p13, p12, p3);
				return 2;
			}
		}
		{
			if (b > 0 && a < 0 && c < 0)
			{
				vec3 p21 = triangside_blend(p2, p1, b / (b - a));
				vec3 p23 = triangside_blend(p2, p3, b / (b - c));
				front[0] = triangle_t(p2, p23, p21);
				back[0] = triangle_t(p23, p3, p21);
				back[1] = triangle_t(p21, p3, p1);
				return 1;
			}
			if (b < 0 && a > 0 && c > 0)
			{
				vec3 p21 = triangside_blend(p2, p1, b / (b - a));
				vec3 p23 = triangside_blend(p2, p3, b / (b - c));
				back[0] = triangle_t(p2, p23, p21);
				front[0] = triangle_t(p23, p3, p21);
				front[1] = triangle_t(p21, p3, p1);
				return 2;
			}
		}
		{
			if (c > 0 && a < 0 && b < 0)
			{
				vec3 p31 = triangside_blend(p3, p1, c / (c - a));
				vec3 p32 = triangside_blend(p3, p2, c / (c - b));
				front[0] = triangle_t(p3, p31, p32);
				back[0] = triangle_t(p31, p1, p2);
				back[1] = triangle_t(p31, p2, p32);
				return 1;
			}
			if (c < 0 && a > 0 && b > 0)
			{
				vec3 p31 = triangside_blend(p3, p1, c / (c - a));
				vec3 p32 = triangside_blend(p3, p2, c / (c - b));
				back[0] = triangle_t(p3, p31, p32);
				front[0] = triangle_t(p31, p1, p2);
				front[1] = triangle_t(p31, p2, p32);
				return 2;
			}
		}

		PRINT("a=" << a << " b=" << b << " c=" << c);
		PRINT("!error");
		return 0;
	}
	void buildbsp(bsp_t* bsp, int depth)
	{
		if (depth > 1500)
		{
			PRINT("depth > 1500");
			return;
		}
		if (bsp->triangles.empty())
			return;

		bsp_t* left = new bsp_t();
		bsp->left = left;
		left->parent = bsp;

		bsp_t* right = new bsp_t();
		bsp->right = right;
		right->parent = bsp;

		vec3 o, n;
		{
			const vec3& p1 = bsp->triangles[0].p[0];
			const vec3& p2 = bsp->triangles[0].p[1];
			const vec3& p3 = bsp->triangles[0].p[2];

			o = p1;
			n = (p2 - p1).cross(p3 - p1).normcopy();

			bsp->cuttri.d = o.dot(n);
			bsp->cuttri.n = n;
			bsp->cuttri.p = o;
		}

		// div triangles
		for (int i = 1; i < bsp->triangles.size(); i++)
		{
			const vec3& p1 = bsp->triangles[i].p[0];
			const vec3& p2 = bsp->triangles[i].p[1];
			const vec3& p3 = bsp->triangles[i].p[2];

			/*if ((p1 - p2).len() < c_precision || (p3 - p2).len() < c_precision || (p1 - p3).len() < c_precision)
			{
				continue;
			}*/

			triangle_t front[2], back[2];
			triangle_t tri(p1, p2, p3);

			int ret = cuttriangle(tri, o, n, front, back);

			if (ret == 0)
			{
				continue;
			}


			if (ret == -1)
			{
				right->triangles.push_back(triangle_t(p1, p2, p3));
			}
			else if (ret == -2)
			{
				left->triangles.push_back(triangle_t(p1, p2, p3));
			}
			else
			{
				if (ret == 1)
				{
					right->triangles.push_back(front[0]);
					left->triangles.push_back(back[0]);
					left->triangles.push_back(back[1]);
				}
				else if (ret == 2)
				{
					right->triangles.push_back(front[0]);
					right->triangles.push_back(front[1]);
					left->triangles.push_back(back[0]);
				}
				else if (ret == 3)
				{
					right->triangles.push_back(front[0]);
					left->triangles.push_back(back[0]);
				}
			}
		}
		if (!left->triangles.empty())
		{
			buildbsp(left, depth + 1);
		}
		if (!right->triangles.empty())
		{
			buildbsp(right, depth + 1);
		}

		//if (left->triangles.empty() && right->triangles.empty())
		//{// leaf
		//	//PRINT("surfacedesc=" << bsp->surfacedesc);
		//	if (bsp->parent && bsp->parent->right == bsp)
		//	{
		//		convex_t convex;
		//		while (bsp->parent)
		//		{
		//			bsp = bsp->parent;
		//			convex.push_back(bsp->cuttri);
		//		}
		//		PRINT("convex=" << convex.size());
		//		convexlist.push_back(convex);
		//	}
		//}
	}
	void walkbsp(submesh& sm, bsp_t* bsp, crvec p, bool bright, real minD, int depth)
	{
		if (!bsp->left && !bsp->right)
		{
			if (bright)
			{
				/*{
					for (int i = 0; i < bsp->cutpathtriangles.size(); i ++)
					{
						triangle_t& it = bsp->cutpathtriangles[i];
						color = blendcor(0xFF00FF00, 0xFFFF0000, i / real(bsp->cutpathtriangles.size()));
						if(i == bsp->cutpathtriangles.size()-1)
							triang0(it.p[0]+vec::UY * 3, it.p[1] + vec::UY * 3, it.p[2] + vec::UY * 3);
					}
				}*/
				color = 0xFFFFFF00;
				//pyramid(p + vec3::UY * 0, 0.25);
				//PRINT("pyramid " << p.x << " " << p.y << " " << p.z << " PATH=" << bsp->path << " minD=" << minD);
			}
			return;
		}
		{
			vec3 o, n;

			const vec3& p1 = bsp->triangles[0].p[0];
			const vec3& p2 = bsp->triangles[0].p[1];
			const vec3& p3 = bsp->triangles[0].p[2];

			o = p1;
			n = (p2 - p1).cross(p3 - p1).normcopy();

			real d = check_front(o, n, p);

			minD = MIN(fabs(d), minD);

			if (fabs(d) <= c_precision)
			{
				PRINT("depth =" << depth << " d=" << d);
				color = 0xFF80aFa0;
				PRINT("!-----------triangle ");
				PRINTVEC3(p1);
				PRINTVEC3(p2);
				PRINTVEC3(p3);
				//triang0(p1 + vec::UY * 3, p2 + vec::UY * 3, p3 + vec::UY * 3);
			}

			if (bsp->right && d > c_precision)
			{
				// color = blendcor(0xFF00FF00, 0xFFFF0000, depth / 5.0f);
				// if (depth > 2)
				 //    triang0(p1 + vec::UY * 3, p2 + vec::UY * 3, p3 + vec::UY * 3);
				walkbsp(sm, bsp->right, p, true, minD, depth + 1);
			}

			if (bsp->left && d < -c_precision)
			{
				//color = blendcor(0xFFFFFFFF, 0xFF0000FF, depth / 5.0f);
			   // if(depth > 2)
			   //     triang0(p2 + vec::UY * 3, p1 + vec::UY * 3, p3 + vec::UY * 3);
				walkbsp(sm, bsp->left, p, false, minD, depth + 1);
			}
		}
	}
	void clear(bsp_t* bsp)
	{
		if (bsp->right)
		{
			clear(bsp->right);
		}
		if (bsp->left)
		{
			clear(bsp->left);
		}
		delete bsp;
	}
	//void checkbsp(submesh& sm, bsp_t* bsp)
	//{
	//	for (int i = 0; i < 5000; i++)
	//	{
	//		vec3 p = vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 20.5;

	//		//walkbsp(sm, bsp, vec(-4.61679,2.29581,2.18805), 0, 1000, 0);
	//		walkbsp(sm, bsp, p, 0, 100000, 0);
	//		// pyramid(vec3::UY * 0.5+ vec3::UX * 3.0 + vec3::UZ * 1.5, 0.1);
	//	}
	//}
	bool checkbsp(bsp_t* bsp, crvec p, bool bright = 0, int depth = 0)
	{
		if (!bsp->left && !bsp->right)
		{
			return bright;
		}
		{
			vec3 o, n;

			const vec3& p1 = bsp->triangles[0].p[0];
			const vec3& p2 = bsp->triangles[0].p[1];
			const vec3& p3 = bsp->triangles[0].p[2];

			o = p1;
			n = (p2 - p1).cross(p3 - p1).normcopy();

			real d = check_front(o, n, p);

			if (bsp->right && d > c_precision)
			{
				if (checkbsp(bsp->right, p, true, depth + 1))
					return true;
			}

			if (bsp->left && d < -c_precision)
			{
				if (checkbsp(bsp->left, p, false, depth + 1))
					return true;
			}
		}
		return false;
	}

	bool checkin(const submesh& sm, bsp_t* bsp, const boundingbox& cube)
	{
		bool ret = checkbsp(bsp, cube.a);
		ret = ret && checkbsp(bsp, cube.b);
		return ret;
	}

	void bsp_by_tris(bsp_t* bsp, std::vector<triangle_t>& triangles)
	{
		for (auto it : triangles)
			bsp->triangles.push_back(it);
		buildbsp(bsp, 0);
	}

	void addmesh(const submesh& sm, bsp_t* bsp)
	{
		for (int i = 1; i < sm.tris.size(); i++)
		{
			int t1, t2, t3;

			t1 = sm.tris[i].vertexIndex[0];
			t2 = sm.tris[i].vertexIndex[1];
			t3 = sm.tris[i].vertexIndex[2];

			const vertex& p1 = sm.vertices[t1];
			const vertex& p2 = sm.vertices[t2];
			const vertex& p3 = sm.vertices[t3];

			bsp->triangles.push_back(BSP::triangle_t(p1.p * c_scale, p2.p * c_scale, p3.p * c_scale));
		}
	}
}
