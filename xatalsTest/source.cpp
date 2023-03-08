/*
example
This example shows how to use the xatlas API to generate a unique set of texture coordinates.
Input: an .obj model file.
Output:
	* an .obj model file (example_output.obj). This is simplistic for example purposes, it doesn't copy materials from the input .obj file.
	* texture coordinates rasterized to images, colored by chart (example_charts*.tga) and by triangle (example_tris*.tga).
*/

#include "stb_image_write.h"
#include "tiny_obj_loader.h"
#include "xatlas.h"

#include <Windows.h>
#include <mutex>
#include <stdarg.h>
#include <assert.h>
#include <time.h>
#include <iostream>

#define ASSERT(_condition) if (!(_condition)) { printf("[FAILED] '%s' %s %d\n", #_condition, __FILE__, __LINE__); assert(_condition); }
#define FOPEN(_file, _filename, _mode) { if (fopen_s(&_file, _filename, _mode) != 0) _file = NULL; }

static void RandomColor(uint8_t * color)
{
	for (int i = 0; i < 3; i++)
		color[i] = uint8_t((rand() % 255 + 192) * 0.5f);
}

static void SetPixel(uint8_t* dest, int destWidth, int x, int y, const uint8_t* color)
{
	uint8_t* pixel = &dest[x * 3 + y * (destWidth * 3)];
	pixel[0] = color[0];
	pixel[1] = color[1];
	pixel[2] = color[2];
}

// https://github.com/miloyip/line/blob/master/line_bresenham.c
// License: public domain.
static void RasterizeLine(uint8_t* dest, int destWidth, const int* p1, const int* p2, const uint8_t* color)
{
	const int dx = abs(p2[0] - p1[0]), sx = p1[0] < p2[0] ? 1 : -1;
	const int dy = abs(p2[1] - p1[1]), sy = p1[1] < p2[1] ? 1 : -1;
	int err = (dx > dy ? dx : -dy) / 2;
	int current[2];
	current[0] = p1[0];
	current[1] = p1[1];
	while (SetPixel(dest, destWidth, current[0], current[1], color), current[0] != p2[0] || current[1] != p2[1])
	{
		const int e2 = err;
		if (e2 > -dx) { err -= dy; current[0] += sx; }
		if (e2 < dy) { err += dx; current[1] += sy; }
	}
}

/*
https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
Copyright Dmitry V. Sokolov
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
static void RasterizeTriangle(uint8_t* dest, int destWidth, const int* t0, const int* t1, const int* t2, const uint8_t* color)
{
	if (t0[1] > t1[1]) std::swap(t0, t1);
	if (t0[1] > t2[1]) std::swap(t0, t2);
	if (t1[1] > t2[1]) std::swap(t1, t2);
	int total_height = t2[1] - t0[1];
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1[1] - t0[1] || t1[1] == t0[1];
		int segment_height = second_half ? t2[1] - t1[1] : t1[1] - t0[1];
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1[1] - t0[1] : 0)) / segment_height;
		int A[2], B[2];
		for (int j = 0; j < 2; j++) {
			A[j] = int(t0[j] + (t2[j] - t0[j]) * alpha);
			B[j] = int(second_half ? t1[j] + (t2[j] - t1[j]) * beta : t0[j] + (t1[j] - t0[j]) * beta);
		}
		if (A[0] > B[0]) std::swap(A, B);
		for (int j = A[0]; j <= B[0]; j++)
			SetPixel(dest, destWidth, j, t0[1] + i, color);
	}
}

int main()
{
	const char* filename = "./models/outDelete.obj";
	bool useUvMesh = false;

	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	unsigned int flags = tinyobj::triangulation;

	if (!tinyobj::LoadObj(shapes, materials, err, filename, NULL, flags)) {
		printf("Error: %s\n", err.c_str());
		return -1;
	}
	if (!shapes.size()) {
		printf("Error: no shapes in obj file\n");
		return -1;
	}

	const clock_t start = clock();
	// Create empty atlas
	xatlas::Atlas* atlas = xatlas::Create();
	// Add meshes to atlas
	int missingUvsCount = 0;
	int totalVertices = 0, totalFaces = 0;
	for (int i = 0; i < (int)shapes.size(); i++) {
		const tinyobj::mesh_t& objMesh = shapes[i].mesh;
		if (useUvMesh) {

		}
		else {
			// not use uv mesh
			xatlas::MeshDecl decl;
			decl.vertexCount = (int)objMesh.positions.size() / 3;
			decl.vertexPositionData = objMesh.positions.data();
			decl.vertexPositionStride = sizeof(float) * 3;
			if (!objMesh.normals.empty()) {
				decl.vertexNormalData = objMesh.normals.data();
				decl.vertexNormalStride = sizeof(float) * 3;
			}
			if (!objMesh.texcoords.empty()) {
				decl.vertexUvData = objMesh.texcoords.data();
				decl.vertexUvStride = sizeof(float) * 2;
			}
			decl.indexCount = (int)objMesh.indices.size();
			decl.indexData = objMesh.indices.data();
			decl.indexFormat = xatlas::IndexFormat::UInt32;
			// add mesh
			xatlas::AddMeshError error = xatlas::AddMesh(atlas, decl);
			if (error != xatlas::AddMeshError::Success) {
				xatlas::Destroy(atlas);
				printf("\rError adding mesh %d '%s': %s\n", i, shapes[i].name.c_str(), xatlas::StringForEnum(error));
				return -1;
			}
			totalVertices += decl.vertexCount;
			if (decl.faceCount > 0)
				totalFaces += decl.faceCount;
			else
				totalFaces += decl.indexCount / 3; // Assume triangles if MeshDecl::faceCount not specified.
		}
	}
	xatlas::AddMeshJoin(atlas); // Not necessary. Only called here so geometry totals are printed after the AddMesh progress indicator.
	printf("   %u total vertices\n", totalVertices);
	printf("   %u total faces\n", totalFaces);

	if (missingUvsCount > 0) {
		printf("   %u/%u meshes missing UVs\n", missingUvsCount, (int)shapes.size());
	}

	// Generate atlas
	printf("Generating atlas\n");
	//xatlas::Generate(atlas);
	xatlas::ChartOptions co;
	co.straightnessWeight = 60.0;
	co.maxCost = 200;
	co.maxIterations = 10;
	co.maxBoundaryLength = 10.55;
	xatlas::ComputeCharts(atlas, co);
	xatlas::PackCharts(atlas);
	const clock_t end = clock();
	printf("   %g ms\n", (end - start) * 1000.0 / (double)CLOCKS_PER_SEC);
	printf("   %d charts\n", atlas->chartCount);
	printf("   %d atlases\n", atlas->atlasCount);
	for (uint32_t i = 0; i < atlas->atlasCount; i++)
		printf("      %d: %0.2f%% utilization\n", i, atlas->utilization[i] * 100.0f);
	printf("   %ux%u resolution\n", atlas->width, atlas->height);
	
	for (int i = 0; i < atlas->meshCount; i++) {
		const xatlas::Mesh& mesh = atlas->meshes[i];
		const tinyobj::mesh_t& objMesh = shapes[i].mesh;
		if (useUvMesh) {

		}
		else {
			// Index count shouldn't change
			ASSERT(mesh.indexCount == objMesh.indices.size());
			// Vertex count should be equal or greater
			ASSERT(mesh.vertexCount >= objMesh.positions.size() / 3);
			// Index order should be preserved
			for (int j = 0; j < mesh.indexCount; j++) {
				const xatlas::Vertex& vertex = mesh.vertexArray[mesh.indexArray[j]];
				ASSERT(vertex.xref == objMesh.indices[j]);
			}
		}
	}

	// Write meshes
	const char* outFilename = "outNewUV2D.obj";
	FILE* file;
	FOPEN(file, outFilename, "w");
	if (file) {
		int firstVertex = 0;
		for (int i = 0; i < atlas->meshCount; i++) {
			const xatlas::Mesh& mesh = atlas->meshes[i];
			for (int v = 0; v < mesh.vertexCount; v++) {
				const xatlas::Vertex& vertex = mesh.vertexArray[v];
				const float* pos = &shapes[i].mesh.positions[vertex.xref * 3];
				fprintf(file, "v %g %g %g\n", pos[0], pos[1], pos[2]);
				if (!shapes[i].mesh.normals.empty()) {
					const float* normal = &shapes[i].mesh.normals[vertex.xref * 3];
					fprintf(file, "vn %g %g %g\n", normal[0], normal[1], normal[2]);
				}
				fprintf(file, "vt %g %g\n", vertex.uv[0] / atlas->width, vertex.uv[1] / atlas->height);
			}
			fprintf(file, "o %s\n", shapes[i].name.c_str());
			fprintf(file, "s off\n");
			
			for (int f = 0; f < mesh.indexCount; f += 3) {
				fprintf(file, "f ");
				for (int j = 0; j < 3; j++) {
					const int index = firstVertex + mesh.indexArray[f + j] + 1;
					fprintf(file, "%d/%d/%d%c", index, index, index, j == 2 ? '\n' : ' ');
				}
			}
			firstVertex += mesh.vertexCount;
		}
		fclose(file);
	}

	if (atlas->width > 0 && atlas->height > 0) {
		printf("Rasterizing result...\n");
		// Dump images.
		std::vector<uint8_t> outputTrisImage, outputChartsImage;
		const uint32_t imageDataSize = atlas->width * atlas->height * 3;
		outputTrisImage.resize(atlas->atlasCount * imageDataSize);
		outputChartsImage.resize(atlas->atlasCount * imageDataSize);
		for (uint32_t i = 0; i < atlas->meshCount; i++) {
			const xatlas::Mesh& mesh = atlas->meshes[i];
			// Rasterize mesh triangles.
			const uint8_t white[] = { 255, 255, 255 };
			const uint32_t faceCount = mesh.indexCount / 3;
			uint32_t faceFirstIndex = 0;
			for (uint32_t f = 0; f < faceCount; f++) {
				int32_t atlasIndex = -1;
				int verts[255][2];
				const uint32_t faceVertexCount = 3;
				for (uint32_t j = 0; j < faceVertexCount; j++) {
					const xatlas::Vertex& v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + j]];
					atlasIndex = v.atlasIndex; // The same for every vertex in the face.
					verts[j][0] = int(v.uv[0]);
					verts[j][1] = int(v.uv[1]);
				}
				if (atlasIndex < 0)
					continue; // Skip faces that weren't atlased.
				uint8_t color[3];
				RandomColor(color);
				uint8_t* imageData = &outputTrisImage[atlasIndex * imageDataSize];

				RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);

				for (uint32_t j = 0; j < faceVertexCount; j++)
					RasterizeLine(imageData, atlas->width, verts[j], verts[(j + 1) % faceVertexCount], white);
				faceFirstIndex += faceVertexCount;
			}
			// Rasterize mesh charts.
			for (uint32_t j = 0; j < mesh.chartCount; j++) {
				const xatlas::Chart* chart = &mesh.chartArray[j];
				uint8_t color[3];
				RandomColor(color);
				for (uint32_t k = 0; k < chart->faceCount; k++) {
					const uint32_t face = chart->faceArray[k];

					const uint32_t faceVertexCount = 3;
					faceFirstIndex = face * 3;

					int verts[255][2];
					for (uint32_t l = 0; l < faceVertexCount; l++) {
						const xatlas::Vertex& v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + l]];
						verts[l][0] = int(v.uv[0]);
						verts[l][1] = int(v.uv[1]);
					}
					uint8_t* imageData = &outputChartsImage[chart->atlasIndex * imageDataSize];

					RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);

					for (uint32_t l = 0; l < faceVertexCount; l++)
						RasterizeLine(imageData, atlas->width, verts[l], verts[(l + 1) % faceVertexCount], white);
				}
			}
		}
		for (uint32_t i = 0; i < atlas->atlasCount; i++) {
			char filename[256];
			snprintf(filename, sizeof(filename), "example_tris%02u.tga", i);
			printf("Writing '%s'...\n", filename);
			stbi_write_tga(filename, atlas->width, atlas->height, 3, &outputTrisImage[i * imageDataSize]);
			snprintf(filename, sizeof(filename), "example_charts%02u.tga", i);
			printf("Writing '%s'...\n", filename);
			stbi_write_tga(filename, atlas->width, atlas->height, 3, &outputChartsImage[i * imageDataSize]);
		}
	}
	xatlas::Destroy(atlas);
	return 0;
}