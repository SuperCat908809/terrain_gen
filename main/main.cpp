#include <format>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <concepts>
#include <algorithm>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define __STDC_LIB_EXT1__
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>


__global__ void _kernel_calc_candidates(int width, int height, int* solids, int* candidates_out) {
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	if (x >= width || y >= height) return;

#define IDX(idx_x, idx_y) (width * (idx_y) + (idx_x))

	candidates_out[IDX(x, y)] = 0;
	if (solids[IDX(x, y)] == 1)return;

	if (x > 0)				if (solids[IDX(x - 1, y)] == 1) candidates_out[IDX(x, y)] = 1;
	if (x <= width - 1)		if (solids[IDX(x + 1, y)] == 1) candidates_out[IDX(x, y)] = 1;
	if (y > 0)				if (solids[IDX(x, y - 1)] == 1) candidates_out[IDX(x, y)] = 1;
	if (y <= height - 1)	if (solids[IDX(x, y + 1)] == 1) candidates_out[IDX(x, y)] = 1;
}

__global__ void _kernel_calc_neighbour_distribs(int width, int height, int* solids, int* candidates, double* weights, double* neighbours_distrib_out) {
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	if (x >= width || y >= height) return;

#define IDX(idx_x, idx_y) (width * (idx_y) + (idx_x))

	neighbours_distrib_out[IDX(x, y)] = 0.0;
	if (solids[IDX(x, y)] == 1) return;
	if (candidates[IDX(x, y)] == 1) return;

	int neighbour_count = 0;

	if (x > 0)				if (solids[IDX(x - 1, y)] == 0) neighbour_count++;
	if (x <= width - 1)		if (solids[IDX(x + 1, y)] == 0) neighbour_count++;
	if (y > 0)				if (solids[IDX(x, y - 1)] == 0) neighbour_count++;
	if (y <= height - 1)	if (solids[IDX(x, y + 1)] == 0) neighbour_count++;

	neighbours_distrib_out[IDX(x, y)] = weights[IDX(x, y)] / neighbour_count;
}

__global__ void _kernel_blur_iteration(
	int width, int height,
	int* solids,
	int* candidates,
	double* neighbour_distribs,
	double* weights_in,
	double* weights_out
) {
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	if (x >= width || y >= height) return;

#define IDX(idx_x, idx_y) (width * (idx_y) + (idx_x))

	if (solids[IDX(x, y)] == 1) return;

	double new_weight = 0.0;

	if (candidates[IDX(x, y)] == 1) new_weight += weights_in[IDX(x, y)];

	if (x > 0)				new_weight += neighbour_distribs[IDX(x-1,y)];
	if (x <= width - 1)		new_weight += neighbour_distribs[IDX(x+1,y)];
	if (y > 0)				new_weight += neighbour_distribs[IDX(x,y-1)];
	if (y <= height - 1)	new_weight += neighbour_distribs[IDX(x,y+1)];

	weights_out[IDX(x, y)] = new_weight;
}

class World {
	int width, height;
	std::vector<double> current{}, next{};
	std::vector<int> solids{};

	int getNeighbourCount(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		int count = 0;
		if (x > 0) count++;
		if (y > 0) count++;
		if (x < width - 1) count++;
		if (y < height - 1) count++;
		return count;
	}

public:

	World(int width, int height) : width(width), height(height) {
		assert(width > 0);
		assert(height > 0);
		ResetCurrent();
		ResetNext();
		ResetSolids();
	}

	void ResetCurrent() {
		current.clear();
		for (int i = 0; i < width * height; i++) {
			current.push_back(0.0);
		}
	}

	void ResetNext() {
		next.clear();
		for (int i = 0; i < width * height; i++) {
			next.push_back(0.0);
		}
	}

	void ResetSolids() {
		solids.clear();
		for (int i = 0; i < width * height; i++) {
			solids.push_back(0);
		}
	}

	int GetWidth() { return width; }
	int GetHeight() { return height; }

	double& GetCellCurrent(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		return current[width * y + x];
	}

	double& GetCellNext(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		return next[width * y + x];
	}

	int& GetCellSolid(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		return solids[width * y + x];
	}

	double GetMaxNonCandidateCell() {
		double best = 0.0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (isCandidate(x, y)) continue;
				double weight = GetCellCurrent(x, y);
				if (weight > best) best = weight;
			}
		}
		return best;
	}

	void IterateBlur() {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {

				if (GetCellSolid(x, y) == 1) continue;

				double current_weight = GetCellCurrent(x, y);
				int neighbours = getNeighbourCount(x, y);

				if (neighboursSolidCell(x, y)) {
					GetCellNext(x, y) += current_weight;
					continue;
				}

				double split_weight = current_weight / neighbours;

#if 1
				for (int dx = std::max(0, x - 1); dx <= std::min(x + 1, width - 1); dx++) {
					if (dx == x) continue;
					GetCellNext(dx, y) += split_weight;
				}
				for (int dy = std::max(0, y - 1); dy <= std::min(y + 1, height - 1); dy++) {
					if (dy == y) continue;
					GetCellNext(x, dy) += split_weight;
				}
#else
				if (x > 0) GetCellNext(x - 1, y) += split_weight;
				if (y > 0) GetCellNext(x, y - 1) += split_weight;
				if (x < width - 1) GetCellNext(x + 1, y) += split_weight;
				if (y < height - 1) GetCellNext(x, y + 1) += split_weight;
#endif
			}
		}

		std::swap(current, next);
		ResetNext();
	}

	void SelectNewSolidCell() {
		double sum = 0.0;
		double selected_weight = 0.0;
		std::pair<int, int> selected_coord = std::make_pair(-1, -1);

		int candidates = 0;

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (!isCandidate(x, y)) continue;

				candidates++;

				double candidate_weight = GetCellCurrent(x, y);
				sum += candidate_weight;
				double discard_selected_probability = sum == 0.0 ? 1.0 : candidate_weight / sum;
				double decision = rand() / (double)RAND_MAX;
				if (decision < discard_selected_probability) {
					selected_weight = candidate_weight;
					selected_coord = std::make_pair(x, y);
				}
			}
		}

		assert(std::get<0>(selected_coord) != -1);

		//std::cout << candidates << " candidates found\n";

		GetCellSolid(std::get<0>(selected_coord), std::get<1>(selected_coord)) = 1;
	}

	bool isCandidate(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		// cannot be a candidate if its already solid
		if (GetCellSolid(x, y) == 1) return false;

		if (neighboursSolidCell(x, y)) return true;

		return false;
	}

	bool neighboursSolidCell(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		if (x > 0 && GetCellSolid(x - 1, y) == 1) return true;
		if (y > 0 && GetCellSolid(x, y - 1) == 1) return true;
		if (x < width - 1 && GetCellSolid(x + 1, y) == 1) return true;
		if (y < height - 1 && GetCellSolid(x, y + 1) == 1) return true;

		return false;
	}
};

std::pair<int, int> getValidStart(World& world) {
	for (int i = 0; i < 50000; i++) {
		int x = rand() % world.GetWidth();
		int y = rand() % world.GetHeight();

		if (world.GetCellSolid(x, y) == 1) continue;
		if (world.neighboursSolidCell(x, y)) continue;

		return { x,y };
	}
	return { -1,-1 };
}

void saveWorldToImage(World& world, std::string path) {
	std::vector<uint8_t> image_data{};
	image_data.reserve(world.GetWidth() * world.GetHeight() * 3);

	double max = world.GetMaxNonCandidateCell();

	for (int y = 0; y < world.GetHeight(); y++) {
		for (int x = 0; x < world.GetWidth(); x++) {
			
			if (world.GetCellSolid(x, y) == 1) {
				image_data.push_back(240);
				image_data.push_back(240);
				image_data.push_back(240);
				continue;
			}

			float val = world.GetCellCurrent(x, y);
			//val = log10f(val);
			//val = std::clamp((val + 60) / (0 - -60), 0.0f, 1.0f);
			val /= max;
			val = std::clamp(val, 0.0f, 1.0f);

			float col = (1 - val) * 0.02f + val * 0.95f; // val c [0, 1], lerp between almost black and almost white
			col = sqrtf(col); // tone mapping

			image_data.push_back(static_cast<uint8_t>(col * 255.99f));
			image_data.push_back(static_cast<uint8_t>(col * 63.99f));
			image_data.push_back(static_cast<uint8_t>(col * 255.99f));
		}
	}

	stbi_flip_vertically_on_write(true);
	stbi_write_jpg(path.c_str(), world.GetWidth(), world.GetHeight(), 3, image_data.data(), 95);
}



int main() {

	int world_width = 80;
	int world_height = world_width;
	int blur_iterations = std::max(world_width, world_height) * 2;
	int particle_additions = std::min(world_width, world_height) * 10;

	World world(world_width, world_height);

	// growth seed
	world.GetCellSolid(world_width / 2, world_height / 2) = 1;

	srand(0);

	for (int particle_iter = 0; particle_iter < particle_additions; particle_iter++) {

		auto start = getValidStart(world);
		if (std::get<0>(start) == -1) {
			std::cout << "\nCould not find valid starting point\n\n";
			exit(-1);
		}
		world.GetCellCurrent(std::get<0>(start), std::get<1>(start)) = 1.0;

		for (int i = 0; i < blur_iterations; i++) {

			world.IterateBlur();

			//std::cout << "Particle iter: " << particle_iter + 1 << ". Completed iter " << i << "\n";
		}

		world.SelectNewSolidCell();

		std::stringstream ss1{};
		ss1 << std::setw(6) << std::setfill('0') << particle_iter;
		std::string path1 = "sim\\waves " + ss1.str() + ".jpg";
		saveWorldToImage(world, path1);

		std::cout << "Particle iter " << particle_iter + 1 << "\n";

		world.ResetCurrent();

		//std::stringstream ss{};
		//ss << std::setw(6) << std::setfill('0') << particle_iter;
		//std::string path = std::format("sim\\iter {}.jpg", ss.str());
		//saveWorldToImage(world, path);
	}

	saveWorldToImage(world, "final_render.jpg");

	std::cout << "\n\nFinished\n";
	return 0;
}