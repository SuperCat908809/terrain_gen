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


const int world_width = 200;
const int world_height = 200;
const int blur_iterations = 2000;
const int particle_additions = 6000;

#define IDX(idx_x, idx_y) (world_width * (idx_y) + (idx_x))
#define CEIL_DIV(num, denom) (((num) + (denom) - 1) / (denom))

const dim3 threads = dim3(8, 8);
const dim3 blocks = dim3(CEIL_DIV(world_width, threads.x), CEIL_DIV(world_height, threads.y));


#if 0
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
#else

__global__ void _kernel_calc_candidates(
	int width, int height,
	int* solids,
	int* candidates_out
);
std::pair<int, int> findValidWeightCoord(std::vector<int>solids, std::vector<int>candidates);
__global__ void _kernel_calc_neighbour_distribs(
	int width, int height,
	int* solids,
	int* candidates,
	double* neighbours_counts
);
__global__ void _kernel_blur_iteration(
	const int width, int height,
	const int* solids,
	const int* candidates,
	const double* neighbour_distribs,
	const double* weights_in,
	double* weights_out
);

void write_image_to_file(int width, int height, int channels, const std::vector<uint8_t>& image, std::string path);
void write_image_to_file(int width, int height, int channels, const std::vector<float>& image, std::string path);

int main() {

	int* d_solids{};
	int* d_candidates{};
	double* d_neighbour_distribs{};
	double* d_weights_src{};
	double* d_weights_dst{};

	assert(cudaSuccess == cudaMalloc((void**)&d_solids, sizeof(int) * world_width * world_height));
	assert(cudaSuccess == cudaMalloc((void**)&d_candidates, sizeof(int) * world_width * world_height));
	assert(cudaSuccess == cudaMalloc((void**)&d_neighbour_distribs, sizeof(double) * world_width * world_height));
	assert(cudaSuccess == cudaMalloc((void**)&d_weights_src, sizeof(double) * world_width * world_height));
	assert(cudaSuccess == cudaMalloc((void**)&d_weights_dst, sizeof(double) * world_width * world_height));

	{
		std::vector<int> solids(world_width * world_height);
		std::fill(solids.begin(), solids.end(), 0);
		solids[IDX(world_width / 2, world_height / 2)] = 1;

		assert(cudaSuccess == cudaMemcpy(d_solids, solids.data(), sizeof(int) * solids.size(), cudaMemcpyHostToDevice));
	}

	srand(0);

	for (int particle_iter = 0; particle_iter < particle_additions; particle_iter++)
	{
		_kernel_calc_candidates<<<blocks, threads>>>(world_width, world_height, d_solids, d_candidates);
		assert(cudaSuccess == cudaGetLastError());
		assert(cudaSuccess == cudaDeviceSynchronize());

		_kernel_calc_neighbour_distribs<<<blocks, threads>>>(world_width, world_height, d_solids, d_candidates, d_neighbour_distribs);
		assert(cudaSuccess == cudaGetLastError());
		assert(cudaSuccess == cudaDeviceSynchronize());

		{
			std::vector<int> solids(world_width * world_height);
			assert(cudaSuccess == cudaMemcpy(solids.data(), d_solids, sizeof(int) * solids.size(), cudaMemcpyDeviceToHost));

			std::vector<int> candidates(world_width * world_height);
			assert(cudaSuccess == cudaMemcpy(candidates.data(), d_candidates, sizeof(int) * candidates.size(), cudaMemcpyDeviceToHost));

			std::pair<int, int> weight_coord = findValidWeightCoord(solids, candidates);

			if (std::get<0>(weight_coord) == -1) {
				std::cout << "Couldn't find valid weight starting point\n";
				exit(-1);
			}

			std::vector<double> weights(world_width * world_height);
			std::fill(weights.begin(), weights.end(), 0.0);
			weights[IDX(std::get<0>(weight_coord), std::get<1>(weight_coord))] = 1.0;

			assert(cudaSuccess == cudaMemcpy(d_weights_src, weights.data(), sizeof(double) * weights.size(), cudaMemcpyHostToDevice));
		}

		for (int blur_iter = 0; blur_iter < blur_iterations; blur_iter++) {
			_kernel_blur_iteration<<<blocks, threads>>>(world_width, world_height, d_solids, d_candidates, d_neighbour_distribs, d_weights_src, d_weights_dst);

			std::swap(d_weights_src, d_weights_dst);
		}

		assert(cudaSuccess == cudaGetLastError());
		assert(cudaSuccess == cudaDeviceSynchronize());

		{
			std::vector<double> weights(world_width * world_height);
			assert(cudaSuccess == cudaMemcpy(weights.data(), d_weights_src, sizeof(double) * weights.size(), cudaMemcpyDeviceToHost));

			std::vector<int> candidates(world_width * world_height);
			assert(cudaSuccess == cudaMemcpy(candidates.data(), d_candidates, sizeof(int) * candidates.size(), cudaMemcpyDeviceToHost));

			double sum = 0.0;
			double selected_weight = 0.0;
			auto selected_coord = std::make_pair(-1, -1);

			for (int y = 0; y < world_height; y++) {
				for (int x = 0; x < world_width; x++) {
					if (candidates[IDX(x, y)] != 1) continue;

					double candidate_weight = weights[IDX(x, y)];
					sum += candidate_weight;
					
					double select_candidate_probability = candidate_weight / sum;
					double decision = rand() / (double)RAND_MAX;

					if (decision < select_candidate_probability) {
						selected_weight = candidate_weight;
						selected_coord = std::make_pair(x, y);
					}
				}
			}

			assert(std::get<0>(selected_coord) != -1);

			std::vector<int> solids(world_width * world_height);
			assert(cudaSuccess == cudaMemcpy(solids.data(), d_solids, sizeof(int) * solids.size(), cudaMemcpyDeviceToHost));

			solids[IDX(std::get<0>(selected_coord), std::get<1>(selected_coord))] = 1;


			assert(cudaSuccess == cudaMemcpy(weights.data(), d_weights_src, sizeof(double) * world_width * world_height, cudaMemcpyDeviceToHost));
			std::vector<float> byte_image{};
			double max = *std::max_element(weights.begin(), weights.end());
			for (int i = 0; i < world_width * world_height; i++) {
				double weight = weights[i];
				int solid = solids[i];

				if (solid == 1) {
					byte_image.push_back(1.0);
					byte_image.push_back(1.0);
					byte_image.push_back(1.0);
					continue;
				}

				double factor = std::clamp(weight / max, 0.0, 1.0);
				factor = (1 - factor) * 0.00 + factor * 1.00;

				byte_image.push_back(factor * 0.95);
				byte_image.push_back(factor * 0.25);
				byte_image.push_back(factor * 0.95);
			}
			std::stringstream ss1{};
			ss1 << "sim\\iter " << std::setw(6) << std::setfill('0') << particle_iter + 1 << ".png";
			write_image_to_file(world_width, world_height, 3, byte_image, ss1.str());

			assert(cudaSuccess == cudaMemcpy(d_solids, solids.data(), sizeof(int) * solids.size(), cudaMemcpyHostToDevice));

			std::cout << "Finished particle " << particle_iter + 1 << "\n";
		}
	}

	assert(cudaSuccess == cudaFree(d_solids));
	assert(cudaSuccess == cudaFree(d_candidates));
	assert(cudaSuccess == cudaFree(d_neighbour_distribs));
	assert(cudaSuccess == cudaFree(d_weights_src));
	assert(cudaSuccess == cudaFree(d_weights_dst));
}

__global__ void _kernel_calc_candidates(
	int width, int height,
	int* solids,
	int* candidates_out
) {
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	if (x >= width || y >= height) return;

#define IDX(idx_x, idx_y) (width * (idx_y) + (idx_x))

	candidates_out[IDX(x, y)] = 0;
	if (solids[IDX(x, y)] == 1) return;

	if (x > 0)				if (solids[IDX(x - 1, y)] == 1) candidates_out[IDX(x, y)] = 1;
	if (x < width - 1)		if (solids[IDX(x + 1, y)] == 1) candidates_out[IDX(x, y)] = 1;
	if (y > 0)				if (solids[IDX(x, y - 1)] == 1) candidates_out[IDX(x, y)] = 1;
	if (y < height - 1)		if (solids[IDX(x, y + 1)] == 1) candidates_out[IDX(x, y)] = 1;
}

std::pair<int, int> findValidWeightCoord(std::vector<int>solids, std::vector<int>candidates) {
#define IDX(idx_x, idx_y) (world_width * (idx_y) + (idx_x))

	for (int i = 0; i < 50000; i++) {
		auto coord = std::make_pair(rand() % world_width, rand() % world_height);
		if (solids[IDX(std::get<0>(coord), std::get<1>(coord))] == 1) continue;
		if (candidates[IDX(std::get<0>(coord), std::get<1>(coord))] == 1) continue;
		return coord;
	}

	return{ -1,-1 };
}

__global__ void _kernel_calc_neighbour_distribs(
	int width, int height,
	int* solids,
	int* candidates,
	double* neighbours_distrib_out
) {
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	if (x >= width || y >= height) return;

#define IDX(idx_x, idx_y) (width * (idx_y) + (idx_x))

	neighbours_distrib_out[IDX(x, y)] = 0.0;
	if (solids[IDX(x, y)] == 1) return;
	if (candidates[IDX(x, y)] == 1) return;

	int neighbour_count = 0;

	if (x > 0)			neighbour_count++;
	if (x < width - 1)	neighbour_count++;
	if (y > 0)			neighbour_count++;
	if (y < height - 1)	neighbour_count++;

	neighbours_distrib_out[IDX(x, y)] = 1.0 / (double)neighbour_count;
}

__global__ void _kernel_blur_iteration(
	const int width, int height,
	const int* solids,
	const int* candidates,
	const double* neighbour_distribs,
	const double* weights_in,
	double* weights_out
) {
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	if (x >= width || y >= height) return;

#define IDX(idx_x, idx_y) (width * (idx_y) + (idx_x))

	if (solids[IDX(x, y)] == 1) {
		weights_out[IDX(x, y)] = 0.0;
		return;
	}

	double new_weight = 0.0;

	if (candidates[IDX(x, y)] == 1) new_weight += weights_in[IDX(x, y)];

	if (x > 0)				new_weight += weights_in[IDX(x - 1, y)] * neighbour_distribs[IDX(x - 1, y)];
	if (x < width - 1)		new_weight += weights_in[IDX(x + 1, y)] * neighbour_distribs[IDX(x + 1, y)];
	if (y > 0)				new_weight += weights_in[IDX(x, y - 1)] * neighbour_distribs[IDX(x, y - 1)];
	if (y < height - 1)		new_weight += weights_in[IDX(x, y + 1)] * neighbour_distribs[IDX(x, y + 1)];

	weights_out[IDX(x, y)] = new_weight;
}

void write_image_to_file(int width, int height, int channels, const std::vector<float>& image, std::string path) {
	std::vector<uint8_t> byte_image(width * height * channels);
	std::transform(image.begin(), image.end(), byte_image.begin(), [](float f) { return static_cast<uint8_t>(f * 255.99f); });
	write_image_to_file(width, height, channels, byte_image, path);
}
void write_image_to_file(int width, int height, int channels, const std::vector<uint8_t>& image, std::string path) {
	stbi_flip_vertically_on_write(true);
	stbi_write_png(path.c_str(), width, height, channels, image.data(), width * channels);
}

#endif