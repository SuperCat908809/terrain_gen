#include <format>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <concepts>
#include <algorithm>

#define __STDC_LIB_EXT1__
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>

class World {
	int width, height;
	std::vector<double> current{}, next{};
	std::vector<int> solids{};

#if 0
	template <typename Func, typename... Args> requires std::invocable<Func, Args...>
	void iterateOverNeighbours(int x, int y, Args... args) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		for (int dx = std::max(0, x - 1); dx <= std::min(x + 1, width - 1); dx++) {
			if (dx == x) continue;
			Func(args...);
		}
		for (int dy = std::max(0, y - 1); dy <= std::min(y + 1, height - 1); dy++) {
			if (dy == y) continue;
			Func(args...);
		}
	}
#endif

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

#if 0
	std::vector<std::pair<int, int>> getNeighbourCoords(int x, int y) {
		assert(0 <= x && x < width);
		assert(0 <= y && y < height);

		std::vector<std::pair<int, int>> list{};

		for (int dx = std::max(0, x - 1); dx <= std::min(x + 1, width - 1); dx++) {
			if (dx == x) continue;
			list.push_back(std::make_pair(dx, y));
		}
		for (int dy = std::max(0, y - 1); dy <= std::min(y + 1, height - 1); dy++) {
			if (dy == y) continue;
			list.push_back(std::make_pair(x, dy));
		}

		return list;
	}
#endif

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

#if 0
				bool skip = false;
				for (auto& coord : neighbours) {
					int neighbour_solid = GetCellSolid(std::get<0>(coord), std::get<1>(coord));
					if (neighbour_solid == 1) {
						GetCellNext(x, y) += current_weight;
						skip = true;
						break;
					}
				}
				if (skip) continue;
#else
				if (neighboursSolidCell(x, y)) {
					GetCellNext(x, y) += current_weight;
					continue;
				}
#endif

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
#if 0
		std::vector<std::pair<std::pair<int, int>, double>> candidates{};

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (isCandidate(x, y)) {
					candidates.push_back({ {x, y}, GetCellCurrent(x, y) });
				}
			}
		}

		// select candidate using reservoir sampling
		double sum = 0.0;
		double selected_weight = 0.0;
		int selected_index = -1;

		int i = 0;
		for (auto candidate : candidates) {
			sum += std::get<1>(candidate);

			if (rand() / (double)RAND_MAX > selected_weight / sum) {
				selected_index = i;
				selected_weight = std::get<1>(candidate);
			}

			i++;
		}

		auto selected_coord = std::get<0>(candidates[selected_index]);

		GetCellSolid(std::get<0>(selected_coord), std::get<1>(selected_coord)) = 1;

		ResetCurrent();
#else
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
#endif
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

#if 0
		auto neighbour_coords = getNeighbourCoords(x, y);
		for (auto coord : neighbour_coords) {
			if (GetCellSolid(std::get<0>(coord), std::get<1>(coord)) == 1) return true;
		}
#else
		if (x > 0 && GetCellSolid(x - 1, y) == 1) return true;
		if (y > 0 && GetCellSolid(x, y - 1) == 1) return true;
		if (x < width - 1 && GetCellSolid(x + 1, y) == 1) return true;
		if (y < height - 1 && GetCellSolid(x, y + 1) == 1) return true;
#endif

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