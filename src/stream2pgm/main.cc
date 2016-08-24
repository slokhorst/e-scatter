/**
 * @file src/viewer/stream2pgm.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include <algorithm>
#include <cinttypes>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

int main(const int argc, char* argv[]) {
    bool auto_gray_scale = false;
    bool secondary = true;
    bool backscatter = false;

    std::map<int,std::map<int,int>> image_map;
    int min_y, max_y;
    int min_x, max_x;
    int min_c, max_c;

    min_x = max_x = 0;
    min_y = max_y = 0;
    min_c = max_c = 0;
    while(true) {
        struct {
            float rx, ry, rz;
            float dx, dy, dz;
            float K;
            int px, py;
        } record;

        std::cin.read(reinterpret_cast<char*>(&record), sizeof(record));
        if(std::cin.eof())
            break;

        if((secondary && (record.K <= 50)) || (backscatter && (record.K > 50))) {
            const int pixel_cnt = ++(image_map[record.py][record.px]);
            if(max_c == 0) {
                min_y = max_y = record.py;
                min_x = max_x = record.px;
                min_c = max_c = pixel_cnt;
            } else {
                min_y = std::min(min_y, record.py);
                max_y = std::max(max_y, record.py);
                min_x = std::min(min_x, record.px);
                max_x = std::max(max_x, record.px);
                min_c = std::min(min_c, pixel_cnt);
                max_c = std::max(max_c, pixel_cnt);
            }
        }
    }

    if(image_map.empty())
        return EXIT_SUCCESS;

    if(!auto_gray_scale) {
        min_c = 0;
        max_c = 65535;
    }

    const int width = (1+max_x-min_x);
    const int height = (1+max_y-min_y);
    std::cout << "P5" << " " << width << " " << height << " " << 65535 << std::endl;
    for(int y = min_y; y <= max_y; y++) {
        std::vector<uint16_t> row_vec(width, 0);
        if(image_map.count(y)) {
            auto row_map = image_map[y];
            for(int x = min_x; x <= max_x; x++) {
                auto cit = row_map.find(x);
                if(cit != row_map.cend()) {
                    uint16_t z = 65535.0*(cit->second-min_c)/(max_c-min_c);
                    row_vec[x-min_x] = (z<<8) | (z>>8);
                }
            }
        }
        std::cout.write(reinterpret_cast<const char*>(row_vec.data()), sizeof(uint16_t)*row_vec.size());
    }

    return EXIT_SUCCESS;
}
