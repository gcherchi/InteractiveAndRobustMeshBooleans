/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2022 G. Cherchi, F. Pellacini, M. Attene and M. Livesu                  *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://www.gianmarcocherchi.com                                                 *
 *                                                                                       *
 *      Fabio Pellacini (fabio.pellacini@uniroma1.it)                                    *
 *      https://pellacini.di.uniroma1.it                                                 *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 * ***************************************************************************************/

#ifdef _MSC_VER // Workaround for known bugs and issues on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#define NOMINMAX // https://stackoverflow.com/questions/1825904/error-c2589-on-stdnumeric-limitsdoublemin
#endif

#include "intersect_custom.h"
#include "numerics.h"
#include <vector>
#include <random>

int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> num_dist(-10000000000, 10000000000);
    std::uniform_int_distribution<int64_t> denom_dist(1024, 1048576); // Un numero di denominatore abbastanza grande


    for (int i = 0; i < 1000; ++i) {
        std::vector<bigrational> ray_v0 = {
            bigrational(num_dist(gen), denom_dist(gen), (i % 2 == 0) ? 1 : -1),
            bigrational(num_dist(gen), denom_dist(gen), 1),
            bigrational(num_dist(gen), denom_dist(gen), 1)
        };

        std::vector<bigrational> ray_v1 = {
            bigrational(num_dist(gen)+12, denom_dist(gen), 1),
            bigrational(num_dist(gen)+1, denom_dist(gen), -1),
            bigrational(num_dist(gen)+1, denom_dist(gen), 1)
        };

        std::vector<bigrational> tv0 = {
            bigrational(num_dist(gen)+9, denom_dist(gen), (i % 2 == 0) ? -1 : 1),
            bigrational(num_dist(gen), denom_dist(gen), 1),
            bigrational(0, 1, 0) // Zero as a bigrational
        };

        std::vector<bigrational> tv1 = {
            bigrational(num_dist(gen)+3, 1, -1),
            bigrational(num_dist(gen), denom_dist(gen), -1),
            bigrational(num_dist(gen), 1, 0) // Zero as a bigrational
        };

        std::vector<bigrational> tv2 = {
            bigrational(num_dist(gen)+5, denom_dist(gen), -1),
            bigrational(num_dist(gen), denom_dist(gen), -1),
            bigrational(num_dist(gen), denom_dist(gen), 1)
        };

        // Chiamata alla funzione di intersezione
        int intersection = segment_triangle_intersect_3d(&ray_v0[0], &ray_v1[0], &tv0[0], &tv1[0], &tv2[0]);

        std::cout << "Esempio " << i + 1 << " - Intersezione: " << intersection << std::endl;
    }

    /*std::vector<bigrational> ray_v0 = {
        bigrational(1404007700713777, 1048576, -1),
        bigrational(18861893487848237, 402653184,1),
        bigrational(1655961844303259, 8388608, 1)
    };

    std::vector<bigrational> ray_v1 = {
        bigrational(8442799307783729, 4194304, 1),
        bigrational(18861893487848237, 402653184, 1),
        bigrational(1655961844303259, 8388608, 1)
    };

    std::vector<bigrational> tv0 = {
        bigrational(699137320401585, 524288, -1),
        bigrational(2370548114715617, 8388608, 1),
        bigrational(0, 1, 0) // Zero as a bigrational
    };

    std::vector<bigrational> tv1 = {
        bigrational(1338875904, 1, -1),
        bigrational(6313650243502081, 67108864, -1),
        bigrational(0, 1, 0) // Zero as a bigrational
    };

    std::vector<bigrational> tv2 = {
        bigrational(1409835321425457, 1048576, -1),
        bigrational(6439575860597473, 134217728, -1),
        bigrational(4967885532909777, 8388608,1)
    };


    int intersection = segment_triangle_intersect_3d(&ray_v0[0], &ray_v1[0], &tv0[0], &tv1[0], &tv2[0]);

    std::cout << "intersection: " << intersection << std::endl;
    */

    return 0;
}