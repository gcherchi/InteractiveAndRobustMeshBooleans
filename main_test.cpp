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

    const bigrational s0[] = {bigrational(bignatural(static_cast<uint32_t>(28954682368)), bignatural(static_cast<uint32_t> (3)), 1),
                        bigrational(bignatural(static_cast<uint32_t>(18845212751614903)), bignatural(static_cast<uint32_t> (36291456)), -1),
                        bigrational(bignatural(static_cast<uint32_t>(1242883924321459)), bignatural(static_cast<uint32_t> (2097152)), -1)
    };

    const bigrational s1[] = {bigrational(bignatural(static_cast<uint32_t>(28954718211)), bignatural(static_cast<uint32_t> (2)), 1),
                        bigrational(bignatural(static_cast<uint32_t>(18845212751614903)), bignatural(static_cast<uint32_t> (36291456)), -1),
                        bigrational(bignatural(static_cast<uint32_t>(1242883924321459)), bignatural(static_cast<uint32_t> (2097152)), -1)
    };


    std::cout << "Values" << std::endl;
    std::cout << "s0 x y z coords: " << s0[0] << " " << s0[1] << " " << s0[2] << std::endl;
    std::cout << "s1 x y z coords: " << s1[0] << " " << s1[1] << " " << s1[2] << std::endl;
    std::cout << std::endl;

    std::cout << "Are equals? " << std::endl;
    const char *x_equals = s0[0] == s1[0] ? "yes" : "no";
    std::cout << "s0_x == s1_x: " << x_equals<< std::endl;

    const char *y_equals = s0[1] == s1[1] ? "yes" : "no";
    std::cout << "s0_y == s1_y: " << y_equals<< std::endl;

    const char *z_equals = s0[2] == s1[2] ? "yes" : "no";

    std::cout << "s0_z == s1_z: " << z_equals<< std::endl;

    std::cout << std::endl;
    const bigrational distance = (s1[0] - s0[0]) * (s1[0] - s0[0]) +
                              (s1[1] - s0[1]) * (s1[1] - s0[1]) +
                              (s1[2] - s0[2]) * (s1[2] - s0[2]);
    std::cout << std::endl;
    std::cout << "Distance: " << distance << std::endl;

    return 0;
}