#include <iostream>
#include <sstream>

#include "gaussianQuadrature.hpp"
#include "RadModel.hpp"

using namespace std;

namespace myLib
{
    void getWeights(RadModel &radModel)
    {
        const int nQuad = 2 * radModel.params.nQuad;

        switch(nQuad)
        {
            case 2:
            {
                const double inv_sqrt_3 = 1.0 / sqrt(3.0);
                radModel.quadMu = { inv_sqrt_3 };
                radModel.quadW = { 1.0 };
                break;
            }
            case 4:
            {
                radModel.quadMu = { 0.8611363115940526, 0.3399810435848563 };
                radModel.quadW = { 0.3478548451374538, 0.6521451548625461 };
                break;
            }
            case 6:
            {
                radModel.quadMu = { 0.9324695142031521, 0.6612093864662645, 0.2386191860831969 };
                radModel.quadW = { 0.1713244923791704, 0.3607615730481386, 0.4679139345726910 };
                break;
            }
            case 8:
            {
                radModel.quadMu = { 0.9602898564975363, 0.7966664774136267, 0.5255324099163290, 0.1834346424956498 };
                radModel.quadW = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620 };
                break;
            }
            case 32:
            {   // point datat taken from https://pomax.github.io/bezierinfo/legendre-gauss.html
                radModel.quadMu = { 0.997263862, 0.985611512, 0.964762256, 0.934906076, 0.896321156, 0.849367614, 0.794483796, 0.732182119,
                                    0.663044267, 0.587715757, 0.506899909, 0.421351276, 0.331868602, 0.239287362, 0.144471962, 0.048307666 };
                radModel.quadW = { 0.007018610, 0.016274395, 0.025392065, 0.034273863, 0.042835898, 0.050998059, 0.058684093, 0.065822223,
                                   0.072345794, 0.078193896, 0.083311924, 0.087652093, 0.091173879, 0.093844399, 0.095638720, 0.096540089 };
                break;
            }
            default:
                ostringstream output;
                output << "nQuad = " << nQuad << " is not valid";
                radModel.log(output);
        }
    }
}